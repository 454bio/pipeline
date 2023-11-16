#include <stdio.h>
#include <stdlib.h>
#include "Model4.h"
#include <string>
#include <map>
#include "EditDist.h"
#include <math.h>
#include <sstream>

int verbose = 0;

struct PhaseParams {
    double ie;
    double cf;
    double dr;
    double err;
};

struct SpotData {
    std::vector<Signal> vals;
    std::string spotName;
};

// default grid search params
double drmin = 0.005;
double drmax = 0.025;
int    drnum = 5;
double iemin = 0.12;
double iemax = 0.21;
int    ienum = 19;
double cfmin = 0.07;
double cfmax = 0.12;
int    cfnum = 10;

void template2bases(char *dnaTemplate, char *basecalls)
{
    char bases[] = {'G', 'C', 'A', 'T'};
    int len = strlen(dnaTemplate);
    int i = 0;
    for(i=0;i<len;i++)
        basecalls[i] = bases[dnaTemplate[i]-1];
    basecalls[i] = 0;
}

std::vector<std::string> Parse(char *str, char delim)
{
    std::vector<std::string> tokens;
    char *ptr = str;
    while (ptr && *ptr != 0) {
        char *delimPtr = strchr(ptr, delim);
        if (delimPtr) {
            *delimPtr = 0;
            delimPtr++;
        }
        tokens.push_back(ptr);
        ptr = delimPtr;
    }

    return tokens;
}

void LoadSpotData(const char *filename, std::vector<SpotData> &spotData)
{
    std::vector<Signal> vals;
    int spotId;
    char name[256];
    int cycle;
    Signal s;
    int lastSpotId = 0;
    std::string spotName;

    FILE *fp = fopen(filename, "r");
    if (fp) {
        char line[1024];
        // read and ignore the header
        fgets(line, sizeof(line), fp);
        while (fgets(line, sizeof(line), fp)) {
            std::vector<std::string> tokens = Parse(line, ',');
            spotId = atoi(tokens[0].c_str());
            // if spot id is different, save the last spot data and clear things out
            if (spotId != lastSpotId) {
                if (lastSpotId != 0) {
                    SpotData d;
                    d.vals = vals;
                    d.spotName = spotName;
                    spotData.push_back(d);
                }
                lastSpotId = spotId;
                vals.clear();
            }
            // now read the new spot data and store
            spotName = tokens[1];
            for(int i=0;i<4;i++)
                s.v[i] = atof(tokens[3+i].c_str());
            vals.push_back(s);
        }
        fclose(fp);

        // save the last spot and vals
        SpotData d;
        d.vals = vals;
        d.spotName = spotName;
        spotData.push_back(d);
    }
}

void NormalizeSpotData(std::vector<SpotData> &spotData)
{
    int numSpots = spotData.size();
    for(int i=0;i<numSpots;i++) {
        int numCycles = spotData[i].vals.size();
        Signal smin;
        Signal smax;
        for (int c=0;c<numCycles;c++) {
            for(int b=0;b<4;b++) {
                if (c==0 || spotData[i].vals[c].v[b] < smin.v[b])
                    smin.v[b] = spotData[i].vals[c].v[b];
                if (c==0 || spotData[i].vals[c].v[b] > smax.v[b])
                    smax.v[b] = spotData[i].vals[c].v[b];
            }
        }
        for (int c=0;c<numCycles;c++) {
            for(int b=0;b<4;b++) {
                spotData[i].vals[c].v[b] -= smin.v[b];
                if ((smax.v[b] - smin.v[b]) > 0)
                    spotData[i].vals[c].v[b] /= (smax.v[b] - smin.v[b]);
            }
        }
    }
}

double CallBases(char *dnaTemplate, std::vector<Signal> &measuredSignal, std::vector<double> &errorPerCycle, double ie, double cf, double dr, int maxNumCycles=0)
{
    int numCycles = measuredSignal.size();
    if (maxNumCycles > 0 and maxNumCycles < numCycles)
        numCycles = maxNumCycles;

    Model4 m;
    m.Init(100);
    m.SetParams(ie, cf, dr);

    std::vector<Signal> dyeIntensities(numCycles);
    std::vector<double> totalSignal(numCycles, 0);

    memset(dnaTemplate, 0, 1024); // TODO really need to include size
    char bases[4] = {'A', 'C', 'G', 'T'};
    double cumulativeError;

    int numIterations = 3;
    for(int iteration=0;iteration<numIterations;iteration++) {
        m.Reset();
        cumulativeError = 0.0;
        for(int cycle=0;cycle<numCycles;cycle++) {
            double best_error = 0.0;
            int best_base = -1;
            Signal best_signal;
            for(int base=0;base<4;base++) {
                // insert a "what-if" base at the current position, and predict what the signal looks like
                dnaTemplate[cycle] = base+1;
                Signal signal = m.GetSignal(dnaTemplate);
                double signalSum = 0.0;
                for(int i=0;i<4;i++)
                    signalSum += signal.v[i];

                double error = 0.0;
                for(int i=0;i<4;i++) {
                    double delta = (measuredSignal[cycle].v[i] - signal.v[i])/signalSum;
                    error += delta*delta;
                }

                // keep track of the lowest error, this is the best predition
                if (error < best_error || best_base == -1) {
                    best_base = base;
                    best_error = error;
                    best_signal = signal;
                }
            }

            // append/replace with best base at current position (cycle)
            dnaTemplate[cycle] = best_base+1;
            for(int i=0;i<5;i++)
                dyeIntensities[cycle].v[i] = best_signal.v[i];
            totalSignal[cycle] = 0.0;
            for(int i=0;i<4;i++)
                totalSignal[cycle] += best_signal.v[i];
            errorPerCycle[cycle] = best_error;

            // update the model - note that we do this after getting the measured signals, because this matches the physical
            // system where the first base is exposed to nucleotides prior to UV cleavage
            m.ApplyUV(dnaTemplate, numCycles);

            cumulativeError += errorPerCycle[cycle];
        }
    }

    return cumulativeError/numCycles;
}

PhaseParams GridSearch(char *dnaTemplate, std::vector<Signal> &measuredSignal, std::vector<double> &errorPerCycle, int maxNumCycles=0)
{
    // grid-search first, then call based on lowest error

    int numCycles = measuredSignal.size();
    PhaseParams params;
    double minerr = 99999.0;

    for(int dri=0;dri<drnum;dri++) {
        double drtest = drmin + (dri/(double)(drnum-1.0)) * (drmax-drmin);
        for(int cfi=0;cfi<cfnum;cfi++) {
            double cftest = cfmin + (cfi/(double)(cfnum-1.0)) * (cfmax-cfmin);
            for(int iei=0;iei<ienum;iei++) {
                double ietest = iemin + (iei/(double)(ienum-1.0)) * (iemax-iemin);
                double err = CallBases(dnaTemplate, measuredSignal, errorPerCycle, ietest, cftest, drtest, maxNumCycles);
                if (err < minerr) {
                    minerr = err;
                    params.ie = ietest;
                    params.cf = cftest;
                    params.dr = drtest;
                }
            }
        }
    }
    params.err = CallBases(dnaTemplate, measuredSignal, errorPerCycle, params.ie, params.cf, params.dr);
    return params;
}

int main(int argc, char *argv[])
{
    const char *spotFile = "color_transformed_spots.csv";
    bool gridsearch = false;
    double ie = 0.0, cf = 0.0, dr = 0.0;
    const char *fastQFileName = "out.fastq";
    bool wantNormalize = false;
    int maxNumCycles = 0;

    int argcc = 1;
    while (argcc < argc) {
        switch (argv[argcc][1]) {
            case 'i':
                spotFile = argv[++argcc];
            break;

            case 'o':
                fastQFileName = argv[++argcc];
            break;

            case 'g':
                gridsearch = true;
            break;

            case 'G': // set the grid search params (iemin iemax ienum cfmin cfmax cfnum drmin drmax drnum)
                iemin = atof(argv[++argcc]);
                iemax = atof(argv[++argcc]);
                ienum = atoi(argv[++argcc]);

                cfmin = atof(argv[++argcc]);
                cfmax = atof(argv[++argcc]);
                cfnum = atoi(argv[++argcc]);

                drmin = atof(argv[++argcc]);
                drmax = atof(argv[++argcc]);
                drnum = atoi(argv[++argcc]);
            break;

            case 'p': // phase params (ie cf dr)
                ie = atof(argv[++argcc]);
                cf = atof(argv[++argcc]);
                dr = atof(argv[++argcc]);
            break;

            case 'c':
                maxNumCycles = atoi(argv[++argcc]);
            break;

            case 'n':
                wantNormalize = true;
            break;

            case 'v':
                verbose++;
        }
        argcc++;
    }

    // load spots
    std::vector<SpotData> spotData;
    LoadSpotData(spotFile, spotData);
    int numSpots = spotData.size();
    printf("Loaded %d spots\n", numSpots);
    int numCycles = spotData[0].vals.size();

    // normalize
    if (wantNormalize)
        NormalizeSpotData(spotData);

    // grid-search phase-correct each spot
    int numReadsAll = 0; // needs to be min 4 bases correct
    double qualScoreAll = 0.0;

    int numReadsHQ = 0;
    double qualScoreHQ = 0.0;
    double readLenHQ = 0;

    int num6Q7 = 0;

    FILE *fastQFile = fopen(fastQFileName, "w");

    char dnaTemplate[1024];
    char basecalls[1024];
    std::vector<double> errorPerCycle(numCycles, 0);
    std::vector<double> qualScoreList;
    int progress = -1;
    for(int i=0;i<numSpots;i++) {
        PhaseParams params;
        if (gridsearch)
            params = GridSearch(dnaTemplate, spotData[i].vals, errorPerCycle, maxNumCycles);
        else {
            params.ie = ie;
            params.cf = cf;
            params.dr = dr;
            params.err = CallBases(dnaTemplate, spotData[i].vals, errorPerCycle, params.ie, params.cf, params.dr);
        }
        template2bases(dnaTemplate, basecalls);

        // assign a phred score to each base, using the error in measured vs predicted at each position
        std::string qscores;
        for(int c=0;c<numCycles;c++) {
            double err = errorPerCycle[c] * 2.0; // range is 0.0 to 1.0, 0.0 is very high quality Q30, 1.0 is Q3 (50% prob)
            if (err < 0.0) err = 0.0;
            if (err > 1.0) err = 1.0;
            double prob = 0.001 + err*(0.5-0.001);
            int QScore = int(-10.0 * log10(prob)); // truncate down to be fair
            char QChar = QScore + 33; // phred-33 encoding
            qscores += QChar;
        }

        if (fastQFile) {
            std::stringstream ss;
            ss << " {";
            if (gridsearch) {
                ss << "\"ie\":" << params.ie << ",\"cf\":" << params.cf << ",\"dr\":" << params.dr << ",";
            }
            ss << "\"err\":" << params.err << "}";
            fprintf(fastQFile, "@spot_%06d%s\n%s\n+\n%s\n", i, ss.str().c_str(), basecalls, qscores.c_str());
        }

        if (gridsearch && verbose>0) {
            double percent = (double)i/(numSpots-1.0);
            int newProgress = (int)(percent*100.0);
            if (newProgress != progress) {
                progress = newProgress;
                printf("\rprogress: %02d%%", progress);
                fflush(stdout);
            }
        }
    }
    if (gridsearch && verbose>0)
        printf("\ndone.\n");

    if (fastQFile)
        fclose(fastQFile);
}

