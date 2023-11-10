
#include "Model4.h"
#include <memory.h>

Model4::Model4()
{
}

Model4::~Model4()
{
}

void Model4::Init(int _strandLen)
{
    strandLen = _strandLen;
    Reset();

    dr = 0.0;
    ie = 0.0;
    cfDefault = 0.0;
    memset(cf, 0, 4 * sizeof(double));
}

void Model4::Reset()
{
    state.resize(strandLen);
    prevState.resize(strandLen);

    memset(state.data(), 0, strandLen * sizeof(double));
    state[0] = 1.0;
}

void Model4::SetParams(double _ie, double _cf, double _dr)
{
    ie = _ie;
    dr = _dr;
    for(int i=0;i<4;i++)
        cf[i] = _cf;
    cfDefault = _cf;
}

void Model4::SetParams(double _ie, double _cf[4], double _dr)
{
    ie = _ie;
    dr = _dr;
    for(int i=0;i<4;i++)
        cf[i] = _cf[i];
    cfDefault = _cf[0]; // maybe needs to be the avg
}

int Model4::Base2Index(char base)
{
    // converts a char to an index - this is ugly and needs to be fixed
    if (base == 'A')
        return 0;
    else if (base == 'C')
        return 1;
    else if (base == 'G')
        return 2;
    else if (base == 'T')
        return 3;
    return 0;
}

void Model4::ApplyUV(char *dnaTemplate, int maxLen)
{
    prevState = state;
    int numExtensions = 3; // technically this goes on forever, but after 3 rounds there is not much left
    std::vector<double> extendAmount(numExtensions, 0);
    std::vector<double> totalChange(strandLen+numExtensions, 0);
    int templateLen = strlen(dnaTemplate);
    double cfDye;

    for(int i=0;i<strandLen;i++) {
        // amount of product available to extend
        extendAmount[0] = state[i] * (1.0 - ie);
        double changeAmount = extendAmount[0];
        if (extendAmount[0] == 0 || (maxLen > 0 and i >= (maxLen-1)))
            continue;
        // continue to extend (carry-forward)
        for(int s=1;s<numExtensions;s++) {
            int templateIndex = i + s - 1;
            if (templateIndex < templateLen)
                cfDye = cf[dnaTemplate[templateIndex]-1];
            else
                cfDye = cfDefault;

            double cfAmount = extendAmount[s-1] * cfDye;
            extendAmount[s] = cfAmount;
            extendAmount[s-1] -= extendAmount[s];
        }

        // update the total change from this strand's position
        for (int s=0;s<numExtensions;s++)
            totalChange[i+1+s] += extendAmount[s];
        totalChange[i] -= changeAmount;
    }

    // apply the change to the current state
    for(int i=0;i<strandLen;i++) {
        state[i] += totalChange[i];
        // droop is applied across all states
        state[i] *= (1.0 - dr);
    }
}

Signal Model4::GetSignal(char *dnaTemplate)
{
    Signal signal;
    int templateLen = strlen(dnaTemplate);
    // each position within our state, sum the signal, binned by known DNA bases
    for(int i=0;i<strandLen;i++) {
        if (i < templateLen)
            signal.v[dnaTemplate[i]-1] += state[i];
        else
            signal.v[4] += state[i];
    }
    return signal;
}

