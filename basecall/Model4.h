#pragma once

#include <vector>

struct Signal {
    double v[5];
    Signal() {
        memset(v, 0, sizeof(double) * 5);
    }
};

class Model4 {
    public:
        Model4();
        virtual ~Model4();
        void Init(int strandLen);

        void SetParams(double ie, double cf, double dr = 0);
        void SetParams(double ie, double cf[4], double dr = 0);
        void ApplyUV(char *dnaTemplate, int maxLen = 0);
        Signal GetSignal(char *dnaTemplate);
        void Reset();
        int Base2Index(char base);

        int strandLen;
        std::vector<double> state;
        std::vector<double> prevState;

        double dr;
        double ie;
        double cf[4];
        double cfDefault;

};

