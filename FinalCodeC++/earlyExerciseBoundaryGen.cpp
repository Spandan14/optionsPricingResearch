#include "functions.h"

using namespace std;

double* earlyExerciseBoundaryGen(int N, double T, double r, double S_0, double K, double v, double q, int binarySteps, double (*function)(int, double, double, double, double, double, double, callPut), callPut callOrPut) {
    static double ans[MAXN + 1];
    double mult = 2 * ((static_cast<int>(callOrPut)) - 1.5);
    for (int n = 0; n < N; ++n) {
        
        double tau = T * (N - n) / N;
        double dt = tau / N;
        double u = exp(v * sqrt(dt));

        double maxBound = S_0 * pow(u, n) * 10;
        double minBound = -maxBound;

        
        double value = S_0;
        
        for (int binaryStep = 0; binaryStep < binarySteps; ++binaryStep) {
            double intrinsic = max(0.0, mult * (value - K));
            double optionValue = function(N, tau, r, value, K, v, q, callOrPut);
            if (mult == 1) {
                if (optionValue - 0.001 < intrinsic && intrinsic != 0) {
                    maxBound = value;
                    value = (value + minBound) / 2;
                } else {
                    minBound = value;
                    value = (value + maxBound) / 2;
                }
            } else {
                if (optionValue - 0.001 > intrinsic || intrinsic == 0) {
                    maxBound = value;
                    value = (value + minBound) / 2;
                } else {
                    minBound = value;
                    value = (value + maxBound) / 2;
                }
            }
        }
        ans[n] = value;
    }
    
    ans[N] = K;
    return ans;
}