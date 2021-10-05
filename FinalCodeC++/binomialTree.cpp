#include "functions.h"

using namespace std;

double binomialTree(int N, double T, double r, double S_0, double K, double v, double q, callPut callOrPut) {
    double mult = 2 * ((static_cast<int>(callOrPut)) - 1.5);

    double dt = T/N;
    double u = exp(v * sqrt(dt));
    double d = 1 / u;
    double p = ((exp((r - q) * dt)) - d) / (u - d);
    double dsct = exp(-r * dt);
    
    double ans[MAXN];

    for (int j = 0; j < N + 1; ++j) {
        ans[j] = max(0.0, mult * (S_0 * (pow(u, N - 2 * j))- K));
    }

    for (int n = N - 1; n > -1; --n) {
        for (int j = 0; j < n + 1; ++j) {
            ans[j] = max(mult * (S_0 * pow(u, n - 2 * j) - K), dsct * (p * ans[j] + (1 - p) * ans[j + 1]));
        }
    }

    return ans[0];
}
