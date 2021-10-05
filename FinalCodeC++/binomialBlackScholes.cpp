#include "functions.h"

using namespace std;

double binomialBlackScholes(int N, double T, double r, double S_0, double K, double v, double q, callPut callOrPut) {
    double mult = 2 * ((static_cast<int>(callOrPut)) - 1.5);
    
    double dt = T / N;
    double u = exp(v * sqrt(dt));
    double d = 1 / u;
    double p = ((exp((r - q) * dt)) - d) / (u - d);
    double dsct = exp(-r * dt);

    double ans[MAXN];

    for (int j = 0; j < N; ++j) {
        ans[j] = S_0 * pow(u, (N - 1 - 2 * j));
    }

    for (int j = 0; j < N; ++j) {
        ans[j] = max(mult * (ans[j] - K), blackScholes(dt, r, ans[j], K, v, q, callOrPut));
    }

    for (int n = N - 2; n > -1; --n) {
        for (int j = 0; j < n + 1; ++j) {
            ans[j] = max(mult * (S_0 * pow(u, (n - 2 * j)) - K), dsct * (p * ans[j] + (1 - p) * ans[j + 1]));
        }
    }

    return ans[0];
}