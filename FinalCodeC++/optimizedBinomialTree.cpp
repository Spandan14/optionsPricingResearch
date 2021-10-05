#include "functions.h"

using namespace std;

double optimizedBinomialTree(int N, double T, double r, double S_0, double K, double v, double q, callPut callOrPut) {
    double mult = 2 * ((static_cast<int>(callOrPut)) - 1.5);

    double dt = T / N;
    double dsct = exp(-r * dt);

    double u = exp(v * sqrt(dt));
    double d = 1 / u;
    double p = ((exp((r - q) * dt)) - d)/(u - d);
    
    double stockPrice[MAXN * 2 + 1];
    double optionPrice[MAXN * 2 + 1];

    stockPrice[N] = S_0;

    for (int j = 1; j < N + 1; ++j) {
        stockPrice[j + N] = stockPrice[N + j - 1] * u;
        stockPrice[N - j] = stockPrice[N - j + 1] * d;
    }

    for (int j = -N; j < N + 1; j += 2) {
        optionPrice[j + N] = max(mult * (stockPrice[j + N] - K), 0.0);
    }

    for (int n = N - 1; n > -1; --n) {
        for (int j = -n; j < n + 1; j+=2) {
            optionPrice[j + N] = max(dsct * (p * optionPrice[j + N + 1] + (1 - p) * optionPrice[j + N - 1]), mult * (stockPrice[j + N] - K));
        }
    }

    return optionPrice[N];
}
