#include "functions.h"

using namespace std;

double trinomialTree(int N, double T, double r, double S_0, double K, double v, double q, callPut callOrPut) {
    double mult = 2 * ((static_cast<int>(callOrPut)) - 1.5);
    double dt = T/N;
    double u = exp(v * sqrt(2 * dt));
    double d = 1 / u;
    double m = 1;

    double R = exp((r-q) * dt);
    double S = pow(v, 2) * dt;

    double pu = (pow(R, 2) + S - R - (R/u) + (1/u))/(pow(u, 2) - u - 1 + (1 / u));
    double pd = (pow(R, 2) + S - R - u * R + u)/((1 / pow(u, 2)) - (1/u) - 1 + u);
    double pm = 1 - pu - pd;

    double dsct = exp(-r * dt);

    double ans[2 * MAXN + 1];

    for (int j = 0; j < 2 * N + 1; ++j) {
        ans[j] = max(0.0, mult * (S_0 * pow(u, (N - j)) - K));
    }
    
    for (int n = N - 1; n > -1; --n) {
        for (int j = 0; j < 2 * n + 1; ++j) {
            ans[j] = max(mult * (S_0 * pow(u, (n - j)) - K), dsct * (pu * ans[j] + pm * ans[j + 1] + pd * ans[j + 2]));
        }
    }

    return ans[0];
}