#include "functions.h"
#include <random>

using namespace std;

double* brownianMotionMultiplePaths(int N, double T, double r, double S_0[MAXDIM], double v[MAXDIM], double q[MAXDIM], double p[MAXDIM][MAXDIM], double w[MAXDIM], int numAssets) {
    double dt = T / N;
    double A[MAXDIM][MAXDIM];
    static double S_basket[MAXN+1];
    double sig[MAXDIM][MAXDIM];

    double count[MAXDIM];
    double brownianPaths[MAXN+1][MAXDIM];

    unsigned seed = rand();
    default_random_engine gen(seed);
    normal_distribution<double> dist(0.0, 1.0);

    for (int i = 0; i < numAssets; ++i) {
        brownianPaths[i][0] = S_0[i];
        for (int j = 0; j < i + 1; ++j) {
            if (i == j) {
                sig[i][j] = pow(v[i], 2);
            } else {
                sig[i][j] = p[i][j] * v[i] * v[j];
                sig[j][i] = sig[i][j];
            }

            double sum = 0;
            for (int k = 0; k < j; ++k) {
                sum += A[i][k] * A[j][k];
            }
            if (i == j) {
                A[i][j] = sqrt(sig[i][i] - sum);
            } else {
                A[i][j] = (1.0 / A[j][j] * (sig[i][j] - sum));
            }
        }
    }

    for (int i = 0; i < N; ++i) {
        double rand[MAXDIM];
        for (int j = 0; j < numAssets; ++j) {
            rand[j] = dist(gen);
        }
        for (int j = 0; j < numAssets; ++j) {
            double sum = 0;
            for (int k = 0; k < numAssets; ++k) {
                sum += sig[j][k] * rand[k];
            }
            count[j] = sum;
        }
        for (int j = 0; j < numAssets; ++j) {
            brownianPaths[j][i + 1] = brownianPaths[j][i] * exp((r - q[j] - .5 * pow(v[j], 2)) * (dt) + count[j] * sqrt(dt));
        }
    }

    for (int i = 0; i < N + 1; ++i) {
        double sum = 0;
        for (int j = 0; j < numAssets; ++j) {
            sum += brownianPaths[j][i] * w[j];
        }
        S_basket[i] = sum;
    }

    return S_basket;
}