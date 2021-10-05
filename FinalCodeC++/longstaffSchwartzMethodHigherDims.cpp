#include "functions.h"
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

using namespace std;

double longstaffSchwartzMethodHigherDims(int N, double T, double r, double S_0[MAXDIM], double K, double v[MAXDIM], double q[MAXDIM], int M, int k, double p[MAXDIM][MAXDIM], double w[MAXDIM], int numAssets, callPut callOrPut) {
    double mult = 2 * ((static_cast<int>(callOrPut)) - 1.5);
    
    double S[LSM_MAXM][LSM_MAXN + 1];
    double P[LSM_MAXM];

    double dt = T / N;

    vector<int> inMoneyPaths;

    double* brownianPath;
    for (int i = 0; i < M; ++i) {
        brownianPath = brownianMotionMultiplePaths(N, T, r, S_0, v, q, p, w, numAssets);
        for (int j = 0; j < N + 1; ++j) {
            S[i][j] = brownianPath[j];
        }
        P[i] = max(mult * (S[i][N] - K), 0.0);
    }   

    for (int time = N - 1; time > 0; --time) {
        double contValues[LSM_MAXM];

        for (int i = 0; i < M; ++i) {
            if (mult * (S[i][time] - K) > 0) {
                inMoneyPaths.push_back(i);
            }
        }

        if (inMoneyPaths.size() == 0) {
            for (int i = 0; i < M; ++i) {
                P[i] *= exp(-r * dt);
            }
            continue;
        }

        double xi[LSM_MAXM];
        double yi[LSM_MAXM];

        for (int i = 0; i < inMoneyPaths.size(); ++i) {
            xi[i] = S[inMoneyPaths[i]][time];
            yi[i] = exp(-r * dt) * P[inMoneyPaths[i]];
        }

        Eigen::MatrixXd beta = leastSquares(xi, yi, k, inMoneyPaths.size());
     
        
        for (int i = 0; i < inMoneyPaths.size(); ++i) {
            double sum = 0;
            for (int j = 0; j < k + 1; ++j) {
                sum += beta(j) * laguerre(j, xi[i]);
            }
            contValues[inMoneyPaths[i]] = sum;
            sum = 0;
        }

        for (int i = 0; i < M; ++i) {
            P[i] = find(inMoneyPaths.begin(), inMoneyPaths.end(), i) != inMoneyPaths.end() && max(mult * (S[i][time] - K), 0.0) > contValues[i] ? max(mult * (S[i][time] - K), 0.0) : exp(-r * dt) * P[i];
        }

        inMoneyPaths.clear();
    }

    double price = 0;
    for (int i = 0; i < M; ++i) {
        price += exp(-r * dt) * P[i];
    }
    return price / M;
}