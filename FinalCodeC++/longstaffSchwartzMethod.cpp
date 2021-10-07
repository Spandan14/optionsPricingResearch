#include "functions.h"
#include <vector>
#include <iostream>

using namespace std;

double longstaffSchwartzMethod(int N, double T, double r, double S_0, double K, double v, double q, int M, int k, callPut callOrPut) {
    double mult = 2 * ((static_cast<int>(callOrPut)) - 1.5);
    
    double S[LSM_MAXM][LSM_MAXN + 1];
    double P[LSM_MAXM];

    double dt = T / N;

    vector<int> inMoneyPaths;

    double* brownianPath;
    for (int i = 0; i < M; ++i) {
        brownianPath = brownianMotion(N, T, r, S_0, v, q);
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

        for (int i = 0; i < inMoneyPaths.size(); ++i) {
            int temp = inMoneyPaths[i];
            P[temp] = max(mult * (S[temp][time] - K), 0.0) > contValues[temp] ? max(mult * (S[temp][time] - K), 0.0) : exp(-r * dt) * P[temp];
        }

        inMoneyPaths.clear();
    }

    double price = 0;
    for (int i = 0; i < M; ++i) {
        price += exp(-r * dt) * P[i];
    }
    return price / M;
}