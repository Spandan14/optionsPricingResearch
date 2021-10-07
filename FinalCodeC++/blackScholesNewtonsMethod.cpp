#include "functions.h"

using namespace std;

double blackScholesNewtonsMethod(double tau, double r, double S_0, double K, double v_0, double q, callPut callOrPut, int newtonSteps, double optionPrice) {
    double v_old = v_0;
    double v_new = 0;
    for (int step = 0; step < newtonSteps; ++step) {
        v_new = v_old - ((blackScholes(tau, r, S_0, K, v_old, q, callOrPut) - optionPrice)/blackScholesPartialWRTVolatility(tau, r, S_0, K, v_old, q, callOrPut));
    }

    return v_new;
}