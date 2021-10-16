#include "functions.h"

using namespace std;

double blackScholesPartialWRTVolatility(double tau, double r, double S_0, double K, double v, double q, callPut callOrPut) {
    double mult = 2 * ((static_cast<int>(callOrPut)) - 1.5);
    
    double d1 = ((log(S_0 / K) + tau * (pow(v, 2) / 2 + r)) / (v * sqrt(tau)));
    double d2 = d1 - v * sqrt(tau);

    double d1Partial = mult * (sqrt(tau) - (tau * (pow(v, 2) / 2 + r) + log(S_0 / K)) / (pow(v, 2) * sqrt(tau)));
    double d2Partial = mult * (d1Partial - sqrt(tau));

    return mult * (S_0 * exp(-q * tau) * normalPDF(mult * d1) * d1Partial - K * exp(-r * tau) * normalPDF(mult * d2) * d2Partial);
}
