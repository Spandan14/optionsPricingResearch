#include "functions.h"

using namespace std;

double blackScholesPartialWRTVolatility(double tau, double r, double S_0, double K, double v, double q, callPut callOrPut) {
    double mult = 2 * ((static_cast<int>(callOrPut)) - 1.5);
    
    double d1 = mult * ((log(S_0 / K) + tau * (pow(v, 2) / 2 + r)) / (v * sqrt(tau)));
    double d2 = mult * (d1 - v * sqrt(tau));

    double d1Partial = sqrt(tau) - (tau * (pow(v, 2) / 2 + r) + log(S_0 / K)) / (pow(v, 2) * sqrt(tau));
    double d2Partial = d1Partial - sqrt(tau);

    return mult * (S_0 * exp(-q * tau) * normalPDF(d1) * d1Partial - K * exp(-r * tau) * normalPDF(d2) * d2Partial);
}