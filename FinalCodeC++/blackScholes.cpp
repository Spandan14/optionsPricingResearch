#include "functions.h"

double blackScholes(double tau, double r, double S_0, double K, double v, double q, callPut callOrPut) {
    double mult = 2 * ((static_cast<int>(callOrPut)) - 1.5);
    
    double d1 = (log(S_0 / K) + (r - q + (0.5 * (pow(v, 2)))) * tau)/(v*sqrt(tau));
    double d2 = d1 - v * sqrt(tau);
    return mult * (S_0 * exp(-q * tau) * normalCDF(mult * d1) - normalCDF(mult * d2) * K * exp(-r * tau));
}


