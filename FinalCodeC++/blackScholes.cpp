#include "functions.h"

double blackScholes(double tau, double r, double S_0, double k, double v, double q, callPut callOrPut) {
    double mult = 2 * ((static_cast<int>(callOrPut)) - 1.5);
    
    double d1 = mult * ((log(S_0 / k) + (r - q + (0.5 * (pow(v, 2)))) * tau)/(v*sqrt(tau)));
    double d2 = mult * (d1 - v * sqrt(tau));
    return mult * (S_0 * exp(-q * tau) * normalCDF(mult * d1) - normalCDF(mult * d2) * k * exp(-r * tau));
}


