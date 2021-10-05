#include "functions.h"

using namespace std;

double binarySearchVolatility(double T, double r, double S_0, double K, double minVol, double maxVol, double q, callPut callOrPut, int binarySteps, int modelSteps, double (*function)(int, double, double, double, double, double, double, callPut), double optionPrice) {
    double currentVol = (minVol + maxVol) / 2;
    for (int step = 0; step < binarySteps; ++step) {
        double currentPrice = function(modelSteps, T, r, S_0, K, currentVol, q, callOrPut);
        if (currentPrice > optionPrice) {
            maxVol = currentVol;
            currentVol = (minVol + maxVol) / 2;
        }
        else if (currentPrice == optionPrice) {
            break;
        } else {
            minVol = currentVol;
            currentVol = (minVol + maxVol) / 2;
        }
    }
    
    return currentVol;
}