#include "functions.h"

using namespace std;

double richardsonExtrapolatedBinomialBlackScholes(int N, double T, double r, double S_0, double K, double v, double q, callPut callOrPut) {
    double doubleStepResult = binomialBlackScholes(2 * N, T, r, S_0, K, v, q, callOrPut);
    double singleStepResult = binomialBlackScholes(N, T, r, S_0, K, v, q, callOrPut);

    return doubleStepResult * 2 - singleStepResult;
}