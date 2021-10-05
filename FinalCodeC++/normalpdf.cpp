#include "functions.h"

double normalPDF(double val) {
    double pi = 3.141592653589793238462643383279502884197169399375105820974944;
    return exp(-0.5 * pow(val, 2)) / sqrt(2 * pi);
}

