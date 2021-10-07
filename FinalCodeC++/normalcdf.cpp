#include "functions.h"

double normalCDF(double val) {
    return 0.5 * (erfc(-val / sqrt(2)));
}

