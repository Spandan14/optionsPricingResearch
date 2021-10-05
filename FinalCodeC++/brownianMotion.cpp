#include "functions.h"
#include <random>

using namespace std;

double* brownianMotion(int N, double T, double r, double S_0, double v, double q) {
    double dt = T / N;
    static double brownianPath[MAXN];
    
    unsigned seed = rand();
    default_random_engine gen(seed);
    normal_distribution<double> dist(0.0, 1.0);

    brownianPath[0] = S_0;

    for (int i = 0; i < N; ++i) {
        double rand = dist(gen);
        brownianPath[i + 1] = brownianPath[i] * exp((r - q - 0.5 * pow(v, 2)) * dt + v * sqrt(dt) * rand);
    }

    return brownianPath;
}
