#ifndef FUNCTIONS
#define FUNCTIONS

#include <math.h>
#include <Eigen/Dense>
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

enum callPut {Put=1, Call=2};

const int MAXN = 75000;
const int LSM_MAXN = 300;
const int LSM_MAXM = 100000;
const int MAXDIM = 4;

double normalCDF(double val);
double blackScholes(double tau, double r, double S_0, double K, double v, double q, callPut callOrPut);
double binomialTree(int N, double T, double r, double S_0, double K, double v, double q, callPut callOrPut);
double optimizedBinomialTree(int N, double T, double r, double S_0, double K, double v, double q, callPut callOrPut);
double binomialBlackScholes(int N, double T, double r, double S_0, double K, double v, double q, callPut callOrPut);
double richardsonExtrapolatedBinomialBlackScholes(int N, double T, double r, double S_0, double K, double v, double q, callPut callOrPut);
double trinomialTree(int N, double T, double r, double S_0, double K, double v, double q, callPut callOrPut);
double* brownianMotion(int N, double T, double r, double S_0, double v, double q);
Eigen::MatrixXd leastSquares(double xi[LSM_MAXM], double yi[LSM_MAXM], int k, int len);
double longstaffSchwartzMethod(int N, double T, double r, double S_0, double K, double v, double q, int M, int k, callPut callOrPut);
double* brownianMotionMultiplePaths(int N, double T, double r, double S_0[MAXDIM], double v[MAXDIM], double q[MAXDIM], double p[MAXDIM][MAXDIM], double w[MAXDIM], int numAssets);
double longstaffSchwartzMethodHigherDims(int N, double T, double r, double S_0[MAXDIM], double K, double v[MAXDIM], double q[MAXDIM], int M, int k, double p[MAXDIM][MAXDIM], double w[MAXDIM], int numAssets, callPut callOrPut);
double normalPDF(double val);
double blackScholesPartialWRTVolatility(double tau, double r, double S_0, double K, double v, double q, callPut callOrPut);
double blackScholesNewtonsMethod(double tau, double r, double S_0, double K, double v_0, double q, callPut callOrPut, int newtonSteps, double optionPrice);
double binarySearchVolatility(double T, double r, double S_0, double K, double minVol, double maxVol, double q, callPut callOrPut, int binarySteps, int modelSteps, double (*function)(int, double, double, double, double, double, double, callPut), double optionPrice);
double* earlyExerciseBoundaryGen(int N, double T, double r, double S_0, double K, double v, double q, int binarySteps, double (*function)(int, double, double, double, double, double, double, callPut), callPut callOrPut);

#endif