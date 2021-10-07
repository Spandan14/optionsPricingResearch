#include "functions.h"

using namespace std;

Eigen::MatrixXd leastSquares(double xi[LSM_MAXM], double yi[LSM_MAXM], int k, int len) {    
    Eigen::MatrixXd L(len, k+1);
    Eigen::VectorXd ymat(len);

    for (int i = 0; i < len; ++i) {
        ymat(i) = yi[i];
    }

    for (int i = 0; i < len; ++i) {
        for (int j = 0; j < k + 1; ++j) {
            L(i, j) = laguerre(j, xi[i]);
        }
    }

    Eigen::MatrixXd temp = L.transpose() * L;
    return (temp.completeOrthogonalDecomposition().pseudoInverse() * L.transpose()) * ymat;
}