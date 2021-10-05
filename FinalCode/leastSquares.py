import numpy as np

def leastSquares(xi, yi, k):
    L = np.zeros((k + 1, len(yi)))
    weight = np.zeros(k + 1)

    for i in range(len(yi)):
        for j in range(k + 1):
            weight[j] = 1
            Z = np.polynomial.Laguerre(weight)
            L[i][j] = Z(xi[i])
            weight[j] = 0

    beta = np.dot(np.dot(np.linalg.pinv(np.dot(np.transpose(L), L)), np.transpose(L)), np.transpose(yi))

    return beta