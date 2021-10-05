import math
import numpy as np

def brownianMotionMultiplePaths(N, T, r, S_0, v, q, p, w, numAssets):
    dt = T / N
    A = np.zeros((numAssets, numAssets))
    S_basket = np.zeros(N+1)
    sig = np.zeros((numAssets, numAssets))
    
    count = np.zeros(numAssets)
    brownianPaths = np.zeros((N + 1, numAssets))
    
    for i in range(numAssets):
        brownianPaths[i][0] = S_0[i]
        for j in range(i + 1):
            if (i == j):
                sig[i][j] = v[i] ** 2
            else:
                sig[i][j] = p[i][j] * v[i] * v[j]
                sig[j][i] = sig[i][j]

        for j in range(i + 1):
            sum = 0
            for k in range(j):
                sum += A[i][k] * A[j][k]
            if (i == j):
                A[i][j] = np.sqrt(sig[i][i] - sum)
            else:
                A[i][j] = (1.0 / A[j][j] * (sig[i][j] - sum))

    for i in range(N):
        rand = np.random.normal(0, 1, numAssets)
        for j in range(numAssets):
            sum = 0
            for k in range(numAssets):
                sum += sig[j][k] * rand[k]
            count[j] = sum
        for j in range(numAssets):
            brownianPaths[j][i + 1] = brownianPaths[j][i] * math.e ** ((r - q[j] - .5 * (v[j] ** 2)) * (dt) + count[j] * math.sqrt(dt))

    for i in range(N + 1):
        sum = 0
        for j in range(numAssets):
            sum += brownianPaths[j][i] * w[j]
        S_basket[i] = sum

    return S_basket