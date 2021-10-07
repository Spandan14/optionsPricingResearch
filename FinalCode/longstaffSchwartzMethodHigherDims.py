import math
import numpy as np
from brownianMotionMultiplePaths import brownianMotionMultiplePaths
from leastSquares import leastSquares

def longstaffSchwartzMethodHigherDims(N, T, r, S_0, K, v, q, M, k, p, w, callOrPut):
    S = np.zeros((N + 1, M))
    P = np.zeros(M)
    dt = T / N
    mult = 2 * (callOrPut.value - 1.5)
    inMoneyPaths = np.array

    for i in range(M):
        S[i] = brownianMotionMultiplePaths(N, T, r, S_0, v, q, p, w, len(S_0))
        P[i] = max(mult * (S[i][N] - K), 0)
    lastRun = True

    for time in range(N - 1, 0, -1):
        laguerreWeights = np.zeros(k + 1)
        continuationValues = np.zeros(M)

        for i in range(M):
            if max(mult * (S[i][time] - K), 0) > 0:
                inMoneyPaths = np.append(inMoneyPaths, i)
        if lastRun:
            inMoneyPaths = np.delete(inMoneyPaths, 0)

        if len(inMoneyPaths) == 0:
            P[i] = math.e ** (-r * dt) * P[i]
            lastRun = False
            continue
        lastRun = True

        xi = [0] * len(inMoneyPaths)
        yi = [0] * len(inMoneyPaths)

        for i in range(0, len(inMoneyPaths)):
            xi[i] = S[inMoneyPaths[i]][time]
            yi[i] = math.e ** (-r * dt) * P[inMoneyPaths[i]]

        beta = leastSquares(xi, yi, k)

        count = 0
        for i in range(len(inMoneyPaths)):
            for j in range(k + 1):
                laguerreWeights[j] = 1
                Z = np.polynomial.Laguerre(laguerreWeights)
                count += beta[j] * Z(xi[i])
                laguerreWeights[j] = 0
            continuationValues[inMoneyPaths[i]] = count
            count = 0

        for i in range(inMoneyPaths):
            if(max(mult*(S[inMoneyPaths[i]][time] - K), 0) > continuationValues[inMoneyPaths[i]]):
                P[i] = (max(mult * (S[inMoneyPaths[i]][time] - K), 0))
            else:
                P[i] = math.e ** (-r * dt) * P[inMoneyPaths[i]]

        inMoneyPaths = np.array

    price = 0
    for i in range(M):
        price += (math.e ** (-r * dt)) * P[i]
    return price / M
