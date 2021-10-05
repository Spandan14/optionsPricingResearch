import math
import numpy as np

def optimizedBinomialTree(N, T, r, S_0, K, v, q, callOrPut):
    dt = T / N
    u = math.e ** (v * math.sqrt(dt)) 
    d = 1 / u
    p = ((math.e ** ((r - q) * dt)) - d)/(u - d)
    dsct = math.e**(-r*dt)
    
    mult = 2 * (callOrPut.value - 1.5)

    stockPrice = np.zeros(2 * N + 1)
    optionPrice = np.zeros(2 * N + 1)

    stockPrice[N] = S_0

    for j in range(1, N + 1):
        stockPrice[j + N] = stockPrice[N + j - 1] * u
        stockPrice[N - j] = stockPrice[N - j + 1] * d

    for j in range(-N, N + 1, 2):
        optionPrice[j + N] = max(mult * (stockPrice[j + N] - K), 0)

    for n in range(N - 1, -1, -1):
        for j in range(-n, n + 1, 2):
            optionPrice[j + N] = max(dsct * (p * optionPrice[j + N + 1] + (1 - p) * optionPrice[j + N - 1]), mult * (stockPrice[j + N] - K))

    return optionPrice[N]