import math
import numpy as np

def trinomialTree(N, T, r, S_0, K, v, q, callOrPut):
    dt = T/N
    u = math.e**(v*math.sqrt(2*dt))
    d = 1/u
    m = 1

    R = math.e ** ((r-q) * dt)
    S = v ** 2 * dt

    pu = (R ** 2 + S - R - (R/u) + (1/u))/(u ** 2 - u - 1 + (1 / u))
    pd = (R ** 2 + S - R - u * R + u)/((1 / u ** 2) - (1/u) - 1 + u)
    pm = 1 - pu - pd
    
    dsct = math.e**(-r*dt)
    mult = 2 * (callOrPut.value - 1.5)

    ans = np.zeros(2*N+1)

    for j in range(0, 2*N + 1):
        ans[j] = max(0, mult * (S_0 * (u ** (N - j)) - K))

    for n in range(N - 1, -1, -1):
        for j in range(0, 2*n + 1):
            ans[j] = max(mult * (S_0 * (u ** (n - j)) - K), dsct * (pu * ans[j] + pm * ans[j + 1] + pd * ans[j + 2]))

    return ans[0]