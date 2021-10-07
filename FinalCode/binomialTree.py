import math
import numpy as np

def binomialTree(N, T, r, S_0, K, v, q, callOrPut):
    dt = T / N
    u = math.e ** (v * math.sqrt(dt))
    d = 1 / u
    p = ((math.e ** ((r - q) * dt)) - d) / (u - d)
    dsct = math.e ** (-r * dt)
    mult = 2 * (callOrPut.value - 1.5)

    ans = np.zeros(N + 1)
    for j in range(0, N + 1):
        ans[j] = max(0, mult * (S_0 * (u ** (N - 2 * j)) - K))

    for n in range(N - 1, -1, -1):
        for j in range(0, n + 1):
            ans[j] = max(mult * (S_0 * (u ** (n - 2 * j)) - K), dsct * (p * ans[j] + (1 - p) * ans[j + 1]))

    return ans[0]
