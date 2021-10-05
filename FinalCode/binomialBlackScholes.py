import math
import numpy as np
from blackScholes import blackScholes

def binomialBlackScholes(N, T, r, S0, K, v, q, callOrPut):
    dt = T / N
    u = math.e ** (v * math.sqrt(dt))
    d = 1 / u
    p = ((math.e ** ((r - q) * dt)) - d) / (u - d)
    dsct = math.e ** (-r * dt)
    mult = 2 * (callOrPut.value - 1.5)

    ans = np.zeros(N)
    for j in range(0, N):
        ans[j] = (S0 * (u ** (N - 1 - 2 * j)))

    for j in range(0, N):
        ans[j] = max(mult * (ans[j] - K), blackScholes(dt, r, ans[j], K, v, q, callOrPut))

    for n in range(N - 2, -1, -1):
        for j in range(0, n + 1):
            ans[j] = max(mult * (S0 * (u ** (n - 2*j)) - K),
                         dsct * (p * ans[j] + (1 - p) * ans[j + 1]))
    return ans[0]