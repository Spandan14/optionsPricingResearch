import math
import numpy as np

def brownianMotion(N, T, r, S_0, v, q):
    dt = T / N
    brownianPath = np.zeros(N + 1)
    brownianPath[0] = S_0
    for i in range(0, N):
        rand = np.random.normal(0, 1)
        brownianPath[i + 1] = brownianPath[i] * math.e ** ((r - q - .5 * (v ** 2)) * (dt) + v * math.sqrt(dt) * rand)
    return brownianPath
