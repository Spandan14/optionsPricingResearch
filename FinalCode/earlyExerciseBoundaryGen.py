import math
import numpy as np
from callPutEnum import CallPut

def earlyExerciseBoundaryGen(N, T, r, S_0, K, v, q, binarySteps, function, callOrPut):
    ans = np.zeros(N + 1)
    mult = 2 * (callOrPut.value - 1.5)

    for n in range(0, N):
        tau = T * (N - n) / N 
        dt = tau / N
        u = math.e ** (v * math.sqrt(dt))

        maxBound = S_0 * u**(n)*10
        minBound = -maxBound

        value = S_0

        for binaryStep in range(0, binarySteps):
            intrinsic = max(0, mult * (value - K))
            optionValue = function(N, tau, r, value, K, v, q, callOrPut)
            if mult == 1:
                if optionValue - 0.001 < intrinsic and intrinsic != 0:
                    maxBound = value
                    value = (value + minBound) / 2
                else:
                    minBound = value
                    value = (value + maxBound) / 2
            if mult == -1:
                if optionValue - 0.001 > intrinsic or intrinsic == 0:
                    maxBound = value
                    value = (value + minBound) / 2
                else:
                    minBound = value
                    value = (value + maxBound) / 2
        
        ans[n] = value

    ans[N] = K
    return ans