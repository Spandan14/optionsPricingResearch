import math
import numpy as np
from scipy.stats import norm

def blackScholesPartialWRTVolatility(tau, r, S_0, K, v, q, callOrPut):
    NPrime = norm.pdf
    mult = 2 * (callOrPut.value - 1.5)

    d1 = mult * ((np.log(S_0 / K) + tau * ((v ** 2) / 2 + r)) / (v * math.sqrt(tau)))
    d2 = mult * (d1 - v * math.sqrt(tau))
    
    d1Partial = math.sqrt(tau) - (tau * ((v ** 2) / 2 + r) + np.log(S_0 / K)) / (v ** 2 * math.sqrt(tau))
    d2Partial = d1Partial - math.sqrt(tau)

    return mult * (S_0 * math.e ** (-q * tau) * NPrime(d1) * d1Partial - K * math.e ** (-r * tau) * NPrime(d2) * d2Partial)


