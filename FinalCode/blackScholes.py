import math
import numpy as np
from scipy.stats import norm

def blackScholes(tau, r, S_0, K, v, q, callOrPut):
    normalCDF = norm.cdf
    mult = 2 * (callOrPut.value - 1.5)

    d1 = (np.log(S_0 / K) + (r - q + (.5 * (v ** 2))) * tau) / (v * math.sqrt(tau))
    d2 = d1 - (v * math.sqrt(tau))
    return mult * (S_0 * (math.e ** (-q * tau)) * normalCDF(mult * d1) - (normalCDF(mult * d2) * K * math.e**(-r * tau)))