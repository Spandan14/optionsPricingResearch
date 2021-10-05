import math
import numpy as np
from blackScholes import blackScholes
from blackScholesPartialWRTVolatility import blackScholesPartialWRTVolatility

def blackScholesNewtonsMethod(tau, r, S_0, K, v_0, q, callOrPut, newtonSteps, optionPrice):
    v_old = v_0
    v_new = 0
    for step in range(0, newtonSteps):
        v_new = v_old - ((blackScholes(tau, r, S_0, K, v_old, q, callOrPut) - optionPrice)/blackScholesPartialWRTVolatility(tau, r, S_0, K, v_old, q, callOrPut))

    return v_new

