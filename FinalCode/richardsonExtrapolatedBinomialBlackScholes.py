import math
import numpy as np
from binomialBlackScholes import binomialBlackScholes

def richardsonExtrapolatedBinomialBlackScholes(N, T, r, S_0, K, v, q, callOrPut):
    doubleStepResult = binomialBlackScholes(2 * N, T, r, S_0, K, v, q, callOrPut)
    singleStepResult = binomialBlackScholes(N, T, r, S_0, K, v, q, callOrPut)

    return doubleStepResult * 2 - singleStepResult