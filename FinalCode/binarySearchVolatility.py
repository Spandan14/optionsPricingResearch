import math
import numpy as np
from blackScholes import blackScholes

def binarySearchVolatilty(T, r, S_0, K, minVol, maxVol, q, callOrPut, binarySteps, modelSteps, function, optionPrice):
    currentVol = (minVol + maxVol) / 2
    for step in range(0, binarySteps):
        currentPrice = function(modelSteps, T, r, S_0, K, currentVol, q, callOrPut)
        if currentPrice > optionPrice:
            maxVol = currentVol
            currentVol = (minVol + maxVol)/2
        elif currentPrice == optionPrice:
            break
        else:
            minVol = currentVol
            currentVol = (minVol + maxVol) / 2

    return currentVol