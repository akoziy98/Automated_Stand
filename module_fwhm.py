import numpy as np
import math

def fwhm(array):
    if type(array) != int:
        maxval = max(array) - min(array)
        minval = min(array)
        fstel = 0
        secel = 0
        n = len(array)
        for i in range(1, n):
            if array[i] > minval + maxval / 2 and array[i - 1] < minval + maxval / 2 and fstel == 0:
                fstel = i
                break

        for i in range(n-1, 0, -1):
            if array[i] < minval + maxval / 2 and array[i - 1] > minval + maxval / 2 and secel == 0:
                secel = i
                break
        #print(fstel, secel)
        return (fstel != 0) * (secel != 0) * (secel - fstel)
    else:
        return 0

def fwtm(array):
    maxval = max(array) - min(array)
    minval = min(array)
    fstel = 0
    secel = 0
    n = len(array)
    for i in range(1, n):
        if array[i] > minval + maxval / 10 and array[i - 1] < minval + maxval / 10:
            fstel = i
            break

    for i in range(n-1, 0, -1):
        if array[i] < minval + maxval / 10 and array[i - 1] > minval + maxval / 10:
            secel = i
            break

    #print(fstel, secel)
    return (fstel != 0) * (secel != 0) * (secel - fstel)
