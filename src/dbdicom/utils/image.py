
import numpy as np


def clip(array, value_range = None):

    array[np.isnan(array)] = 0
    if value_range is None:
        finite = array[np.isfinite(array)]
        value_range = [np.amin(finite), np.amax(finite)]
    return np.clip(array, value_range[0], value_range[1])
    

def scale_to_range(array, bits_allocated, signed=False):
        
    range = 2.0**bits_allocated - 1
    if signed:
        minval = -2.0**(bits_allocated-1)
    else:
        minval = 0
    maximum = np.amax(array)
    minimum = np.amin(array)
    if maximum == minimum:
        slope = 1
    else:
        slope = range / (maximum - minimum)
    intercept = -slope * minimum + minval
    array = array * slope
    array = array + intercept

    if bits_allocated == 8:
        if signed:
            return array.astype(np.int8), slope, intercept
        else:
            return array.astype(np.uint8), slope, intercept
    if bits_allocated == 16:
        if signed:
            return array.astype(np.int16), slope, intercept
        else:
            return array.astype(np.uint16), slope, intercept
    if bits_allocated == 32:
        if signed:
            return array.astype(np.int32), slope, intercept
        else:
            return array.astype(np.uint32), slope, intercept
    if bits_allocated == 64:
        if signed:
            return array.astype(np.int64), slope, intercept
        else:
            return array.astype(np.uint64), slope, intercept

