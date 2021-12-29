import numpy as np

def time_resolution(filename):

    with open (filename) as fn:
        for line in fn:
            line = line.strip().split(';')
            cnt = np.array([int(el) for el in line])

    return cnt
    