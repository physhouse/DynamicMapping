#!/bin/env python

import numpy as np
import scipy.linalg as linalg
import sys
np.set_printoptions(threshold=np.nan)

def generateForce(step, num, f, out, frame):
    for iline in range(0,9):
        out.write(f.readline())

    dt = 0.005
    newframe = np.zeros((numCG, 12))

    for i in range(0, numCG):
        line = f.readline().split()
        if step > 1 :
            newframe[i] = [float(x) for x in line]
            mass = newframe[i][2]
            newframe[i][9] = mass * (newframe[i][6] - frame[i][6]) / dt
            newframe[i][10] = mass * (newframe[i][7] - frame[i][7]) / dt
            newframe[i][11] = mass * (newframe[i][8] - frame[i][8]) / dt
        else:
            newframe[i] = [float(x) for x in line]

        for ind in range(0,12):
            out.write(str(newframe[i][ind]))
            out.write(' ')
        out.write('\n')

    frame = newframe
    return frame

if __name__ == "__main__":
    args = sys.argv

    nFrame = int(args[1])
    numCG  = int(args[2])
    
    f = open('CG_TRJ.lmpstrj', 'r')
    out = open('force.lmpstrj', 'w')
    frame = np.zeros((numCG, 12))

    for step in range(1, nFrame):
        frame = generateForce(step, numCG, f, out, frame)
