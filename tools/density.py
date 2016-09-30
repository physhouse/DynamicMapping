import sys
import math
import numpy as np
np.set_printoptions(threshold=np.nan)

def pbc(x0, x1, L):
    d = x0 - x1

    if d > 0.5*L:
        return d - L
    elif d < -0.5*L:
        return d + L
    else:
        return d

def density(r1, r2, h, L):
    dx = pbc(r1[0], r2[0], L)
    dy = pbc(r1[1], r2[1], L)
    dz = pbc(r1[2], r2[2], L)
    q = math.sqrt(dx*dx + dy*dy + dz*dz) / (2.0*h)
    q2 = q * q
    q3 = q * q2

    w = 8.0 / math.pi
    if q <= 0.5:
        w = w * (1.0 - 6 * q2 + 6 * q3)
    elif q <= 1.0 and q > 0.5:
        w = w * 2.0 * (1 - 3 * q + 3 * q2 - q3)
    else:
        w  = 0.0

    return w

def readFrame(step, numCG, fp, out):
    L = 0.0
    for iline in range(0,9):
        line = fp.readline().split()
        if iline == 5:
            elements = [float(x) for x in line]
            L = elements[1]

    pad = numCG ** (1./3)
    h = L / pad * 1.00  # h should be at least d, giving a safety factor
    R = np.zeros((numCG,3))
    rou = np.zeros((numCG,1))
    
    for i in range(0, numCG):
        line = fp.readline().split()
        R[i] = [float(x) for x in line[3:6]]

    for i in range(0, numCG - 1):
        for j in range(i+1, numCG):
            w = density(R[i], R[j], h, L)
            rou[i] = rou[i] + w
            rou[j] = rou[j] + w

    rou.tofile(out, sep='\n')
    out.write('\n')
    #for i in range(0, numCG):
     #   out.write(str(rou[i]))
      #  out.write('\n')

if __name__ == "__main__":
    args = sys.argv

    nFrame = int(args[1])
    numCG  = int(args[2])

    f = open('CG_TRJ.lmpstrj', 'r')
    out = open('density.dat', 'w')

    for step in range(0, nFrame):
        readFrame(step, numCG, f, out)
