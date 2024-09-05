import scipy
import numpy as np
import random
import gurobipy as gp
from gurobipy import *
import top_improve
import calc_weights
import time
import csv
import glob


def main(argv):
    prefix = argv[1]
    N_axis = int(argv[2])
    number = int(argv[3])
    writer = csv.writer(open("../data/sizefrac%d.csv" % (N_axis), 'a'))
    npzfile = np.load(prefix+"/ad_grb_{:05d}".format(number) + ".log.npz")

    c = npzfile["c"]
    x = npzfile["x"]
    Delta = int(npzfile["Delta"])
    n = round(math.sqrt(len(c)))
    h = npzfile["h"]
    start = time.time()
    m = gp.Model("iptv")

    myzetadic = {}
    N = n*n

    d = {}
    dsign = [1 for i in range(n*n)]
    for i in range(n*n):
        d[i] = m.addVar(name="d"+str(i), vtype="C", lb=0, ub=1)
        if x[i] > 0.99:
            dsign[i] = -1



    zx = {}
    zy = {}
    for i in range(n*(n-1)):
        zx[i] = m.addVar(name="zx"+str(i), vtype="C", lb=0, ub=1.01)
        zy[i] = m.addVar(name="zy"+str(i), vtype="C", lb=0, ub=1.01)

    zz = {}
    zzcount = 0
    opt_with_bnd = npzfile["opt_with_bnd"]

    if opt_with_bnd:
        for i in range(n):
            for j in range(n):
                if i == 0 or j == 0 or i == n-1 or j == n-1:
                    zz[zzcount] = m.addVar(
                        name="zz"+str(zzcount), vtype="C", lb=0)
                    m.addConstr(dsign[i+j*n] * d[i+j*n] +
                                x[i+j*n] <= zz[zzcount])
                    m.addConstr(-dsign[i+j*n] * d[i+j*n] -
                                x[i+j*n] <= zz[zzcount])
                    zzcount += 1

    for j in range(n):
        for i in range(n-1):

            m.addConstr(x[i+1+j*n]+dsign[i+1+j*n] * d[i+1+j*n] -
                        x[i+j*n] - dsign[i+j*n] * d[i+j*n]-zx[j+i*n] <= 0)
            m.addConstr(-x[i+1+j*n]-dsign[i+1+j*n] * d[i+1+j*n] +
                        x[i+j*n] + dsign[i+j*n] * d[i+j*n]-zx[j+i*n] <= 0)
            myzetadic[(i+j*n, i+1+j*n)] = zx[j+i*n]

            m.addConstr(x[(i+1)*n+j] + dsign[(i+1)*n+j] * d[(i+1)*n+j] -
                        x[i*n+j]-dsign[i*n+j] * d[i*n+j]-zy[j+i*n] <= 0)
            m.addConstr(-x[(i+1)*n+j]-dsign[(i+1)*n+j] * d[(i+1)*n+j] +
                        x[i*n+j]+dsign[i*n+j] * d[i*n+j]-zy[j+i*n] <= 0)
            myzetadic[(i*n+j, (i+1)*n+j)] = zy[j+i*n]

    capconstraint = m.addConstr(gp.quicksum(
        d[i] for i in range(n*n))-Delta <= 0)
    m.setObjective(gp.quicksum(h*c[i]*dsign[i]*d[i] for i in range(n*n)) + h*gp.quicksum(zx[i] for i in range(
        n*n-n)) + h*gp.quicksum(zy[i] for i in range(n*n-n)) + h * gp.quicksum(zz[i] for i in range(zzcount)), gp.GRB.MINIMIZE)

    m.Params.OutputFlag = 0
    m.setParam('PreCrush', 1)
    m.optimize()
    end = time.time()
    fracsize=0
    for i in range(N):
        if d[i].X >0.0001 and d[i].X<0.99999:
            fracsize+=1
    writer.writerow([prefix+"/ad_grb_{:05d}".format(number) + ".log.npz",fracsize,Delta])


if __name__ == '__main__':
    main(sys.argv)
