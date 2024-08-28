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


def mycallback(model, where):
    if where == GRB.Callback.MIPSOL:
        """calculates the improved solution"""
        if model._mycounter >= 0 and model._primal == 1:
            n = model._n
            c = model._c
            x = model._x
            myd = [int(model.cbGetSolution(model._d[i])) for i in range(n * n)]
            Delta = model._Delta
            zetaxval = 0
            zetayval = 0
            for j in range(n):
                for i in range(n-1):
                    zetaxval = zetaxval + \
                        abs(model._dsign[i+1+j*n] * myd[i+1+j*n] -
                            model._dsign[i+j*n] * myd[i+j*n] + x[i+1+j*n] - x[i+j*n])
                    zetayval = zetayval + abs(model._dsign[i*n+j] * myd[i*n+j] - model._dsign[(
                        i+1)*n+j] * myd[(i+1)*n+j] + x[i+1+j*n] - x[i+j*n])
            oldval = sum([c[i] * model._dsign[i] * myd[i]
                         for i in range(n * n)]) + zetaxval + zetayval
            if oldval < model._val2:
                bestsol = [int(x[i] + model._dsign[i] * myd[i])
                           for i in range(n * n)]
                result = top_improve.top_improve(
                    n, n * n, model._Delta, 2, c, x, bestsol, range(2))
                zetaxval = 0
                zetayval = 0
                for j in range(n):
                    for i in range(n-1):
                        zetaxval = zetaxval + \
                            abs(result[i+1+j*n] - result[i+j*n])
                        zetayval = zetayval + \
                            abs(result[i*n+j] - result[(i+1)*n+j])
                result = [result[i] - x[i] for i in range(n * n)]
                absresult = [abs(result[i]) for i in range(n * n)]
                mysum = sum([c[i] * result[i]
                            for i in range(n * n)]) + zetaxval + zetayval
                if mysum < min(model._val, oldval):
                    model._val2 = mysum
                    model._posd = absresult
                    model._mycountersec = 1

    if where == GRB.Callback.MIPNODE:
        """add a solution if a better solution was found and not added yet"""
        if model._mycountersec == 1:
            model._mycountersec = 0
            model.cbSetSolution(model._d, model._posd)
            model._val = model.cbUseSolution()

    if where == GRB.Callback.MIPNODE and model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL and model._boxcut == 1:
        model._mycounter4 += 1
        if model._mycounter4 % (5 + 5 * model._mycounter2) == 1:
            model._mycounter4 = 1
            n = model._n
            c = model._c
            x = model._x
            d = model._d
            Delta = model._Delta
            myzetadic = model._myzetadic
            myd = [model.cbGetNodeRel(d[i]) for i in range(n * n)]
            nodes = []
            x1nodes = []
            x0nodes = []
            remcap = Delta
            notnodes = []
            for i in range(n * n):
                conti = True
                """check if integer valued"""
                for j in range(2):
                    if abs(x[i] + model._dsign[i] * myd[i] - j) < 0.0000001:
                        conti = False
                """ if fractional continue"""
                if conti:
                    if x[i] >= -0.000001 and x[i] <= 0.000001:
                        nodes.append(i)
                        x0nodes.append(i)
                    if x[i] >= 1 - 0.000001 and x[i] <= 1 + 0.000001:
                        x1nodes.append(i)
                        nodes.append(i)
            """determine grid"""
            left = n
            right = 0
            bot = n
            top = 0
            for i in nodes:
                left = min(left, i % n)
                right = max(right, i % n)
                top = max(top, math.floor(i / n))
                bot = min(bot, math.floor(i / n))
            gridnodes = []
            x1count = 0
            x0count = 0
            top = top + 1
            right = right + 1
            for i in range(bot, top):
                for j in range(left, right):
                    gridnodes.append(i * n + j)
                    if x[i*n+j] == 1:
                        x1count += 1
                    else:
                        x0count += 1
            # if fractional component is too small compared to grid, do not add cut
            if len(nodes) > 0.75 * len(gridnodes):
                if x1count < x0count:
                    complementcount = x1count
                    takenval = 0
                else:
                    complementcount = x0count
                    takenval = 1
                """remove capacity spent outside fractional component from remcap"""
                for i in range(n * n):
                    if i not in gridnodes:
                        if abs(myd[i]) > 0.99999:
                            notnodes.append(i)
                            remcap -= round(abs(myd[i]))
                gridnodes1 = []
                gridnodes2 = []
                border = math.ceil(math.sqrt(remcap))
                for i in range(bot, top):
                    for j in range(left, right):
                        if i - bot > border and top - i > border and j - left > border and right - j > border:
                            gridnodes2.append(i * n + j)
                        else:
                            gridnodes1.append(i * n + j)
                mykeys = []
                for i in range(bot, top):
                    for j in range(left, right-1):
                        mykeys.append((j+n*i, j+n*i+1))
                for i in range(bot, top-1):
                    for j in range(left, right):
                        mykeys.append((j+n*i, j+n*i+n))

                lowerdel = math.floor(math.sqrt(remcap)) * \
                    math.floor(math.sqrt(remcap))

                factor = min(top-bot, right-left, (math.ceil(math.sqrt(lowerdel)) + round(math.sqrt(lowerdel))),
                             2*(math.sqrt(max(0, (top-bot) * (right-left) - lowerdel-complementcount))))/lowerdel
                for i in range(lowerdel+1, remcap+1):
                    factor = min(factor, min(top-bot, right-left, (math.ceil(math.sqrt(i)) + round(
                        math.sqrt(i))), 2*(math.sqrt(max(0, (top-bot) * (right-left) - i-complementcount))))/i)

                factor2 = 2
                # only add cut if for lp relaxation lhs is significantly smaller than rhs
                if factor * sum(abs(myd[i]) for i in gridnodes) > 1.1**(model._mycounter2+1) * sum(model.cbGetNodeRel(myzetadic[key]) for key in mykeys):
                    model._mycounter2 += 1
                    model.cbLazy(gp.quicksum(factor * d[i] for i in gridnodes if x[i] == takenval) <= (factor+2) * (Delta-remcap) - (
                        factor+factor2) * gp.quicksum(d[i] for i in notnodes) + gp.quicksum(myzetadic[key] for key in mykeys))

    if where == GRB.Callback.MIPNODE and model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL and model._fullyconnected == 1:
        model._mycounter5 += 1
        if model._mycounterfull <= 2 and model._mycounter5 % (10+10*model._mycounterfull) == 1:
            model._mycounter5 = 1
            n = model._n
            c = model._c
            x = model._x
            d = model._d
            Delta = model._Delta
            myzetadic = model._myzetadic
            myd = [model.cbGetNodeRel(d[i]) for i in range(n*n)]
            nodes = []
            x1nodes = []
            x0nodes = []
            remcap = Delta  # capacity in the fractional component
            notnodes = []
            for i in range(n*n):
                conti = True
                """check if integer valued"""
                for j in range(2):
                    if abs(x[i]+model._dsign[i]*myd[i]-j) < 0.0000001:
                        conti = False
                        remcap -= myd[i]
                """ if fractional continue"""
                if conti:
                    if x[i] >= -0.000001 and x[i] <= 0.000001:
                        nodes.append(i)
                        x0nodes.append(i)
                    if x[i] >= 1-0.000001 and x[i] <= 1+0.000001:
                        x1nodes.append(i)
                        nodes.append(i)
            if len(x1nodes) > len(x0nodes):
                biggestlist = x1nodes
            else:
                biggestlist = x0nodes
            if remcap > 0:
                if 0.25*len(biggestlist) < remcap < 0.75*len(biggestlist):
                    addcost = calc_weights.calc_weights(
                        len(biggestlist), biggestlist, n, n*n, remcap, Delta)
                    newset = set()
                    for key in addcost.keys():
                        if key[0] > -1:
                            newset.add(key[0])
                            newset.add(key[1])
                    a1 = max(0, remcap * (len(newset)-remcap))
                    ell = min(len(newset), Delta)
                    factor = a1/remcap
                    if ell-remcap > 0:
                        a2 = max(0, ell * (len(newset)-ell))
                        factor2 = (a1-a2)/(ell-remcap) + factor
                    else:
                        factor2 = 0
                    mysum = sum([-factor * abs(myd[i]) for i in newset]) \
                        + sum([0.5*addcost[key]*model.cbGetNodeRel(myzetadic[key]) for key in addcost.keys() if key[0] > -1]) \
                        - sum([factor2 * abs(myd[i]) for i in notnodes])
                    model.cbLazy(gp.quicksum(factor * d[i] for i in newset) <= gp.quicksum(0.5 * addcost[key] * myzetadic[key]
                                 for key in addcost.keys() if key[0] > -1) - gp.quicksum(factor2 * d[i] for i in notnodes) + factor2*(Delta-remcap))
                    model._mycounterfull += 1


def main(argv):
    prefix = argv[1]
    N_axis = int(argv[2])
    i_alpha = int(argv[3])
    primal = int(argv[4])
    branch = int(argv[5])
    boxcut = int(argv[6])
    number = int(argv[7])
    timelimit = int(argv[8])
    fullyconnected = int(argv[6])
    if N_axis == -1:
        writer = csv.writer(open("../data/hard_instances_for_none_runtime.csv", 'a'))
        npzfile = np.load(prefix)
    elif N_axis == -2: 
        writer = csv.writer(open("../data/hard_instances_for_pbc_runtime.csv", 'a'))
        npzfile = np.load(prefix)
    else:
        writer = csv.writer(open("../data/%d_%d.csv" % (N_axis, i_alpha), 'a'))
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
        d[i] = m.addVar(name="d"+str(i), vtype="I", lb=0, ub=1)
        if x[i] > 0.99:
            dsign[i] = -1

    if branch == 1:
        for i in range(0, n):
            for j in range(0, n):
                prio = 0
                if i % 2 == 0:
                    prio += 8
                d[i*n+j].setAttr("BranchPriority", prio)

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

    m._d = d
    m._x = x
    m._c = c
    m._Delta = Delta
    m._dsign = dsign
    m._n = n
    m._myzetadic = myzetadic

    m.Params.OutputFlag = 0
    m.setParam('PreCrush', 1)
    m.Params.TimeLimit = timelimit
    m._modus = 0
    m._mycountersec = 0
    m._mycounter4 = 0
    m._mycounter5 = 0
    m._val2 = float('inf')
    m._val = float('inf')

    m._primal = primal
    m._fullyconnected = fullyconnected
    m._boxcut = boxcut

    m._mycounter = 0
    m._mycounter2 = 0
    m._mycounterfull = 0
    if boxcut == 1 or fullyconnected == 1:
        m.Params.LazyConstraints = 1
    m.optimize(mycallback)

    end = time.time()
    reg_offset = npzfile["reg_offset"]
    gap = m.MIPGap
    obj = m.getObjective().getValue()
    if N_axis<0:
        writer.writerow([prefix, primal,
                            branch, fullyconnected, boxcut, obj, obj-reg_offset, gap, m.ObjBound, end-start])
    else:
        if N_axis == 96:
            writer.writerow([prefix+"/ad_grb_{:05d}".format(number) + ".log.npz", primal,
                            branch, fullyconnected, boxcut, obj, obj-reg_offset, gap, m.ObjBound, end-start])
        else:
            writer.writerow([prefix+"/ad_grb_{:05d}".format(number) + ".log.npz", primal,
                            branch, fullyconnected, boxcut, obj, obj-reg_offset, gap, end-start])


if __name__ == '__main__':
    main(sys.argv)
