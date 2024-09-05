from pandas import read_csv
import pandas
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy import stats


## Plot Figure5a.
headers = ["Name","Primal","Branch","Cut1","Cut2","Obj","Obj-offset","gap","time"]
df1 =read_csv('../data/64_0.csv',skiprows=[0],names=headers)
df2= read_csv('../data/64_1.csv',skiprows=[0],names=headers)
df3 = read_csv('../data/64_2.csv',skiprows=[0],names=headers)
df4 =read_csv('../data/64_3.csv',skiprows=[0],names=headers)
df5= read_csv('../data/64_4.csv',skiprows=[0],names=headers)
df_merged1 = pandas.concat([df1, df2], ignore_index=True, sort=False)
df_merged2 = pandas.concat([df3, df4], ignore_index=True, sort=False)
df_merged = pandas.concat([df_merged1, df_merged2], ignore_index=True, sort=False)
df = pandas.concat([df_merged, df5], ignore_index=True, sort=False)

dfa = df["time"][(df["Primal"]==0) &  (df["Branch"]==0) &(df["Cut1"]==0) ]
mylist=[]
xlist=[]
total = (dfa[dfa<3700].count())
k=1
l=10
for i in range(0,600*k,l):
    xlist.append(1*i)
    mylist.append(dfa[dfa<1*i].count()/total)

dfb = df["time"][(df["Primal"]==1) &  (df["Branch"]==1) &(df["Cut1"]==1) ]
mylist2=[]
total2 = (dfb[dfb<3700].count())
for i in range(0,600*k,l):
    mylist2.append(dfb[dfb<1*i].count()/total2)

dfb = df["time"][(df["Primal"]==1) &  (df["Branch"]==0) &(df["Cut1"]==0) ]
mylist3=[]
total3 = (dfb[dfb<3700].count())
for i in range(0,600*k,l):
    mylist3.append(dfb[dfb<1*i].count()/total3)

dfb = df["time"][(df["Primal"]==0) &  (df["Branch"]==1) &(df["Cut1"]==0) ]
mylist4=[]
total4 = (dfb[dfb<3700].count())
for i in range(0,600*k,l):
    mylist4.append(dfb[dfb<1*i].count()/total4)

dfb = df["time"][(df["Primal"]==0) &  (df["Branch"]==0) &(df["Cut1"]==1) ]
mylist5=[]
total5 = (dfb[dfb<3700].count())
for i in range(0,600*k,l):
    mylist5.append(dfb[dfb<1*i].count()/total5)
  
dfb = df["time"][(df["Primal"]==1) &  (df["Branch"]==1) &(df["Cut1"]==0) ]
mylist6=[]
total6 = (dfb[dfb<3700].count())
for i in range(0,600*k,l):
    mylist6.append(dfb[dfb<1*i].count()/total6)

dfb = df["time"][(df["Primal"]==1) &  (df["Branch"]==0) &(df["Cut1"]==1) ]
mylist7=[]
total7= (dfb[dfb<3700].count())
for i in range(0,600*k,l):
    mylist7.append(dfb[dfb<1*i].count()/total7)

dfb = df["time"][(df["Primal"]==0) &  (df["Branch"]==1) &(df["Cut1"]==1) ]
mylist8=[]
total8 = (dfb[dfb<3700].count())
for i in range(0,600*k,l):
    mylist8.append(dfb[dfb<1*i].count()/total8)

plt.figure(figsize=(7, 7), dpi=300)
plt.plot(xlist,mylist2, label = "p-b-c")
plt.plot(xlist,mylist8, label = "b-c", marker = ".")
plt.plot(xlist,mylist5, label = "c", marker = "o")
plt.plot(xlist,mylist7, label = "p-c", marker = "v")
plt.plot(xlist,mylist4, label = "b" , marker = "s")
plt.plot(xlist,mylist6, label = "p-b", marker = "D")
plt.plot(xlist,mylist3, label = "p" ,marker = "*")
plt.plot(xlist,mylist, label = "none", marker = "^" )
plt.xlabel("time to reach optimality", fontsize=22)
plt.ylabel("fraction of instances solved", fontsize=22)
plt.tick_params(labelsize=22)
plt.legend(fontsize=22)
plt.tight_layout()
plt.savefig("../figures/figure5a.png", dpi=300)


## Plot Figure5b.
plt.figure(figsize=(7, 7), dpi=300)
headers = ["Name","Primal","Branch","Cut1","Cut2","Obj","Obj-offset","gap","lowerbound","time"]
df1= read_csv('../data/96_1.csv',skiprows=[0],names=headers)
df2 = read_csv('../data/96_2.csv',skiprows=[0],names=headers)
df3 =read_csv('../data/96_3.csv',skiprows=[0],names=headers)
df_merged = pandas.concat([df1, df2], ignore_index=True, sort=False)
df = pandas.concat([df_merged, df3], ignore_index=True, sort=False)

dfa = df["time"][df["Primal"]==0]
mylist=[]
xlist=[]
total = (dfa[dfa<3700].count())
for i in range(0,3600,60):
    xlist.append(1*i)
    mylist.append(dfa[dfa<1*i].count()/total)

dfb = df["time"][df["Primal"]==1]
mylist2=[]
xlist2=[]
total2 = (dfb[dfb<3700].count())
for i in range(0,3600,60):
    xlist2.append(1*i)
    mylist2.append(dfb[dfb<1*i].count()/total2)

plt.plot(xlist,mylist2, color = "red", label = "p-b-c")
plt.plot(xlist,mylist, label = "none", marker = "^" )
plt.xlabel("time to reach optimality", fontsize=22)
plt.ylabel("fraction of instances solved", fontsize=22)
plt.tick_params(labelsize=22)
plt.legend(fontsize=22)
plt.tight_layout()
plt.savefig("../figures/figure5b.png", dpi=300)


## Plot Figure6.
headers = ["Name","Primal","Branch","Cut1","Cut2","Obj","Obj-offset","gap","lowerbound","time"]
df1= read_csv('../data/96_1.csv',skiprows=[0],names=headers)
df2 = read_csv('../data/96_2.csv',skiprows=[0],names=headers)
df3 =read_csv('../data/96_3.csv',skiprows=[0],names=headers)
df_merged = pandas.concat([df1, df2], ignore_index=True, sort=False)
df = pandas.concat([df_merged, df3], ignore_index=True, sort=False)
dgap = df["gap"] [df["gap"] >= 0.001]

fig, axes = plt.subplots(figsize=(10, 9), dpi=300)
ax1=axes.violinplot(dataset = [df["gap"][(df["Primal"] == 0 ) &  (df["gap"] >= -0.001) ].values] , showmeans = True)
ax2=axes.violinplot(dataset = [ df["gap"][(df["Primal"] == 1)   &  (df["gap"] >= -0.001)].values ] , showmeans = True)
labels=[]
labels.append((mpatches.Patch(color=ax1["bodies"][0].get_facecolor().flatten()), "none"))
labels.append((mpatches.Patch(color=ax2["bodies"][0].get_facecolor().flatten()), "p-b-c"))
plt.legend(*zip(*labels), loc=2, fontsize=22)
plt.ylabel("duality gap", fontsize=22)
plt.tick_params(labelsize=22)
plt.tight_layout()
plt.savefig("../figures/figure6.png", dpi=300)


## Plot Figure7a.
headers = ["Name", "size", "Delta"]
df1 =read_csv('../data/sizefrac32.csv',skiprows=[0],names=headers)
headers2 = ["Name","Primal","Branch","Cut1","Cut2","Obj","Obj-offset","gap","time"]
dftemp0= read_csv('../data/32_0.csv',skiprows=[0],names=headers2)
dftemp1= read_csv('../data/32_1.csv',skiprows=[0],names=headers2)
dftemp2=  read_csv('../data/32_2.csv',skiprows=[0],names=headers2)
dftemp3= read_csv('../data/32_3.csv',skiprows=[0],names=headers2)
dftemp4= read_csv('../data/32_4.csv',skiprows=[0],names=headers2)
df_merged1 = pandas.concat([dftemp0, dftemp1], ignore_index=True, sort=False)
df_merged2 = pandas.concat([dftemp2, dftemp3], ignore_index=True, sort=False)
df_merged = pandas.concat([df_merged1, df_merged2], ignore_index=True, sort=False)
df = pandas.concat([df_merged, dftemp4], ignore_index=True, sort=False)
dffinal = df[ (df["Primal"] == 0) & (df["Branch"] == 0) & (df["Cut1"] == 0)]
dffinal.reset_index(drop=True,inplace=True)
df1.insert(3,"time",dffinal['time'])
y = df1["time"][df1["Delta"]==64].values
x = df1["size"][df1["Delta"]==64].values

plt.figure(figsize=(14, 10), dpi=300)
plt.plot(x,y,"bo", label="data points")
res = stats.linregress(x, y)
plt.plot(x, res.intercept + res.slope*x, 'r', label='linear regression')
plt.xlabel("size of the fractional component", fontsize=32)
plt.ylabel("run time", fontsize=32)
plt.legend(fontsize=32)
plt.tick_params(labelsize=32)
plt.tight_layout()
plt.savefig("../figures/figure7a.png", dpi=300)


## Plot Figure7b.
headers = ["Name", "size", "Delta"]
df1 =read_csv('../data/sizefrac64.csv',skiprows=[0],names=headers)
headers2 = ["Name","Primal","Branch","Cut1","Cut2","Obj","Obj-offset","gap","time"]
dftemp0= read_csv('../data/64_0.csv',skiprows=[0],names=headers2)
dftemp1= read_csv('../data/64_1.csv',skiprows=[0],names=headers2)
dftemp2=  read_csv('../data/64_2.csv',skiprows=[0],names=headers2)
dftemp3= read_csv('../data/64_3.csv',skiprows=[0],names=headers2)
dftemp4= read_csv('../data/64_4.csv',skiprows=[0],names=headers2)
df_merged1 = pandas.concat([dftemp0, dftemp1], ignore_index=True, sort=False)
df_merged2 = pandas.concat([dftemp2, dftemp3], ignore_index=True, sort=False)
df_merged = pandas.concat([df_merged1, df_merged2], ignore_index=True, sort=False)
df = pandas.concat([df_merged, dftemp4], ignore_index=True, sort=False)
dffinal = df[ (df["Primal"] == 0) & (df["Branch"] == 0) & (df["Cut1"] == 0)]
dffinal.reset_index(drop=True,inplace=True)
df1.insert(3,"time",dffinal['time'])
y = df1["time"][df1["Delta"]==64*4].values
x = df1["size"][df1["Delta"]==64*4].values

plt.figure(figsize=(14, 10), dpi=300)
plt.plot(x,y,"bo", label="data points")
res = stats.linregress(x, y)
plt.plot(x, res.intercept + res.slope*x, 'r', label='linear regression')
plt.xlabel("size of the fractional component", fontsize=32)
plt.ylabel("run time", fontsize=32)
plt.legend(fontsize=32)
plt.tick_params(labelsize=32)
plt.tight_layout()
plt.savefig("../figures/figure7b.png", dpi=300)


## Plot Figure7c.
headers = ["Name", "size", "Delta"]
df1 =read_csv('../data/sizefrac96.csv',skiprows=[0],names=headers)
headers2 = ["Name","Primal","Branch","Cut1","Cut2","Obj","Obj-offset","gap","lb","time"]
dftemp1= read_csv('../data/96_1.csv',skiprows=[0],names=headers2)
dftemp2=  read_csv('../data/96_2.csv',skiprows=[0],names=headers2)
dftemp3= read_csv('../data/96_3.csv',skiprows=[0],names=headers2)
df_merged1 = pandas.concat([dftemp1, dftemp2], ignore_index=True, sort=False)
df= pandas.concat([df_merged1, dftemp3], ignore_index=True, sort=False)
dffinal = df[ (df["Primal"] == 0) & (df["Branch"] == 0) & (df["Cut1"] == 0)]
dffinal.reset_index(drop=True,inplace=True)
df1.insert(3,"time",dffinal['time'])
y = df1["time"][df1["Delta"]==576].values
x = df1["size"][df1["Delta"]==576].values

plt.figure(figsize=(14, 10), dpi=300)
plt.plot(x,y,"bo", label="data points")
res = stats.linregress(x, y)
plt.plot(x, res.intercept + res.slope*x, 'r', label='linear regression')
plt.xlabel("size of the fractional component", fontsize=32)
plt.ylabel("run time", fontsize=32)
plt.legend(fontsize=32)
plt.tick_params(labelsize=32)
plt.tight_layout()
plt.savefig("../figures/figure7c.png", dpi=300)