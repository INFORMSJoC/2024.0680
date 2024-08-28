import csv
import glob
import os

from pandas import read_csv

loN = [32, 64, 96]
dict_alpha_index = {32:range(5), 64:range(5), 96:range(1,4)}
timelimit = 1

for N in loN:
    for i in dict_alpha_index[N]:
        prefix = "../instances/%d_%d" % (N, i)
        fh = open("../data/%d_%d.csv" % (N, i), 'w')
        writer = csv.writer(fh)
        if N == 96:
            writer.writerow(["Name", "P", "B", "C1", "C2", "Objective without constant", 
                         "Objective", "Gap", "LB", "Time"])
        else:
            writer.writerow(["Name", "P", "B", "C1", "C2", "Objective without constant", 
                         "Objective", "Gap", "Time"])
        fh.close()
        
        myinstlist = glob.glob("*.log.npz", root_dir=prefix +"/")        
        if N == 96:
            combinations = [ [0,0,0], [1,1,1] ]
        else:
            combinations = [
                [0,0,0], [1,0,0], [0,1,0], [0,0,1], [1,1,0], [1,0,1], [0,1,1], [1,1,1] 
            ]
        for combination in combinations:
            for j in range(len(myinstlist)):
                number=j
                os.system("python calc_time_ip.py {} {} {} {} {} {} {} {}".format(
                    prefix, N, i, combination[0], combination[1], combination[2], number, timelimit)
                )  

headers = ["Name"]

df1 = read_csv('../instances/hard_instances_for_none.csv',names=headers).values
df2 = read_csv('../instances/hard_instances_for_pbc.csv',names=headers).values

count1=0
count2=0

fh = open("../data/hard_instances_for_none_runtime.csv", 'w')
writer = csv.writer(fh)
writer.writerow(["Name", "P", "B", "C1", "C2", "Objective without constant", "Objective", "Gap", "LB", "Time"])
fh.close()
fh = open("../data/hard_instances_for_pbc_runtime.csv", 'w')
writer = csv.writer(fh)
writer.writerow(["Name", "P", "B", "C1", "C2", "Objective without constant", "Objective", "Gap", "LB", "Time"])
fh.close()

while count1+count2<len(df1)+len(df2):
    if count1<len(df1):
        item = df1[count1]
        os.system("python calc_time_ip.py {} {} {} {} {} {} {} {}".format(
                    "../"+item[0], -1, -1, 0, 0, 0, -1, 3*timelimit)
                )
        count1+=1
        
    if count2<len(df2):
        item = df2[count2]
        os.system("python calc_time_ip.py {} {} {} {} {} {} {} {}".format(
                    "../"+item[0], -2, -2, 1, 1, 1, -1, 3*timelimit)
                )
        count2+=1

    
