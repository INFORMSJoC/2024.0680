import csv
import glob
import os





loN = [32, 64, 96]
dict_alpha_index = {32:range(5), 64:range(5), 96:range(1,4)}
timelimit = 1

for N in loN:
    fh = open("../data/sizefrac%d.csv" % (N), 'w')
    writer = csv.writer(fh)
    writer.writerow(["Name", "fracsize","Delta"])
    fh.close()
    for i in dict_alpha_index[N]:
        prefix = "../instances/%d_%d" % (N, i)
      
        myinstlist = glob.glob("*.log.npz", root_dir=prefix +"/")        
        for i in range(len(myinstlist)):
            
              
            os.system("python calc_frac_size.py {} {} {}".format(
                    prefix, N,i)
                )
