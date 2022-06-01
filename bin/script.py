import os
for i in range(1,5):
    os.system("nohup ./HMSearchApp %s %s &"%(i,i+1))
