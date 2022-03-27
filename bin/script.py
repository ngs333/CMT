import os
for i in range(1,100,10):
    os.system("nohup ./HMSearchApp %s %s &"%(i,9))
