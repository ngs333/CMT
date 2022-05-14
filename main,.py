import pandas as pd
import numpy as np
import os

data = {}
for i in os.listdir():
    if not i.endswith(".py"):
        df = pd.read_csv(i)
        df = df[df[" activity inference"] != 3]
        uid=i.split("_")[1].split(".")[0]
        data[uid] = df[" activity inference"].value_counts()[0]/len(df)


print(data)
~                                                        