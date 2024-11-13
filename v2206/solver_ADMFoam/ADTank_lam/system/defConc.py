import numpy as np 


names = ['Sac', 'Sh2', 'Gch4']
data = [8.93e-2, 2.505e-7, 1.65]

concDict = dict(zip(names, data))
for name in names:
    print(name, ": ", concDict[name])