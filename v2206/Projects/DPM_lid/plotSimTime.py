import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import re
import os
import gc

fileName = "log_ErgunWenYuDrag"
fileName = "log"
file = open(fileName,"r")
data = file.read()

simTimeStr = re.findall(r'Time = \d+\.?\d+\n',data)
exeTimeStr = re.findall(r'ExecutionTime = \d+\.?\d+\s',data)

# print(len(simTimeStr),len(exeTimeStr))
# print(simTimeStr[1],exeTimeStr[1])

simTime = []
for line in simTimeStr:
    row = re.findall(r'\d+\.?\d*',line)
    simTime.append(float(row[0]))

exeTime = []
for line in exeTimeStr:
    row = re.findall(r'\d+\.?\d*',line)
    exeTime.append(float(row[0]))

sampleLen = int(95*len(simTime)/100)
plt.figure(dpi=800)
plt.plot(simTime[:sampleLen],exeTime[:sampleLen])
# plt.plot(simTime,exeTime)
plt.xlabel("simulation time [s]")
plt.ylabel("execution time [s]")
plt.grid()
plt.savefig('fig.jpeg')

# for x in range(2,4):
#     print(x)