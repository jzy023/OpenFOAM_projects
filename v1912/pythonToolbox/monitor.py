# ---------------------------------------------------------------------
# Python Script 													  -
# -- Reading key numberds from log file 						      -
# ---------------------------------------------------------------------
# import math
import numpy as np
import pandas as pd
from scipy import interpolate as intp
from matplotlib import pyplot as plt
import re
import os
import gc
 

# Load file -----------------------------------------------------------
file = open("Output.log","r")
# file = open("output_test.log","r")
output = file.read()

cpuTime_str = re.findall(r'Time = \d+\.?\d+\n',output)
deltaT_str = re.findall(r'deltaT = \d+\.?\d+',output)
cumuErr_str = re.findall(r'cumulative = [-]?\d+\.?[e0-9]?[-]?\d+',output)
smoothRes_str = re.findall(r'smoothSolver.+Final residual.+,',output)
ubar_str = re.findall(r'uncorrected Ubar = [-]?\d+\.?\d+',output)
# dGrP_str = re.findall(r'pressure gradient = [-]?\d+\.?\d+[e]?[-+]?\d+',output))
dGrP_str = re.findall(r'dGradP_ = [-]?\d+\.?[e0-9]?[-]?\d+',output)

l1 = len(cpuTime_str)
l2 = len(deltaT_str)
l3 = int(len(cumuErr_str)/2)
l4 = int(len(smoothRes_str)/4)
l5 = int(len(ubar_str)/3)

print(l1)
print(l2)
print(l3)
print(l4)
print(l5)
if l1 != l2:
	print('findall criteria does NOT match')

# CPU time log --------------------------------------------------------
file_cpuTime = open('cpuTime.log',"w")
for i in range(0,l1):
	file_cpuTime.write(cpuTime_str[i])
	file_cpuTime.write("\n")
	pass
file_cpuTime.close()
file_cpuTime = open('cpuTime.log',"r")
cpuTime = file_cpuTime.read()
cpuTime = re.findall(r'\d+\.?\d+',cpuTime)
cpuTime = np.array(cpuTime,dtype=np.float32)

# Delta T log ---------------------------------------------------------
file_deltaT = open('deltaT.log',"w")
for i in range(0,l2):
	file_deltaT.write(deltaT_str[i])
	file_deltaT.write("\n")
	pass
file_deltaT.close()
file_deltaT = open('deltaT.log',"r")
deltaT = file_deltaT.read()
deltaT = re.findall(r'\d+\.?\d+',deltaT)
deltaT = np.array(deltaT,dtype=np.float32)

# Cumulative Err log --------------------------------------------------
file_cumuErr = open('cumuErr.log',"w")
for i in range(0,l3):
	file_cumuErr.write(cumuErr_str[i*2])
	file_cumuErr.write("\n")
	pass
file_cumuErr.close()
file_cumuErr = open('cumuErr.log',"r")
cumuErr = file_cumuErr.read()
cumuErr = re.findall(r'[-]?\d+\.?[e0-9]?[-]?\d+',cumuErr)
cumuErr = np.array(cumuErr,dtype=np.float32)

# Residual log --------------------------------------------------------
file_uxRes = open('uxRes.log',"w")
file_uyRes = open('uyRes.log',"w")
file_uzRes = open('uzRes.log',"w")
file_kRes = open('kRes.log',"w") 
for i in range(0,l4):
	file_uxRes.write(smoothRes_str[i*4])
	file_uxRes.write("\n")
	file_uyRes.write(smoothRes_str[i*4+1])
	file_uyRes.write("\n")
	file_uzRes.write(smoothRes_str[i*4+2])
	file_uzRes.write("\n")
	file_kRes.write(smoothRes_str[i*4+3])
	file_kRes.write("\n")
	pass
file_uxRes.close()
file_uyRes.close()
file_uzRes.close()
file_kRes.close()
file_uxRes = open('uxRes.log',"r")
file_uyRes = open('uyRes.log',"r")
file_uzRes = open('uzRes.log',"r")
file_kRes = open('kRes.log',"r")

uxRes = file_uxRes.read()
uxRes = re.findall(r'\d+\.?\d*[a-zA-Z]?[-+]?\d+',uxRes)
uxRes = np.array(uxRes,dtype=np.float32)

uyRes = file_uyRes.read()
uyRes = re.findall(r'\d+\.?\d*[a-zA-Z]?[-+]?\d+',uyRes)
uyRes = np.array(uyRes,dtype=np.float32)

uzRes = file_uzRes.read()
uzRes = re.findall(r'\d+\.?\d*[a-zA-Z]?[-+]?\d+',uzRes) 
uzRes = np.array(uzRes,dtype=np.float32)

kRes = file_kRes.read()
kRes = re.findall(r'\d+\.?\d*[a-zA-Z]?[-+]?\d+',kRes)
kRes = np.array(kRes,dtype=np.float32)

uxRes_0 = np.zeros(l4)
uyRes_0 = np.zeros(l4)
uzRes_0 = np.zeros(l4)
kRes_0 = np.zeros(l4)

for i in range(0,l4):
	uxRes_0[i] = uxRes[i*2+1]
	uyRes_0[i] = uyRes[i*2+1]
	uzRes_0[i] = uzRes[i*2+1]
	kRes_0[i] = kRes[i*2+1]
	pass

# Ubar log -----------------------------------------------------------
file_Ubar = open('Ubar.log',"w")
for i in range(0,l5):
	file_Ubar.write(ubar_str[i*3])
	file_Ubar.write("\n")
	pass
file_Ubar.close()
file_Ubar = open('Ubar.log',"r")
Ubar = file_Ubar.read()
Ubar = re.findall(r'[-]?\d+\.?\d+',Ubar)
Ubar = np.array(Ubar,dtype=np.float32)

# Ubar log -----------------------------------------------------------
file_dGrP = open('dGrP.log',"w")
for i in range(0,l5):
	file_dGrP.write(dGrP_str[i*3])
	file_dGrP.write("\n")
	pass
file_dGrP.close()
file_dGrP = open('dGrP.log',"r")
dGrP = file_dGrP.read()
dGrP = re.findall(r'[-]?\d+\.?\d+',dGrP)
dGrP = np.array(dGrP,dtype=np.float32)

# Plotting ------------------------------------------------------------
gc.collect()

# plt.figure(1)
# plt.plot(cpuTime,deltaT)
# plt.xlabel('cpuTime [s]')
# plt.ylabel('deltaT [s]')
# plt.grid()
# plt.show()

# plt.figure(2)
# plt.plot(cpuTime,abs(cumuErr))
# plt.xlabel('cpuTime [s]')
# plt.ylabel('cumuErr')
# plt.grid()
# plt.show()

plt.figure(3)
plt.plot(cpuTime[0:200],Ubar[0:200])
plt.xlabel('cpuTime [s]')
plt.ylabel('Ubar')
plt.grid()
plt.show()

plt.figure(4)
plt.plot(cpuTime[0:200],dGrP[0:200])
plt.xlabel('cpuTime [s]')
plt.ylabel('gradP')
plt.grid()
plt.show()

# plt.figure(4)
# plt.plot(cpuTime,uxRes_0,'-k',label='Ux')
# plt.plot(cpuTime,uyRes_0,'-r',label='Uy')
# plt.plot(cpuTime,uzRes_0,'-b',label='Uz')
# plt.plot(cpuTime,kRes_0,'-g',label='k')
# plt.xlabel('cpuTime [s]')
# plt.ylabel('Residual')
# plt.legend(loc='upper right')
# plt.grid()
# plt.show()

