# ---------------------------------------------------------------------
# Python Script 													  -
# -- Reading key numberds from log file 						      -
# ---------------------------------------------------------------------
import math
import numpy as np
import pandas as pd
from scipy import interpolate as intp
from matplotlib import pyplot as plt
import re
import os
import gc
import sys

# Log Law ------------------------------------------------------------
# z0 = 0.1		# [m]		suface roughness
# uTau = 1.089	# [m/s]		friction velocity 1.08 for k = 0.41 

Cplus = -6.5
z0 = 0.01		# [m]		suface roughness
tauW = 0.045    # [m2/s2]
uTau = math.sqrt(tauW)

karman = 0.40	# [null]	Von Karman constant
H = 100		    # [m]		boundary layer height
nu = 1.5e-5		# [m^2/s]	kinematic viscousity of air
rho = 1.225		# [kg/m^3]	density of air
res = 64		# [null]	domain resolution

# Velocity profile
Y = np.concatenate((np.linspace(0.5*0.1*H/7,0.1*H-0.5*0.1*H/7,7),
	np.linspace(0.1*H+0.5*0.9*H/57,H-0.5*0.9*H/57,57)))/H
uPlus = 1/karman*np.log((Y*H-z0)/z0) + Cplus
Utheory = uPlus * uTau

# Load file -----------------------------------------------------------
file = pd.read_csv("data_.csv")
y0 = file.loc[:, "Points:1"].values
y0 = y0/max(y0)
uMean0 = file.loc[:, "UMean:0"].values
uMean = np.interp(Y, y0, uMean0)

plt.semilogx(Y,Utheory,c='b',linestyle='-',linewidth=1,marker='none',markerfacecolor='none')
plt.semilogx(Y,uMean,c='r',linestyle='none',linewidth=1,marker='^',markerfacecolor='none')
plt.grid()
plt.savefig('test.png')



# file_Rsgs1 = open("ABL/Rsgs_MT200_LUST","r")
# file_Rsgs2 = open("ABL/Rsgs_MT210_LUST","r")
# data_Rsgs1 = file_Rsgs1.read()
# data_Rsgs2 = file_Rsgs2.read()
# Rsgs1 = re.findall(r'[-]?\d+\.?\d*[e0-9]*[-+]?\d*',data_Rsgs1)
# Rsgs2 = re.findall(r'[-]?\d+\.?\d*[e0-9]*[-+]?\d*',data_Rsgs2)
# Rsgs1 = np.array(Rsgs1,dtype=np.float32)
# Rsgs2 = np.array(Rsgs2,dtype=np.float32)

# file_Rij1 = open("ABL/Rij_MT200_LUST","r")
# file_Rij2 = open("ABL/Rij_MT210_LUST","r")
# data_Rij1 = file_Rij1.read()
# data_Rij2 = file_Rij2.read()
# R1 = re.findall(r'[-]?\d+\.?\d*[e0-9]*[-+]?\d*',data_Rij1)
# R2 = re.findall(r'[-]?\d+\.?\d*[e0-9]*[-+]?\d*',data_Rij2)
# R1 = np.array(R1,dtype=np.float32)
# R2 = np.array(R2,dtype=np.float32)

# file_U1 = open("ABL/UMean_MT200_LUST","r")	# Ufix = 24
# file_U2 = open("ABL/UMean_MT210_LUST","r")
# data_U1 = file_U1.read()
# data_U2 = file_U2.read()
# U1 = re.findall(r'[-]?\d+\.?\d*[e0-9]*[-+]?\d*',data_U1)
# U2 = re.findall(r'[-]?\d+\.?\d*[e0-9]*[-+]?\d*',data_U2)
# U1 = np.array(U1,dtype=np.float32)
# U2 = np.array(U2,dtype=np.float32)

# file_nut1 = open("ABL/nut_MT200_LUST","r")	# Ufix = 24
# file_nut2 = open("ABL/nut_MT210_LUST","r")
# data_nut1 = file_nut1.read()
# data_nut2 = file_nut2.read()
# nut1 = re.findall(r'[-]?\d+\.?\d*[e0-9]*[-+]?\d*',data_nut1)
# nut2 = re.findall(r'[-]?\d+\.?\d*[e0-9]*[-+]?\d*',data_nut2)
# nut1 = np.array(nut1,dtype=np.float32)
# nut2 = np.array(nut2,dtype=np.float32)

# # Compute Uave --------------------------------------------------------
# Ux1 = inputAve.computeU(U1,64,2,7)
# Rij1 = inputAve.computeR(R1,64,2,7)
# Rsgs1 = inputAve.computeR(Rsgs2,64,2,7)
# utau1 = (karman*Ux1[0])/(np.log(7.14/z0))
# print(utau1)

# Ux2 = inputAve.computeU(U2,64,2,7)
# Rij2 = inputAve.computeR(R2,64,2,7)
# Rsgs2 = inputAve.computeR(Rsgs2,64,2,7)
# utau2 = (karman*Ux2[0])/(np.log(7.14/z0))
# print(utau2)

# # Plotting ------------------------------------------------------------
# # gc.collect()

# plt.figure(num=1,figsize=(6, 4.5), dpi=120, facecolor='w', edgecolor='k')
# plt.gcf().subplots_adjust(bottom=0.2)
# plt.gcf().subplots_adjust(left=0.14)
# plt.rcParams.update({'font.size': 12})
# plt.semilogx(Y,Utheory,c='k',linestyle='-',linewidth=1.2,
# 	marker='',markerfacecolor='none',mew=0.8,ms=5)
# plt.semilogx(Y,Ux1/utau1,c='b',linestyle=':',linewidth=0.8,
# 	marker='',markerfacecolor='none',mew=0.8,ms=5)
# plt.semilogx(Y,Ux2/utau2,c='r',linestyle='-.',linewidth=0.8,
# 	marker='',markerfacecolor='none',mew=0.8,ms=5)
# plt.legend([
# 	'$U_{x}$ Log law',
# 	'$U_{x}$ ABL BC0',
# 	'$U_{x}$ ABL BC1'
# 	])
# plt.ylim((10,25))
# plt.xlabel('y/H')
# plt.ylabel('u+')
# plt.grid()


# plt.figure(num=2,figsize=(6, 4.5), dpi=120, facecolor='none', edgecolor='k')
# plt.gcf().subplots_adjust(bottom=0.2)
# plt.gcf().subplots_adjust(left=0.14)
# plt.rcParams.update({'font.size': 12})
# plt.plot(Y,Rij2[:,0]/(utau2**2),c='k',linestyle='-.',linewidth=0.8,
# 	marker='',markerfacecolor='none',mew=0.8,ms=5)
# plt.plot(Y,Rsgs2[:,0]/(utau2**2),c='k',linestyle=':',linewidth=0.8,
# 	marker='',markerfacecolor='none',mew=0.8,ms=5)
# plt.plot(Y,(Rij2[:,0]+Rsgs2[:,0])/(utau2**2),c='k',linestyle='-',linewidth=1,
# 	marker='',markerfacecolor='none',mew=0.8,ms=5)

# plt.plot(Y,Rij2[:,3]/(utau2**2),c='r',linestyle='-.',linewidth=0.8,
# 	marker='',markerfacecolor='none',mew=0.8,ms=5)
# plt.plot(Y,Rij2[:,5]/(utau2**2),c='b',linestyle='-.',linewidth=0.8,
# 	marker='',markerfacecolor='none',mew=0.8,ms=5)
# plt.plot(Y,Rij2[:,1]/(utau2**2),c='g',linestyle='-.',linewidth=0.8,
# 	marker='',markerfacecolor='none',mew=0.8,ms=5)

# plt.plot(Y,Rsgs2[:,3]/(utau2**2),c='r',linestyle=':',linewidth=0.8,
# 	marker='',markerfacecolor='none',mew=0.8,ms=5)
# plt.plot(Y,Rsgs2[:,5]/(utau2**2),c='b',linestyle=':',linewidth=0.8,
# 	marker='',markerfacecolor='none',mew=0.8,ms=5)
# plt.plot(Y,Rsgs2[:,1]/(utau2**2),c='g',linestyle=':',linewidth=0.8,
# 	marker='',markerfacecolor='none',mew=0.8,ms=5)

# plt.plot(Y,(Rij2[:,3]+Rsgs2[:,3])/(utau2**2),c='r',linestyle='-',linewidth=1,
# 	marker='',markerfacecolor='none',mew=0.8,ms=5)
# plt.plot(Y,(Rij2[:,5]+Rsgs2[:,5])/(utau2**2),c='b',linestyle='-',linewidth=1,
# 	marker='',markerfacecolor='none',mew=0.8,ms=5)
# plt.plot(Y,(Rij2[:,1]+Rsgs2[:,1])/(utau2**2),c='g',linestyle='-',linewidth=1,
# 	marker='',markerfacecolor='none',mew=0.8,ms=5)
# plt.legend([
# 	'$R_{resolved}$',
# 	'$R_{subgrid}$',
# 	'$R_{total}$'
# 	])
# plt.ylabel(r'$R_{ij}/U_{\tau}^2$')
# plt.xlabel('y/H')
# plt.grid()

# plt.show()