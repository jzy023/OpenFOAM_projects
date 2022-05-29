# ---------------------------------------------------------------------
# Python Script 													  -
# -- initialize the mean Ux profile in the log law region		      -
# ---------------------------------------------------------------------
# import math
import numpy as np
import pandas as pd
from scipy import interpolate as intp
from matplotlib import pyplot as plt
import csv
import sys
import re

# Parameters ---------------------------------------------------------
H = 1    		# [m]		boundary layer height
z0 = 0.1		# [m]		suface roughness
utau = 5.5   	# [m/s]		friction velocity 1.0635
karman = 0.41	# [null]	Von Karman constant
nu = 1e-5		# [m^2/s]	kinematic viscousity of air
rho = 1.225		# [kg/m^3]	density of air
W = np.pi*H
L = 2*np.pi*H

# Grid generation ----------------------------------------------------
C = 2
res = 181
Z = np.linspace(0,1,res)
Z = H*np.flip(1-np.tanh(C*Z)/np.tanh(C))
n_bad = np.sum(Z[:]<z0)
Z = Z[n_bad+1:len(Z)]
# Velocity profile
Ux = utau/karman*np.log(Z/z0)
# validation
# yplus = (utau*Z)/nu
# uplus = Ux/utau
# plt.figure(1)
# plt.semilogx(yplus,uplus,'*')
# plt.xlabel('y+')
# plt.ylabel('u+')
# plt.grid() 
# plt.show()

# Read Files --------------------------------------------------------
file_L0 = open("L.orig","r")
L0 = file_L0.read()
L0 = re.findall(r'\d+\.?[e0-9]*[-+]?\d*',L0)
L0 = np.array(L0,dtype=np.float32)
print(L0)

file_R0 = open("R.orig","r")
R0 = file_R0.read()
R0 = re.findall(r'\d+\.?[e0-9]*[-+]?\d*',R0)
R0 = np.array(R0,dtype=np.float32)
# R0 = R0.astype(np.float)
R0 = np.reshape(R0,(len(L0),int(len(R0)/len(L0))))
print(R0)

# Modify R ----------------------------------------------------------
R = R0
for i in range(0,len(R0)):
	for j in range(0,len(R0[0])):
		R[i][j] = R[i][j]*0.01
	pass
pass

# Write Files -------------------------------------------------------
file_Z = open('points',"w")
file_Z.write("(\n")
for i in range(0,len(Z)):
	file_Z.write("(%i %.8e %i)\n" %(0,Z[i],0))
	pass
file_Z.write(")")
file_Z.close()

file_U = open('U',"w")
file_U.write("(\n")
for i in range(0,len(Ux)):
	file_U.write("(%.8e %i %i)\n" %(Ux[i],0,0))
	pass
file_U.write(")")
file_U.close()

file_R = open('R',"w")
file_R.write("(\n")
for i in range(0,len(Z)):
	file_R.write("(")
	for j in range(0,len(R0[0])):	
		file_R.write("%.8e " %R[i][j])
	pass
	file_R.write(")\n")
pass
file_R.write(")")
file_R.close()

# Bulk Velocity
U_ = (Ux[0:len(Ux)-1]+Ux[1:len(Ux)])/2
Z_ = Z[1:len(Z)]-Z[0:len(Z)-1]
Ub = np.sum(np.multiply(U_,Z_))/H
print("Bulk velocity = %.4f" % Ub)

