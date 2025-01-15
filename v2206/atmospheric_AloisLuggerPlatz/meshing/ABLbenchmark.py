# ---------------------------------------------------------------------
# Python Script 													  -
# -- Reading key numberds from log file 						      -
# ---------------------------------------------------------------------
import math
import numpy as np
import pandas as pd
from scipy import interpolate as intp
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
# import re
# import os
# import gc
# import sys


# Log Law ------------------------------------------------------------
Cplus = -6.5
z0 = 0.01		# [m]		suface roughness
tauW = 0.045    # [m2/s2]
uTau = math.sqrt(tauW)

karman = 0.40	# [null]	Von Karman constant
H = 100		    # [m]		boundary layer height
nu = 1.5e-5		# [m^2/s]	kinematic viscousity of air
rho = 1.225		# [kg/m^3]	density of air
res = 64		# [null]	domain resolution


# Log data -----------------------------------------------------------
# Velocity profile
Y = np.concatenate((np.linspace(0.5*0.1*H/7,0.1*H-0.5*0.1*H/7,7),
	np.linspace(0.1*H+0.5*0.9*H/57,H-0.5*0.9*H/57,57)))/H
uPlus = 1/karman*np.log((Y*H-z0)/z0) + Cplus
Utheory = uPlus * uTau

file = pd.read_csv("data_uProfile.csv")
y0 = file.loc[:, "Points:1"].values
y0 = y0/max(y0)
uMean0 = file.loc[:, "UMean:2"].values
uMean = np.interp(Y, y0, uMean0)


# Plot ---------------------------------------------------------------
# Setting Formats
fontSize_axis_label = 20
fontSize_axis_ticks = 18
fontSize_title = 16
fontSize_legend = 14

format = mtick.ScalarFormatter(useMathText=True)
format.set_useOffset(True)
format.set_scientific(True) 
format.set_powerlimits((-1,1))
position = [0.145,0.14,0.825,0.80]

# Plotting
plt.figure(figsize=(8, 5))
plt.grid()
# logging
plt.semilogx(Y,Utheory,'k-',linewidth=1,markersize='8',markerfacecolor='none')
plt.semilogx(Y,uMean,'r^',linewidth=1,markersize='8',markerfacecolor='none')
# formating
plt.xlabel(r"$y/H$", fontsize = fontSize_axis_label)
plt.ylabel(r"$u^{+}$", fontsize = fontSize_axis_label)
plt.legend(['log-law theory', 'simulation'], 
            fontsize=fontSize_legend,
            loc="lower right")  # Set legend font size
plt.tick_params(axis='both', which='major', labelsize=16)  # Set tick font size

plt.subplots_adjust(left=0.12, bottom=0.16, right=0.95, top=0.9) 
plt.savefig('u_profile_semi.png', dpi = 200)



plt.figure(figsize=(8, 5))
plt.grid()
# logging
plt.plot(Y,Utheory,'k-',linewidth=1,markersize='8',markerfacecolor='none')
plt.plot(Y,uMean,'r^',linewidth=1,markersize='8',markerfacecolor='none')
# formating
plt.xlabel(r"$y/H$", fontsize = fontSize_axis_label)
plt.ylabel(r"$u^{+}$", fontsize = fontSize_axis_label)
plt.legend(['log-law theory', 'simulation'], 
            fontsize=fontSize_legend,
            loc="lower right")  # Set legend font size
plt.tick_params(axis='both', which='major', labelsize=16)  # Set tick font size

plt.subplots_adjust(left=0.12, bottom=0.16, right=0.95, top=0.9)
plt.savefig('u_profile.png', dpi = 200)


