import math
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from fractions import Fraction

## =================================================================
## Setting Formats
fontSize_axis_label = 16
fontSize_axis_ticks = 14
fontSize_title = 16
fontSize_legend = 12

## =================================================================
# number of cells on each direction
Nx = 200
Ny = 600

data = pd.read_csv("t40.csv")
# data = pd.read_csv("t60.csv")
# extract the location data from the columns
x = data['Points:0'].to_numpy()[:((Nx+1)*(Ny+1))]
y = data['Points:1'].to_numpy()[:((Nx+1)*(Ny+1))]
# extract the velocity data from the columns
rho = data['T'].to_numpy()[:((Nx+1)*(Ny+1))]
u = data['U:0'].to_numpy()[:((Nx+1)*(Ny+1))]
v = data['U:1'].to_numpy()[:((Nx+1)*(Ny+1))]

rhoWall = np.zeros((Ny+1,Nx+1))
rhoc    = np.zeros((Ny,Nx))
XWall   = np.zeros((Ny+1,Nx+1))
YWall   = np.zeros((Ny+1,Nx+1))
Xc      = np.zeros((Ny,Nx))
Yc      = np.zeros((Ny,Nx))
UWall   = np.zeros((Ny+1,Nx+1))
VWall   = np.zeros((Ny+1,Nx+1))
Uc      = np.zeros((Ny,Nx))
Vc      = np.zeros((Ny,Nx))


# build (U, V) and (X, Y) for on the walls
for i in range(Ny+1):
    for j in range(Nx+1):
        rhoWall[i][j] = rho[(Nx+1)*i + j]
        UWall[i][j] = u[(Nx+1)*i + j]
        VWall[i][j] = v[(Nx+1)*i + j]
        XWall[i][j] = x[(Nx+1)*i + j]
        YWall[i][j] = y[(Nx+1)*i + j]
# build (U, V) and (X, Y) for the cell center
for i in range(Ny):
    for j in range(Nx):
        rhoc[i][j] = (rhoWall[i][j] + rhoWall[i+1][j] + rhoWall[i][j+1] + rhoWall[i+1][j+1])/4
        Uc[i][j] = (UWall[i][j] + UWall[i+1][j] + UWall[i][j+1] + UWall[i+1][j+1])/4
        Vc[i][j] = (VWall[i][j] + VWall[i+1][j] + VWall[i][j+1] + VWall[i+1][j+1])/4
        Xc[i][j] = (XWall[i][j] + XWall[i+1][j] + XWall[i][j+1] + XWall[i+1][j+1])/4
        Yc[i][j] = (YWall[i][j] + YWall[i+1][j] + YWall[i][j+1] + YWall[i+1][j+1])/4

# geo information
c0X = np.min(Xc)
c0Y = np.min(Yc)
cmX = np.max(Xc)
cmY = np.max(Yc)

dcX = cmX - c0X
dcY = cmY - c0Y
cellX = dcX/Nx
cellY = dcY/Ny
Xcf = Xc.flatten()
Ycf = Yc.flatten()

# find wave number kappa
kappaNorm = 2*np.pi*np.max([cellX,cellY])
kappa = np.zeros(np.max([Nx,Ny]))
for i in range(len(kappa)):
    kappa[i] = kappaNorm*i

# find energy Ek
rhock = np.fft.fft2(rhoc,axes=(0, 1))/(Nx*Ny)
Uck = np.fft.fft2(Uc,axes=(0, 1))/(Nx*Ny)
Vck = np.fft.fft2(Vc,axes=(0, 1))/(Nx*Ny)

Ekrho = (np.abs(rhock))
Ekx = (np.abs(Uck**2))
Eky = (np.abs(Vck**2))

# calculating Ek along x axis
Erho = np.zeros(len(kappa))
Ex = np.zeros(len(kappa))
Ey = np.zeros(len(kappa))
# for x in range(Nx):
#     for y in range(Ny):
#         Erho[x] += Ekrho[y][x]
#         Ex[x] += Ekx[y][x]
#         Ey[x] += Eky[y][x]
# Erho /= Ny
# Ex /= Ny
# Ey /= Ny
for y in range(Ny):
    for x in range(Nx):
        Erho[y] += Ekrho[y][x]
        Ex[y] += Ekx[y][x]
        Ey[y] += Eky[y][x]
Erho /= Nx
Ex /= Nx
Ey /= Nx

# for c in range(len(Xcf)):
    # xx = ((Xcf[c] - c0X) - 0.5*dcX)/dcX*Nx
    # yy = ((Ycf[c] - c0Y) - 0.5*dcY)/dcY*Ny
    # ki = int(np.round(np.sqrt(xx*xx + yy*yy)))
    # E[ki] += Ek[c]
# E /= kappaNorm


cf = int(Ny/2)
fig,ax = plt.subplots(figsize=(6,6))
# plot Kolmogorov
KolSlope1= -5/3
KolSlope1 = Fraction(KolSlope1).limit_denominator()
KolLegend1 = KolSlope1
f1 = np.logspace(np.log10(kappa[1]), np.log10(kappa[-cf]), 100)
s1 = (f1**(KolSlope1)) * 1e-9
KolSlope2= -11/3 #-11/3
KolSlope2 = Fraction(KolSlope2).limit_denominator()
KolLegend2 = KolSlope2
f2 = np.logspace(np.log10(kappa[1]), np.log10(kappa[-cf]), 100)
s2 = (f2**(KolSlope2)) * 1e-13
ax.loglog(f1, s1, '--', f2, s2, '--', color=(0.5, 0.5, 0.5), label='_nolegend_')
ax.text(1.1e-1, 6e-8, KolLegend1, fontsize=12, color=(0.5, 0.5, 0.5))
ax.text(1.1e-1, 6e-10, KolLegend2, fontsize=12, color=(0.5, 0.5, 0.5))
# plot the energy spectrum
ax.loglog(kappa[:-cf], Ex[:-cf])
ax.loglog(kappa[:-cf], Ey[:-cf])
# ax.loglog(kappa[:-cf], Erho[:-cf])
ax.set_ylim(1e-14,1e-4)
ax.set_xlabel("Wavenumber ($\kappa$)",fontsize = fontSize_axis_label)
ax.set_ylabel("Energy Spectrum (E)",fontsize = fontSize_axis_label)
ax.tick_params(axis="x",labelsize=fontSize_axis_ticks,colors="black")
ax.tick_params(axis="y",labelsize=fontSize_axis_ticks,colors="black")
ax.legend(["E($\kappa$) on u", "E($\kappa$) on v","E($\kappa$) of rho"],
           loc='lower left',prop={'size': fontSize_legend})
ax.grid()
plt.savefig("Ek.jpeg",dpi=800,bbox_inches="tight")



## =================================================================
# sampling_rate = 1000
# file0 = open("postProcessing/probe0/0/U","r")
# file1 = open("postProcessing/probe1/0/U","r")
# file2 = open("postProcessing/probe2/0/U","r")

# p0 = file0.read()
# p1 = file1.read()
# p2 = file2.read()
# p0 = re.findall(r'[-]?\d+\.?\d*[e0-9]*[-+]?\d*',p0)
# p1 = re.findall(r'[-]?\d+\.?\d*[e0-9]*[-+]?\d*',p1)
# p2 = re.findall(r'[-]?\d+\.?\d*[e0-9]*[-+]?\d*',p2)
# p0 = p0[5:]
# p1 = p1[5:]
# p2 = p2[5:]
# p0 = np.array(p0,dtype=np.float32)
# p1 = np.array(p1,dtype=np.float32)
# p2 = np.array(p2,dtype=np.float32)

# length = int(len(p0)/4)
# window = np.hamming(length)
# time0 = np.zeros(length)

# vel0 = np.zeros((3,length))
# vel1 = np.zeros((3,length))
# vel2 = np.zeros((3,length))

# i = 0           
# cutoff0 = 0.25 # lower cutoff time
# cutoff1 = 1.5 # upper cutoff time 
# for d in range(length):
#     if (cutoff0 < p0[4*d+0] and p0[4*d+0] < cutoff1):
#         # probe 1
#         time0[i]  = p0[4*d+0]  # time
#         vel0[0,i] = p0[4*d+1]  # vx
#         vel0[1,i] = p0[4*d+2]  # vy
#         vel0[2,i] = p0[4*d+3]  # vz
#         # probe 1
#         # time1[i]  = p1[4*d+0]  # time
#         vel1[0,i] = p1[4*d+1]  # vx
#         vel1[1,i] = p1[4*d+2]  # vy
#         vel1[2,i] = p1[4*d+3]  # vz
#         # probe 2
#         # time2[i]  = p2[4*d+0]  # time
#         vel2[0,i] = p2[4*d+1]  # vx
#         vel2[1,i] = p2[4*d+2]  # vy
#         vel2[2,i] = p2[4*d+3]  # vz
#     i += 1
# # print(vel0[0])
# # plt.plot(vel0[0])
# # plt.savefig("test.jpeg")


# ## =================================================================
# u0 = vel0[0]*window
# v0 = vel0[1]*window
# u1 = vel1[0]*window
# v1 = vel1[1]*window
# u2 = vel2[0]*window
# v2 = vel2[1]*window

# uk0 = np.fft.fft(u0)
# vk0 = np.fft.fft(v0)
# uk1 = np.fft.fft(u1)
# vk1 = np.fft.fft(v1)
# uk2 = np.fft.fft(u2)
# vk2 = np.fft.fft(v2)

# power_spectrum0 = 0.5* (np.abs(uk0)**2 + np.abs(vk0)**2) / (length * np.sum(window**2))
# power_spectrum1 = 0.5* (np.abs(uk1)**2 + np.abs(vk1)**2) / (length * np.sum(window**2))
# power_spectrum2 = 0.5* (np.abs(uk2)**2 + np.abs(vk2)**2) / (length * np.sum(window**2))

# freq_vector = np.fft.fftfreq(length, 1/sampling_rate)
# positive_freqs = np.where(freq_vector > 0)

# f = np.logspace(np.log10(freq_vector[positive_freqs][1]), np.log10(freq_vector[positive_freqs][-1]), 100)
# s = (f**(-5/3)) * 2e-4 # amplitude scaling factor

## =================================================================
# uk0 = abs(np.fft.fft(vel0[0])/length)**2
# vk0 = abs(np.fft.fft(vel0[1])/length)**2
# EK = 0.5*(uk0 + vk0)

# plt.loglog(freq_vector[positive_freqs], power_spectrum0[positive_freqs])
# plt.loglog(freq_vector[positive_freqs], power_spectrum1[positive_freqs])
# plt.loglog(freq_vector[positive_freqs], power_spectrum2[positive_freqs])
# plt.loglog(f, s, 'r--')
# plt.savefig("test.jpeg")