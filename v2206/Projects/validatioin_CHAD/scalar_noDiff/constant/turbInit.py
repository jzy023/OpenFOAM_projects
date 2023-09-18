import math

targetRe = 1e3
Vref = 10
L = 1
nu = 1e-3
Re = Vref*L/nu

Iturb = 0.05
# Iturb = 0.16*math.pow(Re,-1/8)
Cmu = 0.09
Lref = 0.01*L

k = 1.5*(Iturb*Vref)*(Iturb*Vref)
epsilon = math.pow(Cmu, 1)*math.pow(k, 1.5)/Lref
omega = epsilon/(k * math.pow(Cmu, 7/4))
# omega = math.sqrt(k)/(math.pow(Cmu,0.25)*Lref)
# omega = math.pow(Cmu,0.75)*math.sqrt(k)/(Lref)
nut = k/omega

print("Re:\t",Re)
print("I:\t",Iturb)
print("k:\t", k)
print("epsilon:", epsilon)
print("omega:\t",omega)
print("nut:\t",nut)
