import numpy as np

Cp = 4182
kappa = 0.598
g = 9.81
H = 0.12
nu = 2.0e-5
rho = 1000
mu = nu*rho # 0.05
Re = rho*((H/(2*np.pi))*np.sqrt(g*(H/(2*np.pi))))/mu
# Re = rho*(H/2)*np.sqrt(g*(H/2))/mu
Pr = Cp*mu/kappa
print("mu is: ", mu)
print("Reynolds number is: ", Re)
print("Prandtl number is: ", Pr)
print("alpha is: ", nu/Pr)