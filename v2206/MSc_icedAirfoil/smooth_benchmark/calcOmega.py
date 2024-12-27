
nu = 1.5e-5
y0 = 0.0137
beta1 = 0.072

I = 0.01
Uref = 15.0
Lref = 1.0
Cmu = 0.09

k0 = 1.5*((I*abs(Uref))**2)
omega0 =  (k0**0.5) / ((Cmu**0.25)*Lref)
omega_wall = (10*6*nu) / (beta1*(y0**2))

print("k_ref: ", k0)
print("omega_ref: ", omega0)
print("omega_wall: ", omega_wall)