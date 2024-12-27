
# constant parameters
beta1 = 0.072
Cmu = 0.09

# case setups
I = 0.005
y0 = 0.002
nu = 1.5e-5
Uref = 15.0
Lref = 1.0
Lturb = 0.07*Lref

# intermediates
Re = Uref*Lref/nu
# Iturb = 0.16*(Re**(-1/8))
Iturb = I

# results
k0 = 1.5*((Iturb*abs(Uref))**2)
epsilon0 = (Cmu**0.75)*(k0**1.5)/Lturb
omega0 = (k0**0.5) / ((Cmu**0.25)*Lturb)
omega_wall = (10*6*nu) / (beta1*(y0**2))

# outputs
print("turbulent intensity (I): ", Iturb)
print("k_ref: ", k0)
print("epsilon_ref: ", epsilon0)
print("omega_ref: ", omega0)
print("omega_wall: ", omega_wall)