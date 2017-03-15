
from scipy import *


# fundamental constants
epsilon0 = 8.85419e-12  # C**2/Jm = [surface charge density / Efield]
me = 9.10939e-31 # kg
e = 1.60218e-19 # C
hbar = 1.05457e-34 # Js
c = 2.99792e8 # m/s
eV = 1.60218e-19 # J


# derived constants
alpha = e**2 / (4*pi*epsilon0*hbar*c) # Fine structure const = 1/137 (unitless)
a0 = hbar / (alpha * me * c)  # Bohr radius = 0.529e-10m
E1 = - alpha**2 * me * c**2 / 2 # Ground state energy of Hydrogen atom = -13.6 eV
Rydberg = E1  # -13.6 eV
# En = E1/n**2



