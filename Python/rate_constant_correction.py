import math

# constants
GAS_CONSTANT = 8.3144621  # J / K / mol
PLANCK_CONSTANT = 6.62606957e-34  # J * s
BOLTZMANN_CONSTANT = 1.3806488e-23  # J / K
SPEED_OF_LIGHT = 2.99792458e8  # m / s

# input parameters
temperature = 298.15 # K
activation_barrier = 100.0 # kJ / mol
frequency = -500.0 # cm-1 : imaginary frequency

# Eyring's equation (H. EYRING, J. Chem. Phys., 1935, 3, 107–115.)
E = math.exp(- activation_barrier * 1000.0 /(GAS_CONSTANT * temperature))
Z = (BOLTZMANN_CONSTANT * temperature / PLANCK_CONSTANT)
k = Z * E

# Wigner's tunneling correction (E. Wigner, Phys. Rev., 1932, 40, 749–759.)
# 100.0 is the scaling factor to convert cm-1 into m-1
wigner_factor = (PLANCK_CONSTANT * (abs(frequency) * 100.0 * SPEED_OF_LIGHT) / (BOLTZMANN_CONSTANT * temperature))
# the transmission coefficient
kappa = 1 + 1/24 * wigner_factor ** 2
Z_corr = kappa * Z
k_corr = Z_corr * E

print("rate constant from Eyring's equation :      {:e}".format(k))
print("corrected rate constant by Wigner's model : {:e}".format(k_corr), "(κ = {:e})".format(kappa))

