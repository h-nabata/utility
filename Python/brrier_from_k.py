import math

# constants
GAS_CONSTANT = 8.3144621  # J / K / mol
PLANCK_CONSTANT = 6.62606957e-34  # J * s
BOLTZMANN_CONSTANT = 1.3806488e-23  # J / K
SPEED_OF_LIGHT = 2.99792458e8  # m / s

# input parameters
temperature = 300.15  # K
rate_constant = 1.00e-03  # kJ / mol
frequency = -300.0  # cm-1 : imaginary frequency

# Wigner's tunneling correction (E. Wigner, Phys. Rev., 1932, 40, 749–759.)
# 100.0 is the scaling factor to convert cm-1 into m-1
wigner_factor = (PLANCK_CONSTANT * abs(frequency) * 100.0 * SPEED_OF_LIGHT / (BOLTZMANN_CONSTANT * temperature))
# the transmission coefficient
kappa = 1 + 1/24 * wigner_factor ** 2

# Eyring's equation (H. EYRING, J. Chem. Phys., 1935, 3, 107–115.)
Coef = (GAS_CONSTANT * temperature) / 1000.0
Anti = rate_constant * PLANCK_CONSTANT / (BOLTZMANN_CONSTANT * temperature)

dE = -1 * Coef * math.log(Anti)
dE_corr = -1 * Coef * math.log(Anti/kappa)

print("反応温度 =", temperature, "(K) / 反応速度定数 =", rate_constant, "(1/s) / 虚振動数 =", frequency, "(cm-1)")
print("Wigner補正項 κ = {:e}".format(kappa))
print("活性化障壁       : {:12}".format(dE))
print("    補正後       : {:12}".format(dE_corr))
