function gamma()
    # a simple barrier estimation in TST
    
    # constants
    GAS_CONSTANT = 8.3144621            # J / K / mol
    PLANCK_CONSTANT = 6.62606957e-34    # J * s
    BOLTZMANN_CONSTANT = 1.3806488e-23  # J / K

    # parameters (input)
    time = 10 * 24 * 3600   # s (time scale)
    temp = 273.15 + 150     # K (reaction temperature)
  
    gamma = - GAS_CONSTANT * temp * log(PLANCK_CONSTANT/(time * BOLTZMANN_CONSTANT * temp))

    println(gamma)           # J/(mol*K)
    println(gamma / 1000)    # kJ/(mol*K)

end
