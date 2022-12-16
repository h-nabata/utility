function main()
    GAS_CONSTANT = 8.3144621            # J / K / mol
    PLANCK_CONSTANT = 6.62606957e-34    # J * s
    BOLTZMANN_CONSTANT = 1.3806488e-23  # J / K
    SPEED_OF_LIGHT = 2.99792458e10      # cm / s
    AVOGADRO_CONSTANT = 6.0221415e23    # 1 / mol
    AMU_to_KG = 1.66053886e-27          # UNIT CONVERSION
    temperature = 298.15                # K

    ### References
    # [1] M. Mammen, E. Shakhnovich, J. Deutch, and G. Whitesides, J. Org. Chem. 1998, 63, 12, 3821–3830. https://pubs.acs.org/doi/10.1021/jo970944f
    # [2] "GoodVibes" by patonlab ( https://github.com/patonlab/GoodVibes )

    solv_molecular_vol = 125.412  # Ang**3  (DMF ... Cavity volume = 125.412 ~ 126.089 Ang**3) (const.))
    solv_molarity = 10.763        # mol / L (const.)
    molecular_mass = 73.09        # g / mol   73.09 ... DMF;  119.38 ... CHCl3
    conc = 3.1                    # mol / L (the sum of concentration of solutes (const.))

    # "For cubic molecules in a 3D cubic array, Cfree = 8; for spherical molecules, Cfree = 6.3" from ref. [1]
    # (*) 1 Ang**3 = 1e-27 L
    v_free = 8 * ((1e27 / (solv_molarity * AVOGADRO_CONSTANT)) ^ (1/3) - solv_molecular_vol ^ (1/3)) ^ 3    # Ang**3
    freespace = (v_free * 1e-27) * solv_molarity * AVOGADRO_CONSTANT

    # eq. (6a) from ref. [1]
    lmda = ((2.0 * π * molecular_mass * AMU_to_KG * BOLTZMANN_CONSTANT * temperature) ^ 0.5) / PLANCK_CONSTANT
    ndens = (conc * 1000) * AVOGADRO_CONSTANT / freespace
    entropy = GAS_CONSTANT * (2.5 + log(lmda ^ 3 / ndens))  # J / (mol*K)

    println(entropy / (1000 * 2625.50227))   # hartree
    println(298.15*entropy/1000)   # J/(mol*K)
end
