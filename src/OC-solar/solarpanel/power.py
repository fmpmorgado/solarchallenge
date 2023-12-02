import numpy as np

def compute_IV_PindadoCubas(G, T, solar_panel, method = 'constant', V = 0):
    #https://oa.upm.es/67369/1/IEEETran_Aerospece_2021.pdf
    #https://oa.upm.es/43747/1/PindadoRE20016.pdf
    #https://core.ac.uk/download/pdf/211559749.pdf
    #https://research-information.bris.ac.uk/ws/portalfiles/portal/216857427/SSC19_WP1_16_Tom_Etchells.pdf

    if method != 'constant' and method != 'optimal':
        raise NameError("Options for method: constant, optimal")

    a = solar_panel.design_data['a']
    N_s = solar_panel.design_data['N_s']
    I_sc = solar_panel.design_data['I_sc']
    V_oc = solar_panel.design_data['V_oc']
    I_mp = solar_panel.design_data['I_mp']
    V_mp = solar_panel.design_data['V_mp']
    alpha_Isc = solar_panel.design_data['alpha_Isc']
    alpha_Voc = solar_panel.design_data['alpha_Voc']
    alpha_Imp = solar_panel.design_data['alpha_Imp']
    alpha_Vmp= solar_panel.design_data['alpha_Vmp']
    Gref = solar_panel.design_data['Gref']
    Tref = solar_panel.design_data['Tref']

    #Boltzmann constant in eV/K, 8.617332478e-05
    k = 8.617332478e-05
    VT = N_s * T * k

    #Account for Temperature:
    I_sc = G / Gref * (I_sc + alpha_Isc * (T - Tref))
    I_mp = G / Gref * (I_mp + alpha_Imp * (T - Tref))
    V_oc = V_oc + a * VT * np.log(G / Gref) + alpha_Voc * (T - Tref)
    V_mp = V_mp + a * VT * np.log(G / Gref) + alpha_Vmp * (T - Tref)

    if method == 'optimal': 
        V = V_mp

    phi = I_sc / I_mp * (I_sc / (I_sc-I_mp)) * ((V_oc - V_mp) / V_oc)
    
    if V <= V_mp:
        I = I_sc * (1.0 - (1.0 - I_mp / I_sc) * (V / V_mp) ** (I_mp/(I_sc - I_mp)))
    else:
        I = I_mp * V_mp / V * (1 - ((V - V_mp) / (V_oc - V_mp)) ** phi)

    return 0

"""
def compute_diode_parameters(eff_irradiance, temp_cell, alpha_sc, a_ref,
                             I_L_ref, I_o_ref, R_sh_ref, R_s, EgRef=1.121,
                             dEgdT=-0.0002677, irrad_ref = 1000.0, temp_ref = 25.0):
    #https://pvpmc.sandia.gov/modeling-guide/2-dc-module-iv/single-diode-equivalent-circuit-models/de-soto-five-parameter-module-model/
    #https://github.com/pvlib/pvlib-python/blob/main/pvlib/pvsystem.py#L1653

    #No spectral effect
    #De Soto model (De Soto et al. 2006)

    #TODO
    #Boltzmann constant in eV/K, 8.617332478e-05
    k = 8.617332478e-05
    
    temp_cell_K = temp_cell + 273.15
    temp_ref_K = temp_ref + 273.15

    E_g = EgRef * (1 + dEgdT*(temp_cell_K - temp_ref_K))

    IL = eff_irradiance / irrad_ref * (I_L_ref + alpha_sc * (temp_cell_K - temp_ref_K))
    I0 = (I_o_ref * (((temp_cell_K)/ temp_ref_K) ** 3) * (np.exp(EgRef / (k*(temp_ref_K)) - (E_g / (k*(temp_cell_K))))))

    Rsh = R_sh_ref * (irrad_ref / eff_irradiance)
    Rs = R_s

    nNsVth = a_ref * (temp_cell_K / temp_ref_K)

    return (IL, I0, Rs, Rsh, nNsVth)
"""
