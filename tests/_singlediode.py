"""
testing single-diode methods with PVlib library
"""

import pytest
import pvlib
import numpy as np
from src.solarpanel.power import 

RADIANCE = 888
TCELL = 55

@pytest.fixture
def cec_module_ZT320P():

    """
    Define Zytech Solar ZT320P module parameters for testing.
    """

    parameters = {
        'Name': 'Zytech Solar ZT320P',
        'Manufacturer': 'Zytech Solar',
        'Technology': 'Multi-c-Si',
        'Bifacial': 0,
        'STC': 320.42,
        'PTC': 289.8,
        'A_c': 1.931,
        'Length': 1.95,
        'Width': 0.99,
        'N_s': 72,
        'I_sc_ref': 9.12,
        'V_oc_ref': 46.6,
        'I_mp_ref': 8.66,
        'V_mp_ref': 37,
        'alpha_sc': 0.00440496, 
        'beta_oc': -0.149073,
        'T_NOCT': 46.4,
        'a_ref': 1.87379,
        'I_L_ref': 9.21845,
        'I_o_ref': 1.44659e-10,
        'R_s': 0.475581,
        'R_sh_ref': 604.222,
        'Adjust': 5.83834,
        'gamma_r': -0.4308,
        'BIPV': 'N',
        'Version': 'SAM 2021.12.02',
        'Date': '11/16/2022'
    }

    return parameters
@pytest.skip
def test_diode_desoto(cec_module_ZT320P):

    EgRef = 1.121
    dEgdT = -0.0002677

    pvlib_sol = pvlib.pvsystem.calcparams_desoto(
                effective_irradiance=RADIANCE, temp_cell=TCELL,
                alpha_sc=cec_module_ZT320P['alpha_sc'], a_ref=cec_module_ZT320P['a_ref'],
                I_L_ref=cec_module_ZT320P['I_L_ref'], I_o_ref=cec_module_ZT320P['I_o_ref'],
                R_sh_ref=cec_module_ZT320P['R_sh_ref'], R_s=cec_module_ZT320P['R_s'],
                EgRef=EgRef, dEgdT=dEgdT)
    
    sol = compute_diode_parameters(eff_irradiance=RADIANCE, temp_cell=TCELL,
                                   alpha_sc=cec_module_ZT320P['alpha_sc'], a_ref=cec_module_ZT320P['a_ref'],
                                   I_L_ref=cec_module_ZT320P['I_L_ref'], I_o_ref=cec_module_ZT320P['I_o_ref'],
                                   R_sh_ref=cec_module_ZT320P['R_sh_ref'], R_s=cec_module_ZT320P['R_s'],
                                   EgRef=EgRef, dEgdT=dEgdT)
    
    assert np.isclose(sol, pvlib_sol).all()