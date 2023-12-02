import numpy as np

# PV module data from a typical datasheet (e.g. Kyocera Solar KD225GX LPB)
module_data = {'celltype': 'multiSi', # technology
               'STC': 224.99, # STC power
               'PTC': 203.3, # PTC power
               'v_mp': 29.8, # Maximum power voltage
               'i_mp': 7.55, # Maximum power current
               'v_oc': 36.9, # Open-circuit voltage
               'i_sc': 8.18, # Short-circuit current
               'alpha_sc': 0.001636, # Temperature Coeff. Short Circuit Current [A/C]
               'beta_voc': -0.12177, # Temperature Coeff. Open Circuit Voltage [V/C]
               'gamma_pmp': -0.43, # Temperature coefficient of power at maximum point [%/C]
               'cells_in_series': 60, # Number of cells in series
               'temp_ref': 25}  # Reference temperature conditions

# Import the pvlib library
import pvlib

# 1st step: Estimating the parameters for the CEC single diode model
""" WARNING - This function relies on NREL's SAM tool. So PySAM, its Python API, needs to be installed 
in the same computer. Otherwise, you can expect the following error: 'ImportError if NREL-PySAM is not installed.'
"""
cec_fit_params = pvlib.ivtools.sdm.fit_cec_sam(module_data['celltype'], module_data['v_mp'], module_data['i_mp'],
                                  module_data['v_oc'], module_data['i_sc'], module_data['alpha_sc'],
                                  module_data['beta_voc'], module_data['gamma_pmp'], 
                                  module_data['cells_in_series'], module_data['temp_ref'])

#I_L_ref, I_o_ref, R_s, R_sh_ref, a_ref and Adjust.

# Let's have a look to the output
print(cec_fit_params, 60*8.61733326E-5*1*(25+273.15))

# Effective irradiance values (W/m2)
irrad = np.array([200,400,600,800,1000])
# Average cell temperature (degrees Celsius)
temp_cell = np.array([40, 40, 40, 40, 40])

# 2nd step: Apply model to estimate the 5 parameters of the single diode equation using the CEC model
diode_params = pvlib.pvsystem.calcparams_cec(irrad, temp_cell, module_data['alpha_sc'], cec_fit_params[4], 
                                            cec_fit_params[0], cec_fit_params[1], cec_fit_params[3], 
                                            cec_fit_params[2], cec_fit_params[5])

# The result of the function returns a Tuple of 5 parameters to be used in the single diode equation
print('Number of elements returned: ', len(diode_params))

"""
class Earth:
    orbital_elements = {"a": ,
                        "e": ,
                        "i": ,
                        ""}
"""
class Satellite():
    def __init__(self):
        self.panel = {}
        self.position = {}
        self.velocity = {}
        self.attitude = {}

class SolarPanel():
    def __init__(self):
        self.data = module_data
        self.teste = 1

solar = SolarPanel()
print(solar.teste)


def incident_power(date):

    #Stephen Boltzmann law
    E = 4*np.pi*2



    pass

def compute_angle_incidence():
    pass

from astropy.time import Time
from astropy.coordinates import SkyCoord, GCRS

from time import time
a = time()
c = SkyCoord(0, 0, frame ='icrs', unit='deg', obstime=Time([2000, 2010], format='jyear'))#.transform_to('gcrs').distance.astronomical_unit

#sc = SkyCoord(ra=[15], dec=[-70], unit='deg', obstime=Time([2000, 2010], format='jyear'))
print(time()-a)

print(c)

from src.solarpanel.module import SolarModule
m = SolarModule.from_database("Zytech Solar ZT320P")

parameters = {'celltype': 'multiSi', # technology
               'STC': 224.99, # STC power
               'PTC': 203.3, # PTC power
               'v_mp': 29.8, # Maximum power voltage
               'i_mp': 7.55, # Maximum power current
               'v_oc': 36.9, # Open-circuit voltage
               'i_sc': 8.18, # Short-circuit current
               'alpha_sc': 0.001636, # Temperature Coeff. Short Circuit Current [A/C]
               'beta_voc': -0.12177, # Temperature Coeff. Open Circuit Voltage [V/C]
               'gamma_pmp': -0.43, # Temperature coefficient of power at maximum point [%/C]
               'cells_in_series': 60, # Number of cells in series
               'temp_ref': 25}  # Reference temperature conditions


n = SolarModule(parameters)
print(m.parameters)

from src.solarpanel.power import compute_diode_parameters, compute_IV_PindadoCubas


diode_params = pvlib.pvsystem.calcparams_desoto(irrad, temp_cell, module_data['alpha_sc'], cec_fit_params[4], 
                                            cec_fit_params[0], cec_fit_params[1], cec_fit_params[3], 
                                            cec_fit_params[2])

p = compute_diode_parameters(irrad, temp_cell, module_data['alpha_sc'], cec_fit_params[4], cec_fit_params[0], cec_fit_params[1], cec_fit_params[3], 
                                            cec_fit_params[2])


from src.solarpanel.solarpanel import SolarPanel, SolarPanel2


from datetime import date

j = SolarPanel(0,0,date(2022,1,1),0,0)
j.power(0,0,date(2023,1,1))



solar_panel = {
    'a': 1,
    'N_s': 3,
    'I_sc': 0.506,
    'V_oc': 2.667,
    'I_mp': 0.487,
    'V_mp': 2.371,
    'alpha_Isc': 0.32/1000.0,
    'alpha_Voc': -6.0/1000.0,
    'alpha_Vmp': -6.1/1000.0,
    'alpha_Imp': 0.28/1000.0,
    'Tref': 28,
    'Gref': 1367,
}


m = SolarPanel2(solar_panel, 1, 1)

TCELL = 20
RADIANCE = 1300

sol = compute_IV_PindadoCubas(G=RADIANCE, T=TCELL, solar_panel=m, method='constant', V = 2.6)