import numpy as np
from solarchallenge.constants import BOLTZMANN
#TODO Missing SolarPanel import (solve circular import)


def compute_IV_PindadoCubas(solar_irr: float, temperature: float, solar_panel,
                            voltage: None|float=None) -> tuple[float, float]:
    """
    Method to compute the IV of a solar panel.

    Method is an alternative to the single diode model. It uses the most common parameters
    found in solar panel datasheets.

    More information on the method can be found at:
    - https://oa.upm.es/67369/1/IEEETran_Aerospece_2021.pdf
    - https://oa.upm.es/43747/1/PindadoRE20016.pdf

    Parameters
    ----------
    solar_irr : float
        Effective Solar irradiance [W/m**2] (inclination already accounted for)
    temperature : float
        Temperature of the solar panel [ºC or K]
    solar_panel : SolarPanel
        Object of class SolarPanel
    voltage : float | None
        Operating voltage. If None, the maximum power voltage is used


    Returns
    -------
    current : float
        Current of the solar panel [A]
    voltage : float
        Voltage of the solar panel [V]
    """

    # Solar Panel characteristics
    a = solar_panel.a
    n_s = solar_panel.n_s
    i_sc = solar_panel.i_sc
    v_oc = solar_panel.v_oc
    i_mp = solar_panel.i_mp
    v_mp = solar_panel.v_mp
    alpha_i_sc = solar_panel.alpha_i_sc
    alpha_v_oc = solar_panel.alpha_v_oc
    alpha_i_mp = solar_panel.alpha_i_mp
    alpha_v_mp= solar_panel.alpha_v_mp
    solar_irr_ref = solar_panel.solar_irr_ref
    tref = solar_panel.tref

    # Thermal voltage
    v_t = n_s * temperature * BOLTZMANN

    #Check for negative solar irradiances (due to incidence angle theta > 90 deg or theta < -90 deg)
    try:
        if solar_irr <= 0.001: solar_irr = 0.001    
    except:
        solar_irr[solar_irr <= 0.001] = 0.001

    #Account for temperature and solar irradiance:
    i_sc = solar_irr / solar_irr_ref * (i_sc + alpha_i_sc * (temperature - tref))
    i_mp = solar_irr / solar_irr_ref * (i_mp + alpha_i_mp * (temperature - tref))
    v_oc = v_oc + a * v_t * np.log(solar_irr / solar_irr_ref) + alpha_v_oc * (temperature - tref)
    v_mp = v_mp + a * v_t * np.log(solar_irr / solar_irr_ref) + alpha_v_mp * (temperature - tref)

    if voltage is None: 
        voltage = v_mp

    phi = i_sc / i_mp * (i_sc / (i_sc-i_mp)) * ((v_oc - v_mp) / v_oc)
    
    #This check is to see if any np.array was used as an argument of a function,
    #as the function allows to vectorize the computation of power for arrays of
    #angle incidence and solar irradiance at constant temperature.
    try:
        if voltage <= v_mp:
            current = i_sc * (1.0 - (1.0 - i_mp / i_sc) * (voltage / v_mp) ** (i_mp / (i_sc - i_mp)))
        else:
            current = i_mp * v_mp / voltage * (1 - ((voltage - v_mp) / (v_oc - v_mp)) ** phi)
    except:
        current = np.zeros(v_mp.shape)

        if type(voltage) is int or type(voltage) is float:
            voltage = np.ones(v_mp.shape) * voltage

        index = (voltage - v_mp <= 0)
        current[index] = i_sc[index] * (1.0 - (1.0 - i_mp[index] / i_sc[index]) * (voltage[index] / v_mp[index]) ** (i_mp[index] /(i_sc[index] - i_mp[index])))

        index = (voltage - v_mp > 0)
        current[index] = i_mp[index] * v_mp[index] / voltage[index] * (1 - ((voltage[index] - v_mp[index]) / (v_oc[index] - v_mp[index])) ** phi[index])

    return current, voltage