import numpy as np
from solarchallenge.constants import IR_EARTH, ALBEDO_EARTH, STEFAN_BOLTZMANN


def heat_earth_IR(emittance_back: float, ff: float, area: float) -> float:
    """
    Compute Earth Infrared (IR) heat onto the solar panel

    The solar panel is assumed to be pointing at zenith, so the back of the solar panel
    is directly exposed to Earth IR at all times.

    Parameters
    ----------
    emittance_back : float
        Emittance value of the solar panel back part.
    ff : float
        Form Factor (from Earth)
    area : float
        Area of the solar panel [m**2]
    
    Returns
    -------
    q_earth_ir : float
        Heat impinging the back of the solar panel due to Earth IR [W]
    """

    q_earth_ir = IR_EARTH * ff * emittance_back * area

    return q_earth_ir


def heat_earth_albedo(absorptance_back: float, area: float, ff: float,
                      solar_irr: float, eta: float) -> float:
    """
    Compute Earth albedo heat onto the solar panel

    The solar panel is assumed to be pointing at zenith, so the back of the solar panel
    is directly exposed to the solar reflection from Earth (albedo).

    Parameters
    ----------
    absorptance_back : float
        Absorptance value of the solar panel back part.
    area : float
        Area of the solar panel [m**2]
    ff : float
        Form Factor (from Earth)
    solar_irr : float
        Solar irradiance [W/m**2]
    eta : float
        Zenith angle [deg]
        
    Returns
    -------
    q_earth_alb : float
        Heat impinging the back of the solar panel due to albedo [W]
    """

    q_earth_alb = solar_irr * ALBEDO_EARTH* ff * np.cos(np.radians(eta)) * area * absorptance_back

    #q_earth_alb is positive only
    q_earth_alb[q_earth_alb < 0] = 0

    return q_earth_alb


def heat_solar_radiation(absorptance_back: float, absorptance_front: float, 
                         area: float, solar_irr: float, theta: float, shadow: bool) -> float:
    """
    Compute the solar radiation onto the solar panel.

    Parameters
    ----------
    absorptance_back : float
        Absorptance value of the solar panel back part.
    absorptance_front : float
        Absorptance value of the solar panel front part (Photovoltaic panels).
    area : float
        Area of the solar panel [m**2]
    solar_irr : float
        Solar irradiance [W/m**2]
    theta : float
        Incidence angle between solar rays and front part of solar panel [deg]
    shadow : bool
        Boolean value to indicate if solar panel is shadowed by Earth (Eclipse)    
    
    Returns
    -------
    q_solar_back : float
        Heat impinging the back of the solar panel due to solar irradiation [W]
    q_solar_front : float
        Heat impinging the front of the solar panel due to solar irradiation [W]
    """
 
    q_solar_front = absorptance_front * solar_irr * area * np.cos(np.radians(theta)) * ~shadow
    q_solar_back = absorptance_back * solar_irr * area * np.cos(np.radians(theta + 180.0)) * ~shadow

    q_solar_back[q_solar_back < 0] = 0
    q_solar_front[q_solar_front < 0] = 0

    return q_solar_back, q_solar_front


def heat_panel_radiation(emittance_back: float, emittance_front: float,
                         area: float, temperature: float) -> float:
    """
    Compute the radiation from the solar panel to space.

    Parameters
    ----------
    emittance_back : float
        Emittance value of the solar panel back part.
    emittance_front : float
        Emittance value of the solar panel front part (Photovoltaic panels).
    area : float
        Area of the solar panel [m**2]
    temperature : float
        temperature of the solar panel [ºC]

    Returns
    -------
    q_rad : float
        Radiation from the solar panel to space [W]    
    """

    q_rad = STEFAN_BOLTZMANN * (emittance_back + emittance_front) * area * (temperature + 273.15) ** 4
    return q_rad
