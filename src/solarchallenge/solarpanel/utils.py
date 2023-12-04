import numpy as np
import numpy.typing as npt
from solarchallenge.constants import A_SURF_SUN, R_MEAN_EARTH, STEFAN_BOLTZMANN, T_SURF_SUN


def compute_angle_solar_incidence(r_planet: npt.ArrayLike, r_panel: npt.ArrayLike) -> npt.ArrayLike:
    """
    Computes the angle between the solar rays and the orientation of the solar panel

    This function assumes the solar panel is pointing to zenith

    Parameters
    ----------
    r_planet : np.array[floats]
        List of position vectors of the planet in relation to the sun (EME2000 frame) [km]
    r_panel : np.array[floats]
        List of position vectors of the solar panel in relation to the planet (J2000 frame) [km]
    
    Returns
    -------
    theta : np.array[float]
        List of incidence angles between solar panel orientation and solar rays [deg]
    """
    normal_panel = r_panel / np.linalg.norm(r_panel, axis = -1).reshape(-1,1)

    r_sun_panel = r_planet + r_panel
    r_sun_panel /= np.linalg.norm(r_sun_panel, axis = -1).reshape(-1,1)

    theta = np.arccos(np.clip(np.sum(-normal_panel * r_sun_panel, axis = -1), -1.0, 1.0))
    theta = np.degrees(theta)

    return theta


def compute_solar_zenith_angle(r_planet: npt.ArrayLike, r_panel: npt.ArrayLike) -> npt.ArrayLike:
    """
    Computes the Solar zenith angle with respect to the panel

    This function assumes the solar panel is pointing to zenith 

    Parameters
    ----------
    r_planet : np.array[floats]
        List of position vectors of the planet in relation to the sun (EME2000 frame) [km]
    r_panel : np.array[floats]
        List of position vectors of the solar panel in relation to the planet (J2000 frame) [km]
    
    Returns
    -------
    eta : np.array[float]
        List of solar zenith angles [deg]
    """

    #Sun position with respect to the planet
    r_sun = -r_planet

    normal_sun = r_sun / np.linalg.norm(r_sun, axis = -1).reshape(-1,1)
    normal_panel = r_panel / np.linalg.norm(r_panel, axis = -1).reshape(-1,1)

    eta = np.arccos(np.clip(np.sum(normal_panel * normal_sun, axis = -1), -1.0, 1.0))
    eta = np.degrees(eta)
    
    return eta


def compute_solar_radiance(r_planet: npt.ArrayLike, r_panel: npt.ArrayLike) -> npt.ArrayLike:
    """
    Computes the solar irradiance as a function of the distance between sun and panel

    Parameters
    ----------
    r_planet : np.array[floats]
        List of position vectors of the planet in relation to the sun (EME2000 frame) [km]
    r_panel : np.array[floats]
        List of position vectors of the solar panel in relation to the planet (J2000 frame) [km]
    
    Returns
    -------
    solar_irr : np.array[float]
        List of solar irradiance [W/m**2]
    """

    r_sun_panel = r_planet + r_panel
    distance_sun_panel = np.linalg.norm(r_sun_panel, axis = -1).reshape(-1,1)

    solar_irr = (STEFAN_BOLTZMANN * A_SURF_SUN * 1E6 * T_SURF_SUN ** 4) / (4 * np.pi * (distance_sun_panel * 1E3) ** 2)
    
    return solar_irr.reshape(-1)


def compute_form_factor(r_panel: npt.ArrayLike):
    """
    Computes the local form factor with respect to the planet

    Parameters
    ----------
    r_panel : np.array[floats]
        List of position vectors of the solar panel in relation to the planet (J2000 frame) [km]
    
    Returns
    -------
    ff : np.array[float]
        List of form factors
    """
    
    ff = (R_MEAN_EARTH / np.linalg.norm(r_panel, axis = -1).reshape(-1)) ** 2

    return ff


def check_is_shadowed(r_planet: npt.ArrayLike, r_panel: npt.ArrayLike, eta: npt.ArrayLike):

    #Two conditions need to be verified
    
    #1 - Minimum distance of solar to Eun-Earth has to be less that Earth radius
    #2 - Zenith angle has to be > 90 degrees

    r_sun_planet = r_planet
    r_sun_panel = r_planet + r_panel
    d = np.linalg.norm(np.cross(r_sun_planet, -r_sun_panel ,axis=-1), axis = -1) / np.linalg.norm(r_sun_planet, axis = -1)

    shadow = (eta > 90) * (d < R_MEAN_EARTH)
    
    return shadow