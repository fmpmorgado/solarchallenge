import numpy as np
import numpy.typing as npt

def constant_accel_wrapper(accel: float):
    """
    Function wrapper to use in scipy.solve_ivp.

    Applies constant acceleration during orbit propagation.

    Parameters
    ----------
    accel : float
        Acceleration in velocity direction [km/s**2]
    
    Returns
    -------
    constant_accel : function
        Function that computes the acceleration in x, y, z direction

    Additional Info (Nested function)
    ---------------------------------
    t : float
        Current time of the solve_ivp step
    state : np.array[float]
        State vector [x, y, z, vx, vy, vz] (km, km/s).
    mu : float
        Standard Gravitational parameter [km**3/s**2]
    """

    def constant_accel(t: float, state: npt.ArrayLike , mu: float) -> npt.ArrayLike:

        v_vec = state[3:]
        v_direction = v_vec/np.linalg.norm(v_vec)

        a_x = accel * v_direction[0]
        a_y = accel * v_direction[1]
        a_z = accel * v_direction[2]
        
        return np.array([0, 0, 0, a_x, a_y, a_z])
    
    return constant_accel
    
def J2_perturbation_wrapper(J2: float, R: float):
    """
    Function wrapper to use in scipy.solve_ivp.

    Applies J2 perturbation during orbit propagation.

    Parameters
    ----------
    J2 : float
        J2 factor of the orbit attractor
    R : float
        Attractor radius [km]
    
    Returns
    -------
    J2_perturbation : function
        Function that computes the J2 perturbation x, y, z direction

    Additional Info (Nested function)
    ---------------------------------
    t : float
        Current time of the solve_ivp step
    state : np.array[float]
        State vector [x, y, z, vx, vy, vz] (km, km/s).
    mu : float
        Standard Gravitational parameter [km**3/s**2]
    """

    def J2_perturbation(t: float, state: npt.ArrayLike , mu: float) -> npt.ArrayLike:

        r = state[:3]

        #Distance from center of the attractor
        r_norm = np.linalg.norm(r)

        factor = (3.0 / 2.0) * mu * J2 * (R**2) / (r_norm**5)

        a_x = r[0] * factor * (5.0 * r[2] ** 2 / r_norm**2 - 1)
        a_y = r[1] * factor * (5.0 * r[2] ** 2 / r_norm**2 - 1)
        a_z = r[2] * factor * (5.0 * r[2] ** 2 / r_norm**2 - 3)

        return np.array([0, 0, 0, a_x, a_y, a_z])
    
    return J2_perturbation

def J3_perturbation_wrapper(J3: float, R: float):
    """
    Function wrapper to use in scipy.solve_ivp.

    Applies J3 perturbation during orbit propagation.

    Parameters
    ----------
    J3 : float
        J3 factor of the orbit attractor
    R : float
        Attractor radius [km]
    
    Returns
    -------
    J3_perturbation : function
        Function that computes the J3 perturbation x, y, z direction

    Additional Info (Nested function)
    ---------------------------------
    t : float
        Current time of the solve_ivp step
    state : np.array[float]
        State vector [x, y, z, vx, vy, vz] (km, km/s).
    mu : float
        Standard Gravitational parameter [km**3/s**2]
    """

    def J3_perturbation(t: float, state: npt.ArrayLike , mu: float) -> npt.ArrayLike:
        r = state[:3]

        #Distance from center of the attractor
        r_norm = np.linalg.norm(r)

        factor = (1.0 / 2.0) * mu * J3 * (R**3) / (r_norm**5)

        cos_phi = r[2] / r_norm

        a_x = 5.0 * r[0] / r_norm * (7.0 * cos_phi**3 - 3.0 * cos_phi)
        a_y = 5.0 * r[1] / r_norm * (7.0 * cos_phi**3 - 3.0 * cos_phi)
        a_z = 3.0 * (35.0 / 3.0 * cos_phi**4 - 10.0 * cos_phi**2 + 1)
        
        return np.array([0, 0, 0, a_x, a_y, a_z]) * factor

    return J3_perturbation

"""
def atmospheric_drag_wrapper(R, CD, A_over_m, H0, rho0):
    def atmospheric_drag_exponential(t, state, mu):

        r_vec = state[:3]
        v_vec = state[3:]
        r = np.linarg.norm(r_vec)
        v = np.linarg.norm(v_vec)

        B = CD * A_over_m
        rho = rho0 * np.exp(-(r - R) / H0)

        return -(1.0 / 2.0) * rho * B * v * v_vec
"""