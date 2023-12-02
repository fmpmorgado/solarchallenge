import numpy as np
from datetime import datetime
from .orbit import coe2rv, rv2coe, Orbit

#http://murison.alpheratz.net/dynamics/twobody/KeplerIterations_summary.pdf
def compute_kepler(orbit: Orbit, epoch: datetime, time_perihelion: datetime,
                   tolerance: float=1e-6, max_iter: int=100, 
                   e_anom_i: float | None=None) -> tuple[float, float, float]:

    """Computes the position and velocity for a given epoch
       using Kepler equation.

    Parameters
    ----------
    orbit : Orbit
        Object of class Orbit
    epoch : datetime 
        Datetime of target epoch
    time_perihelion : datetime
        Datetime of perihelion passage
    tolerance : float 
        Tolerance for the Newton-Raphson solve
    max_iter : int
        Maximum number of iterations for the Newton-Raphson solve
    e_anom_i : float
        Initial eccentric anomaly value for the Newton-Raphson solve
    
    Returns
    -------
    r_target : list[float]
        Position vector for computed true anomaly [km]
    v_target : list[float]
        Velocity for computed true anomaly [km/s]
    nu : float
        Computed true anomaly [deg]
    """

    a, ecc, inc, raan, argp, __ = rv2coe(r=orbit.r, v=orbit.v, mu=orbit.attractor.mu)
    mu = orbit.attractor.mu

    #Computation of the orbit period
    period = np.sqrt(4 * np.pi ** 2 / mu * a ** 3)

    #Mean anomaly
    m_anom = 2 * np.pi / period * (epoch - time_perihelion).total_seconds()
    m_anom %= (2 * np.pi)

    #Eccentric anomaly
    e_anom = newton_raphson(m_anom=m_anom, ecc=ecc, tolerance=tolerance,
                            max_iter=max_iter, e_anom_i=e_anom_i)

    #True anomaly
    nu = 2.0 * np.arctan(np.sqrt((1.0 + ecc)/(1.0 - ecc))*np.tan(0.5 * e_anom))

    #Convert from rad to deg
    nu *= 180.0/np.pi

    #Compute position and velocity with computed true anomaly
    target_r, target_v = coe2rv(a=a, ecc=ecc, inc=inc, raan=raan,
                                argp=argp, nu=nu, mu=mu)

    return target_r, target_v, nu

def newton_raphson(m_anom: float, ecc: float, tolerance: float=1e-12,
                   max_iter: float=100, e_anom_i: float | None=None) -> float:
    """Newton-Raphson solver to compute the eccentric anomaly
        for the Kepler equation

    Parameters
    ----------
    m_anom : float
        Mean anomaly [rad]
    ecc : float
        Eccentricity
    tolerance : float 
        Tolerance for the Newton-Raphson solve
    max_iter : int
        Maximum number of iterations for the Newton-Raphson solve
    e_anom_i : float
        Initial eccentric anomaly value for the Newton-Raphson solve
    
    Returns
    -------
    e_anom : list[float]
        Eccentric anomaly [rad]
    """

    if not e_anom_i:
        e_anom = m_anom
    else:
        e_anom = e_anom_i

    n_iter = 0

    while (n_iter < max_iter):
        e_anom -= (e_anom - (ecc * np.sin(e_anom)) - m_anom) / (1.0 - (ecc * np.cos(e_anom)))
        if abs((e_anom - (ecc * np.sin(e_anom)) - m_anom) / (1.0 - (ecc * np.cos(e_anom)))) <= tolerance:
            break
        n_iter += 1

    if n_iter == max_iter:
        raise Exception("Newton Raphson fo the computation of eccentric anomaly did not converge.")

    return e_anom