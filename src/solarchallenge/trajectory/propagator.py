"""
    The trajectory.propagator modul contains functions to propagate
    the orbit in time
"""
import numpy as np
import numpy.typing as npt
from scipy.integrate import solve_ivp
from .orbit import coe2rv, rv2coe, Orbit
from datetime import datetime

def propagate_cowell(orbit: Orbit, tf: float=86400.0, rtol: float=1E-11,
                     atol: float=1E-12, method: str = 'DOP853',
                     t_eval: None | npt.ArrayLike=None, perturbations: None | list=None
                     ) -> tuple[npt.ArrayLike, npt.ArrayLike]:
    """Propagates orbit using Cowell's formulation
    using scipy.solve_ivp module

    Parameters
    ----------
    orbit : Orbit
        orbit object to propagate
    tf : float
        Time of propagation [sec]
    rtol : float
        Relative tolerance for scipy.solve_ivp
    atol : float
        Absolute tolerance for scipy.solve_ivp
    method : str
        Name of the time propagator method to be used
        by scipy.solve_ivp: (ex: DOP853, RK43, ...)
    t_eval : np.array[float]
        Array of time epochs to store the solution. If t_eval
        set to None, the function returns the last stored solution
    perturbations : list[function]
        List of orbit perturbation functions to be applied

    Returns
    -------
    r : np.array[float]
        List of position vectors [km]
    v : np.array[float]
        List of Velocity vectors [km/s]
    """

    # Define function f to account for the perturbations
    def f(t0, state, k):
        du_kep = dsat_dt(t0, state, k)
        du_pert = np.zeros(du_kep.shape)
        
        if perturbations != None:
            for pert in perturbations:
                du_pert += pert(t0, state, k)
        
        return du_kep + du_pert

    r = orbit.r
    v = orbit.v

    # Initial state to solve the IVP
    u0 = np.array([r[0], r[1], r[2], v[0], v[1], v[2]])

    mu = orbit.attractor.mu

    result = solve_ivp(
        f,
        (0, tf),
        u0,
        args=(mu,),
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
        method=method,
        dense_output=False,
    )

    if t_eval is None:
        r = result.y[:3,-1]
        v = result.y[3:,-1]
    else:
        r = result.y[:3,:]
        v = result.y[3:,:]  
    
    return r, v


#https://www.sciencedirect.com/science/article/pii/S1110016821000016
def dsat_dt(t0: float, u: list[float], mu: float) -> npt.ArrayLike:
    """Differential equation for the initial value two body problem.
       Based on the approach used in poliastro library

    Parameters
    ----------
    t0 : float
        Time.
    u : list[float]
        Satellite state vector [x, y, z, vx, vy, vz] (km, km/s).
    mu : float
        Standard gravitational parameter.

    Returns
    -------
    du: list[float]
        Derivatives list to solve the IVP problem
    """

    x, y, z, vx, vy, vz = u
    r3 = (x ** 2 + y ** 2 + z ** 2) ** 1.5

    du = np.array([vx, vy, vz, -mu * x / r3, -mu * y / r3, -mu * z / r3])

    return du


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
        raise Exception("Newton Raphson for the computation of eccentric anomaly did not converge.")

    return e_anom