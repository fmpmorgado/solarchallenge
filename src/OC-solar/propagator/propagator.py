
import numpy as np
from scipy.integrate import solve_ivp
from .perturbations import *

accel = 2E-2

def constant_accel_factory(accel):
    def constant_accel(t0, u, k):
        v = u[3:]
        norm_v = (v[0] ** 2 + v[1] ** 2 + v[2] ** 2) ** 0.5
        return accel * v / norm_v

    return constant_accel

def f(t0, state, k):
    du_kep = dsat_dt(t0, state, k)
    ax, ay, az = constant_accel_factory(accel)(t0, state, k)
    du_ad = np.array([0, 0, 0, ax, ay, az])

    return du_kep + du_ad


def propagate_cowell(orbit, tf=86400, rtol = 1E-11, atol = 1E-12, method = 'DOP853', t_eval = None, perturbations = None):

    def f(t0, state, k):
        
        du_kep = dsat_dt(t0, state, k)
        du_per = np.zeros(du_kep.shape)
        if perturbations != None:
            for per in perturbations:
                du_per += per(t0, state, k)
        #if perturbations is not None and len(perturbations) == 2:
        #    print(state, du_per, du_kep); exit()
        return du_kep + du_per

  #  f = dsat_dt
    r = orbit.r
    v = orbit.v

    u0 = [r[0], r[1], r[2], v[0], v[1], v[2]]
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
        dense_output=True,
    )

    if t_eval is None:
        r = result.y[:3,-1]
        v = result.y[3:,-1]
    else:
        r = result.y[:3,:]
        v = result.y[3:,:]  
    
    return r, v




#https://www.sciencedirect.com/science/article/pii/S1110016821000016
def dsat_dt(t0, u, mu):
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
    """

    x, y, z, vx, vy, vz = u
    r3 = (x**2 + y**2 + z**2) ** 1.5

    #dsat_dt = d[r v]_dt = [v a]
    du = np.array([vx, vy, vz, -mu * x / r3, -mu * y / r3, -mu * z / r3])

    return du