import numpy as np
import pytest
from solarchallenge.trajectory.orbit import Orbit, rv2coe
from solarchallenge.trajectory.kepler import compute_kepler
from solarchallenge.bodies.bodies import Earth
from solarchallenge.propagator.propagator import propagate_cowell
from datetime import datetime, timedelta

@pytest.fixture
def iss_orbit():
        
    r = [859.07256, -4137.20368, 5295.56871] #km
    v = [7.37289205, 2.08223573, 0.43999979] #km/s

    a, ecc, inc, raan, argp, nu = rv2coe(r=r, v=v, mu=Earth.mu)

    orbit = Orbit.from_coe(a=a, ecc=ecc, inc=inc, raan=raan, argp=argp,
                           nu=0, epoch=datetime(1,1,1), attractor=Earth)

    return orbit

def test_kepler_one_orbit_period(iss_orbit):
    
    expected_r = iss_orbit.r 
    expected_v = iss_orbit.v

    a = rv2coe(r=iss_orbit.r, v=iss_orbit.v, mu=iss_orbit.attractor.mu)[0]
    
    period = np.sqrt(4 * np.pi ** 2 / iss_orbit.attractor.mu * a ** 3)

    r_kep, v_kep, __ = compute_kepler(orbit=iss_orbit, epoch=datetime(1,1,1)+timedelta(seconds=1*period),
                                      time_perihelion=datetime(1,1,1))

    assert np.isclose(r_kep, expected_r, rtol = 1E-5).all()
    assert np.isclose(v_kep, expected_v, rtol = 1E-5).all()


def test_cowell_one_orbit_period(iss_orbit):

    expected_r = iss_orbit.r 
    expected_v = iss_orbit.v

    a = rv2coe(r=iss_orbit.r, v=iss_orbit.v, mu=iss_orbit.attractor.mu)[0]
    period = np.sqrt(4 * np.pi ** 2 / iss_orbit.attractor.mu * a ** 3)

    r_cow, v_cow = propagate_cowell(orbit=iss_orbit, tf=period)

    assert np.isclose(r_cow, expected_r, rtol = 1E-5).all()
    assert np.isclose(v_cow, expected_v, rtol = 1E-5).all()

def test_cowell_is_kepler(iss_orbit):

    time_after_epoch = 123456

    r_kep, v_kep, __ = compute_kepler(orbit=iss_orbit, epoch=datetime(1,1,1)+timedelta(seconds=time_after_epoch),
                                    time_perihelion=datetime(1,1,1))
    
    r_cow, v_cow = propagate_cowell(orbit=iss_orbit, tf=time_after_epoch)

    assert np.isclose(r_cow, r_kep, rtol = 1E-5).all()
    assert np.isclose(v_cow, v_kep, rtol = 1E-5).all()
