from solarchallenge.solarpanel.solarpanel import SolarPanel
from solarchallenge.visualization.plotting import plot_temperature, plot_heat, plot_power, plot_orbit, plot_coe
from solarchallenge.model import Model
from solarchallenge.constants import (J2_EARTH, J3_EARTH, R_MEAN_EARTH, J2000,
                                      A_EARTH, ECC_EARTH, INC_EARTH, RAAN_EARTH,
                                      ARGP_EARTH, NU_EARTH, D)
from solarchallenge.solarpanel.utils import compute_angle_solar_incidence, compute_solar_radiance
from solarchallenge.trajectory.orbit import Orbit
from solarchallenge.trajectory.perturbations import constant_accel_wrapper, J2_perturbation_wrapper, J3_perturbation_wrapper
from solarchallenge.bodies.bodies import Sun, Earth

### Setting orbits

panel_orbit = Orbit.from_vector(r = [-2384.46, 4729.01, 3050.46],
                          v = [-7.36138, -2.98997, 1.64354],
                          epoch = J2000,
                          attractor = Earth)

panel_perturbations = [constant_accel_wrapper(1E-6),
                 J2_perturbation_wrapper(J2=J2_EARTH, R=R_MEAN_EARTH)]

earth_orbit = Orbit.from_coe(a=A_EARTH, ecc=ECC_EARTH, inc=INC_EARTH, raan=RAAN_EARTH, argp=ARGP_EARTH, nu=NU_EARTH,
                             epoch=J2000,
                             attractor=Sun)

earth_perturbations = []

### Setting solar panel

panel_information = {
    'a': 1,     # ideality factor of the diode
    'n_s': 1,   # Number of series connected cells
    'i_sc': 0.506, #short circuit current
    'v_oc': 2.667, #open circuit voltage
    'i_mp': 0.487, #current at maximum power
    'v_mp': 2.371, #voltage at maximum power
    'alpha_i_sc': 0.32/1000.0, #thermal coefficients
    'alpha_v_oc': -6.0/1000.0,
    'alpha_v_mp': -6.1/1000.0,
    'alpha_i_mp': 0.28/1000.0,
    'tref': 28,  # temperature of reference
    'solar_irr_ref': 1367, #solar irradiance of reference
    'eff': 1 # Panel efficiency (To account for other losses)
}

# Only used if the thermal model is activated, but required to initialize the object
additional_panel_information = {
     "area": 1,  # Area
     "mass": 1, # Mass
     "c": 1000, # Specific heat capacity
     "emittance_back": 1.0,
     "emittance_front": 1.0,
     "absorptance_back": 1.0,
     "absorptance_front": 1.0,
}

panel = SolarPanel.from_dict(panel_information | additional_panel_information)

### Simulation driver
simulator = Model()

simulator.set_orbit_body(earth_orbit)
simulator.set_orbit_solar_panel(panel_orbit)
simulator.set_solar_panel(panel)

#Propagate the orbit for a given interval
simulator.propagate_orbit(start=J2000 + 28 * D, # The datetime for start argument needs to be later than the epochs of the orbits
                          end=J2000 + 29 * D,
                          npoints=10000, 
                          perturbations_panel = panel_perturbations,
                          perturbations_planet = earth_perturbations)

#Compute power
simulator.compute_power(temperature=28, voltage=2, maximum_power=True, thermal_model=True)

### Plots
plot_orbit(simulator.r_panel, R_MEAN_EARTH * 0.5)

plot_coe(simulator.r_panel, simulator.v_panel, simulator.time, simulator.orbit_panel.attractor.mu)

plot_temperature(time=simulator.time, temperature=simulator.temperature)

plot_heat(time=simulator.time, q_earth_ir=simulator.q_earth_ir, q_earth_alb=simulator.q_earth_alb,
          q_solar_front=simulator.q_solar_front, q_solar_back=simulator.q_solar_back,
          q_sat_rad=simulator.q_sat_rad)

solar_irr = compute_solar_radiance(simulator.r_planet, simulator.r_panel)
theta = compute_angle_solar_incidence(simulator.r_planet, simulator.r_panel)
plot_power(time=simulator.time, power=simulator.power, solar_irr=solar_irr, theta=theta)