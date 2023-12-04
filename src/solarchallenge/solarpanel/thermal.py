# https://2018.american.conference.modelica.org/Papers/5_Posielek.pdf
# https://www.azurspace.com/images/0005979-01-01_DB_4G32C_Advanced.pdf
# https://www.mdpi.com/2226-4310/10/2/108

# https://www.researchgate.net/publication/252978798_Analytical_and_numerical_approaches_of_a_solar_array_thermal_analysis_in_a_low-earth_orbit_satellite
# https://oa.upm.es/67369/1/IEEETran_Aerospece_2021.pdf
# https://ieeexplore.ieee.org/document/7763870
# https://tfaws.nasa.gov/wp-content/uploads/On-Orbit_Thermal_Environments_TFAWS_2014.pdf

from .utils import compute_solar_zenith_angle, compute_form_factor, check_is_shadowed
from .heatflux import heat_earth_albedo, heat_earth_IR, heat_panel_radiation, heat_solar_radiation
import numpy as np
from .solarpanel import SolarPanel

def steady_state(q_in, panel: SolarPanel, solar_irr, theta, emittance_front, emittance_back, area, temp0, voltage = None):

    # Using a mid point algorithm to compute the initial temperature (q_in - q_out = 0)

    temp = temp0
    temp_min = -150 # [ºC]
    temp_max = 150 # [ºC]

    q_out = 0

    while abs(q_out - q_in) > 1E-6:
        
        power, __, __ = panel.compute_power(solar_irr=solar_irr, angle_incidence=theta, temperature=temp, voltage=voltage)
        q_out = power + heat_panel_radiation(emittance_back=emittance_back, emittance_front=emittance_front,
                                     area=area, temperature=temp)

        if q_out < q_in:
            temp_min = temp
            temp = (temp + temp_max)/2.0

        elif q_out > q_in:
            temp_max = temp
            temp = (temp + temp_min)/2.0

    return temp


def energy_solver(panel: SolarPanel, time, r_panel, r_planet, solar_irr, theta, voltage = None):

    # Compute required parameters for heatflux computation
    # This parameters are not dependent on temperature

    eta = compute_solar_zenith_angle(r_planet=r_planet, r_panel=r_panel)
    ff = compute_form_factor(r_panel=r_panel)
    shadow = check_is_shadowed(r_planet=r_planet, r_panel=r_panel, eta=eta)

    C = panel.mass * panel.c # heat capacitance
    area = panel.area
    absorptance_front = panel.absorptance_front
    emittance_front = panel.emittance_front
    absorptance_back = panel.absorptance_back
    emittance_back = panel.absorptance_front

    # Compute heat flux that are not dependent on temperature:    
    q_earth_ir = heat_earth_IR(emittance_back=emittance_back, ff = ff, area = area)

    q_earth_alb = heat_earth_albedo(absorptance_back=absorptance_back, area=area, ff = ff,
                                 solar_irr=solar_irr, eta=eta)
    
    q_solar_back, q_solar_front = heat_solar_radiation(absorptance_back=absorptance_back,
                                                       absorptance_front=absorptance_front,
                                                       area=area,
                                                       solar_irr=solar_irr,
                                                       theta=theta,
                                                       shadow=shadow)
    
    q_in = q_earth_ir + q_earth_alb + q_solar_back + q_solar_front
    
    #Initialize numpy arrays:
    v_array = np.zeros(time.shape)
    i_array = np.zeros(time.shape)
    power_array = np.zeros(time.shape)
    temp_array = np.zeros(time.shape)
    q_sat_rad = np.zeros(time.shape)

    # Find steady state at t0 (q_in - q_out = 0)
    temp_array[0] = steady_state(q_in[0], panel, solar_irr[0], theta[0],
                       emittance_front, emittance_back, area,
                       temp0 = 25, voltage = voltage)

    # Compute time-steps
    dt_array = time[1:] - time[:-1]

    # Propagate in time using forward Euler method
    for i, dt in enumerate(dt_array):
        power_array[i], i_array[i], v_array[i] = panel.compute_power(solar_irr=solar_irr[i],
                                                 angle_incidence=theta[i], temperature=temp_array[i],
                                                 voltage=voltage)
        
        q_sat_rad[i] = heat_panel_radiation(emittance_back=emittance_back,
                                            emittance_front=emittance_front,
                                            area=area, temperature=temp_array[i])

        temp_array[i+1] = temp_array[i] + (q_in[i] - power_array[i] - q_sat_rad[i]) / C * dt

    # Computation of the power at the end time
    power_array[-1], i_array[-1], v_array[-1] = panel.compute_power(solar_irr=solar_irr[-1],
                                                angle_incidence=theta[-1],
                                                temperature=temp_array[-1],
                                                voltage=voltage)
    
    q_sat_rad[-1] = heat_panel_radiation(emittance_back=emittance_back,
                                         emittance_front=emittance_front,
                                         area=area, temperature=temp_array[-1])

    #Store the heats
    q_total = {"q_earth_ir" : q_earth_ir,
               "q_earth_alb" : q_earth_alb,
               "q_solar_front": q_solar_front,
               "q_solar_back" : q_solar_back,
               "q_sat_rad" : q_sat_rad}      

    return power_array, i_array, v_array, temp_array, q_total