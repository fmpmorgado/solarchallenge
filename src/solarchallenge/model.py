import numpy as np
import numpy.typing as npt
from solarchallenge.trajectory.propagator import propagate_cowell
from solarchallenge.solarpanel.solarpanel import SolarPanel
from solarchallenge.trajectory.orbit import Orbit
from solarchallenge.solarpanel.thermal import energy_solver
from solarchallenge.solarpanel.utils import compute_solar_radiance, compute_angle_solar_incidence
from dataclasses import dataclass
from datetime import datetime


@dataclass
class OrbitData:
    """
    Class to store the Orbit data after propagation (inherited by the Model class)

    Parameters:
    ----------
    time : np.array[float]
        Correspondent time of stored position and velocity in UTC format
    r_panel : np.array[float]
        Solar panel position array in x, y and z direction - J2000 frame [km]
    v_panel : np.array[float]
        Solar panel velocity array in x, y and z direction - J2000 frame [km/s]
    r_planet : np.array[float]
        Earth position array in x, y and z direction - EME2000 frame [km]
    v_planet : np.array[float]
        Earth velocity array in x, y and z direction - EME2000 frame [km/s]    
    """

    time: npt.ArrayLike
    r_panel: npt.ArrayLike
    v_panel: npt.ArrayLike
    r_planet: npt.ArrayLike
    v_planet: npt.ArrayLike

    def store_orbit_data(self, time: npt.ArrayLike, r_panel: npt.ArrayLike,
                         v_panel: npt.ArrayLike, r_planet:npt.ArrayLike, v_planet: npt.ArrayLike) -> None:
        """
        Stores the orbit data in the respective object variables and checks for consistency.

        Parameters:
        ----------
        time : np.array[float]
            Correspondent time of stored position and velocity in UTC format
        r_panel : np.array[float]
            Solar panel position array in x, y and z direction - J2000 frame [km]
        v_panel : np.array[float]
            Solar panel velocity array in x, y and z direction - J2000 frame [km/s]
        r_planet : np.array[float]
            Earth position array in x, y and z direction - EME2000 frame [km]
        v_planet : np.array[float]
            Earth velocity array in x, y and z direction - EME2000 frame [km/s] 
        """

        # Check if lists are same size
        if len(time) != len(r_panel) or \
           len(time) != len(v_panel) or \
           len(time) != len(r_planet) or \
           len(time) != len(v_planet):
            raise ValueError("Lists need to have the same size.")

        # Store the orbit data
        self.time = time
        self.r_panel = r_panel
        self.v_panel = v_panel
        self.r_planet = r_planet
        self.v_planet = v_planet


@dataclass
class HeatData:
    """
    Class to store the Heat data acting on the panel during simulated time
    (inherited by the Model class)

    Parameters:
    ----------
    q_earth_ir : np.array[float]
        Heat due to Earth IR [W]
    q_earth_alb : np.array[float]
        Heat due to albedo effect [W]
    q_solar_front : np.array[float]
        Heat due to solar irradiance on the photovoltaic side [W]
    q_solar_back : np.array[float]
        Heat due to solar irradiance on the back side [W]
    q_sat_rad : np.array[float]
        Heat leaving the satellite through radiation (Black body radiation) [W]    
    """

    q_earth_ir: npt.ArrayLike
    q_earth_alb: npt.ArrayLike
    q_solar_front: npt.ArrayLike
    q_solar_back: npt.ArrayLike
    q_sat_rad: npt.ArrayLike

    def store_heat_data(self, q_earth_ir: npt.ArrayLike, q_earth_alb: npt.ArrayLike,
                        q_solar_front: npt.ArrayLike, q_solar_back: npt.ArrayLike,
                        q_sat_rad: npt.ArrayLike):
        """
        Stores the heat data in the respective object variables.

        Parameters:
        ----------
        q_earth_ir : np.array[float]
            Heat due to Earth IR [W]
        q_earth_alb : np.array[float]
            Heat due to albedo effect [W]
        q_solar_front : np.array[float]
            Heat due to solar irradiance on the photovoltaic side [W]
        q_solar_back : np.array[float]
            Heat due to solar irradiance on the back side [W]
        q_sat_rad : np.array[float]
            Heat leaving the satellite through radiation (Black body radiation) [W]    
        """

        self.q_earth_ir = q_earth_ir
        self.q_earth_alb = q_earth_alb
        self.q_solar_front = q_solar_front
        self.q_solar_back = q_solar_back
        self.q_sat_rad = q_sat_rad
        

@dataclass
class PowerData:
    """
    Class to store the computed power and temperature data during simulated time
    (inherited by the Model class)

    Parameters:
    ----------
    voltage : np.array[float]
        Voltage of the solar panel [V]
    current : np.array[float]
        Current of the solar panel [A]
    power : np.array[float]
        Power produced by the solar panel [W]
    temperature : np.array[float]
        Temperature of the solar panel [ºC]
    """

    voltage: npt.ArrayLike
    current: npt.ArrayLike
    power: npt.ArrayLike
    temperature: npt.ArrayLike

    def store_power_data(self, voltage: npt.ArrayLike, current: npt.ArrayLike,
                         power: npt.ArrayLike, temperature: npt.ArrayLike):
        """
        Stores the power data in the respective object variables.

        Parameters:
        ----------
        voltage : np.array[float]
            Voltage of the solar panel [V]
        current : np.array[float]
            Current of the solar panel [A]
        power : np.array[float]
            Power produced by the solar panel [W]
        temperature : np.array[float]
            Temperature of the solar panel [ºC]
        """

        self.voltage = voltage
        self.current = current
        self.power = power
        self.temperature = temperature


class Model(OrbitData, PowerData, HeatData):
    """
    The Model class is used as an interface with the submodules
    to propagate the planet and solar panel orbits, compute the
    produced power and perform thermal analysis. The class also
    stores the computed data.

    Parameters:
    ----------
    orbit_panel : Orbit
        Object of Orbit class containing information on the panel orbit.
    orbit_body : Orbit
        Object of Orbit class containing information on the planet orbit.
    panel : SolarPanel
        Object of SolarPanel class containing information on the solar
        panel characteristics.
    """

    def __init__ (self):
        self.orbit_panel: Orbit
        self.orbit_body: Orbit
        self.panel: SolarPanel

    def set_orbit_solar_panel(self, orbit_panel: Orbit):
        """
        Set the solar panel orbit

        Parameters
        ----------
        orbit_panel : Orbit
            Object of Orbit class containing information of the panel orbit.
        """
        self.orbit_panel = orbit_panel

    def set_orbit_body(self, orbit_body: Orbit):
        """
        Set the planet orbit

        Parameters
        ----------
        orbit_panel : Orbit
            Object of Orbit class containing information of the panel orbit.
        """
        self.orbit_body = orbit_body

    def set_solar_panel(self, panel: SolarPanel):
        """
        Set the solar panel

        Parameters
        ----------
        panel : SolarPanel
            Object of SolarPanel class containing information on the solar
            panel characteristics.
        """
        self.panel = panel

    def propagate_orbit(self, start: datetime, end: datetime, npoints: int=10000,
                        ivp_method: str='DOP853',
                        rtol: float=1E-11, atol: float=1E-12,
                        perturbations_panel: list|None=None,
                        perturbations_planet: list|None=None):
        """
        Propagate the orbit of the planet and the solar panel by solving the
        orbit initial value problem, accounting for perturbations.

        The orbits are propagated until the specified end time. The solution
        between the start time and end time is stored using an array size
        of npoints evenly distributed in time.

        The propagation of the solar panel and planet start at the orbit state
        epoch. After propagation, the solution is stored.

        Parameters
        ----------
        start : Datetime
            Start time to store results from propagation.
        end : Datetime
            End time of propagation.
        npoints : int
            Number of points to store propagation solution. The points are
            distributed between start and end time.
        ivp_method : str
            Method used to solve the ivp problem (DOP853, RK45). The method
            is passed to the scipy.solve_ivp function.
        rtol : float
            Relative tolerance for orbit propagation.
        atol : float
            Absolute tolerance for orbit propagation.
        perturbations_panel : list or None
            Perturbations for the solar panel orbit, given as a list of
            functions that can be called by scipy.solve_ivp
            (See solarchallenge.trajectory.perturbations).
            If None, no perturbations are assumed.
        perturbations_planet : list or None
            Perturbations for the planet orbit, given as a list of functions
            that can be called by scipy.solve_ivp
            (See solarchallenge.trajectory.perturbations).
            If None, no perturbations are assumed.
        """

        #Verify if orbits are stored:
        if not self.orbit_panel:
            raise Exception("Need to set the orbit of the solar panel.")
        if not self.orbit_body:
            raise Exception("Need to set the orbit of the planet (ex: Earth orbit).")

        #Verify if end date is after start date:
        if end < start:
            raise ValueError("Start date needs to be before end date")
        
        #Verify if start date is not before orbit epochs
        if start < self.orbit_panel.epoch or start < self.orbit_body.epoch:
            raise ValueError("Start date needs to be after orbit epoch")

        #Compute time to start and stop storing propagation solution
        t0_body = (start - self.orbit_body.epoch).total_seconds()
        tf_body = (end - self.orbit_body.epoch).total_seconds()

        t0_panel = (start - self.orbit_panel.epoch).total_seconds()
        tf_panel = (end - self.orbit_panel.epoch).total_seconds()

        #Compute time array to store solution.
        #Passed as t_eval argument to scipy.solve_ivp()
        t_list_body = np.linspace(t0_body, tf_body, npoints)
        t_list_panel = np.linspace(t0_panel, tf_panel, npoints)

        #Propagate orbit
        r_planet, v_planet = propagate_cowell(orbit=self.orbit_body, tf=tf_body,
                                              rtol=rtol, atol=atol, method=ivp_method,
                                              t_eval=t_list_body, perturbations=perturbations_planet)
        r_panel, v_panel = propagate_cowell(orbit=self.orbit_panel, tf=tf_panel,
                                            rtol=rtol, atol=atol, method=ivp_method,
                                            t_eval=t_list_panel, perturbations=perturbations_panel)

        t_list_UTC = np.linspace(start.timestamp(), end.timestamp(), npoints)

        r_planet = np.swapaxes(r_planet, 0, 1)
        v_planet = np.swapaxes(v_planet, 0, 1)
        r_panel = np.swapaxes(r_panel, 0, 1)
        v_panel = np.swapaxes(v_panel, 0, 1)

        #Store orbit data
        self.store_orbit_data(time=t_list_UTC, r_panel=r_panel, v_panel=v_panel,
                              r_planet=r_planet, v_planet=v_planet) 

    def compute_power(self, temperature=28, voltage=2.6, method_power='pindado_cubas',
                      maximum_power=False, thermal_model = False):
        """
        Compute the power generated by the solar panel in the propagated orbit
        for a given temperature and voltage of operation. The power is calculated
        using the Pindado_Cubas method. 

        If maximum_power is set to True, the voltage of operation is set to be the voltage
        at maximum power.

        If thermal_model is set to True, the temperature of the solar panel for
        the computed orbit propagation is evaluated using a lumped mass thermal model.
        The use of a thermal model is useful as the power generated by the solar panel
        is affected by the temperature.

        Parameters
        ----------
        temperature : float
            Temperature of the solar panel [ºC]
        voltage : float
            Voltage of operation [ºC]
        method_power : str
            Method to compute the solar panel produced power (pindado_cubas)
        maximum_power : bool
            Boolean flag. If True, the voltage at maximum power is used.
            If False, the voltage passed as an argument is used.
        thermal_model : bool
            Boolean flag. If True, a lumped mass thermal model is applied to
            compute the temperature of the solar panel for the propagated orbit.
            The initial temperature value is retrieved by solving the steady-state
            solution. If False, the temperature passed as an argument is used
        """
        if not self.panel:
            raise Exception("Need to set the solar panel.")

        if thermal_model == True:
            try:
                self.panel.area #Check if variable is defined
            except:
                raise Exception("Parameters for the thermal analysis need to be specified.")

        if maximum_power is True:
            voltage = None

        #Compute angle of incidence
        theta = compute_angle_solar_incidence(self.r_planet, self.r_panel)
        
        #Compute solar irradiance
        solar_irr = compute_solar_radiance(self.r_planet, self.r_panel)

        if thermal_model is False:
            power, cur, vol = self.panel.compute_power(method=method_power,
                                                       solar_irr=solar_irr,
                                                       angle_incidence=theta,
                                                       temperature=temperature,
                                                       voltage=voltage)

            self.store_power_data(voltage=vol, current=cur, power=power, temperature=np.ones(power.shape) * temperature)
    
        else: #Thermal model is active
            power, cur, vol, temp, heat = energy_solver(panel=self.panel, time=self.time,
                                                                       r_panel=self.r_panel,
                                                                       r_planet=self.r_planet,
                                                                       solar_irr=solar_irr,
                                                                       theta=theta,
                                                                       voltage=voltage)
    
            self.store_power_data(voltage=vol, current=cur, power=power, temperature=temp)
            self.store_heat_data(**heat)
        """
        import plotly.express as px
        import pandas as pd

        df = pd.DataFrame(dict(
        x = self.time,
        y = T
        ))

        df = pd.DataFrame(Q)
        df['x'] = self.time

        fig = px.line(df, x = 'x', y = df.columns[:-1], title='Life expectancy in Canada')
        fig.show()

        """

"""
import plotly.graph_objs as go

def ms(x, y, z, radius, resolution=20):
    u, v = np.mgrid[0:2*np.pi:resolution*2j, 0:np.pi:resolution*1j]
    X = radius * np.cos(u)*np.sin(v) + x
    Y = radius * np.sin(u)*np.sin(v) + y
    Z = radius * np.cos(v) + z
    return (X, Y, Z)

def plot_orbit(r):
    fig = go.Figure()
    fig.add_trace(
        go.Scatter3d(
            x=r[:,0],
            y=r[:,1],
            z=r[:,2],
            mode='lines',
        ))
    
    (X, Y, Z) = ms(0,0,0,6000)
    fig.add_trace(
        go.Surface(
            x=X,
            y=Y,
            z=Z,
        ))


    #fig.show(renderer="png")
#from plotly.offline import iplot, init_notebook_mode
#init_notebook_mode()
#iplot(fig)
"""