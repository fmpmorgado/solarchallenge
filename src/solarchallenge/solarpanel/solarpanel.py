import numpy as np
import numpy.typing as npt
from dataclasses import dataclass
from solarchallenge.solarpanel.power import compute_IV_PindadoCubas

@dataclass
class SolarPanel:
    """
    Class to store solar panel characteristics, based on the standard
    data sheets.
    
    Parameters
    ----------
    i_sc : float
        Short circuit current [A]
    v_oc : float
        Open circuit voltage [V]
    i_mp : float
        Current at maximum power [A]
    v_mp : float
        Voltage at maximum power [V]
    alpha_i_sc : float
        Short circuit current thermal coefficient [A/ºC]
    alpha_v_oc : float
        Open circuit voltage thermal coefficient [V/ºC]
    alpha_i_mp : float
        Current at maximum power thermal coefficient [A/ºC]
    alpha_v_mp : float
        Voltage at maximum power thermal coefficient [V/ºC]
    tref : float
        Temperature of reference [ºC]
    solar_irr_ref : float
        Solar irradiance of reference [W/m**2]
    a : float
        Ideality factor of the diode
    n_s : int
        Number of series connected cells
    eff : float
        Efficiency factor
    """

    i_sc: float
    v_oc: float
    i_mp: float
    v_mp: float
    alpha_i_sc: float
    alpha_v_oc: float
    alpha_i_mp: float
    alpha_v_mp: float
    tref: float
    solar_irr_ref: float
    a: float = 1.0
    n_s: int = 1
    eff: float = 1.0

    @classmethod
    def from_dict(cls, data):
        """
        Create an Object of class SolarPanel by reading a dict

        Parameters
        ----------
        data : dict
            Dictionary with Solar Panel characteristics

        Returns
        -------
        Contructor initialization of class SolarPanel 
        """

        return cls(**data)
    
    def compute_power(self, method: str='pindado_cubas', solar_irr: float=1367.0,
                      angle_incidence: float=0.0, temperature: float=28,
                      voltage: float|None=None) -> tuple[float, float, float] :
        """
        Compute the power produced by the solar panel.

        Currently, the power is computed using a similar method to the single diode model.

        Parameters:
        -----------
        method : str
            Method to compute the current at given conditions (I-V). Currently only pindado_cubas is
            implemented.
        solar_irr : float
            Solar irradiance [W/m**2]
        angle_incidence : float
            Angle between solar rays and panel orientation [deg]
        temperature : float
            Temperature of the solar panel [ºC]
        voltage : float | None
            Operating voltage. If None, the maximum power voltage is used

        Returns:
        --------
        power : float
            Power produced by the solar panel at given conditions [W]
        current : float
            Current of the solar panel [A]
        voltage : float
            Voltage of the solar panel [V]
        """

        if method != 'pindado_cubas':
            raise NotImplementedError("Select an available method: pindado_cubas.")
        
        current, voltage = compute_IV_PindadoCubas(solar_irr=solar_irr*np.cos(np.radians(angle_incidence)),
                                       temperature=temperature, voltage=voltage, solar_panel=self)

        power = current * voltage * self.eff

        return power, current, voltage