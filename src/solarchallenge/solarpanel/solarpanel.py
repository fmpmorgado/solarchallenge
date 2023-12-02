import numpy as np
from datetime import date
from dateutil.relativedelta import relativedelta

class SolarPanel:
    #https://research-information.bris.ac.uk/ws/portalfiles/portal/216857427/SSC19_WP1_16_Tom_Etchells.pdf

    def __init__(self, efficiency: float, area: float, mission_start: date, I_d: float = 0.77, degrad: float = 0):
        
        self.eff = efficiency
        self.area = area
        self.I_d = I_d
        self.mission_start = mission_start
        self.degrad = degrad

    def power(self, irradiation: float, inclination: float, time_instant: date) -> float:
        
        life_years =(time_instant-self.mission_start).total_seconds()/3.1689E-8
        p = irradiation * self.area * self.eff * self.I_d * (1 - self.degrad)**life_years * np.cos(inclination)
        
        return p
    
class SolarPanel2:
    
    def __init__(self, design_data, area, mass):
        self.design_data = design_data
        self.area = area
        self.mass = mass