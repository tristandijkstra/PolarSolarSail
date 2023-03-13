'''
mass.py

This file calculates the mass of the passive elements of the thermal control system.
This does not include heaters and coolers, which require greater localization to 
estimate
'''
import numpy as np

class Properties():
    def __init__(self):
        radius = 1.5                                                                # radius in [m]
        length = 5                                                                  # length in [m]
        self.area_facing_sun = np.pi * radius ** 2                                  # area in [m^2]
        self.area_total = (2 * self.area_facing_sun) + (2 * np.pi * radius * length)# total surface area in [m^2]
    
    def insulation(self):
        # Multi-layer insulation, assuming aluminized Kapton properties
        self.kap_thickness = 0.0127                                                 # thickness in [cm]
        self.kap_density = 0.019                                                    # density in [g/cm^3]
        self.kap_n = 20                                                             # number of layers
        mli_mass = self.kap_n * self.kap_thickness * self.kap_density * self.area_total * (100**2 / 1000)
        return mli_mass
    
    def heat_shield(self, layers): 
        # Heat shield estimations
        mli_mass = self.kap_density * self.kap_thickness * self.area_facing_sun * (100**2 / 1000) * self.kap_n * layers
        ceramic_density = 96                                                        # ceramic cloth density
        ceramic_thickness = 0.013                                                   # associated thickness
        ceramic_mass = ceramic_density * ceramic_thickness * self.area_facing_sun
        total_mass = mli_mass + ceramic_mass
        return total_mass

    def passive_thermal_mass(self, layers): 
        self.passive_thermal_mass = self.insulation() + self.heat_shield(layers)    # mass of passive thermal elements
        return self.passive_thermal_mass
    
    def cost(self):
        unit_cost = self.passive_thermal_mass * 22800                               # FY10 unit cost in $K/kg
        avg_inflation = 0.0228                                                      # average inflation rate
        unit_cost_2023 = unit_cost * (1 + avg_inflation) ** 13                      # adjusting for inflation
        materials_cost = unit_cost_2023 * 0.94 * 10 ** (-6)                         # exchange rate adjustment
        return materials_cost




