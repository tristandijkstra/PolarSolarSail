"""
materials.py

This file includes properties for all material options, to ease
tradeoff analysis.

 
"""

### Sail Materials

standard = {
    "name": "standard",
    "emissivity": 0.30,
    "reflectivity": 0.92,
    "absorptivity": 0.08,
    "conductivity": 0.155,
    "density": 1540,
    "specific_heat": 1172.3
}

aluminum = {
    "name": "aluminum",
    "emissivity": 0.10,
    "reflectivity": 0.92,
    "absorptivity": 0.08,
    "conductivity": 237,
    "density": 2810,
    "specific_heat": 900,
}
# Emissivity source: https://www.engineeringtoolbox.com/radiation-heat-emissivity-aluminum-d_433.html
# Reflectivity source: https://laserbeamproducts.wordpress.com/2014/06/19/reflectivity-of-aluminium-uv-visible-and-infrared/
# Conductivity source: https://www.periodic-table.org/aluminium-thermal-conductivity/
# Absorptivity taken as 1 - Reflectivity.

chromium = {
    "name": "chromium",
    "emissivity": 0.30,
    "reflectivity": 0.90,
    "absorptivity": 0.10,
    "conductivity": 93.7,
    "density": 7150,
    "specific_heat": 450,
}
# Emissivity source: https://neutrium.net/heat-transfer/total-normal-emissivities-of-selected-materials/
# Conductivity source: https://material-properties.org/chromium-thermal-properties-melting-point-thermal-conductivity-expansion/

cp1 = {
    "name": "cp-1",
    "emissivity": 0.00,
    "reflectivity": 0.00,
    "absorptivity": 0.00,
    "conductivity": 0.155,
    "density": 1540,
    "specific_heat": 1172.3,
}

### Spacecraft Bus Materials
aluminum_alloy = {
    "name": "aluminum alloy",
    "emissivity": 0.08,
    "reflectivity": 0.98,
    "absorptivity": 0.02,
    "conductivity": 130,
    "density": 2810,
    "specific_heat": 960,
}
# Emissivity source: https://www.engineeringtoolbox.com/radiation-heat-emissivity-aluminum-d_433.html
# Conductivity source: https://en.wikipedia.org/wiki/7075_aluminium_alloy


titanium = {
    "name": "titanium",
    "emissivity": 0.19,
    "reflectivity": 0.70,
    "absorptivity": 0.30,
    "conductivity": 21.9,
    "density": 4540,
    "specific_heat": 520,
}
# Emissivity source: https://www.engineeringtoolbox.com/emissivity-coefficients-d_447.html
# Reflectivity source: https://refractiveindex.info/?shelf=3d&book=metals&page=titanium

cfrp = {
    "name": "carbon fibre",
    "emissivity": 0.77,
    "reflectivity": 0.17,
    "absorptivity": 0.83,
    "conductivity": 13.7,
    "density": 1800,
    "specific_heat": 800,
}
# Emissivity source: https://www.engineeringtoolbox.com/emissivity-coefficients-d_447.html
# Reflectivity source: https://refractiveindex.info/?shelf=main&book=C&page=Phillip

### Coating Materials
mli = {
    "name": "multi-layer insulation",
    "emissivity": 0.03,
    "reflectivity": 0.98,
    "absorptivity": 0.02,
    "conductivity": 1e-5,
    "density": 19,
    "specific_heat": 1000,
}
# Emissivity source: http://www.thermalengineer.com/library/effective_emittance.htm
# Reflectivity: this assumes the outer surface is aluminum for reflection

solar_black_shield = {
    "name": "solar black shield",
    "emissivity": 0.03,
    "reflectivity": 0.04,
    "absorptivity": 0.96,
    "conductivity": 0,
    "density": 3200,
    "specific_heat": 800,
}
# Emissivity and absorptivity sources: https://enbio.eu/wp-content/uploads/2018/10/SolarBlack.pdf

solar_black = {
    "name": "solar black",
    "emissivity": 0.78,
    "reflectivity": 0.04,
    "absorptivity": 0.96,
    "conductivity": 0,
    "density": 3200,
    "specific_heat": 800,
}
# Emissivity and absorptivity sources: https://enbio.eu/wp-content/uploads/2018/10/SolarBlack.pdf

az93 = {
    "name": "az93 white paint",
    "emissivity": 0.99,
    "reflectivity": 0.87,
    "absorptivity": 0.13,
    "conductivity": 0,
    "density": 3600,
    "specific_heat": 800,
}
# Density: Assumed same of magnesium oxide due to lack of available data.
# Emissivity, Reflectivity, Absorptivity: https://www.aztechnology.com/product/1/az-93

mag_oxide = {
    "name": "magnesium oxide white paint",
    "emissivity": 0.90,
    "reflectivity": 0.91,
    "absorptivity": 0.09,
    "conductivity": 0,
    "density": 3600,
    "specific_heat": 800,
}
# Thermal properties: http://solarmirror.com/fom/fom-serve/cache/43.html
# Density: https://en.wikipedia.org/wiki/Magnesium_oxide

osr = {
    "name": "optical solar reflectors",
    "emissivity": 0.93,
    "reflectivity": 0.98,
    "absorptivity": 0.02,
    "conductivity": 3,
    "density": 2600,
    "specific_heat": 741,
}
# Reflectivity source: (assuming quartz outer layer over metal)
# Density: search "optical solar reflectors density", the source is a downloadable PDF from Excelitas.

### Solar Panel Materials
gallium_arsenide = {
    "name": "gallium arsenide",
    "emissivity": 0.75,
    "reflectivity": 0.38,
    "absorptivity": 0.62,
    "conductivity": 56,
    "density": 1.76,
    "specific_heat": 327,
}
# Emissivity: http://eprints.gla.ac.uk/150163/
# Reflectivity: https://phys.org/news/2017-08-solar-cells-optics-micro-nanoscale.html
# Density: https://www.spectrolab.com/DataSheets/Panel/panels.pdf
# Note that the density is in kg/m^2, not kg/m^3 as is true for the other values.


### Heat Shield Materials
ceramic_cloth = {
    "name": "ceramic cloth",
    "emissivity": 0.67,
    "reflectivity": 0.81,
    "absorptivity": 0.67,
    "conductivity": 0.35,
    "density": 960,
    "specific_heat": 1130,
}
# Emissivity: https://www.coleparmer.com/tech-article/emissivity-of-specific-materials#anchor12
# Reflectivity: https://refractiveindex.info/?shelf=main&book=Bi12SiO20&page=Gospodinov
# Reflectivity and absorptivity taken from another member of the silicate family.
# Density: https://www.ceceramicfiber.com/Article/Ceramicfiberblanketd_1.html
# Conductivity: https://www.researchgate.net/publication/320443315_High-Temperature_Thermal_Conductivity_of_Ceramic_Fibers/figures?lo=1&utm_source=bing&utm_medium=organic

### Other Materials
batteries = {
    "name": "batteries",
    "emissivity": 0.03,
    "reflectivity": 0.98,
    "absorptivity": 0.02,
    "conductivity": 0.25,
    "density": 8760,
    "specific_heat": 820,
}

hydrazine = {
    "name": "hydrazine",
    "emissivity": 0.10,
    "reflectivity": 0.92,
    "absorptivity": 0.08,
    "conductivity": 237,
    "density": 1020,
    "specific_heat": 100,
}

metis = {
    "name": "metis",
    "emissivity": 0.19,
    "reflectivity": 0.70,
    "absorptivity": 0.30,
    "conductivity": 0.25,
    "density": 24.55/0.192,
    "specific_heat": 200,
}

cdm = {
    "name": "cdm",
    "emissivity": 0.10,
    "reflectivity": 0.92,
    "absorptivity": 0.08,
    "conductivity": 2,
    "density": 16/0.015,
    "specific_heat": 500,
}

def sail_material(choice):
    if choice == "aluminum":
        return aluminum
    elif choice == "chromium":
        return chromium
    elif choice == "cp-1":
        return cp1
    elif choice == "standard":
        return standard
    else:
        raise Exception(
            "The only currently supported options are 'aluminum', 'chromium', and 'cp-1'."
        )


def bus_material(choice):
    if choice == "aluminum":
        return aluminum_alloy
    elif choice == "titanium":
        return titanium
    elif choice == "carbon fibre":
        return cfrp
    elif choice == "multi-layer insulation":
        return mli
    elif choice == "solar black":
        return solar_black
    elif choice == "az93 white paint":
        return az93
    elif choice == "magnesium oxide paint":
        return mag_oxide
    elif choice == "optical solar reflectors":
        return osr
    else:
        raise Exception(
            "The only currently supported options are 'aluminum', 'titanium', 'carbon fibre', 'multi-layer insulation', 'solarblack', 'az93 white paint', 'magnesium oxide paint', or 'osr'."
        )


def panel_material(choice):
    if choice == "gallium arsenide":
        return gallium_arsenide
    elif choice == "osr":
        return osr
    else:
        raise Exception(
            "The only currently supported option is 'gallium arsenide' and 'osr'."
        )


def shield_material(choice):
    if choice == "ceramic cloth":
        return ceramic_cloth
    elif choice == "multi-layer insulation":
        return mli
    elif choice == "solar black shield":
        return solar_black_shield
    else:
        raise Exception("The only currently supported option is 'ceramic cloth'.")


def boom_material(choice):
    if choice == "carbon fibre":
        return cfrp
    elif choice == "multi-layer insulation":
        return mli
    elif choice == "aluminum":
        return aluminum
    else:
        raise Exception("The only currently supported option is 'carbon fibre'.")
    

def antenna_material(choice):
    if choice == "solar black":
        return solar_black
    elif choice == "titanium":
        return titanium
    else:
        raise Exception("The only currently supported option is 'carbon fibre'.")
    
def internal_material(choice):
    if choice == "battery":
        return batteries
    elif choice == "hydrazine":
        return hydrazine
    elif choice == "metis":
        return metis
    elif choice == "cdm":
        return cdm
    else:
        raise Exception("The only currently supported options are 'batteries', 'hydrazine', 'metis', and 'cdm'.")
