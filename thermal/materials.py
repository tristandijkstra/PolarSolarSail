'''
materials.py

This file includes properties for all material options, to ease
tradeoff analysis.

'''

### Sail Materials

aluminum = {'name': 'aluminum', 'emissivity': 0.10, 'reflectivity': 0.98, 'absorptivity': 0.02}
# Emissivity source: https://www.engineeringtoolbox.com/radiation-heat-emissivity-aluminum-d_433.html
# Reflectivity source: https://laserbeamproducts.wordpress.com/2014/06/19/reflectivity-of-aluminium-uv-visible-and-infrared/
# Absorptivity taken as 1 - Reflectivity.

chromium = {'name': 'chromium', 'emissivity': 0.36, 'reflectivity': 0.90, 'absorptivity': 0.10}
# Emissivity source: https://neutrium.net/heat-transfer/total-normal-emissivities-of-selected-materials/

### Spacecraft Bus Materials
aluminum_alloy = {'name': 'aluminum alloy', 'emissivity': 0.08, 'reflectivity': 0.98, 'absorptivity': 0.02}
# Emissivity source: https://www.engineeringtoolbox.com/radiation-heat-emissivity-aluminum-d_433.html

titanium = {'name': 'titanium', 'emissivity': 0.19, 'reflectivity': 0.70, 'absorptivity': 0.30}
# Emissivity source: https://www.engineeringtoolbox.com/emissivity-coefficients-d_447.html
# Reflectivity source: https://refractiveindex.info/?shelf=3d&book=metals&page=titanium

cfrp = {'name': 'carbon fibre', 'emissivity': , 'reflectivity': , 'absorptivity': }


### Coating Materials
mli = {'name': 'multi-layer insulation', 'emissivity': , 'reflectivity': , 'absorptivity': }
solar_black = {'name': 'solar black', 'emissivity': , 'reflectivity': , 'absorptivity': 0.96}
az93 = {'name': 'az93 white paint', 'emissivity': 0.92, 'reflectivity': , 'absorptivity': 0.13}
mag_oxide = {'name': 'magnesium oxide white paint', 'emissivity': 0.90, 'reflectivity': , 'absorptivity': 0.09}
osr = {'name': 'optical solar reflectors', 'emissivity': , 'reflectivity': , 'absorptivity': }


### Solar Panel Materials
gallium_arsenide = {'name': 'gallium arsenide', 'emissivity': , 'reflectivity': , 'absorptivity': }

### Heat Shield Materials
ceramic_cloth = {'name': 'ceramic cloth', 'emissivity': , 'reflectivity': , 'absorptivity': }

### Boom Materials
carbon_fibre = {'name': 'carbon fibre', 'emissivity': , 'reflectivity': , 'absorptivity': }


def sail_material(choice):
    if choice == 'standard':
        return aluminum, chromium
    else:
        raise Exception("The only currently supported option is 'standard'.")


def structural_material(choice):
    if choice == 'aluminum':
        return aluminum_alloy
    elif choice == 'titanium':
        return titanium
    elif choice == 'carbon fibre':
        return cfrp
    else:
        raise Exception("The only currently supported options are 'aluminum', 'titanium', and 'carbon fibre'.")
    
def coating_material(choice):
    if choice == 'multi-layer insulation':
        return mli
    elif choice == 'solar black':
        return solar_black
    elif choice == 'az93 white paint':
        return az93
    elif choice == 'magnesium oxide paint':
        return mag_oxide
    elif choice == 'optical solar reflectors':
        return osr
    elif choice == None:
        return
    else:
        raise Exception("The only currently supported options are 'multi-layer insulation' and 'solar black'.")
    
def panel_material(choice):
    if choice == 'gallium arsenide':
        return osr, gallium_arsenide
    else:
        raise Exception("The only currently supported option is 'gallium arsenide'.")
    
def shield_material(choice):
    if choice == 'ceramic cloth':
        return ceramic_cloth, mli
    else:
        raise Exception("The only currently supported option is 'ceramic cloth'.")
    
def boom_material(choice):
    if choice == 'carbon fibre':
        return carbon_fibre
    else:
        raise Exception("The only currently supported option is 'carbon fibre'.")