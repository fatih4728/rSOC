# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 17:07:10 2025

Kirchhoffs law

@author: smfadurm
"""

import numpy as np
import cantera as ct
from CoolProp.CoolProp import PropsSI

T = 800 + 273.15
Species = ['H2', 'O2', 'H2O']

def formationEnthalpy(T, Species, fH_T0=[0, 0, -241.818]):
    dH_T = []
    P = ct.one_atm
    T0 = 298.
    fH_T0 = np.array(fH_T0)
    for species in Species:
        Cp_T0 = PropsSI('Hmolar', 'P', P, 'T', T0, species)
        Cp_T = PropsSI('Hmolar', 'P', P, 'T', T, species)
        dCp = (Cp_T - Cp_T0)
        dH_T.append(dCp*1e-3)
    dH_T = np.array(dH_T)
    return fH_T0 + dH_T

fH_T = formationEnthalpy(T, Species)
print(fH_T)



