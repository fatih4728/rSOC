# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 10:13:47 2025

@author: smfadurm
"""

import cantera as ct
import numpy as np

def get_thermo_properties_dict(species_list, T, P, mechanism='gri30.yaml'):
    """
    Returns enthalpy, entropy, and Gibbs free energy for a list of pure species at given T and P.
    
    Parameters:
        species_list (list of str): Species names, e.g., ['H2', 'O2'].
        T (float): Temperature in K.
        P (float): Pressure in Pa.
        mechanism (str): Cantera mechanism file (default: 'gri30.yaml').
        
    Returns:
        dict: {species: {'H': value_kJ_per_mol, 'S': value_kJ_per_molK, 'G': value_kJ_per_mol}, ...}
    """
    gas = ct.Solution(mechanism)
    results = {}
    
    for sp in species_list:
        gas.TPX = T, P, f'{sp}:1.0'
        results[sp] = {
            'H': gas.enthalpy_mole * 1e-3,   # kJ/mol
            'S': gas.entropy_mole * 1e-3,    # kJ/(mol·K)
            'G': gas.gibbs_mole * 1e-3       # kJ/mol
        }
    
    return results



if __name__=='__main__':
    # Example usage
    T = 800 + 273.15  # K
    P = ct.one_atm
    species_list = ['H2', 'O2', 'H2O']  
    TD = get_thermo_properties_dict(species_list, T, P)
    
    # # Example: get entropy of O2 directly
    # entropy_O2 = thermo_dict['O2']['S']
    # print(f'\nEntropy of O2: {entropy_O2:.2f} kJ/(mol·K)')
    
    # G_array = np.array([thermo_dict[sp]['G'] for sp in species_list])
    # dG = G_array[2] - (G_array[0] + G_array[1]*0.5)
    # print(dG)

    dG = TD['H2O']['G'] - TD['H2']['G'] - TD['O2']['G'] * 0.5


