# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 13:11:46 2025

@author: smfadurm
"""

from ase.build import fcc111, add_adsorbate
from ase.calculators.emt import EMT



# Build Ni(111) surface
ni_surf = fcc111('YSZ', size=(3,3,3), vacuum=10.0)

# Add H adsorbate
add_adsorbate(ni_surf, 'H', height=1.0, position='ontop')

ni_surf.calc = EMT()
energy = ni_surf.get_potential_energy()

# Compute total energy
energy = ni_surf.get_potential_energy()
