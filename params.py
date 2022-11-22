# -*- coding: utf-8 -*-
"""
Created on 4th of Nov 2022
Thibaut Houdy

Last edited 18th of Nov 2022
Giorgi Kistauri

"""

import numpy as np

theta12 = np.radians(33.44)
theta13 = np.radians(8.57)
theta23 = np.radians(49.2)
delta_CP = np.radians(215)
delta_m21 = 7.42e-5 #eV2
delta_m31 = 2.515*pow(10,-3) #eV2
delta_m32 = -2.498*pow(10,-3) #in eV^2 in NO
a_matter = 1/3500  # GfNe/sqrt(2) for neutrino

def neutrinos(b_neutrinos):

    if b_neutrinos:
        a_matter = 1/3500  # GfNe/sqrt(2) for neutrino    
    else:
        a_matter = - 1/3500  # -GfNe/sqrt(2) for antineutrinos
        delta_CP = -1*delta_CP

def normal_hierarchy(b_normal_hierarchy):

    #For muon neutrino to electron neutrino oscillation probability
    #oscillation parameters from NuFit 5.1    
    # normal order parameters
    if b_normal_hierarchy:
        delta_CP = np.radians(194)
        delta_m31 = 2.515*pow(10,-3) #in eV^2
    else: 
    # inverted order parameters
        delta_CP = np.radians(287)
        delta_m32 = -2.498*pow(10,-3)#in eV^2
    
def change_mixing_angles(l_theta12, l_theta13, l_theta23):
    theta12 = l_theta12
    theta23 = l_theta23
    theta13 = l_theta13