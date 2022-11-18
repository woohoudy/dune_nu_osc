# -*- coding: utf-8 -*-
"""
Created on 4th of Nov 2022
Thibaut Houdy

Last edited 18th of Nov 2022
Giorgi Kistauri

"""

import numpy as np
import random as rand
import params as osc_par

#cumulative function for measured energy spectrum 
def cumulative(y):
    cumul, f_cumul = 0, []
    f_cumul.append(0)
    for yy in y:
        cumul += yy
        f_cumul.append(cumul)
    return f_cumul

#generate energy spectrum for near or far detector
def mc_function(x_ene, f_y, Ntrials, sigma): 
    new_spectrum =[]
#    uncertainty = 0.15
    for i in range(int(Ntrials)):
        t_rand = rand.uniform(0,1)
        max_bin = max(np.where(t_rand>f_y)[0])
        try : 
            true_ene = x_ene[max_bin]
            #rint(true_ene)
            rand_energie = rand.gauss(true_ene,sigma)
            if rand_energie > 0.5:
                new_energie = rand.gauss(rand_energie,0.15*rand_energie)
                #print(new_energie)
                #new_spectrum.append(rand_energie)
                if new_energie > 0.5:
                    new_spectrum.append(new_energie)
        except : 
            continue
    #print(new_spectrum)
    return new_spectrum

#simulate neutrino spectrum using measured energy spectrum
def model_flux(x,y,N,sigma):
    #spectrum = np.loadtxt(txtfile)
    #x, y = [i[0] for i in spectrum], [i[1] for i in spectrum]
    y_cumul = cumulative(y)
    y_cumul = y_cumul/max(y_cumul)
    flux = np.array(mc_function(x,y_cumul,N,sigma))
    return flux



#  calculate oscillation probability (Normal Ordering)  [E]-GeV, [L]-km, deltaCP & a
def probability_oscillation(E=[2.5], dCP=par.delta_CP, b_neutrino=True, b_normal_hierarchy=True):   

    c1 = np.sin(osc_par.theta23)**2*np.sin(2*osc_par.theta13)**2
    c2 = np.sin(2*osc_par.theta23)*np.sin(2*osc_par.theta13)*np.sin(2*osc_par.theta12)
    c3 = np.cos(osc_par.theta23)**2*np.sin(2*osc_par.theta12)**2

    osc_par.normal_hierarchy(b_normal_hierarchy)
    osc_par.neutrinos(b_neutrino)

    if b_normal_hierarchy:
        delta_m3l = par.delta_m31
    else:   
        delta_m3l = par.delta_m32

    L = 1285 # Length in Km

    def D3l(E):
      return(1.267*par.delta_m3l*L/E)

    def D21(E):
      return(1.267*par.delta_m21*L/E)

    if dCP<-180 or dCP>180 : 
        dCP = par.delta_CP
   
    Proba = c1 * np.sin(D3l(E)-par.a_matter*L)**2/(D3l(E)-par.a_matter*L)**2 * D3l(E)**2 \
                  + c2 *np.sin(D3l(E)-par.a_matter*L)/(D3l(E)-par.a_matter*L)*D3l(E)*np.sin(par.a_matter*L)/par.a_matter/L \
                    * D21(E)*np.cos(D3l(E)+np.radians(dCP)) + c3*D21(E)**2*np.sin(par.a_matter*L)**2/(par.a_matter*L)**2
    
    return Proba

 