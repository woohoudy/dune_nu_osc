# -*- coding: utf-8 -*-
"""
Created on 4th of Nov 2022
Thibaut Houdy

Last edited 15th of Nov 2022
Giorgi Kistauri

"""
import math
import numpy as np
import matplotlib.pyplot as plt
import os, sys, re
from scipy.signal import savgol_filter, find_peaks
#import ROOT
import argparse
import time
import random as rand
#import uncertainties as u
from uncertainties import unumpy as unp
from uncertainties import ufloat as ufl 
b_debug = False
b_save = False

Usage='''
Usage:
py toy_mc.py -i                                                           #To get info about the soft usage
'''

def defineParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infos', dest="info", action="store_true", default=False, help='get informations')
    args = parser.parse_args()
    return args


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
                if new_energie > 0.5:
                    new_spectrum.append(new_energie)
        except : 
            continue
    #print(new_spectrum)
    return new_spectrum



def load_parameters(b_normal_hierarchy):

    #For muon neutrino to electron neutrino oscillation probability
    #oscillation parameters from NuFit 5.1    
    # normal order parameters
    if b_normal_hierarchy:
        theta12 = np.radians(33.44)
        theta13 = np.radians(8.57)
        theta23 = np.radians(49.2)
        bestfit_delta_CP = np.radians(194)
        delta_m3l = 2.515*pow(10,-3) #in eV^2
    else: 
    # inverted order parameters
        theta12 = np.radians(33.45)
        theta13 = np.radians(8.6)
        theta23 = np.radians(49.5)
        bestfit_delta_CP = np.radians(287)
        delta_m3l = -2.498*pow(10,-3)
    c1 = np.sin(theta23)**2*np.sin(2*theta13)**2
    c2 = np.sin(2*theta23)*np.sin(2*theta13)*np.sin(2*theta12)
    c3 = np.cos(theta23)**2*np.sin(2*theta12)**2


    return c1, c2, c3, bestfit_delta_CP, delta_m3l

#  calculate oscillation probability (Normal Ordering)  [E]-GeV, [L]-km, deltaCP & a

def probability_oscillation(E, dCP, b_neutrino, b_normal_hierarchy):   

    c1, c2, c3, bestfit_delta_CP , delta_m3l = load_parameters(b_normal_hierarchy)
   
    L = 1000 # Length in Km
    delta_m21= 7.42*pow(10,-5) #in eV^2

    def D3l(E):
      return(1.267*delta_m3l*L/E)

    def D21(E):
      return(1.267*delta_m21*L/E)

    if b_neutrino :
        matter_effect = 1/3500  # GfNe/sqrt(2) for neutrino    
    else: 
        matter_effect = - 1/3500  # -GfNe/sqrt(2) for antineutrinos
        dCP = -dCP

   
    Proba = c1 * np.sin(D3l(E)-matter_effect*L)**2/(D3l(E)-matter_effect*L)**2 * D3l(E)**2 \
                  + c2 *np.sin(D3l(E)-matter_effect*L)/(D3l(E)-matter_effect*L)*D3l(E)*np.sin(matter_effect*L)/matter_effect/L \
                    * D21(E)*np.cos(D3l(E)+np.radians(dCP)) + c3*D21(E)**2*np.sin(matter_effect*L)**2/(matter_effect*L)**2
    
    return Proba


#simulate neutrino flux using measured energy spectrum
def model_flux(x,y,N,sigma):
    #spectrum = np.loadtxt(txtfile)
    #x, y = [i[0] for i in spectrum], [i[1] for i in spectrum]
    y_cumul = cumulative(y)
    y_cumul = y_cumul/max(y_cumul)
    flux = np.array(mc_function(x,y_cumul,N,sigma))
    return flux

 

if __name__ == '__main__':

    N_ND = 1e6
    N_FD = 1e3

    near_detector_spectrum = np.loadtxt('../input.txt')
    far_detector_spectrum = np.loadtxt('../output.txt')

    x_input, y_input = [i[0] for i in near_detector_spectrum], [i[1]*1e6 for i in near_detector_spectrum]
    x_output, y_output = [i[0] for i in far_detector_spectrum], [i[1] for i in far_detector_spectrum]

    fig, axs = plt.subplots(1, 1, constrained_layout=True, figsize=(10, 4))
    scale_factor = np.sum(y_input)/np.sum(y_output)
    axs.errorbar(x_input, y_input, np.sqrt(y_input), [(x_input[1]-x_input[0])/4 for _ in x_input],label='input')
    axs.errorbar(x_output, np.array(y_output)*scale_factor, np.sqrt(np.array(y_output))*scale_factor, [(x_output[1]-x_output[0])/4 for _ in x_output], label='output')



 

    y_spec, x_spec = np.histogram(model_flux(x_input,y_input,N_ND,0.3), bins=50, range=[0.5,6])
    #y_s, x_s = np.histogram(model_flux(x_output,y_output,N_FD,0.3),bins = 25)

    x_after = x_spec[:-1]
    y_after = 1.9*y_spec

    # compare input and simulated energies
    #fig, axs = plt.subplots(1, 1, constrained_layout=True, figsize=(10, 4))
    #axs.plot(x_input, y_input/max(y_input),label='input')
    #axs.plot(x_after, y_after/max(y_after),label='simulated')
    #axs.legend(loc='best')


    # spectrum at FD
    models = [  y_after*probability_oscillation(x_after,0, b_neutrino=True, b_normal_hierarchy=True),
                y_after*probability_oscillation(x_after,90, b_neutrino=True, b_normal_hierarchy=True),
                y_after*probability_oscillation(x_after,-90, b_neutrino=True, b_normal_hierarchy=True),
                y_after*probability_oscillation(x_after,180, b_neutrino=True, b_normal_hierarchy=True)]

    model_labels = ["$\delta_{CP}=0$", "$\delta_{CP}=\pi/2$", "$\delta_{CP}=_\pi/2$", "$\delta_{CP}=\pi$"]
    models_colors = ['g', 'b', 'r', 'orange']

    fig_2, axs_2 = plt.subplots(1, 2, constrained_layout=True, figsize=(10, 4))
    fig_2.suptitle("Spectrum at FD")
    axs_2[0].errorbar(x_output, y_output, np.sqrt(y_output), [(x_output[1]-x_output[0])/4 for _ in x_output], linestyle=None, label="measured")
    for it, model in enumerate(models):
        alpha = np.max(y_output[1:])/np.max(model[1:])
        #alpha = np.sum(y_output[1:])/np.sum(model[1:])
        axs_2[0].plot(x_after, model*alpha, color =models_colors[it] ,label=model_labels[it])
    axs_2[0].legend(loc='best')
    axs_2[0].set_title("Neutrinos in NO")


    models = [  y_after*probability_oscillation(x_after,0, b_neutrino=True, b_normal_hierarchy=False),
                y_after*probability_oscillation(x_after,90, b_neutrino=True, b_normal_hierarchy=False),
                y_after*probability_oscillation(x_after,-90, b_neutrino=True, b_normal_hierarchy=False),
                y_after*probability_oscillation(x_after,180, b_neutrino=True, b_normal_hierarchy=False)]
    axs_2[1].errorbar(x_output, y_output, np.sqrt(y_output),[(x_output[1]-x_output[0])/4 for _ in x_output], label="measured")
    for it, model in enumerate(models):
        alpha = np.max(y_output[1:])/np.max(model[1:])
        #alpha = np.sum(y_output[1:])/np.sum(model[1:])
        axs_2[1].plot(x_after, model*alpha, color =models_colors[it] ,label=model_labels[it])
    axs_2[1].legend(loc='best')
    axs_2[1].set_title("Neutrinos in IO")
    
    plt.show()
'''
a=[]
b=[]
c=[]
d=[]
E= 2.5
a = [   probability_oscillation(E,0, b_neutrino=True, b_normal_hierarchy=True),
        probability_oscillation(E,90, b_neutrino=True, b_normal_hierarchy=True),
        probability_oscillation(E,180, b_neutrino=True, b_normal_hierarchy=True),
        probability_oscillation(E,-90, b_neutrino=True, b_normal_hierarchy=True) ]
b = [   probability_oscillation(E,0, b_neutrino=True, b_normal_hierarchy=False),
        probability_oscillation(E,-90, b_neutrino=True, b_normal_hierarchy=False),
        probability_oscillation(E,-180, b_neutrino=True, b_normal_hierarchy=False),
        probability_oscillation(E,90, b_neutrino=True, b_normal_hierarchy=False) ]
c = [   probability_oscillation(E,0, b_neutrino=False, b_normal_hierarchy=True),
        probability_oscillation(E,90, b_neutrino=False, b_normal_hierarchy=True),
        probability_oscillation(E,180, b_neutrino=False,b_normal_hierarchy=True),
        probability_oscillation(E,-90, b_neutrino=False, b_normal_hierarchy=True) ]        
d = [   probability_oscillation(E,0, b_neutrino=False, b_normal_hierarchy=False),
        probability_oscillation(E,-90, b_neutrino=False, b_normal_hierarchy=False),
        probability_oscillation(E,-180, b_neutrino=False, b_normal_hierarchy=False),
        probability_oscillation(E,90, b_neutrino=False, b_normal_hierarchy=False) ]

print(a,c,b,d)
plt.plot(a,c, linestyle = '--',marker='o',color = 'g', label = 'Normal order')        
plt.plot(b,d, linestyle = '-.',marker='+', color = 'b', label = 'Inverted order')
plt.xlabel("Neutrinos")
plt.ylabel("Antineutrinos")
plt.legend(loc='best')
plt.show()


#e = np.arange(0.1, 10., 0.001) #Neutrino energy uniform 0-10 GeV 
model_incoming_flux = np.array(model_incoming_flux)
e = np.sort(model_incoming_flux)
#print(type(e))
#print(p)


#visualize
plt.figure(1) #Ntrials
#plt.subplot(221)
#dCP & a positive for neutrinos
plt.plot(e,p_NO(e,0, a1),color='g', label='d_CP = 0')
plt.plot(e,p_NO(e,90,a1),color='r', label='d_CP = pi/2')
plt.plot(e,p_NO(e,-90,a1), color='b',label='d_CP = -pi/2')
plt.plot(e,p_NO(e,180,a1), color='orange',label='d_CP = pi') 
plt.xlim([0, 6])
plt.ylim([0, 0.2])
#plt.xlabel('Energy')
plt.ylabel('Probability')
#plt.xscale('log')
plt.title("Neutrino_NO")
plt.legend()

plt.subplot(222)
#dCP & a negative for antineutrinos
plt.plot(e,p_NO(e,0,a2), color='g',label='d_CP = 0')
plt.plot(e,p_NO(e,-90, a2), color='r', label='d_CP = pi/2')
plt.plot(e,p_NO(e,90,a2), color='b',label='d_CP = -pi/2')
plt.plot(e,p_NO(e,-180, a2), color='orange', label='d_CP = pi')
plt.xlim([0, 6])
plt.ylim([0, 0.2])
#plt.xlabel('Energy')
plt.ylabel('Probability')
#plt.xscale('log')
plt.title("Antineutrino_NO")
plt.legend()

plt.subplot(223)
#dCP & a positive for neutrinos
plt.plot(e,p_IO(e,0, a1),color='g', label='d_CP = 0')
plt.plot(e,p_IO(e,90,a1),color='r', label='d_CP = pi/2') 
plt.plot(e,p_IO(e,-90, a1),color='b', label='d_CP = -pi/2')
plt.plot(e,p_IO(e,180,a1),color='orange', label='d_CP = pi') 
plt.xlim([0, 6])
plt.ylim([0, 0.2])
plt.xlabel('Energy')
plt.ylabel('Probability')
#plt.xscale('log')
plt.title("Neutrino_IO")
plt.legend()

plt.subplot(224)
#dCP & a negative for antineutrinos
plt.plot(e,p_IO(e,0,a2), color='g',label='d_CP = 0')
plt.plot(e,p_IO(e,-90, a2), color='r', label='d_CP = pi/2')
plt.plot(e,p_IO(e,90,a2), color='b',label='d_CP = -pi/2')
plt.plot(e,p_IO(e,-180, a2), color='orange', label='d_CP = pi')
plt.xlim([0, 6])
plt.ylim([0, 0.2])
plt.xlabel('Energy')
plt.ylabel('Probability')
#plt.xscale('log')
plt.title("Antineutrino_IO")
plt.legend()

plt.show()
plt.savefig('osc_pic.png', bbox_inches='tight')
'''