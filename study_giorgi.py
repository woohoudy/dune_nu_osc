# -*- coding: utf-8 -*-
"""
Created on 4th of Nov 2022
Thibaut Houdy

Last edited 18th of Nov 2022
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
import scipy.stats
import lib_oscil as osc
b_debug = False
b_save = False

Usage='''
Usage:
py toy_mc.py -i                                                           #To get info about the soft usage
'''


if __name__ == '__main__':

    N_ND = 1e6
    
    near_detector_spectrum = np.loadtxt('../input.txt')
    far_detector_spectrum = np.loadtxt('../output.txt')

    x_input, y_input = [i[0] for i in near_detector_spectrum], [i[1] for i in near_detector_spectrum]
    x_output, y_output = [i[0] for i in far_detector_spectrum], [i[1] for i in far_detector_spectrum]

    #fig, axs = plt.subplots(1, 1, constrained_layout=True, figsize=(10, 4))
    #scale_factor = np.sum(y_input)/np.sum(y_output)
    #axs.errorbar(x_input, y_input, np.sqrt(y_input), [(x_input[1]-x_input[0])/4 for _ in x_input],label='input')
    #axs.errorbar(x_output, np.array(y_output)*scale_factor, np.sqrt(np.array(y_output))*scale_factor, [(x_output[1]-x_output[0])/4 for _ in x_output], label='output')



 

    y_spec, x_spec = np.histogram(osc.model_flux(x_input,y_input,N_ND,0.25), bins=50, range=[0.5,6])
    #y_s, x_s = np.histogram(model_flux(x_output,y_output,N_FD,0.3),bins = 25)

    x_after = x_spec[:-1]
    y_after = 1.9*y_spec*1e6
    #x_s = x_s[:-1]

    # compare input and simulated energies
    fig, axs = plt.subplots(1, 1, constrained_layout=True, figsize=(10, 4))
    axs.plot(x_input, y_input/max(y_input),label='input')
    axs.plot(x_after, y_after/max(y_after),label='simulated')
    axs.legend(loc='best')
 #   plt.show()

    # spectrum at FD
    models = [  y_after*osc.probability_oscillation(x_after,0, b_neutrino=True, b_normal_hierarchy=True),
                y_after*osc.probability_oscillation(x_after,90, b_neutrino=True, b_normal_hierarchy=True),
                y_after*osc.probability_oscillation(x_after,-90, b_neutrino=True, b_normal_hierarchy=True),
                y_after*osc.probability_oscillation(x_after,180, b_neutrino=True, b_normal_hierarchy=True)]

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


    models = [  y_after*osc.probability_oscillation(x_after,0, b_neutrino=True, b_normal_hierarchy=False),
                y_after*osc.probability_oscillation(x_after,90, b_neutrino=True, b_normal_hierarchy=False),
                y_after*osc.probability_oscillation(x_after,-90, b_neutrino=True, b_normal_hierarchy=False),
                y_after*osc.probability_oscillation(x_after,180, b_neutrino=True, b_normal_hierarchy=False)]
    axs_2[1].errorbar(x_output, y_output, np.sqrt(y_output),[(x_output[1]-x_output[0])/4 for _ in x_output], label="measured")
    for it, model in enumerate(models):
        alpha = np.max(y_output[1:])/np.max(model[1:])
        #alpha = np.sum(y_output[1:])/np.sum(model[1:])
        axs_2[1].plot(x_after, model*alpha, color =models_colors[it] ,label=model_labels[it])
    axs_2[1].legend(loc='best')
    axs_2[1].set_title("Neutrinos in IO")

    plt.show()


    #fig_3, axs_3 = plt.subplots(1, 3, constrained_layout=True, figsize=(10, 4))
    #data = y_s
    model = y_after*probability_oscillation(x_after,0, b_neutrino=True, b_normal_hierarchy=True)
    print(len(model))
    bins = np.linspace(0,10,50)
    plt.hist(model, bins =bins )
    plt.show()
'''
    chi2 = sum ( [ ((model - data)/np.sqrt(model))])  
    fig_3[1].suptitle("chi2")
    axs_3[1].hist(all_chi2, bins = bins, color='b')
    #axs_3.plot(bins, scipy.stats.chi2.pdf(bins,))
    axs_3[1].legend(loc='best')
    axs_3[1].legend("chi2_neutrino_NO")
'''
    

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