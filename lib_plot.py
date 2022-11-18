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

def compare_input_simulation():
    # compare input and simulated energies
    fig, axs = plt.subplots(1, 1, constrained_layout=True, figsize=(10, 4))
    axs.plot(x_input, y_input/max(y_input),label='input')
    axs.plot(x_after, y_after/max(y_after),label='simulated')
    axs.legend(loc='best')
 #   plt.show()

def show_spectrum_FD():
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


 def compare_deltaCP():
 	blablabla