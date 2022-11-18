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
import random as rand
import lib_oscil as osc_mod
import lib_plot as osc_plot

Usage='''
Usage:
py study_giorgi.py -i                                                           #To get info about the soft usage
'''

def show_spectrum_Laura():
    # spectrum at FD
    models = [  y_after*osc_mod.probability_oscillation(x_after,0, b_neutrino=True, b_normal_hierarchy=True),
                y_after*osc_mod.probability_oscillation(x_after,90, b_neutrino=True, b_normal_hierarchy=True),
                y_after*osc_mod.probability_oscillation(x_after,-90, b_neutrino=True, b_normal_hierarchy=True),
                y_after*osc_mod.probability_oscillation(x_after,180, b_neutrino=True, b_normal_hierarchy=True)]

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

    
    models = [  y_after*osc_mod.probability_oscillation(x_after,0, b_neutrino=True, b_normal_hierarchy=False),
                y_after*osc_mod.probability_oscillation(x_after,90, b_neutrino=True, b_normal_hierarchy=False),
                y_after*osc_mod.probability_oscillation(x_after,-90, b_neutrino=True, b_normal_hierarchy=False),
                y_after*osc_mod.probability_oscillation(x_after,180, b_neutrino=True, b_normal_hierarchy=False)]
    axs_2[1].errorbar(x_output, y_output, np.sqrt(y_output),[(x_output[1]-x_output[0])/4 for _ in x_output], label="measured")
    for it, model in enumerate(models):
        alpha = np.max(y_output[1:])/np.max(model[1:])
        #alpha = np.sum(y_output[1:])/np.sum(model[1:])
        axs_2[1].plot(x_after, model*alpha, color =models_colors[it] ,label=model_labels[it])
    axs_2[1].legend(loc='best')
    axs_2[1].set_title("Neutrinos in IO")

def check_chi2(model, data):

    chi2 = sum ( [ ((model - data)/np.sqrt(model))])  
    fig_3[1].suptitle("chi2")
    axs_3[1].hist(all_chi2, bins = bins, color='b')
    #axs_3.plot(bins, scipy.stats.chi2.pdf(bins,))
    axs_3[1].legend(loc='best')
    axs_3[1].legend("chi2_neutrino_NO")

if __name__ == '__main__':


    N_ND = 1e6    

    #Muonic neutrinos at ND from arXiv:2109.01304v1
    near_detector_spectrum = np.loadtxt('data/input.txt')
    ene_ND_article, spec_ND_article = [i[0] for i in near_detector_spectrum], [i[1]*1e6 for i in near_detector_spectrum]

    #Electronic neutrinos at FD from arXiv:2109.01304v1
    far_detector_spectrum = np.loadtxt('data/output.txt')
    ene_FD_article, spec_FD_article = [i[0] for i in far_detector_spectrum], [i[1] for i in far_detector_spectrum]

    spec_ND_mc = osc_mod.model_flux(ene_FD_article, spec_FD_article, N_ND, 0.25)
    spec_nue_mc = spec_ND_mc*osc_mod.probability_oscillation(spec_ND_mc, -500, b_neutrino=True, b_normal_hierarchy=True)
    spec_FD_mc = rand.gauss(spec_nue_mc, 0.5*spec_nue_mc)
    plt.hist(spec_nue_mc, bins=50, range=[0,6])
    plt.hist(spec_FD_mc, bins=50, range=[0,6])
    plt.show()
    
    print("done")
    
    y_spec, x_spec = np.histogram(spec_ND_mc, bins=50, range=[0.5,6])
    ene_mc_in = x_spec[:-1]
    spec_mc_in = y_spec/max(y_spec)*max(spec_ND_mc)

    #Generate cross check plots:
    #osc_plot.check_input_with_mc(x_input, y_input, ene_mc_in, spec_mc_in)

    model = spec_mc_in*osc_mod.probability_oscillation(ene_mc_in,-500, b_neutrino=True, b_normal_hierarchy=True)

    model = model*max(spec_FD_article)/max(model)
    osc_plot.compare_data_model("Comparison data/MC expected FD", ene_FD_article, spec_FD_article, "data nu_e FD", ene_mc_in, model, "FD nu_mu oscillated spectrum")
    plt.show()
