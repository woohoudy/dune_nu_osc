# -*- coding: utf-8 -*-
"""
Created on 4th of Nov 2022
Thibaut Houdy

Last edited 22th of Nov 2022
Giorgi Kistauri

"""
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.markers as mrks
from scipy import interpolate 
import os, sys, re
import random as rand
import lib_oscil as osc_mod
import lib_plot as osc_plot
import params as osc_par
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve

Usage='''
Usage:
py study_giorgi.py -i                                                           #To get info about the soft usage
'''


N_ND = 1e6    
Baseline = 1285 # Length in Km


def show_spectrum_Laura(Energy):

    #deltaCP = [-90, 0, 90, 180, -90]
    deltaCP = np.linspace(-180,180,360)
    deltaCP_spec = [89, 179, 269, 359]

    a = [ osc_mod.probability_oscillation(E=Energy, L=1000, dCP=_dCP,b_normal_hierarchy=True , b_neutrino=True) for _dCP in deltaCP]
    b = [ osc_mod.probability_oscillation(E=Energy, L=1000, dCP=_dCP,b_normal_hierarchy=False ,b_neutrino=True ) for _dCP in deltaCP]
    c = [ osc_mod.probability_oscillation(E=Energy, L=1000, dCP=_dCP,b_normal_hierarchy=True , b_neutrino=False) for _dCP in deltaCP]        
    d = [ osc_mod.probability_oscillation(E=Energy, L=1000, dCP=_dCP,b_normal_hierarchy=False , b_neutrino=False) for _dCP in deltaCP] 
    osc_par.theta23 = np.radians(44)
    a2 = [ osc_mod.probability_oscillation(E=Energy, L=1000, dCP=_dCP,b_normal_hierarchy=True , b_neutrino=True) for _dCP in deltaCP]
    b2 = [ osc_mod.probability_oscillation(E=Energy, L=1000, dCP=_dCP,b_normal_hierarchy=False ,b_neutrino=True ) for _dCP in deltaCP]
    c2 = [ osc_mod.probability_oscillation(E=Energy, L=1000, dCP=_dCP,b_normal_hierarchy=True , b_neutrino=False) for _dCP in deltaCP]        
    d2 = [ osc_mod.probability_oscillation(E=Energy, L=1000, dCP=_dCP,b_normal_hierarchy=False , b_neutrino=False) for _dCP in deltaCP] 
    

    markers = ['P','o','X','s','P']
    linestyles = ['--','--','-.','-.']
    colors = ['y','black','y','black',]
    labels = ['NO','IO']
    fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=(6, 4))
    
    #plt.plot(a,c, ls ="--",label ='NO')
    #plt.scatter(a,c, marker = "+",color = ["red", "black","red", "black","red"])
    p_nu = [a, b, a2, b2]
    p_antinu = [c, d, c2, d2]

    for i in range(len(p_nu)):
        #print(i, p_nu[i], p_antinu[i])
        ax.plot(np.array(p_nu[i]), np.array(p_antinu[i]), linestyle = linestyles[i], color = colors[i]) 
        for j in range(len(deltaCP_spec)):
            ax.scatter(p_nu[i][deltaCP_spec[j]], p_antinu[i][deltaCP_spec[j]], marker = markers[j], s = 50, facecolor = colors[i], linewidth = 1, edgecolors = 'red') 

    
    ax.set_xlabel(r'$P ( \nu_\mu \to \nu_e )$', size = 12)
    ax.set_ylabel(r'$P ( \~\nu_\mu \to \~\nu_e )$', size = 12)
    #ax.set_xlim(0.005,0.08)
    #ax.set_ylim(0.005,0.08)
    ax.set_title("L = 1000 km, E = 2.5 GeV", size = 13)
    ax.legend()
    plt.xlim(0.005,0.08)
    plt.ylim(0.005,0.08)

#show_spectrum_Laura(2.5)
#plt.show()



    y_lines = mlines.Line2D([],[], color = 'y', label = 'Normal order')
    b_lines = mlines.Line2D([],[], color = 'black', label = 'Inverted order')
    dashed_lines = mlines.Line2D([],[], color = 'black', linestyle = '--', label = '$ \Theta > 45^\circ$')
    dotdashed_lines = mlines.Line2D([],[], color = 'black', linestyle = '-.', label = '$ \Theta < 45^\circ$')
    marker1 = mlines.Line2D([],[], color = 'black', linestyle = 'none', marker = markers[0], markersize = 6, label = '$ \delta_{CP} = -90^\circ $')
    marker2 = mlines.Line2D([],[], color = 'black', linestyle = 'none',marker = markers[1], markersize = 6, label = '$ \delta_{CP} = 0^\circ $')
    marker3 = mlines.Line2D([],[], color = 'black', linestyle = 'none',marker = markers[2], markersize = 6, label = '$ \delta_{CP} = 90^\circ $')
    marker4 = mlines.Line2D([],[], color = 'black', linestyle = 'none',marker = markers[3], markersize = 6, label = '$ \delta_{CP} = 180^\circ $')
    legend1 = ax.legend(handles = [marker1, marker2, marker3, marker4], loc = 'upper right')
    legend2 = ax.legend(handles = [y_lines, b_lines, dashed_lines, dotdashed_lines], loc='lower left')
    ax.add_artist(legend1)
    ax.add_artist(legend2)
    #ax.legend(handles = [b_lines], loc='lower left')



def show_spectrum_mystery(ene_in, spec_in):

    deltaCP = [-90, 0, 90, 180, -90]
    # spectrum at FD
    models = [ spec_in*osc_mod.probability_oscillation(E=ene_in,L=Baseline,dCP=_dCP,b_normal_hierarchy=True ,b_neutrino=True) for _dCP in deltaCP]

    model_labels = ["$\delta_{CP}=%i$"%_dCP for _dCP in deltaCP]
    models_colors = ['g', 'b', 'r', 'orange', 'yellow', '']

    fig_2, axs_2 = plt.subplots(1, 2, constrained_layout=True, figsize=(10, 4))
    fig_2.suptitle("Spectrum at FD")
    axs_2[0].errorbar(x_output, y_output, np.sqrt(y_output), [(x_output[1]-x_output[0])/4 for _ in x_output], linestyle=None, label="measured")
    for it, model in enumeratlse(models):
        alpha = np.max(y_output[1:])/np.max(model[1:])
        axs_2[0].plot(x_after, model*alpha, color =models_colors[it] ,label=model_labels[it])
    axs_2[0].legend(loc='best')
    axs_2[0].set_title("Neutrinos in NO")
    
    models = [  spec_in*osc_mod.probability_oscillation(E=x_after, L=Baseline, dCP=_dCP,b_normal_hierarchy=False, b_neutrino=True ) for _dCP in deltaCP]
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


    #Muonic neutrinos at ND from arXiv:2109.01304v1
    near_detector_spectrum = np.loadtxt('data/input.txt')
    ene_ND_article, spec_ND_article = [i[0] for i in near_detector_spectrum], [i[1]*1e6 for i in near_detector_spectrum]

    #Electronic neutrinos at FD from arXiv:2109.01304v1
    far_detector_spectrum = np.loadtxt('data/output.txt')
    ene_FD_article, spec_FD_article = [i[0] for i in far_detector_spectrum], [i[1] for i in far_detector_spectrum]

    list_ND_mc = osc_mod.model_flux(ene_ND_article, spec_ND_article, N_ND, 0.25)

    spec_ND_mc, ene_ND_mc = np.histogram(list_ND_mc, bins=50, range=[0.2,6])
    ene_ND_mc = ene_ND_mc[:-1]
    #print(ene_ND_mc)

    spec_nue_mc = spec_ND_mc*osc_mod.probability_oscillation(E=ene_ND_mc, L=Baseline, dCP=-500, b_neutrino=True, b_normal_hierarchy=True)
    resolution = Gaussian1DKernel(stddev=2.5)
    spec_FD_mc = convolve(spec_nue_mc, resolution, normalize_kernel=True, boundary="extend")

    model = spec_FD_mc*max(spec_FD_article)/max(spec_FD_mc)
    osc_plot.compare_data_model("Comparison data/MC expected FD", ene_FD_article, spec_FD_article, "data nu_e FD", ene_ND_mc, model, "FD nu_mu oscillated spectrum")

    show_spectrum_mystery(ene_ND_mc,spec_ND_mc)
    plt.show()

