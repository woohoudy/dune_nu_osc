# -*- coding: utf-8 -*-
"""
Created on 4th of Nov 2022
Thibaut Houdy

Last edited 24th of Nov 2022
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


def show_spectrum_Laura(Energy=2.5, l_baseline=1000):

    #deltaCP = [-90, 0, 90, 180, -90]
    deltaCP = np.linspace(-180,180,360)
    deltaCP_spec = [89, 179, 269, 359]
    osc_par.theta23 = np.radians(49.2)

    a = [ osc_mod.probability_oscillation(E=Energy, L=l_baseline, dCP=_dCP,b_normal_hierarchy=True , b_neutrino=True) for _dCP in deltaCP]
    b = [ osc_mod.probability_oscillation(E=Energy, L=l_baseline, dCP=_dCP,b_normal_hierarchy=False ,b_neutrino=True ) for _dCP in deltaCP]
    c = [ osc_mod.probability_oscillation(E=Energy, L=l_baseline, dCP=_dCP,b_normal_hierarchy=True , b_neutrino=False) for _dCP in deltaCP]        
    d = [ osc_mod.probability_oscillation(E=Energy, L=l_baseline, dCP=_dCP,b_normal_hierarchy=False , b_neutrino=False) for _dCP in deltaCP] 
    osc_par.theta23 = np.radians(42)
    a2 = [ osc_mod.probability_oscillation(E=Energy, L=l_baseline, dCP=_dCP,b_normal_hierarchy=True , b_neutrino=True) for _dCP in deltaCP]
    b2 = [ osc_mod.probability_oscillation(E=Energy, L=l_baseline, dCP=_dCP,b_normal_hierarchy=False ,b_neutrino=True ) for _dCP in deltaCP]
    c2 = [ osc_mod.probability_oscillation(E=Energy, L=l_baseline, dCP=_dCP,b_normal_hierarchy=True , b_neutrino=False) for _dCP in deltaCP]        
    d2 = [ osc_mod.probability_oscillation(E=Energy, L=l_baseline, dCP=_dCP,b_normal_hierarchy=False , b_neutrino=False) for _dCP in deltaCP] 
    

    markers = ['P','o','X','s','P']
    linestyles = ['--','--','-.','-.']
    colors = ['y','black','y','black',]
    labels = ['NO','IO']
    fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=(10, 8))
    
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
    ax.set_title("L = %i km, E = %.1f GeV"%(l_baseline, Energy), size = 13)
    ax.legend()
    plt.xlim(0.005,0.08)
    plt.ylim(0.005,0.08)



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
    plt.show()


def show_spectrum_mystery():

    #Muonic neutrinos at ND from arXiv:2109.01304v1
    near_detector_spectrum = np.loadtxt('data/input.txt')
    ene_ND_article, spec_ND_article = [i[0] for i in near_detector_spectrum], [i[1]*1e6 for i in near_detector_spectrum]

    #Electronic neutrinos at FD from arXiv:2109.01304v1
    far_detector_spectrum = np.loadtxt('data/output.txt')
    ene_FD_article, spec_FD_article = [i[0] for i in far_detector_spectrum], [i[1] for i in far_detector_spectrum]
 


    list_ND_mc = osc_mod.model_flux(ene_ND_article, spec_ND_article, N_ND, 0.25)

    spec_ND_mc, ene_ND_mc = np.histogram(list_ND_mc, bins=50, range=[0.2,6])
    ene_ND_mc = ene_ND_mc[:-1]

    deltaCP = [-90, 0, 90, 180, -90]
    # spectrum at FD
    models_NO, models_IO, model_labels = [],[], []
    resolution = Gaussian1DKernel(stddev=2.5)
    for _dCP in deltaCP:
        spec_nue_mc_NO = spec_ND_mc*osc_mod.probability_oscillation( E=ene_ND_mc, L=Baseline,dCP=_dCP, b_normal_hierarchy=True, b_neutrino=True)
        models_NO.append(convolve(spec_nue_mc_NO, resolution, normalize_kernel=True, boundary="extend"))
        spec_nue_mc_IO = spec_ND_mc*osc_mod.probability_oscillation( E=ene_ND_mc, L=Baseline,dCP=_dCP, b_normal_hierarchy=False, b_neutrino=True)
        models_IO.append(convolve(spec_nue_mc_IO, resolution, normalize_kernel=True, boundary="extend"))
        model_labels.append("$\delta_{CP}=%i$"%_dCP)
    
    models_colors = ['green', 'blue', 'red', 'orange', 'yellow', 'grey']
    fig_2, axs_2 = plt.subplots(1, 2, constrained_layout=True, figsize=(10, 4))
    fig_2.suptitle("Spectrum at FD")
    axs_2[0].errorbar(ene_FD_article, spec_FD_article, np.sqrt(spec_FD_article), [(ene_FD_article[1]-ene_FD_article[0])/4 for _ in ene_FD_article], color='red',marker='.', alpha=0.8, linestyle='None', label="data")
    axs_2[1].errorbar(ene_FD_article, spec_FD_article, np.sqrt(spec_FD_article),[(ene_FD_article[1]-ene_FD_article[0])/4 for _ in ene_FD_article], color='red', marker='.', alpha=0.8,linestyle='None', label="data")
    axs_2[0].bar(ene_FD_article, spec_FD_article, width=0.25, color='gray', alpha=0.5, edgecolor='black', linewidth=1)
    axs_2[1].bar(ene_FD_article, spec_FD_article, width=0.25, color='gray', alpha=0.5, edgecolor='black', linewidth=1)
    
    for it, model in enumerate(models_NO):
        alpha_NO = np.max(spec_FD_article[1:])/np.max(model[1:])
        axs_2[0].plot(ene_ND_mc, model*alpha_NO, color =models_colors[it] ,label=model_labels[it])
        alpha_IO = np.max(spec_FD_article[1:])/np.max(models_IO[it][1:])
        axs_2[1].plot(ene_ND_mc, models_IO[it]*alpha_IO, color =models_colors[it] ,label=model_labels[it])
    axs_2[0].legend(loc='best')
    axs_2[0].set_title("Neutrinos in NO")
    axs_2[0].set_xlabel('Energy [GeV]', size = 12)
    axs_2[0].set_ylabel('# events/ GeV', size = 12)
    axs_2[1].set_xlabel('Energy [GeV]', size = 12)
    axs_2[1].set_ylabel('# events/ GeV', size = 12)
    axs_2[1].legend(loc='best')
    axs_2[1].set_title("Neutrinos in IO")
    plt.show()


def plot_input_mc():
    #Muonic neutrinos at ND from arXiv:2109.01304v1
    near_detector_spectrum = np.loadtxt('data/input.txt')
    ene_ND_article, spec_ND_article = [i[0] for i in near_detector_spectrum], [i[1]*1e6 for i in near_detector_spectrum]

    #Electronic neutrinos at FD from arXiv:2109.01304v1
    far_detector_spectrum = np.loadtxt('data/output.txt')
    ene_FD_article, spec_FD_article = [i[0] for i in far_detector_spectrum], [i[1] for i in far_detector_spectrum]

    list_ND_mc = osc_mod.model_flux(ene_ND_article, spec_ND_article, N_ND, 0.25)

    spec_ND_mc, ene_ND_mc = np.histogram(list_ND_mc, bins=100, range=[0.2,6])
    spec_ND_mc= spec_ND_mc*max(spec_ND_article)/max(spec_ND_mc)
    ene_ND_mc = ene_ND_mc[:-1]
    #print(ene_ND_mc)
    osc_plot.check_input_with_mc(ene_ND_article, spec_ND_article, ene_ND_mc, spec_ND_mc )
    plt.show()

def plot_FD_data():


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
    #show_spectrum_mystery(ene_ND_mc,spec_ND_mc)
    plt.show()


def dont_need_at_the_moment():
        # spectrum at FD
    models_NO, models_IO, model_labels = [],[], []
    resolution = Gaussian1DKernel(stddev=2.5)
    deltaCP = np.linspace(-180,180, 360)
    for _dCP in deltaCP:
        spec_nue_mc_NO = spec_ND_mc*osc_mod.probability_oscillation( E=ene_ND_mc, L=Baseline,dCP=_dCP, b_normal_hierarchy=True, b_neutrino=True)
        models_NO.append(convolve(spec_nue_mc_NO, resolution, normalize_kernel=True, boundary="extend"))
        spec_nue_mc_IO = spec_ND_mc*osc_mod.probability_oscillation( E=ene_ND_mc, L=Baseline,dCP=_dCP, b_normal_hierarchy=False, b_neutrino=True)
        models_IO.append(convolve(spec_nue_mc_IO, resolution, normalize_kernel=True, boundary="extend"))
        model_labels.append("$\delta_{CP}=%i$"%_dCP)

    for it, model in enumerate(models_NO):
        alpha_NO = np.max(spec_FD_article[1:])/np.max(model[1:])
        axs_2[0].plot(ene_ND_mc, model*alpha_NO, color =models_colors[it] ,label=model_labels[it])
        alpha_IO = np.max(spec_FD_article[1:])/np.max(models_IO[it][1:])
        axs_2[1].plot(ene_ND_mc, models_IO[it]*alpha_IO, color =models_colors[it] ,label=model_labels[it])

    fig_3[1].suptitle("chi2")
    axs_3[1].hist(chi2s, bins = bins, color='b')
    #axs_3.plot(bins, scipy.stats.chi2.pdf(bins,))
    axs_3[1].legend(loc='best')
    axs_3[1].legend("chi2_neutrino_NO")

def chi2(data, model):
    chi2 = 0 
    for i in range(len(model)):   
        chi2 +=  ((model[i] - data[i])/np.sqrt(data[i]))**2

    return chi2 

def minimize_chi2(dCPs):
    chi2s = []
    for dCP in dCPs :  #dcp, theta23, hierarchy
        model = neutrino(param)
        chi2s.append(chi2(model, data))
    return params[min(chi2s)]


def show_chi2():

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
    t = interpolate.splrep(ene_FD_article, spec_FD_article, s=0)
    new_ene_FD_article = np.linspace(min(ene_FD_article),max(ene_FD_article),50)
    new_spec_FD_article = interpolate.splev(new_ene_FD_article,t,der =0)

    deltaCP = [-90, 0, 90, 180, -90]
    for dCP in deltaCP:
        chi2s =[]
        spec_nue_mc = spec_ND_mc*osc_mod.probability_oscillation(E=ene_ND_mc, L=Baseline, dCP=dCP, b_neutrino=True, b_normal_hierarchy=True)
        resolution = Gaussian1DKernel(stddev=2.5)
        spec_FD_mc = convolve(spec_nue_mc, resolution, normalize_kernel=True, boundary="extend")
        model = spec_FD_mc*max(spec_FD_article)/max(spec_FD_mc)

        chi2s.append(chi2(new_spec_FD_article,model))

    plt.plot(dCP, chi2s, linestyle='blue')
    plt.show()
#    plt.plot(new_ene_FD_article,new_spec_FD_article, linestyle='none', marker='o')
#   plt.plot(ene_FD_article,spec_FD_article,linestyle='none', marker='+')
    

    
#    minimize_chi2(deltaCP)

#    plt.show()



if __name__ == '__main__':
    
    #show_spectrum_Laura(Energy=2.5, l_baseline=1285)
    #plot_input_mc()
    #plot_FD_data()
    #show_spectrum_mystery()
    #show_chi2()
   
    


