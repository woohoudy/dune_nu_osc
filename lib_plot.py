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
import lib_oscil as osc_mod
import params as osc_par


def check_input_with_mc(ene_article, spec_article, ene_mc, spec_mc):
    # compare input and simulated energies
    fig, axs = plt.subplots(1, 1, constrained_layout=True, figsize=(10, 4))
    axs.plot(ene_article, spec_article,label='article')
    axs.plot(ene_mc, spec_mc,label='simulated')
    axs.set_xlabel("Energy [GeV]")
    axs.set_ylabel("#/GeV")
    fig.suptitle("Generating our own MC input data from article")    
    axs.legend(loc='best')

def compare_data_model(fig_title, ene_data, spec_data, label_data, ene_model, spec_model, label_model):

    fig, axs = plt.subplots(1, 1, constrained_layout=True, figsize=(10, 4))
    axs.errorbar(ene_data, spec_data, np.sqrt(spec_data), [(ene_data[1]-ene_data[0])/4 for _ in ene_data], linestyle = 'None', label=label_data)
    axs.plot(ene_model, np.array(spec_model), label=label_model)
    axs.set_xlabel("Energy [GeV]")
    axs.set_ylabel("#/GeV")
    fig.suptitle(fig_title)    
    axs.legend(loc='best')

def compare_deltaCP():
 	blablabla