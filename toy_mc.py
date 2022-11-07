# -*- coding: utf-8 -*-
"""
Created on 4th of Nov 2022
Thibaut Houdy



"""
import math
import numpy as np
import matplotlib.pyplot as plt
import os, sys, re
from scipy.signal import savgol_filter, find_peaks
import ROOT
import argparse
import time
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

def LoadParameters(track_file):
  
  POT = 6e20
  InitalErrorOnFlux = 0.1
  InitalErrorOnSigma = 0.2
  fluxMean = 1.92e13/1e21
  sigmaMean = 2.5e-9
  ErrorOnFlux =  InitalErrorOnFlux*fluxMean
  ErrorOnSigma = InitalErrorOnSigma*sigmaMean
  bckNear=3000
  bckFar=200
  backgroundActivated=false
  silence=false
  correlation=false
  correlationFactor=0.
  andthisvalue=5

def p(E,L):  #  nm to ne oscillation probability  [E]-GeV, [L]-km

    delta_m21=7.42*pow(10,-5) #in eV^2
    a = 1/3500  # GfNe/v2 
    # normal order parameters
    theta12_no=np.radians(33.44)
    theta13_no=np.radians(8.57)
    theta23_no=np.radians(49.2)
    deltaCP_no=np.radians(0)
    delta_m31=2.515*pow(10,-3) #in eV^2

    # inverted order parameters

    theta12_io=np.radians(33.45)
    theta13_io=np.radians(8.6)
    theta23_io=np.radians(49.5)
    deltaCP_io=np.radians(287)
    delta_m32=-2.498*pow(10,-3) #in eV^2

    
    P_NO = pow(np.sin(theta23_no),2)*pow(np.sin(2*theta13_no),2) \
              * pow(np.sin(delta_m31*L*1.267/E-a*L),2)/(pow(delta_m31*L*1.267/E-a*L,2)) \
                * pow(delta_m31*L*1.267/E,2) \
                  + np.sin(2*theta23_no)*np.sin(2*theta13_no)*np.sin(2*theta12_no)*np.sin(delta_m31*L*1.267/E-a*L)/(delta_m31*L*1.267/E-a*L)*delta_m31*L*1.267/E*np.sin(a*L)/a/L \
                    * delta_m21*L*1.267/E*np.cos(delta_m31*L*1.267/E+deltaCP_no) \
                      + pow(np.cos(theta23_no),2)*pow(np.sin(2*theta12_no),2)*pow(np.sin(a*L),2)/pow(a*L,2)*pow(delta_m21*L*1.267/E,2) 
    P_IO = pow(np.sin(theta23_io),2)*pow(np.sin(2*theta13_io),2) \
              * pow(np.sin(delta_m31*L*1.267/E-a*L),2)/(pow(delta_m31*L*1.267/E-a*L,2)) \
                * pow(delta_m31*L*1.267/E,2) \
                  + np.sin(2*theta23_io)*np.sin(2*theta13_io)*np.sin(2*theta12_io)*np.sin(delta_m31*L*1.267/E-a*L)/(delta_m31*L*1.267/E-a*L)*delta_m31*L*1.267/E*np.sin(a*L)/a/L \
                    * delta_m21*L*1.267/E*np.cos(delta_m31*L*1.267/E+deltaCP_io) \
                      + pow(np.cos(theta23_io),2)*pow(np.sin(2*theta12_io),2)*pow(np.sin(a*L),2)/pow(a*L,2)*pow(delta_m21*L*1.267/E,2) 
    P_aNO = pow(np.sin(theta23_no),2)*pow(np.sin(2*theta13_no),2) \
              * pow(np.sin(delta_m31*L*1.267/E+a*L),2)/(pow(delta_m31*L*1.267/E+a*L,2)) \
                * pow(delta_m31*L*1.267/E,2) \
                  - np.sin(2*theta23_no)*np.sin(2*theta13_no)*np.sin(2*theta12_no)*np.sin(delta_m31*L*1.267/E+a*L*1.267)/(delta_m31*L*1.267/E+a*L)*delta_m31*L*1.267/E*np.sin(a*L)/a/L \
                    * delta_m21*L*1.267/E*np.cos(delta_m31*L*1.267/E-deltaCP_no) \
                      + pow(np.cos(theta23_no),2)*pow(np.sin(2*theta12_no),2)*pow(np.sin(-a*L),2)/pow(a*L,2)*pow(delta_m21*L*1.267/E,2)                   
    P_aIO = pow(np.sin(theta23_io),2)*pow(np.sin(2*theta13_io),2) \
              * pow(np.sin(delta_m31*L*1.267/E+a*L),2)/(pow(delta_m31*L*1.267/E+a*L,2)) \
                * pow(delta_m31*L*1.267/E,2) \
                  - np.sin(2*theta23_io)*np.sin(2*theta13_io)*np.sin(2*theta12_io)*np.sin(delta_m31*L*1.267/E+a*L*1.267)/(delta_m31*L*1.267/E+a*L)*delta_m31*L*1.267/E*np.sin(-a*L)/a/L \
                    * delta_m21*L*1.267/E*np.cos(delta_m31*L*1.267/E-deltaCP_io) \
                      + pow(np.cos(theta23_io),2)*pow(np.sin(2*theta12_io),2)*pow(np.sin(-a*L),2)/pow(a*L,2)*pow(delta_m21*L*1.267/E,2) 

    return P_NO, P_IO, P_aNO, P_aIO

e = np.arange(0.1, 10., 0.001) #Neutrino energy 0-10GeV
L = 1300 # Length in Km
#print(p)
plt.figure(1)
plt.subplot(211)
plt.plot(e,p(e,L)[0],color='g', label='NO')
plt.plot(e,p(e,L)[1],color='r', label='IO') 
#plt.xlim([0, 10])
plt.ylim([0, 0.5])
plt.xlabel('Energy')
plt.ylabel('Probability')
plt.xscale('log')
plt.title("Neutrino Oscillation")
plt.legend()
plt.subplot(212)
plt.plot(e,p(e,L)[2], color='b',label='NO')
plt.plot(e,p(e,L)[3], color='orange', label='IO')
#plt.xlim([0, 10])
plt.ylim([0, 0.5])
plt.xlabel('Energy')
plt.ylabel('Probability')
plt.xscale('log')
plt.title("Antineutrino Oscillation")
plt.legend()
plt.show()

