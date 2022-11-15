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

def cumulative(y):
    cumul, f_cumul = 0, []
    f_cumul.append(0)
    for yy in y:
        cumul += yy
        f_cumul.append(cumul)
    return f_cumul

def mc_function(x_ene, f_y, Ntrials,sigma): #generate energy spectrum for near or far detector
    new_spectrum = []
    for i in range(int(Ntrials)):
        t_rand = rand.uniform(0,1)
        #min_bin = min(np.where(t_rand<f_y)[0])
        max_bin = max(np.where(t_rand>f_y)[0])
        try : 
            true_ene1 = x_ene[max_bin]
            #true_ene2 = x_ene[min_bin]
            rand_energie = rand.gauss(true_ene1,sigma)
            #rand_energie = rand.gauss(true_ene2,0.25)        
            #rand_energie = rand.uniform(true_ene1,true_ene2)
            if rand_energie > 0:
                new_spectrum.append(rand_energie)
        except : continue
            #true_ene = x_ene[max_bin]           
    return new_spectrum

def D31(E):
  return(1.267*delta_m31*L/E)

def D21(E):
  return(1.267*delta_m21*L/E)

def D32(E):
  return(1.267*delta_m32*L/E)  

#For muon neutrino to electron neutrino oscillation probability

L = 1285 # Length in Km

#oscillation parameters from NuFit 5.1

delta_m21=7.42*pow(10,-5) #in eV^2

# normal order parameters
theta12_no =np.radians(33.44)
theta13_no =np.radians(8.57)
theta23_no =np.radians(49.2)
#deltaCP_no =np.radians(194)
delta_m31 = 2.515*pow(10,-3) #in eV^2

# inverted order parameters

theta12_io =np.radians(33.45)
theta13_io =np.radians(8.6)
theta23_io =np.radians(49.5)
#deltaCP_io =np.radians(287)
delta_m32 =-2.498*pow(10,-3) #in eV^2  

a1 = 1/3500  # GfNe/sqrt(2) for neutrinos
a2 = - 1/3500  # -GfNe/sqrt(2) for antineutrinos

#for NO
c1 = np.sin(theta23_no)**2*np.sin(2*theta13_no)**2
c2 = np.sin(2*theta23_no)*np.sin(2*theta13_no)*np.sin(2*theta12_no)
c3 = np.cos(theta23_no)**2*np.sin(2*theta12_no)**2

#for IO
c4 = np.sin(theta23_io)**2*np.sin(2*theta13_io)**2
c5 = np.sin(2*theta23_io)*np.sin(2*theta13_io)*np.sin(2*theta12_io)
c6 = np.cos(theta23_io)**2*np.sin(2*theta12_io)**2

#  calculate oscillation probability (Normal Ordering)  [E]-GeV, [L]-km, deltaCP & a
def p_NO(E,dCP,a):   
    deltaCP_no = np.radians(dCP)

    # p for NO
    P_NO = c1 * np.sin(D31(E)-a*L)**2/(D31(E)-a*L)**2 * D31(E)**2 \
                  + c2 *np.sin(D31(E)-a*L)/(D31(E)-a*L)*D31(E)*np.sin(a*L)/a/L \
                    * D21(E)*np.cos(D31(E)+deltaCP_no) + c3*D21(E)**2*np.sin(a*L)**2/(a*L)**2
    return P_NO

# Same for Inverted Ordering    
def p_IO(E,dCP,a):  
    deltaCP_io=np.radians(dCP)
    # p for IO               
    P_IO = c4 * np.sin(D32(E)-a*L)**2/(D32(E)-a*L)**2 * D32(E)**2 \
                  + c5 *np.sin(D32(E)-a*L)/(D32(E)-a*L)*D32(E)*np.sin(a*L)/a/L \
                    * D21(E)*np.cos(D32(E)+deltaCP_io) + c6*(D21(E))**2*np.sin(a*L)**2/(a*L)**2                 
    return P_IO

#simulate neutrino flux using measured energy spectrum
def model_flux(x,y,N,sigma):
    #spectrum = np.loadtxt(txtfile)
    #x, y = [i[0] for i in spectrum], [i[1] for i in spectrum]
    y_cumul = cumulative(y)
    y_cumul = y_cumul/max(y_cumul)
    flux = np.array(mc_function(x,y_cumul,N,sigma))
    return flux





if __name__ == '__main__':
    alpha1 = 0.0165   #including efficiency, flux and sigma
    alpha2 = 0.02
    alpha3 = 0.0135
    alpha4 = 0.0155
    N_ND = 1.5*1e6
    N_FD = 1e3
    near_detector_spectrum = np.loadtxt('../input.txt')
    far_detector_spectrum = np.loadtxt('../output.txt')
    x_input, y_input = [i[0] for i in near_detector_spectrum], [i[1] for i in near_detector_spectrum]
    x_output, y_output = [i[0] for i in far_detector_spectrum], [i[1] for i in far_detector_spectrum]

    y_spec, x_spec = np.histogram(model_flux(x_input,y_input,N_ND,0.3), bins=50)
    #y_s, x_s = np.histogram(model_flux(x_output,y_output,N_FD,0.3),bins = 25)

    x_after = x_spec[:-1]
    y_after = 1.9*y_spec
    #x_final = x_s[:-1]
    #y_final = y_s
    #print(y_after)
    
    # spectrum at FD
    plt.plot(x_output, y_output, label="measured")
    plt.plot(x_after, y_after*p_NO(x_after,0, a1)*alpha1, color ='g' ,label="model for d_CP=0")
    plt.plot(x_after, y_after*p_NO(x_after,90, a1)*alpha2, color ='r', label="d_CP=pi/2")
    plt.plot(x_after, y_after*p_NO(x_after,-90, a1)*alpha3, color = 'b',label="d_CP=-pi/2")
    plt.plot(x_after, y_after*p_NO(x_after,180, a1)*alpha4, color = 'orange',label="d_CP=pi")

    #plt.plot(x_input , y_input/sum(y_input), label="before")
    #plt.plot(x_after, y_after, label="after")

    #plt.plot(x_output , y_output/sum(y_output), label="before")
    #plt.plot(x_final, y_final/sum(y_s), label="after")

    plt.xlim([0.2,7])
    plt.ylim([0,300])
    plt.legend(loc='best')
    plt.title("Spectrum at FD")
    plt.show()






'''
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