# -*- coding: utf-8 -*-
"""
Created on 4th of Nov 2022
Thibaut Houdy



"""
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
