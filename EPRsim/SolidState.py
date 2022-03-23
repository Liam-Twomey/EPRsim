#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 09:58:57 2017

@author: Stephan Rein

Copyright (c) 2017, Stephan Rein, M.Sc., University of Freiburg

2017-06-10

"""

from copy import deepcopy
import numpy as np
from EPRsim.Presettings import convert_user_input_and_Set_up_defaults
from EPRsim.resfield_full import stick_spectrum_calculation
from EPRsim.Hamiltonian_Eig import ZFS_Hamiltonian, HF_Eig
from EPRsim.spectral_processing import create_conv_spectrum
from EPRsim.Convolutions import pseudo_modulation
import EPRsim.Tools as tool
import time as time

def solid_state_kernel(Par1, SimPar1):
    Par = deepcopy(Par1)
    SimPar = deepcopy(SimPar1)
    #st= time.time()
    Par, SimPar = convert_user_input_and_Set_up_defaults(Par, SimPar)
    if Par.warning == 1:
        return (np.linspace(Par.Range[0], Par.Range[1], int(Par.Points)),
                np.zeros(int(Par.Points)), Par.warning)
    ZFS_Hamiltonian(Par)
    HF_Eig(Par)
    #print(time.time()-st)
    #st= time.time()
    intensity, resonance, Par = stick_spectrum_calculation(Par)
    #print(time.time()-st)
    #st= time.time()
    magnetic_field, spectrum = create_conv_spectrum(Par, intensity, resonance)
    #print(time.time()-st)
    #st= time.time()
    """Do pseudo-field modulation if necessary"""
    if Par.Harmonic == 1:
        spectrum = tool.pseudo_field_modulation(0.001, magnetic_field, spectrum)
        #magnetic_field, spectrum = pseudo_modulation(Exp,Opt,Sys, magnetic_field,spectrum)
    if Par.verbosity:
        print_info(Par)
    return magnetic_field, spectrum, Par.warning


def print_info(Par):
    print("\n*******************************************")
    print("*************RUN SOLID STATE*************")
    print("*******************************************\n")
    print("Hamiltonian point group: " + Par.Point_Group)
    print("Octants for powder average: " + str(Par.nOctants))
    print("Interpolated theta/phi: " + str(Par._ntheta) + "/" + str(Par._nphi))
    print("Hilbert space dimension: "+str(Par.Hilbert_dim))
    print("Electron spin dimension: "+str(Par.e_dimension))
    print("Nuclear spin dimension: " + str(Par.dim_nuc_tot))
    print("Possible transitions: "+str(Par.all_trans_dim))
    print("\"Allowed\" EPR transitions: " + str(Par.allowed_EPR_trans))
    print("Treated transitions: " + str(Par.Transdim))
    return
