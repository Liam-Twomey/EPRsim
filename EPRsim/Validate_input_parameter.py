#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:14:47 2019

@author: stephan
"""

import numpy as np
import EPRsim.Tools as tool


con = tool.physical_constants()

def validate_spin_system(Param):
    """ get_spin_system(Sys)

    in:  Sys             (one class with definitions of spin system parameters)

    out: g               (g-tensor. Given as 3x3 matrix)
         A               (vector of hyperfine tensors each given as 3x3 matrix)
         Nucs            (string of couling nuclei)
         lw              (Lorentzian linewidth given in mT as FWHM)
         lwG             (Gaussian linewwidth given in mT as FWHM)
         tcorr           (rotational correlation time given in s)
         equiv           (vector with numbers of equivalent coupling nuclei)
         Iso             (boolen. True = isotropic regime. False = fast-motion)

    Returns the information of the Sys class. Sets default parameter if they
    are not given in the Sys class.
    """
    Param.iso = True
    if isinstance(Param.g, int) or isinstance(Param.g, float):
        Param.g = np.ones(3)*Param.g
    elif len(Param.g) == 3:
        Param.g = Param.g
    else:
        Param.g = np.array([Param.g[0], Param.g[0], Param.g[1]])
    Param.g = np.asarray(Param.g)
    if Param.A is not None:
        get_hyperfines_from_Sys(Param)
        Param.A = Param.A*1e6
        Param._equiv = np.ones(Param._nA)*Param.n
        Param._nofA = int(np.sum(np.ones(Param._nA)*Param.n))
    else:
        Param.A = 1.0*np.zeros(3)
        Param._nA = 1
        Param._equiv = np.ones(Param._nA)
    if isinstance(Param.lw, float) or isinstance(Param.lw, int):
        Param.lw = [Param.lw, 0]
    if len(Param.lw) == 1:
        Param.lw.append(0.0)
    if Param.lw[1] < 0.01:
        Param.lw[1] = 0.01
    # Swap Gaussian and Lorentzian
    Param.lw = [Param.lw[1],  Param.lw[0]]
    if Param.tcorr is not None or Param.logtcorr is not None:
        Param.iso = False
    if Param.logtcorr is not None:
        Param.tcorr = np.power(10, float(Param.logtcorr))
    Param.Bfield = np.linspace(Param.Range[0], Param.Range[1], Param.Points)
    Param.mwFreq *= con.GHz2Hz
    return


def _redefine_tcorr(Param):
    if Param.tcorr is None:
        return
    if isinstance(Param.tcorr, int) or isinstance(Param.tcorr, float):
        Param._tcorriso = Param.tcorr
    elif len(Param.tcorr) == 1:
        Param._tcorriso = Param.tcorr[0]
    elif len(Param.tcorr) == 2:
        Param._tcorriso = (1/3)*(Param.tcorr[0]*2+Param.tcorr[1])
    return


def get_hyperfines_from_Sys(Param):
    """ get_hyperfines_from_Sys(Sys)

    in:  Sys             (one class with definitions of spin system parameters)

    out: A               (vector of hyperfine tensors each given as 3x3 matrix)
         Nucs            (string of couling nuclei)
         nA              (number of hyperfine tensors)

    Takes the Sys class and returns information about the hyperfine couplings.
    Is used as sub-function of get_spin_system().
    """
    if Param.A is None or Param.Nucs is None:
        Param._nA = 0
        return
    try:
        dim = len(Param.Nucs.split(','))
        Param._nA = int(dim)
    except:
        Param._nA = 1
    if Param._nA == 1:
        if isinstance(Param.A, int) or isinstance(Param.A, float):
            Param.A = np.ones(3)*Param.A
        elif len(Param.A) == 1:
            Param.A = np.ones(3)*Param.A
        elif len(Param.A) == 3:
            Param.A = Param.A
        else:
            Param.A = np.array([Param.A[0], Param.A[0], Param.A[1]])
    else:
        for i in range(0, Param._nA):
            if isinstance(Param.A[i], int) or isinstance(Param.A[i], float):
                Param.A[i] = np.ones(3)*Param.A[i]
            elif len(Param.A[i]) == 1:
                Param.A[i] = np.ones(3)*Param.A[i]
            elif len(Param.A[i]) == 3:
                Param.A[i] = Param.A[i]
            else:
                A = Param.A[i]
                A = np.array([A[0], A[0], A[1]])
                Param.A[i] = A
    Param.A = np.asarray(Param.A)
    return

def get_D_from_Sys(Param):
    if Param.D is None:
        return