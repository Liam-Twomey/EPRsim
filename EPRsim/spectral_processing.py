# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 14:35:34 2017

Algorithms to create either Lorenzian or Gaussian spectrum from a large
interpolated number of orientations.

@author: Stephan Rein, University of Freiburg, 2017
"""

import numpy as np
from EPRsim.Interpolation_lib import (spline_interpolation_angle_grid,
                                      field_interpol)
from EPRsim.Direct_conversion_to_Field import (create_Lorentzian,
                                               create_Gaussian)


def create_conv_spectrum(Par, intensity, resonance):
    Par.Gaussian = True
    Par._nphi = None
    Par._ntheta = None
    determine_which_broadening(Par)
    interpolation_number(Par, intensity, resonance, True)
    inten, res, theta, phi = spline_interpolation_angle_grid(Par, intensity,
                                                             resonance)
    if Par.Gaussian:
        field, signal = create_Gaussian(Par, res, inten)
    else:
        field, signal = create_Lorentzian(Par, res, inten)
    field, signal = field_interpol(field, signal, Par)
    return field, signal


def determine_which_broadening(Par):
    if Par.lw[1] > 2*Par.lw[0]:
        Par.Gaussian = False
    return


def interpolation_number(Par, intensity, resonance, Intref=False):
    thetamax_ges = 0
    phimax_ges = 0
    # Determine angle grid for spline projections
    thetamax_ges = 0
    phimax_ges = 0
    lenthetha = len(resonance[0, :, 0])
    nTransitions = len(resonance[:, 0, 0])
    for i in range(0, nTransitions):
        phimax = np.max(np.abs(np.gradient(resonance[i, lenthetha-1, :])))
        thetamax = np.max(np.abs(np.gradient(resonance[i, :, 0])))
        if thetamax_ges < thetamax:
            thetamax_ges = thetamax
        thetamax = np.max(np.abs(np.gradient(resonance[i, :, lenthetha-1])))
        if thetamax_ges < thetamax:
            thetamax_ges = thetamax
        if phimax_ges < phimax:
            phimax_ges = phimax
    thetamax_ges = 2*Par.nKnots*thetamax_ges
    phimax_ges = 2*Par.nKnots*phimax_ges
    Range = np.amax(resonance)-np.amin(resonance)
    maxlw = np.sqrt(np.max(Par.lw))
    numges = np.sqrt(thetamax_ges**2+phimax_ges**2)
    num_theta = 2*(thetamax_ges*numges)/maxlw
    num_phi = 2*(numges*phimax_ges)/maxlw
    num_phi = int(round(num_phi/Range))
    num_theta = int(round(num_theta/Range))
    # Minimal number of points
    if Par.Point_Group == "O3":
        num_phi = 4
        num_theta = 4
    elif Par.Point_Group == "Dhinfty":
        if num_phi > num_theta:
            num_theta = 4
        else:
            num_phi = 4
    else:
        if num_phi < 15:
            num_phi = 15
        if num_theta < 15:
            num_theta = 15
    if Par.Point_Group == "Dhinfty":
        maxnangle = 700
    else:
        maxnangle = 500
        if num_phi > maxnangle:
            num_phi = maxnangle
        if num_theta > maxnangle:
            num_theta = maxnangle
    # scale the number of phi due to the number of octants
    num_phi = num_phi*Par.nOctants
    if hasattr(Par, 'Interpolative_Refinement') and Intref:
        num_theta = int(round(Par.Interpolative_Refinement*num_theta))
        num_phi = int(round(Par.Interpolative_Refinement*num_phi))
    if num_phi < 4:
        num_phi = 4
    if num_theta < 4:
        num_theta = 4
    Par._nphi = num_phi
    Par._ntheta = num_theta
    return
