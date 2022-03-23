# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 14:35:34 2017

Algorithms to create either Lorenzian or Gaussian spectrum from a large
interpolated number of orientations. This is the alternative to the
projective method.

Main functions:
create_Gaussian()
create_Lorentzian()

Internal subfunctions:
field_grid_segmentation()

External subfunctions from my own library:
voigt_convolution_Gauss()
voigt_convolution_Lorentz()

@author: Stephan Rein, University of Freiburg, 2017
"""

import numpy as np
import math as math
from EPRsim.Convolutions import (voigt_convolution_Gauss,
                                 voigt_convolution_Lorentz)
try:
    from numba import float64, jit
    Numba = 1
except ImportError:
    Numba = 0
from scipy import stats


def jit_dec():
    """
    Creates a decorator for the Gauß/Lorentz functions. If Numba
    is available, the decorator is a Numba jit (just-in-time compilation),
    otherwise the identity decorator is called (this decorator does nothing).
    The jit decorator is called with defined static types to improve the
    performance. Cashing is enabled to avoid avoid compilation each time
    the Python program is invoked.
    """
    if Numba == 1:
        return(jit(float64[:](float64[:], float64[:], float64[:], float64,
                   float64), nopython=True, cache=True))
    else:
        return dec_identity


def dec_identity(ob):
    """
    Identity decorator. If Numba is not available, this decorator is used. The
    decorator does nothing.
    """
    return ob


def create_Gaussian(Par, resonance,intensity):
    """ create_Gaussian(Par,Par,Par,resonance,intensity,Par._nphi,Par._ntheta)

     input: Par,Par,Par,resonance,intensity,Par._nphi,Par._ntheta

     output: -field
             -signal

     Algorithm:
     Takes the interpolated input resonance fields and corresponding
     intensities and runs in a loop over every orientation summing up the
     Gaussian sub-spectra.
     s(B) = I*exp(-0.5*(B-r)^2/si^2)   with si = FWHM/(2*(2*np.log(2))^(1/2))
     This is done in a vectorized fashion iterating
     over phi for each transition and evaluating all theta values in one
     matrix operation. This matrix becomes large if there are a lot of
     theta steps.
     Then the algorithm uses a local Gaussian approach where
     the Gaussian functions are only evaluation in the range between -1.5 FWHM
     to +1.5 FWHM aorund the center. This time over the field is iterated and
     only the Gaussians are evaluated which have there center in the above
     described local region (-3 to + 3 FWHM).
     The classical sin(theta) weighting is used as the same number of phi
     Knots are evaluated for each theta.
     At the end a convolution with a Lorentzian function is carried out if
     wished.

     (c) Stephan Rein, 29.11.2017
     """
    field, signal, w, w2, intensity = preparation(Par, intensity)
    lw = Par.lw[0]/(2*np.sqrt(2*np.log(2)))
    lw_sq = lw*lw
    FWHM = Par.lw[0]
    res = resonance.flatten(order='C')
    ints = intensity.flatten(order='C')
    length = len(res)
    if length > 4096 and length < 8192:
        ints, res, bin_num = stats.binned_statistic(res, ints, statistic='sum',
                                                    bins=4096)
    elif length > 8192:
        ints, res, bin_num = stats.binned_statistic(res, ints, statistic='sum',
                                                    bins=8192)

    signal = eval_Gauss(res, ints, field, FWHM, lw_sq)
    if Par.lw[1] >= 0.005:
        signal = voigt_convolution_Lorentz(Par, field, signal)
    return field, signal


@jit_dec()
def eval_Gauss(res, inte, field, FWHM, lw_sq):
    signal = np.zeros(len(field))
    # Main Gaussian Loop
    minres = np.min(res)-3.0*FWHM
    maxres = np.max(res)+3.0*FWHM
    for i in range(0, len(field)):
        if field[i] > maxres or field[i] < minres:
            pass
        else:
            abc = field[i]-res
            ind = np.where(np.abs(abc) < 3*FWHM)[0]
            b = inte[ind]
            a = abc[ind]
            signal[i] = np.dot(b, np.exp((-0.5*(a**2))/(lw_sq)))
    return signal


def create_Lorentzian(Par, resonance, intensity):
    """create_Lorentzian(Par,Par,Par,resonance,intensity,Par._nphi,Par._ntheta)

     input: Par,Par,Par,resonance,intensity,Par._nphi,Par._ntheta

     output: -field
             -signal

     Algorithm:
     Takes the interpolated input resonance fields r and corresponding
     intensities I and runs in a loop over every orientation summing up the
     Lorentzian sub-spectra.
     s(B) = I/((B-r)^2+si^2))   with si = 0.5*FWHM
     This is done in a vectorized fashion iterating
     over phi for each transition and evaluating all theta values in one
     matrix operation. This matrix becomes large if there are a lot of
     theta steps.
     Then the algorithm uses a local Lorentzian approach where
     the Lorentzian functions are only evaluation in the range between -50 FWHM
     to +50 FWHM aorund the center. This time over the field is iterated and
     only the Lorentzians are evaluated which have there center in the above
     described local region (-50 to + 50 FWHM).
     The classical sin(theta) weighting is used as the same number of phi
     Knots are evaluated for each theta.
     At the end a convolution with a Gaussian function is carried out if
     wished.

     (c) Stephan Rein, 29.11.2017
     """
    # Allocations
    field, signal, w, w2, intensity = preparation(Par, intensity)
    lw_sq = Par.lw[1]/2.0
    FWHM = Par.lw[1]
    fac = lw_sq
    lw_sq_2 = lw_sq**2
    intensity = intensity*fac
    res = resonance.flatten(order='C')
    ints = intensity.flatten(order='C')
    length = len(res)
    if length > 4096 and length < 8192:
        ints, res, bin_num = stats.binned_statistic(res, ints, statistic='sum',
                                                    bins=4096)
    elif length > 8192:
        ints, res, bin_num = stats.binned_statistic(res, ints, statistic='sum',
                                                    bins=8192)

    signal = eval_Lorentz(res, ints, field, FWHM, lw_sq_2)
    if Par.lw[0] >= 0.005:
        signal = voigt_convolution_Gauss(Par, field, signal)
    return field, signal


@jit_dec()
def eval_Lorentz(res, inte, field, FWHM, lw_sq):
    signal = np.zeros(len(field))
    # Main Gaussian Loop
    minres = np.min(res)-30*FWHM
    maxres = np.max(res)+30*FWHM
    for i in range(0, len(field)):
        if field[i] > maxres or field[i] < minres:
            pass
        else:
            abc = field[i]-res
            ind = np.where(np.abs(abc) < 30*FWHM)[0]
            b = inte[ind]
            a = abc[ind]
            signal[i] = np.sum(np.divide(b, ((a**2)+lw_sq)))
    return signal


def preparation(Par, intensity):
    """
    Preparation for evaluating the Gaussian and Lorentzian functions.

    """
    # Allocations
    field = np.linspace(Par.Range[0], Par.Range[1], Par.Points, endpoint=True)
    field = field_grid_segmentation(field, Par)
    signal = np.zeros(len(field))
    w = np.ones(Par._ntheta)
    w2 = np.ones(Par._nphi)
    w[Par._ntheta-1] = 0.5
    w2[Par._nphi-1] = 0.5
    w2[0] = 0.5
    # Weightig factors for the border of octant and sin(theta) weighting
    for k in range(0, Par._ntheta):
        theta = (math.pi/2)*((k)/(Par._ntheta-1))
        intensity[:, k, :] = (intensity[:, k, :] *
                              np.asscalar(w[k])*math.sin(theta))
    for k in range(0, Par._nphi):
        intensity[:, :, k] = intensity[:, :, k]*w2[k]
    return field, signal, w, w2, intensity


def field_grid_segmentation(field, Par):
    """
     Algorithm:
     Checks if all field points are really necessary to evaluate the Gaussian
     or Lorentzian functions. If the stepsize between 2 steps is smaller
     as 4 times the FWHM it rounds:
     n = round((0.2*FWHM)/(stepsize))
     to the integer n. The magnetic field Points are than divided by n and
     rounded to an integer. The result is a magnetic field with less points.
     (E.g. n = 2 and Par.Points = 1024 leads to Fieldpoints  = 512).

     (c) Stephan Rein, 29.11.2017
     """
    Points = len(field)
    stepsize = (max(field)-min(field))/Points
    linwidth = max(Par.lw)
    if Points > 256 and 0.2*linwidth > stepsize:
        n = int(round((0.2*linwidth)/(stepsize)))
        if n > 5:
            n = 5
        inter_dim = int(round(Points/n))
    else:
        inter_dim = Points
    field = np.linspace(min(field), max(field), inter_dim, endpoint=True)
    return field
