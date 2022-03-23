# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 14:35:34 2017

Algorithms to convolute a spectrum with either a Gaussian Function,
a Lorentzian Function, or the Fouriertransformed of a Besses function
of first kind (pseuod-field modulation).
For all functions extrapolation algorithms are used if the signal
is above a given treshold at the edges of the magnetic field range (see
individual functions)

Main functions:
voigt_convolution_Gauss()
voigt_convolution_Lorentz()
pseudo_modulation()


@author: Stephan Rein, University of Freiburg, 2017
"""

import numpy as np
import math as math
from scipy import signal, special, fftpack, interpolate


def voigt_convolution_Gauss(Sys, field,spectrum):
    """voigt_convolution_Gauss(Sys, field,spectrum)

     input: Sys, field,spectrum

     output:  spectrum

     Algorithm:
     Gives back the voigt-convoluted of an initial Lorentzian Function.
     Uses fast convoultion via FFT:
     s(g)*s(l) = ifft(fft(s(g).fft(s(l)))) with . as pointwise multiplication.
     With the  Gaussian function g:
     g(B) = exp(-0.5*(B-c)^2/(lw^2))
     With the Gaussian linewidth (FWHM) lw and the center of field c.
     If the convoluted signal has not decayed yet at the maximum/minimum field
     position it is linearly extrapolated for 5 points assuming a local
     linearisation of the initial Lorentzian spectrum. Nevertheless, wrong
     intensities are probable on the edges of the spectrum.
     At the end the spectrum is transformed back to the initial range.
     The treshold for switching on the extrapolation is 1e-3 of the maximum
     of the signal at the minimum/maximum field.

     (c) Stephan Rein, 30.11.2017
     """

    # Input field information
    inputlen = len(field)
    inputlowfield = min(field)
    inputmaxfield = max(field)

    # Linear extrapolation for not decayed signals
    if spectrum[0] > max(spectrum)*1e-3:
        stepsize = (max(field)-min(field))/len(field)
        field_tmp_tmp  =  np.arange(min(field)-5*stepsize,max(field)+stepsize,stepsize)
        spl = interpolate.splrep(field,spectrum,k=1)
        field = field_tmp_tmp
        spectrum  = interpolate.splev(field,spl,ext=0)

    if spectrum[len(spectrum)-1] > max(spectrum)*1e-3:
        stepsize = (max(field)-min(field))/len(field)
        field_tmp_tmp  =  np.arange(min(field)-stepsize,max(field)+stepsize*5,stepsize)
        spl = interpolate.splrep(field,spectrum,k=1)
        field = field_tmp_tmp
        stepsize = (max(field)-min(field))/len(field)
        spectrum  = interpolate.splev(field,spl,ext=0)

    # Direct vector convolution
    npoints = len(spectrum)
    spectrum_conv = np.zeros(npoints)
    step = field[1]-field[0]
    std2 = (Sys.lw[0]/(2*np.sqrt(2*np.log(2))))/step
    window = signal.gaussian(npoints, std=std2)
    spectrum_conv = signal.fftconvolve(spectrum, window, mode='same')

    # Spline the full function and give back the user defined field range
    field_tmp = np.linspace(inputlowfield, inputmaxfield,inputlen , endpoint = True)
    spl = interpolate.splrep(field,spectrum_conv,k=3)
    spectrum_conv  = interpolate.splev(field_tmp,spl)
    spectrum = spectrum_conv

    return spectrum_conv


def voigt_convolution_Lorentz(Sys, field,spectrum):
    """ voigt_convolution_Lorentz(Sys, field,spectrum)

     input: Sys, field,spectrum

     output:  spectrum

     Algorithm:
     Gives back the voigt-convoluted of an initial Gaussian Function.
     Uses fast convoultion via FFT:
     s(g)*s(l) = ifft(fft(s(g).fft(s(l)))) with . as pointwise multiplication.
     With the Lorentzian function l:
     l(B) = (lw/math.pi)/((B-c)^2+(0.5*lw)**2)
     With the Lorentzian linewidth (FWHM) lw and the center of field c.
     If the convoluted signal has not decayed yet at the maximum/minimum field
     position it is linearly extrapolated for 5 points assuming a local
     linearisation of the initial Gaussian spectrum. Nevertheless, wrong
     intensities are probable on the edges of the spectrum.
     At the end the spectrum is transformed back to the initial range.
     The treshold for switching on the extrapolation is 1e-3 of the maximum
     of the signal at the minimum/maximum field.

     (c) Stephan Rein, 30.11.2017
     """

    # Input field information
    inputlen = len(field)
    inputlowfield = min(field)
    inputmaxfield = max(field)

    # Linear extrapolation for not decayed signals
    if spectrum[0] > max(spectrum)*1e-3:
        stepsize = (max(field)-min(field))/len(field)
        field_tmp_tmp  =  np.arange(min(field)-2*stepsize,max(field)+stepsize,stepsize)
        spl = interpolate.splrep(field,spectrum,k=1)
        field = field_tmp_tmp
        spectrum  = interpolate.splev(field,spl,ext=0)

    if spectrum[len(spectrum)-1] > max(spectrum)*1e-3:
        stepsize = (max(field)-min(field))/len(field)
        field_tmp_tmp  =  np.arange(min(field)-stepsize,max(field)+stepsize*2,stepsize)
        spl = interpolate.splrep(field,spectrum,k=1)
        field = field_tmp_tmp
        stepsize = (max(field)-min(field))/len(field)
        spectrum  = interpolate.splev(field,spl,ext=0)

    # Direct vector convolution
    middle = field[0]+(field[len(field)-1]-field[0])/2.0
    std = Sys.lw[1]
    Lorentzian = ((0.5*std)/math.pi)/((field-middle)**2+(0.5*std)**2)
    spectrum_conv = signal.fftconvolve(spectrum, Lorentzian, mode='same')

    # Spline the full function and give back the user defined field range
    field_tmp = np.linspace(inputlowfield,inputmaxfield,inputlen , endpoint = True)
    spl = interpolate.splrep(field,spectrum_conv,k=3)
    spectrum_conv  = interpolate.splev(field_tmp,spl)
    spectrum = spectrum_conv
    return spectrum_conv



def pseudo_modulation(Exp,Opt,Sys, field,spectrum) :
    """ pseudo_modulation(Exp,Opt,Sys, field,spectrum)

     input: Exp,Opt,Sys, field,spectrum

     output:  field,spectrum

     Algorithm:
     Pseudo-field-modulation is carried out in inverse field domain by
     multipliying the Fourier transformed of the signal with a Bessel function
     of first kind (bessel(1,)):
     S_mod = ifft(fft(S)*1j*bessel) with bessel = bessel(1, A*B^-1)
     with the modulation amplitude A chosen to be 1/4 of the FWHM of the
     maximum linewidth.
     This can lead to large artifacts if the signal has not yet decayed at the
     edges. Therefore the signal is first linearly extrapolated for 512 points.
     This should hold for nearly all cases. Additionally the extrapolated
     spectrum is again extrapolated with 1024 points to reach zero at
     the first/last point of the extrapolation. The pseudo-modulation is
     carried out for this extrapolated spectrum. Afterwards the modulated
     spectrum is reduced to the initial field range.
     The treshold for switching on the extrapolation is 1e-3 of the maximum
     of the signal at the minimum/maximum field.

     (c) Stephan Rein, 30.11.2017
     """

    #Define field constants for field to frequency conversion
    hbar = 6.62606957*1e-34
    beta = 9.27400097*1e-24
    Hz2mT = (hbar*1000)/beta
    MHz2mT = 1/(2.0023/(Hz2mT*1e6))
    value =0.25*max(Sys.lw)/MHz2mT #Modulation amplitude in inverse space
    np.set_printoptions(threshold=np.nan)
    # Linear extrapolation to zero if the signal has not decayed at low field
    if spectrum[0] > max(abs(spectrum))*1e-3:
        stepsize = (max(field)-min(field))/len(field)
        field_tmp_tmp  =  np.arange(min(field)-512*stepsize,max(field)+stepsize,stepsize)
        spl = interpolate.splrep(field,spectrum,k=1)
        field = field_tmp_tmp
        stepsize = (max(field)-min(field))/len(field)
        spectrum  = interpolate.splev(field_tmp_tmp,spl,ext=0)
        field_tmp1 = np.arange(min(field)-stepsize*1024,max(field)+stepsize*1024,stepsize)
        spectrum = np.append(0,spectrum)
        field_tmp2 = np.append(min(field)-stepsize*1024,field)
        spl = interpolate.splrep(field_tmp2,spectrum,k=1)
        spectrum  = interpolate.splev(field_tmp1,spl,ext=0)
        field = field_tmp1
        #print(field)

    # Linear extrapolation to zero if the signal has not decayed at high field
    if spectrum[len(spectrum)-1] > max(abs(spectrum))*1e-3:
        stepsize = (max(field)-min(field))/len(field)
        field_tmp_tmp  =  np.arange(min(field)-stepsize,max(field)+stepsize*512,stepsize)
        spl = interpolate.splrep(field,spectrum,k=1)
        field = field_tmp_tmp
        stepsize = (max(field)-min(field))/len(field)
        spectrum  = interpolate.splev(field,spl,ext=0)
        field_tmp1 = np.arange(min(field)-stepsize*1024,max(field)+stepsize*1024,stepsize)
        spectrum = np.append(spectrum,0)
        field_tmp2 = np.append(field,max(field)+stepsize*1024)
        spl = interpolate.splrep(field_tmp2,spectrum,k=1)
        spectrum  = interpolate.splev(field_tmp1,spl,ext=0)
        field = field_tmp1

    # Set up x-axis in inverse field domain
    fieldstep = (field[1]-field[0])*(2.0023/(Hz2mT*1e6))
    fourier_x = np.linspace(0, 1/fieldstep, len(spectrum), endpoint=True)
    level =value
    datafft = fftpack.fft(spectrum)
    for i in range(int(len(datafft)/2)+1,len(datafft)):
        datafft[i] = 0
    # Calculate the Bessesl function of first kind and make the convolution
    bessel_fft = special.jv(1,fourier_x*level)
    pseudomod = 1j*bessel_fft* datafft
    # Inverse Fourier transform of the convoluted signal
    signal = np.real(fftpack.ifft(pseudomod))

    # Spline the full function and give back the user defined field range
    field_tmp = np.linspace(Exp.Range[0],Exp.Range[1],Opt.nPoints, endpoint = True)
    spl = interpolate.splrep(field,signal,k=3)
    signal  = interpolate.splev(field_tmp,spl)
    field = field_tmp
    spectrum = signal

    return  field,spectrum