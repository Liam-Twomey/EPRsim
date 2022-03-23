# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 14:35:34 2017

Different spline interpolations (linear, cubic, bilinear, bicubic)
for field or angle grids. Mostly resonance fields and corresponding
intensities are interpolated over theta and phi to obtain a finer grid without
the need of expensive quantum chemical calculations.
As every calculation runs through the function
 spline_interpolation_angle_grid() this is used to weight the intensities
 with a probality function either about theta or ph (or both).

Main functions:
spline_interpolation_angle_grid()
spline_interpolation_analytical_resInt()
spline_interpol_phi_proj()
field_interpol()

Internal subfunctions:
make_orientationdepency_phi()
make_orientationdepency_theta()

@author: Stephan Rein, University of Freiburg, 2017
"""

import numpy as np
import math as math
from scipy import interpolate










def spline_interpolation_angle_grid(Par,intensity, res):
    """ spline_interpolation_angle_grid(Par,Par, intensity, res,Par._nphi,Par._ntheta):

     input: Par,Par, intensity, res,Par._nphi,Par._ntheta

     output: intensity, res, theta, phi

     Algorithm:
     The number of theta and phi values are taken to make a rectangular
     homoneneous grid over theta and phi.
     The intensities and resonance frequencies are interpolated for the given
     number of theta and phi steps. This yields a total number of:
     n = Par._ntheta*Par._nphi
     points on the unit hemisphere (or a number of octants).
     The spline interpolation is always bivariant and bicubic for
     the resonance positions. For the intensities it is bicubic as default but
     is set to bilinear if the spectrum exceeds field range and transition
     probabilities are evaluated only on the edges of the field.

     (c) Stephan Rein, University of Freiburg, 31.10.2017
     """

    #Allocation of theta and phi grids
    res_inter = np.zeros((Par.Transdim ,Par._ntheta,Par._nphi), dtype=np.float64)
    int_inter = np.zeros((Par.Transdim ,Par._ntheta,Par._nphi), dtype=np.float64)
    k1 = np.linspace(0, math.pi/2, num=Par.nKnots, endpoint=True)
    k2 = np.linspace(0, Par.nOctants*np.pi/2, num=Par.phinKnots, endpoint=True)
    k3 = np.linspace(0, math.pi/2, num=Par._ntheta, endpoint=True)
    k4 = np.linspace(0, Par.nOctants*np.pi/2, num=Par._nphi, endpoint=True)

    #Bicubic spline interpolation (Bilinear for intensities if the spectrum
    #exceeds field range)
    if Par.field_warning == False:
        for k in range(0,Par.Transdim ):
            f1 =interpolate.RectBivariateSpline(k1,k2, intensity[k],kx = 3,ky = 3)
            int_inter[k] = f1(k3, k4)
            f2 =interpolate.RectBivariateSpline(k1,k2, res[k],kx =3,ky =3)
            res_inter[k] = f2(k3, k4)
    else:
        for k in range(0,Par.Transdim ):
            f1 =interpolate.RectBivariateSpline(k1,k2, intensity[k],kx = 1,ky =1)
            int_inter[k] = f1(k3, k4)
            f2 =interpolate.RectBivariateSpline(k1,k2, res[k],kx = 3,ky = 3)
            res_inter[k] = f2(k3, k4)
    #rename it for output
    res = res_inter
    intensity = int_inter
    theta = k3
    phi = k4
    #Make orientation dependent weighting if wished
    if  hasattr(Par, 'oritheta'):
        sigma_theta = Par.oritheta[1]
        pos_theta = Par.oritheta[0]
        intensity =  make_orientationdepency_theta(intensity,theta,pos_theta,sigma_theta)
    if  hasattr(Par, 'oriphi'):
        sigma_phi = Par.oriphi[1]
        pos_phi = Par.oriphi[0]
        intensity =  make_orientationdepency_phi(intensity,phi,pos_phi,sigma_phi)
    return  intensity, res, theta, phi












def make_orientationdepency_theta(intensity,theta,pos_theta,sigma_theta):
    """ make_orientationdepency_theta(intensity,theta,pos_theta,sigma_theta)

    input: intensity,theta,pos_theta,sigma_theta

    output: intensity

    Algorithm:
    Makes a weighting of the unbiased intensities with a Gaussian distribution
    about a given theta angle pos_theta. The standard deviation of the Gaussian
    is given by sigma_theta.
    The output is a intessity distribtion with:
    I = I*P(theta)  and P(theta) = exp(-(theta-pos_theta)^2/sigma_theta^2)
    It is calculated in a vectorized fashion for all phi values and
    all transitions at the same time.

    (c) Stephan Rein, University of Freiburg, 01.11.2017
    """
    Par._ntheta = len(theta)
    #Weightig factors auf Gaussian distributions for theta
    for k in range(0,Par._ntheta):
        intensity[:,k,:] =intensity[:,k,:]*(np.exp(-0.5*(theta[k]-pos_theta)**2/sigma_theta**2)+
        np.exp(-0.5*(theta[Par._ntheta-k-1]-pos_theta-np.pi/2)**2/sigma_theta**2)
        +np.exp(-0.5*(theta[Par._ntheta-k-1]-pos_theta+np.pi/2)**2/sigma_theta**2)+
        np.exp(-0.5*(theta[k]-pos_theta+np.pi)**2/sigma_theta**2)+
        np.exp(-0.5*(theta[k]-pos_theta-np.pi)**2/sigma_theta**2))
    return intensity











def make_orientationdepency_phi(Par, intensity,phi,pos_phi,sigma_phi):
    """ make_orientationdepency_phi(intensity,phi,pos_phi,sigma_phi)

    input: intensity,phi,pos_phi,sigma_phi

    output: intensity

    Algorithm:
    Makes a weighting of the unbiased intensities with a Gaussian distribution
    about a given phi angle pos_phi. The standard deviation of the Gaussian
    is given by sigma_phi.
    The output is a intessity distribtion with:
    I = I*P(phi)  and P(phi) = exp(-(phi-pos_phi)^2/sigma_phi^2).
    It is calculated in a vectorized fashion for all theta values and
    all transitions at the same time.

    (c) Stephan Rein, University of Freiburg, 01.11.2017
    """
    Par._nphi = len(phi)
    maxphi = max(phi)
    #Weightig factors auf Gaussian distributions for phi
    for k in range(0,Par._nphi):
        intensity[:,:,k] =intensity[:,:,k]*(np.exp(-0.5*(phi[k]-pos_phi)**2/sigma_phi**2)+
        np.exp(-0.5*(phi[k]-pos_phi+2*maxphi)**2/sigma_phi**2)+np.exp(-0.5*(phi[k]-pos_phi-2*maxphi)**2/sigma_phi**2))
    return intensity











def spline_interpolation_analytical_resInt(Par, intensity,resonance,ntransitions,phi):
     """ spline_interpolation_analytical_resInt(intensity,resonance,ntransitions,Par._ntheta,Par._nphi,phi)

     input: intensity,resonance,ntransitions,Par._ntheta,Par._nphi,phi

     output: resonance, intensity, phi_proj

     Algorithm:
     Creates spile interpolation for a triangle grid used in the projective
     spectrum reconstruction method. ALlocates first all matrices with
     zeros. Afterwards the resonances, intensities and phi values are
     interpolated where the resulting number of values depends on theta
     while the original values are coming from a square grid.
     Therefore it runs for phi between 0 and pi/2 for all transitions:
     i = 0,1,2 ....., n  and   n = 1, ..... Par._ntheta-1
     phi(i,n) = (pi/2)*i/n
     The same is done with linear splines for the instensity and resonance
     positions.

     (c) Stephan Rein, University of Freiburg, 31.10.2017
     """
     #Allocation of the zero matrices
     phimax = max(phi)
     intensity2 =  np.zeros((ntransitions,Par._nphi,Par._ntheta))
     resonance2=  np.zeros((ntransitions,Par._nphi,Par._ntheta))
     phi_proj = np.zeros((Par._nphi,Par._ntheta))
     #Linear spline interpolation
     k =0
     for i in range(0,Par._ntheta):
         phi_proj[i,0:i+1] = np.linspace(0, phimax, i+1, endpoint=True)
         k = np.append(k,phi_proj[i,0:i+1])
     for tr in range(0,ntransitions):
         for i in range(0,Par._ntheta):

           spl = interpolate.splrep(phi, intensity[tr,i,:],k=1)
           intensity2[tr,i,0:i+1] = interpolate.splev(phi_proj[i,0:i+1], spl)
           spl2 = interpolate.splrep(phi, resonance[tr,i,:],k=1)
           resonance2[tr,i,0:i+1] = interpolate.splev(phi_proj[i,0:i+1], spl2)
     #rename it for output
     resonance =  resonance2
     intensity= intensity2
     return resonance, intensity, phi_proj










def spline_interpol_phi_proj(Par, phimax):
     """ spline_interpol_phi_proj(Par._ntheta, phimax)

     input: Par._ntheta, phimax

     output: phi_proj

     Algorithm:
     Creates spile interpolation for a triangle grid.
     Therefore it runs for phi between 0 and pi/2:
     i = 0,1,2 ....., n  and   n = 1, ..... Par._ntheta-1
     phi(i,n) = (pi/2)*i/n

     (c) Stephan Rein, University of Freiburg, 31.10.2017
     """
     #Allocation
     phi_proj = np.zeros((Par._ntheta,Par._ntheta))
     for i in range(0,Par._ntheta):
         phi_proj[i,0:i+1] = np.linspace(0, phimax, i+1, endpoint=True)
     return phi_proj










def field_interpol(field,signal,Par):
    """ field_interpol(field,signal,Par,Par)

    input: field,signal,Par,Par

    output: field,signal

    Algorithm:
    Takes aribtray spectrum/field functions and returns the
    field and corresponding spectrum for the given user range and
    the user given number of points.

    (c) Stephan Rein, University of Freiburg, 31.10.2017
    """
    #Allocation
    k1 = np.linspace(Par.Range[0], Par.Range[1], Par.Points, endpoint=True)
    spl = interpolate.splrep(field, signal,k=3)
    signal = interpolate.splev(k1, spl)
    field = k1
    return field, signal
