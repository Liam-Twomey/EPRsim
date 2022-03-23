#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 11:03:37 2019

@author: stephan
"""
# Load all external libraries
import numpy as np
from scipy import interpolate
import EPRsim.Validate_input_parameter as Val
import EPRsim.Tools as tool
global Numba
try:
    from numba import jit, float64, int64
    Numba = 1
except ImportError:
    Numba = 0


# *****************************************************************************
# Physical constants and unit conversion factors + global default settings
# *****************************************************************************

# Load physical constans
con = tool.physical_constants()
# Define specific default values for the fast motion program
Breit_Rabi_thresh = 1e-7         # Threshold for the fixpoint iteration in T
defPoints = 1024                 # Default number of points


# *****************************************************************************
# Function decorators (needed for Numba just-in-time compilation)
# *****************************************************************************


def dec_fm_lw_kernel():
    """
    Creates a decorator for the function fast_motion_lw_kernel(). If Numba
    is available the decorator is a Numba jit (just-in-time compilation),
    otherwise the identity decorator is called (this decorator does nothing).
    The jit decorator is called with defined static types to improve the
    performance. Cashing is enabled to avoid avoid compilation times each time
    you invoke a Python program.
    """
    if Numba == 1:
        return (jit(float64[:](float64[:, :, :], float64[:, :, :], float64,
                    float64, float64, float64, float64, float64,
                    float64[:, :], float64[:], int64), nopython=True,
                    cache=True))
    else:
        return dec_identity


def dec_Lorentzian():
    """
    Creates a decorator for the create_Lorentzian(). If Numba
    is available the decorator is a Numba jit (just-in-time compilation),
    otherwise the identity decorator is called (this decorator does nothing).
    The jit decorator is called with defined static types to improve the
    performance. Cashing is enabled to avoid avoid compilation times each time
    you invoke a Python program.
    """
    if Numba == 1:
        return (jit(float64[:](float64[:], float64, float64, int64),
                    nopython=True, cache=True))
    else:
        return dec_identity


def dec_identity(ob):
    """
    Identity decorator. If Numba is not available this decorator is used. The
    decorator does nothing.
    """
    return ob


def fast_motion_kernel(Param, SimPar):
    """ fast_motion_kernel(Sys,Exp)

    in:  Sys           (one class with definitions of spin system parameters)
         Exp           (one class with definitions of experimental parameters)

    out: Bfield        (magnetic field vector)
         Int           (intensity of one cw-EPR spectrum)
         Harmonic      (boolean variable. 1 = first derivative. 0 = absorptive)

    Kernel function for the calculation of one cw-EPR spectrum in the
    fast motion regime or the isotropic regime. The calculation is
    basically the same for isototropic and fast-motion spectra except of the
    linewidt. Further differences between fast-motion and isotropic
    calculations were set up before calling this function (expansion of
    hyperfine vectors in the redefine_nuclear_coupling() function).
    The steps in the calculation are:

        1. Set up the spin system and the experimental parameters
        2. Calculate all resonance fields
        3. Binning of resonance fields to the magnetic field vector
        4. Create the spectrum by setting up the Lorentzian linewidth(s)
        5. Convolute it with a Gaussian lineshape function if required
        6. Normalize the signal to the integral of the absorptive signal
        7. Return the spectrum
    """
    res, indices, DeltaA, Deltag, giso = calculate_resfields(SimPar.A,
                                                             SimPar.g,
                                                             SimPar._Nucsvec,
                                                             SimPar._equiv,
                                                             SimPar.mwFreq,
                                                             SimPar._Iequiv,
                                                             SimPar._g_n)
    if not hasattr(SimPar, 'Bfield'):
        SimPar.Bfield = do_default_range(res, SimPar.lw[0], SimPar.lw[1],
                                         SimPar._tcorriso)
    resonances = res
    if SimPar.iso:
        Int = create_isotropic_spectrum(SimPar.Bfield, resonances, indices,
                                        SimPar._equiv,
                                        SimPar._I, Param.Harmonic,
                                        SimPar.lw[0],
                                        SimPar.lw[1])
    else:
        Int = create_fastmotion_spectrum(SimPar.Bfield, resonances, indices,
                                         SimPar._I,
                                         Deltag, DeltaA, SimPar._tcorriso,
                                         SimPar.lw[0], SimPar.lw[1],
                                         giso, SimPar.mwFreq, Param.Harmonic)
    Int = tool.normalize2area(Int, Param.Harmonic)

    if Param.verbosity:
        print_Info(SimPar)
    warning = 0
    if Param.S > 0.5:
        warning = 2
    return SimPar.Bfield, Int, warning


# *****************************************************************************
# Functions for the fast-motion and isotropic calculations
# *****************************************************************************


def create_fastmotion_spectrum(Bfield, resonances, indices, I, Deltag, DeltaA,
                               tcorr, lw, lwG, giso, mfreq, Harmonic=1):
    """create_fastmotion_spectrum(Bfield,resonance,resfield,indices,equiv,I,
                               Deltag,DeltaA,tcorr,lw,lwG,giso,mfreq,Harmonic)

    in:  Bfield        (magnetic field vector in mT)
         resonances    (vector with resonance indices in Bfield)
         indices       (list with all index vectors of the coupled system)
         I             (vector with nuclear spin quantum number for all nuclei)
         Deltag        (g-giso*1 tensor)
         DeltaA        (A-Aiso*1 tensor)
         tcorr         (rotational correlation time in s)
         lw            (Lorentzian linewidth given in mT as FWHM)
         lwG           (Gaussian linewwidth given in mT as FWHM)
         giso          (isotropic g-value of the system)
         mfreq         (experimental microwave fequency in Hz)
         Harmonic=1    (1 = first derivative. 2 = absorptive)

    out: Int           (intensity vector of the spectrum)

    Kernel function for the calculation of a fast-motion spectrum using
    anistropic linewidths. The anisotropic linewidths are calculated in the
    fast_motion_lw() function.
    If ForBoost is enabled, compiled fortran files are used to carry out
    the critical fast_motion_lw() function.
    """
    Biso = mfreq*(con.h/(con.beta*giso))
    Int = np.zeros(len(Bfield))
    lw_s = fast_motion_lw(Deltag, DeltaA, I, indices, Biso, giso,
                          tcorr, len(resonances))
    lw_s += lw
    resonances = np.digitize(resonances, Bfield)
    for i in range(0, len(resonances)):
        g1 = resonances[i]
        try:
            L = create_Lorentzian(Bfield, Bfield[g1], lw_s[i], Harmonic)
        except:
            L = 0
        Int += L
    if lwG > 0:
        Int = tool.convolution_G(lwG, Bfield, Int)
    return Int


def create_isotropic_spectrum(Bfield, resonances, indices, equiv, I, Harmonic,
                              lw, lwG):
    """create_isotropic_spectrum(Bfield,resonance,indices,equiv,I,Harmonic,lw,
                                 lwG)

    in:  Bfield        (magnetic field vector in mT)
         resonances    (vector with resonance indices in Bfield)
         indices       (list with all index vectors of the coupled system)
         equiv         (vector with number of equivalent nuclei)
         I             (vector with nuclear spin quantum number for all nuclei)
         Harmonic=1    (1 = first derivative. 2 = absorptive)
         lw            (Lorentzian linewidth given in mT as FWHM)
         lwG           (Gaussian linewwidth given in mT as FWHM)

    out: Int           (intensity vector of the spectrum)

    Kernel function for the calculation of a isotropic spectrum using
    identical linewidths. The linewidths are calculated in the
    create_Lorentzian() function.
    """
    Int = np.zeros(len(Bfield))
    fac = 1
    if np.max(equiv) > 1:
        for i in range(0, len(indices)):
            intensitites = tool.generalized_Pascal(int(equiv[i]), I[i])
            fac = np.kron(fac, intensitites)
    else:
        fac = np.ones(len(resonances))
    resonances = np.digitize(resonances, Bfield)
    for i in range(0, len(resonances)):
        g1 = resonances[i]
        try:
            L = fac[i]*create_Lorentzian(Bfield, Bfield[g1], lw, Harmonic)
        except:
            L = 0
        Int += L
    if lwG > 0:
        Int = tool.convolution_G(lwG, Bfield, Int)
    return Int


# *****************************************************************************
# Subfunctions for calculation of resonance fields
# *****************************************************************************


def calculate_resfields(A, g, Nucs, equiv, mfreq, Iequiv, g_n):
    """calculate_resfields(A,g,Nucs,equiv,mfreq)

    in:  A                  (vector of hyperfine tensors given as 3x3 matrix)
         g                  (g-tensor as 3x3 matrix)
         Nucs               (string with coupling nuclei)
         equiv              (vector with numbers of equivalent coupling nuclei)
         mfreq              (experimental microwave frequency in Hz)

    out: resfields          (list with resonance fields)
         indices_list       (list with indices of mI quantum numbers)
         DeltaA_list        (list with A-Aiso*1 tensors)
         Deltag_list        (list with g-giso*1 tensors)
         I_list             (list with nuclear spin quantum numbers)
         giso               (isotropic g-value)

    The function calculates the resonance fields for a given spin spin system.
    The resoance postion are calculated using a fixpoint itertion of an
    implicit expression of the field dependent resonance condition. The
    fixpoint iteration is carried out by the function Breit_Rabi_iteration().
    """
    # Create list to append resonances
    indices_list = []
    DeltaA_list = []
    Deltag_list = []
    resfield_list = []
    hyperfine_dim = len(Iequiv)
    for k in range(0, hyperfine_dim):
        if hyperfine_dim == 1:
            giso, aiso, DeltaA, Deltag = tensor_readout(A, g)
        else:
            giso, aiso, DeltaA, Deltag = tensor_readout(A[k], g)
        # Fix-point iteration (analytical formula for the resfield calculation)
        resfields, indices = Breit_Rabi_iteration(aiso, giso, Iequiv[k],
                                                  g_n[k], mfreq, k)
        resfield_list.append(resfields)
        indices_list.append(indices)
        DeltaA_list.append(DeltaA)
        Deltag_list.append(Deltag)
    if len(indices_list) > 1:
        resfields, indices_list = do_hf_splitting(indices_list, resfield_list)
    return resfields, indices_list, DeltaA_list, Deltag_list, giso


def Breit_Rabi_iteration(aiso, giso, I, g_n, mfreq, k=0):
    """Breit_Rabi_iteration(aiso,giso,I,g_n,mfreq,k=0)

    in:  aiso                  (isotropic hyperfine constant)
         giso                  (isotropic g-value)
         I                     (nuclear spin quantum number)
         g_n                   (nuclear g-value for the used isotope)
         mfreq                 (experimental microwave frequency in Hz)
         k=0                   (k = 0:first nuclei, k > 0 all other splittings)

    out: resfields             (list of resonance fields for the given input)
         indices               (list of mI quantum numbers)

    The function calulates the resonance fields for an electron spin system
    coupled to exactly on nucleus. If k > 0 the Biso resonance position is
    subtracted from the resonance fields. This resonance splittings are
    subsequently introduced to the k = 0 resonance fields via first order
    corrections. The function uses the Breit-Rabi formula for the resonance
    condition and solves the implicit function via fixpoint iteration.
    """
    # Allocation and construction of  mI values
    indices = []
    resfields = []
    B_iso = Biso(giso, mfreq)
    mI = np.zeros(int(2*I)+1)
    for i in range(0, len(mI)):
        mI[i] = -I+i
    aiso = aiso+1e-5
    if aiso < 0:
        fac = -1.0
    else:
        fac = 1.0
    epsilon = (aiso/2.0)/mfreq
    gamma = giso*con.beta+g_n*con.beta_n
    # Fixpoint iteration
    for trans in range(0, len(mI)):
        B0_0 = 0
        for i in range(0, 20):
            B0_k = ((con.h*aiso)/(gamma*(1-epsilon**2))*(-1.0*mI[trans] +
                    fac*np.sqrt(mI[trans]**2+(1-epsilon**2) *
                    (np.power(2*epsilon, -2.0)-(I+1.0/2)**2))))
            if abs(B0_k-B0_0) < Breit_Rabi_thresh:
                break
            else:
                epsilon = (aiso/2.0)/(mfreq+con.beta_n*g_n*B0_k/con.h)
                B0_0 = B0_k
        if k > 0:
            B0_k -= B_iso
        resfields = np.append(resfields, B0_k*con.T2mT)
        indices = np.append(indices, mI[trans])
    return resfields, indices


def Biso(giso, mfreq, unit='mT'):
    """Biso(giso,mfreq,unit='mT')

    in:  giso              (isotropic g-value)
         mfreq             (experimental frequency in Hz)
         unit='mT'         (unit in which the magnetic field value is returned)

    out: B_iso             (magnetic field value for the given giso)

    The function convertes a frequency into a magnetic field value for a given
    isotropic g-value. The default unit is mT. Alternatively it can be returned
    in T (unit='T') or Gauss (unit='G')
    """
    B_iso = (con.h*mfreq)/(giso*con.beta)
    if unit == 'mT':
        pass
    elif unit == 'T':
        B_iso *= 1e-3
    elif unit == 'G':
        B_iso *= 10
    return B_iso


def do_hf_splitting(indices_list, resfield_list):
    """do_hf_splitting(indices_list,resfield_list)

    in:  indices_list            (initial list with mI values)
         resfield_list           (initial list of not coupled resonance fields)

    out: resfields               (list of coupled resonance fields)
         indices_list            (multidimensional list of mI values)

    The function does the first order hyperfine splittings of an initial
    resfield vector (splitting with k = 0, see Breit_Rabi_iteration()). The
    function returns the resonance fields for all coupled hyperfine
    interactions with the corresponding mI values (they are required for
    anisotropic linewidths in the fast-motion regime).
    """
    resfields = 0
    dim = []
    n = len(indices_list)
    for i in range(0, n):
            dim.append(len(indices_list[i]))
    for i in range(0, n):
        for k in range(0, n):
            if k != i:
                if i <= k:
                    indices_list[i] = np.kron(indices_list[i], np.ones(dim[k]))
                    resfield_list[i] = np.kron(resfield_list[i],
                                               np.ones(dim[k]))
                else:
                    indices_list[i] = np.kron(np.ones(dim[k]), indices_list[i])
                    resfield_list[i] = np.kron(np.ones(dim[k]),
                                               resfield_list[i])
        resfields += resfield_list[i]
    return resfields, indices_list


def tensor_readout(A, g):
    """tensor_readout(A,g)

    in:  A                              (hyperfine tensor, given as 3x3 matrix)
         g                              (g-tensor, given as 3x3 matrix)

    out: giso                           (isotropic g-value)
         aiso                           (isotropic hyperfine constant)
         DeltaA                         (A-Aiso*1 tensor)
         Deltag                         (g-giso*1 tensor)

    Take a g-tensor and an A-tensor and returns the isotropic values (the
    trace of the 3x3 matrices) as well as the matrices giving the anisotropies
    A-Aiso*1 and g-giso*1.
    """
    g_ten = np.diag(g)
    giso = 1.0/3*np.trace(g_ten)
    giso_ten = np.eye(3)*giso
    Deltag = g_ten-giso_ten
    try:
        A_ten = np.diag(A)
        aiso = 1.0/3*np.trace(A_ten)
        Aiso_ten = np.eye(3)*aiso
        DeltaA = A_ten-Aiso_ten
    except:
        aiso = 1e-9
        DeltaA = np.zeros((3, 3))
    return giso, aiso, DeltaA, Deltag


# *****************************************************************************
# Subfunctions for calculation linewidths/lineshapes
# *****************************************************************************


def fast_motion_lw(Deltag, DeltaA, I, mI, Biso, giso, tcorr, Nresonances):
    """fast_motion_lw(Deltag,DeltaA,I,mI,Biso,giso,tcorr,Nresonances)

    in:  Deltag_list         (g-giso*1 tensor)
         DeltaA              (vector with with A-Aiso*1 tensors)
         I                   (list with nuclear spin quantum numbers)
         mI                  (list with mI projection quantum numbers)
         Biso                (central magnetic field)
         giso                (isotropic g-value)
         tcorr               (rotational correlation time, given in s)
         Nresonances         (number of resonances)

    out: lw_s                (vector with anisotropic linewidths in mT as FWHM)

    This function sets up the input for the calculation of anistropic
    linewidths in the fast-motion regime. The calculation is than carried out
    by the fast_motion_lw_kernel() function. The function returns a vector
    which contains the anistotropic linewidhts.
    """

    mI = np.asarray(mI)
    lw_s = np.zeros(Nresonances)
    omega_0 = giso*Biso*con.beta/con.h
    MHz2mT = (1e9*con.h)/(giso*con.beta)
    # Do all possible precalculation which are index independent
    j0 = tcorr*2*np.pi
    j1 = j0/(1+(omega_0**2)*(j0**2))
    field_dep = (con.beta*Biso/con.h)
    # Formulas according to N.M.Atherton (p.331-332)
    a11 = (2.0/15)*j0+(j1*1.0)/10
    a22 = (1.0/20)*j0+(j1*7.0)/60
    b11 = (4.0/15)*j0+(1.0/5)*j1
    c11 = (1.0/12)*j0-(1.0/60)*j1
    d11 = (4.0/15)*j0+(1.0/10)*j1
    DeltaA = np.asarray(DeltaA)
    Deltag = np.asarray(Deltag)
    I = np.asarray(I)
    lw_s = fast_motion_lw_kernel(DeltaA, Deltag, a11, a22, b11, c11, d11,
                                 field_dep, mI, I, Nresonances)
    lw_s *= 1e-6*MHz2mT
    return lw_s


@dec_fm_lw_kernel()
def fast_motion_lw_kernel_expired(DeltaA, Deltag, a11, a22,
                                  b11, c11, d11, field_dep,
                                  mI, I, Nresonances):
    """fast_motion_lw_kernel(DeltaA,Deltag,a11,a22,b11,c11,d11,field_dep
                            ,indices,mI,I,Nresonances)

    in:  DeltaA       (vector with with A-Aiso*1 tensors)
         Deltag       (g-Aiso*1 tensor)
         a11          (scalar factor for Kivelson formula)
         a22          (scalar factor for Kivelson formula)
         b11          (scalar factor for Kivelson formula)
         c11          (scalar factor for Kivelson formula)
         d11          (scalar factor for Kivelson formula)
         field_dep    (scalar factor for Kivelson formula)
         mI           (vector with mI projection quantum numbers)
         I            (vector with nuclear spin quantum numbers)
         Nresonances  (number of resonances)

    out: lw_s         (vector with anisotropic linewidths given in Hz and FWHM)

    Kernel function for the calcualtion of anistropic linewidths in the fast-
    motion regime. The anistropic linewidth formula is according to
    N.M.Atherton (p.329-330). It should ne noted that the pseudo-secular
    contributions are not neglected.
    Just-in-time complilation of this function is used if Numba is available.
    """
    lw_s = np.zeros(Nresonances)
    for j in range(0, len(mI)):
        DeltaAA = np.sum(DeltaA[j]*DeltaA[j])
        II1 = I[j]*(I[j]+1)
        if j == 0:
            A = ((field_dep**2)*(np.sum(Deltag[j]*Deltag[j])) *
                 a11+II1*DeltaAA*a22)
        else:
            A = II1*(DeltaAA)*a22
        B = field_dep*(np.sum(Deltag[j]*DeltaA[j]))*b11
        C = DeltaAA*c11
        for i in range(0, Nresonances):
            D = 0
            if len(mI) > 1:
                for k in range(j+1, len(mI)):
                    DT = np.sum(DeltaA[j]*DeltaA[k])*d11
                    D += DT*mI[k, i]*mI[j, i]
            # Set up quadratic formula for anistropic linewidth
            lw_s[i] += (A+B*mI[j, i]+C*(mI[j, i]**2)+D)*2
    return lw_s



def fast_motion_lw_kernel(DeltaA, Deltag, a11, a22, b11, c11, d11, field_dep,
                          mI, I, Nresonances):
    """fast_motion_lw_kernel(DeltaA,Deltag,a11,a22,b11,c11,d11,field_dep
                            ,indices,mI,I,Nresonances)

    in:  DeltaA       (vector with with A-Aiso*1 tensors)
         Deltag       (g-Aiso*1 tensor)
         a11          (scalar factor for Kivelson formula)
         a22          (scalar factor for Kivelson formula)
         b11          (scalar factor for Kivelson formula)
         c11          (scalar factor for Kivelson formula)
         d11          (scalar factor for Kivelson formula)
         field_dep    (scalar factor for Kivelson formula)
         mI           (vector with mI projection quantum numbers)
         I            (vector with nuclear spin quantum numbers)
         Nresonances  (number of resonances)

    out: lw_s         (vector with anisotropic linewidths given in Hz and FWHM)

    Kernel function for the calcualtion of anistropic linewidths in the fast-
    motion regime. The anistropic linewidth formula is according to
    N.M.Atherton (p.329-330). It should ne noted that the pseudo-secular
    contributions are not neglected.
    Just-in-time complilation of this function is used if Numba is available.
    """
    n = len(mI)
    Dij = np.zeros((Nresonances, n), dtype=np.float64)
    DeltaAA = np.zeros(n)
    DeltaAg = np.zeros(n)
    dgg1 = np.sum(Deltag[0]*Deltag[0])
    for j in range(n):
        DeltaAA[j] = np.sum(DeltaA[j]*DeltaA[j])
        DeltaAg[j] = np.sum(DeltaA[j]*Deltag[j])
    if n > 1:
        for j in range(0, n):
            for k in range(j+1, n):
                DAA2 = np.sum(DeltaA[j]*DeltaA[k])*d11
                Dij[:, j] += DAA2*mI[k, :]*mI[j, :]
    lw_s = np.zeros(Nresonances)
    for j in range(0, len(mI)):
        II1 = I[j]*(I[j]+1)
        if j == 0:
            A = field_dep**2*dgg1*a11+II1*DeltaAA[0]*a22
        else:
            A = II1*DeltaAA[j]*a22
        B = field_dep*DeltaAg[j]*b11
        C = DeltaAA[j]*c11
        lw_s += (A+B*mI[j, :]+C*(mI[j, :]**2)+Dij[:, j])*2
    return lw_s


@dec_Lorentzian()
def create_Lorentzian(x, x0, lw, Harmonic=1):
    """create_Lorentzian(x,x0,lw,Harmonic=1)

    in:  x                          (x-vector, here the magnetic field vector)
         x0                         (x0, here the resonance field)
         lw                         (Lorentzian linewidth, given in mT as FWHM)
         Harmonic=1                 (1 = first derivative. 2 = absorptive)

    out: L                          (spectrum with Lorentzian lineshape)
    """
    res = x-x0
    if Harmonic == 0:
        L = ((0.5*lw)/np.pi)/(res**2+(0.5*lw)**2)
    else:
        L = -16.0*((res*lw)/np.pi)/(4.0*(res)**2+lw**2)**2
    return L


def field_interpolation(Int, Points, Bfield):
    """field_interpolation(Int,Points,Bfield)

    in:  Int               (intestity vector of the input signal)
         Points            (whished number of points of the output signal)
         Bfield            (magnetic field vector of the input signal)

    out: field             (interpolated field vector with Points elements)
         signal            (interpolated intensity vector with Points elements)

    Takes aribtray intensity vector and produces a cubic spline
    representation. Subsequently returns the field and corresponding spectrum
    for a user-defined number of points (Points).
    """
    Range = [np.min(Bfield), np.max(Bfield)]
    k1 = np.linspace(Range[0], Range[1], Points, endpoint=True)
    spl = interpolate.splrep(Bfield, Int, k=3)
    signal = interpolate.splev(k1, spl)
    field = k1
    return field, signal


def do_default_range(resfield, lw, lwG, tcorr=None):
    """do_default_range(resfield,lw,lwG,Iso =True, tcorr = None)

    in:  -resfield     (resonance fields)
         -lw           (linewidth parameters Lorentzian)
         -lwG          (linewidth parameter Gaussian)
         -Iso          (boolean variable. If true than it is a isotropic calc.)
         -tcorr        (rotational correlation time. None for isotropic calc.)

    Defines a magnetic field range if not present in the Exp class.
    """
    if isinstance(resfield, float) or isinstance(resfield, int):
        minres = resfield
        maxres = resfield
    else:
        minres = np.min(resfield)
        maxres = np.max(resfield)
    if tcorr is not None:
        factor = 6
        lwtcorr = 5e9*tcorr*((1e-3*maxres)**2)
    else:
        factor = 10
        lwtcorr = 0
    minres -= 5*(lw+lwG)+lwtcorr
    maxres += 5*(lw+lwG)+lwtcorr
    diff = (maxres-minres)/factor
    Bfield = np.linspace(minres-diff, maxres+diff, defPoints)
    return Bfield


def print_Info(SimPar):
    print("\n*******************************************")
    print("************RUN FAST MOTION*************")
    print("*******************************************\n")
    if SimPar.iso:
        print("Isotropic calculation")
    else:
        print("Redfield calculation")
        print("Rotational correlation time: " +str(SimPar._tcorriso)+" ns")
    print("Nuclear spins : "+str(SimPar._Nucsvec))
    print("Number of equivalent nuclei: " + str(SimPar._equiv))
    return