#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A collection of functions, used for EPR data processing

Copyright, Stephan Rein, 2019
"""

import numpy as np
from scipy import fftpack, special
from scipy import signal
try:
    import matplotlib.pyplot as plt
except ImportError:
    pass

class physical_constants:
    """
    Collection of constants used in EPR simulations
    """
    def __init__(self):
        self.kb = 1.38064852e-23
        self.kb_info = "Boltzman constant in J/K"
        self.T = 300.0
        self.t_info = "Standard temperature in K"
        self.h = 6.62607015e-34
        self.h_info = "Planks constant in J*s"
        self.J2Hz = 1.50919e33
        self.J2Hz_info = "Joule to Hz"
        self.beta = 9.27401007e-24
        self.beta_info = "Bohr magenton in J/T"
        self.beta_n = 5.050783699*1e-27
        self.beta_n_info = "Nuclear Bohr magneton in J/T"
        self.g_free = 2.0023193
        self.g_free_info = "g-factor of the free electron"
        self.m0 = 9.10938356e-31
        self.m0_info = "Rest mass of the electron in kg"
        self.MHz2mT = (1e9*self.h)/(self.g_free*self.beta)
        self.MHz2mT_info = "Conversion MHz to mT"
        self.mT2MHz = 1/self.MHz2mT
        self.mt2MHz_info = "Conversion mT to MHz"
        self.hbar = self.h/(2*np.pi)
        self.hbar_info = "Planks constant in angular frequencies"
        self.g_n = 5.5856946803
        self.g_n_info = "Nuclear g-factor for a proton"
        self.eps0 = 8.85417e-12
        self.eps0_info = "Electric constant in (A*s)/(V*m)"
        self.T2mT = 1e3
        self.T2mT_info = "Conversion T to mT"
        self.mT2T = 1e-3
        self.mT2T_info = "Conversion mT to T"
        self.GHz2Hz = 1e9
        self.GHz2Hz_info = "Conversion GHz to Hz"
        self.Hz2GHz = 1e-9
        self.Hz2GHz_info = "Conversion Hz to GHz"
        self.Hz2MHz = 1e-6
        self.MHz2Hz = 1e6

    def print_all_constants(self):
        """
        Prints a list of available constants and their meaning

        Examples
        --------

        Get a constant, for example the mass of a electron

        >>> import cwEPRsim.Tools as tool
        >>> con = tool.physical_constants()
        >>> m = con.m0
        >>> print(m)
        9.10938356e-31
        >>> con.m0_info
        'Rest mass of the electron in kg'

        Example for plotting all available constants and their meaning

        >>> import cwEPRsim.Tools as tool
        >>> con = tool.physical_constants
        >>> con.print_all_constants()
        INFORMATION ABOUT ALL CONSTANTS:
        kb:      Boltzman constant in J/K
        T:       Standard temperature in K
        h:       Plank constant in J*s
        J2Hz:    Joule to Hertz
        beta:    Bohr magenton in J/T
        beta_n:  Nuclear Bohr magneton in J/T
        g_free:  g-factor of the free electron
        hbar:    Reduced Plank constant in angular frequencies in J*s
        g_n:     Nuclear g-factor for a proton
        eps0:    Electric constant in (A*s)/(V*m)
        m0:      Rest mass of the electron in kg
        MHz2mT:  Constant MHz to mT (for g factor from the free electron)
        mt2MHz:  Constant mT to MHz (for g factor from the free electron)
        T2mT:    Constant T to mT
        mT2T:    Constant mT to T
        GHz2Hz:  Constant GHz to Hz
        Hz2GHz:  Constant Hz to GHz
        Information about a specific variable can be obtained
        by asking for a member variable with an appended _info (e.g. kb_info).
        This variable is then a string, containing the information
        about the constant.
        """
        print("INFORMATION ABOUT ALL CONSTANTS:\n")
        print("kb:      Boltzman constant in J/K")
        print("T:       Standard temperature in K")
        print("h:       Plank constant in J*s")
        print("J2Hz:    Joule to Hertz")
        print("beta:    Bohr magenton in J/T")
        print("beta_n:  Nuclear Bohr magneton in J/T")
        print("g_free:  g-factor of the free electron")
        print("hbar:    Reduced Plank constant in angular frequencies in J*s")
        print("g_n:     Nuclear g-factor for a proton ")
        print("eps0:    Electric constant in (A*s)/(V*m)")
        print("m0:      Rest mass of the electron in kg")
        print("MHz2mT:  Constant MHz to mT")
        print("mt2MHz:  Constant mT to MHz")
        print("T2mT:    Constant T to mT")
        print("mT2T:    Constant mT to T")
        print("GHz2Hz:  Constant GHz to Hz")
        print("Hz2GHz:  Constant Hz to GHz")
        print("\nInformation about a specific variable can be obtained  " +
              " by asking for a member variable with an appended _info " +
              "(e.g. kb_info). This variable " +
              "is then a string, containing the information " +
              "about the constant.")
        return


def pseudo_field_modulation(modamp, field, spectrum):
    """
    Pseudo-field modulation of EPR signals

    Parameters
    ----------
    amp :   :class:`float`
            Modulation amplitude in mT.

    field : :class:`numpy.ndarray`
            Magnetic field vector.

    spc :   :class:`numpy.ndarray`
            Signal vector of the EPR signal.


    Returns
    -------
    spcm :  :class:`numpy.ndarray`
            Field modulated intesity vector of the EPR signal.

    Notes
    -------
    The pseudo-field modulation is calculated as pointwise multiplication
    of the Fourier transform of the signal vector :math:`S(B_0)` with a
    Bessel function of first kind :math:`J(B_0^{-1}, a_0)`:

    .. math:: A = \\mathfrak{F}^{-1}\\left\\lbrace \\mathfrak{F}
              \\left\\lbrace S(B_0)
              \\right\\rbrace \\odot J_{1}(B_0^{-1}, a_0) \\right\\rbrace ,

    where :math:`a_0` is the modulation amplitude and the Bessel function
    of first kind is given as:

    .. math:: J_{1}(B_0^{-1},a_0) =
              \\sum\limits_{m = 0}^\\infty
              \\dfrac{(-1)^m}{m! \\Gamma (m+2)}
              \\left(\\dfrac{a_0k}{2}\\right)^{2m+1},

    where :math:`\\Gamma` is the Gamma function.
    Note that non-zero intenities of the signal vector at the edges of the
    magetic field will lead to artifacts.


    References
    ----------
    [1] : J. S. Hyde, M. Pasenkiewicz-Gierula, A. Jesmanowicz, W. E. Antholine,
    Appl. Magn. Reson. 1990, 1, 483-496.


    Examples
    --------
    Example for a data vector spc and a magnetic field vector field.

    >>> import numpy as np
    >>> import cwEPRsim.Tools as tool
    >>> field, spc = np.loadtxt("Absorptive_EPR_spectrum.txt")
    >>> modAmp = 0.1
    >>> spc_mod = tool.pseudo_field_modulation(modAmp, field, spc)
    """
    modamp = max(0.001, modamp)
    # Set up x-axis in inverse field domain
    fieldstep = (field[1]-field[0])
    fourier_x = fftpack.fftfreq(len(spectrum), fieldstep)
    # Modulation amplitude in conjugated Fourier space
    level = modamp*np.pi
    datafft = fftpack.fft(spectrum)
    bessel_fft = special.jv(1, fourier_x*level)
    pseudomod = 1j*bessel_fft*datafft
    # Inverse Fourier transform of the convoluted signal
    spcm = np.real(fftpack.ifft(pseudomod))
    return spcm


def degree_in_rad(angle):
    angle = (np.pi*angle)/180
    return angle


def GHz2mT(freq, giso=None):
    """
    Conversion from GHz to mT

    Parameters
    ----------
    freq :  :class:`numpy.ndarray` or :class:`float`
            Frequency in GHz.
    giso :   :class:`float`, optional
            Isotropic g-value.

    Returns
    -------
    field : :class:`numpy.ndarray` or :class:`float`
            Magnetic field in mT.

    Notes
    -------
    The function convertes a frequency value or vector into
    a magnetic field value or vector, respectively.
    The conversion is calculates using:

    .. math:: B_0 = \\dfrac{h \\cdot \\omega \\cdot 10^{12}}
              {g_{\\mathrm{iso}}\\cdot\\beta} \; ,

    where :math:`\\omega` is the microwave frequency in GHz,
    :math:`\\beta` is the Bohr magneton,
    :math:`h` is the Plank constant, :math:`g_{\\mathrm{iso}}` the g-value and
    :math:`B_0` is the resulting magnetic field value.

    If no giso value is given, the g-value of the free electron
    (:math:`g_{\\mathrm{iso}}` = 2.0023193) is used.


    Examples
    --------

    >>> import cwEPRsim.Tools as tool
    >>> freq = 9.6
    >>> B0 = tool.GHz2mT(freq)
    >>> print(B0)
    >>> 342.5518855564894

    >>> import cwEPRsim.Tools as tool
    >>> freq = [96, 97.5]
    >>> giso = 2.30
    >>> B0 = tool.GHz2mT(freq, giso)
    >>> print(B0)
    >>> [2982.16631174 3028.76266037]
    """
    cons = physical_constants()
    if giso is None:
        giso = cons.g_free
    if not isinstance(freq, float) or not isinstance(freq, int):
        freq = np.asarray(freq)
    field = (cons.h*freq*1e12)/(giso*cons.beta)
    return field


def mT2GHz(field, giso=None):
    """
    Conversion from mT to GHz

    Parameters
    ----------
    field : :class:`numpy.ndarray` or :class:`float`
            Magnetic field in mT.
    giso :   :class:`float`, optional
            Isotropic g-value.

    Returns
    -------
    freq  : :class:`numpy.ndarray` or :class:`float`
            Frequency in GHz.

    Notes
    -------
    The function convertes a frequency value or vector into
    a magnetic field value or vector, respectively.
    The conversion is calculates using:

    .. math:: \\omega = \\dfrac{B_0 \\cdot g_{\\mathrm{iso}}\\cdot\\beta
              \\cdot 10^{-12}}{h} \; ,

    where :math:`B_0` is the resulting magnetic field,
    :math:`\\beta` is the Bohr magneton,
    :math:`h` is the Plank constant, :math:`g_{\\mathrm{iso}}` the g-value and
    :math:`\\omega` is the resulting microwave frequency in GHz.

    If no giso value is given, the g-value of the free electron
    (:math:`g_{\\mathrm{iso}}` = 2.0023193) is used.


    Examples
    --------

    >>> import cwEPRsim.Tools as tool
    >>> B0 = 332.5
    >>> giso = 2.05
    >>> freq = tool.mT2GHz(B0, giso)
    >>> print(freq)
    >>> 9.540190522093772

    >>> import cwEPRsim.Tools as tool
    >>> freq = tool.mT2GHz(B0)
    >>> print(freq)
    >>> 9.31829639417826

    >>> import cwEPRsim.Tools as tool
    >>> B0 = [1120, 1152.5]
    >>> freq = tool.mT2GHz(B0)
    >>> print(B0)
    >>> [31.38794575 32.29875667]
    """
    cons = physical_constants()
    if giso is None:
        giso = cons.g_free
    if not isinstance(field, float) or not isinstance(field, int):
        field = np.asarray(field)
    mwFreq = ((giso*field*1e-12*cons.beta)/cons.h)
    return mwFreq


def modulation_amplitude(modamp, field, spectrum):
    """
    Modulation amplitude to simulate the effect of overmodulation

    Parameters
    ----------
    modamp :  :class:`float`
              Modulation amplitude in mT.

    field :   :class:`numpy.ndarray`
              Magnetic field vector in mT.

    spectrum : :class:`numpy.ndarray`
               Intensity vector of the EPR signal.

    Returns
    -------
    spectrumm : :class:`numpy.ndarray`
                Intensity vector of the EPR signal after taking the modulation
                amplitude into account.

    See Also
    --------
    pseudo_field_modulation : Pseudo-field modulation of EPR spectra

    Notes
    -------

    The function does not calculate an approximation of a derivative
    (for this see pseudo_field_modulation()). It merely takes into account
    the effect of field modulation on the shape/width of the EPR signal.
    A large modulation  amplitude compared to the line-width of the signal
    will lead to significant distortions of the EPR spectrum.

    The function calculates a pseudo-field modulation of the signal
    using the given modulation amplitude and subsequently reintegrates the
    signal.


    Examples
    --------

    >>> import numpy as np
    >>> import cwEPRsim.Tools as tool
    >>> field, spc = np.loadtxt("EPR_spectrum.txt")
    >>> ModAmp = 0.5
    >>> spcm = tool.modulation_amplitude(ModAmp, field, spc)
    """
    if modamp != 0:
        spectrum = pseudo_field_modulation(modamp, field, spectrum)
        spectrum = np.cumsum(spectrum)
    return spectrum


def phase_offset(gamma, spc, unit="rad"):
    """
    Phase offset of a real valued EPR intensity vector


    Parameters
    ----------
    gamma : :class:`float`
            Phase angle for the microwave phase offset.
    spc :   :class:`numpy.ndarray`
            Real-valued EPR signal vector.
    unit :  :class:`str`, optional
            Defines the unit: 'degree' for degree, 'rad' for radian (default).

    Returns
    -------
    spcp  : :class:`numpy.ndarray`
            Real-valued EPR signal vector with phase offset.

    Notes
    -------

    The function takes the signal vector and applies a Hilbert
    transform to reconstruct the imaginary part. The signal
    vector :math:`S` has now a real and imaginary part.
    The signal with phase offset is calculated by
    rotating :math:`S_{\\mathrm{off}}` in the complex plane about the
    angle :math:`\\gamma`:

    .. math:: S_{\\mathrm{off}} =
              \\mathfrak{Re} \\left\\lbrace
              \\mathrm{e}^{-\\mathrm{i} \\cdot \\gamma} \\cdot S
              \\right\\rbrace


    Only the real part of the rotated spectrum is returned.

    Examples
    --------

    >>> import numpy as np
    >>> import cwEPRsim.Tools as tool
    >>> field, spc = np.loadtxt("EPR_spectrum.txt")
    >>> gamma = 20.5
    >>> spcp = tool.pseudo_field_modulation(gamma, spc, 'degree')
    """
    if unit == "degree":
        gamma = (np.pi*gamma)/180
    spc_im = signal.hilbert(spc)
    spc_im = np.exp(-1j*gamma)*spc_im
    return np.real(spc_im)


def convolution_L(width, field, spectrum):
    """
    Convolution of a signal with a Lorentzian function

    Parameters
    ----------
    width :   :class:`float`
              Line-width in mT, given as full width at half maximum (FWHM)

    field :   :class:`numpy.ndarray`
              Magnetic field vector in mT.

    spectrum : :class:`numpy.ndarray`
               Intensity vector of the EPR signal.

    Returns
    -------
    spectrumc : :class:`numpy.ndarray`
                Intensity vector of the EPR signal convoluted with a Lorentzian

    Notes
    -------
    The function convolutes the EPR signal with a Lorentzian function of a
    given line-width. The function uses a FFT approach to calculate the
    signal convolution.


    Examples
    --------

    >>> import numpy as np
    >>> import cwEPRsim.Tools as tool
    >>> field, spc = np.loadtxt("EPR_spectrum.txt")
    >>> FWHM = 0.5
    >>> spcc = tool.convolution_L(FWHM, field, spc)
    """

    npoints = len(spectrum)
    mid = npoints//2
    Bcenter = field[mid]
    std = width
    res = field-Bcenter
    L = ((0.5*std)/np.pi)/(res**2+(0.5*std)**2)
    spectrum_conv = np.convolve(spectrum, L, mode='same')
    spectrum = spectrum_conv
    return spectrum_conv


def convolution_G(width, field, spectrum):
    """
    Convolution of a signal with a Gaussian function

    Parameters
    ----------
    width :   :class:`float`
              Line-width in mT, given as full width at half maximum (FWHM)

    field :   :class:`numpy.ndarray`
              Magnetic field vector in mT.

    spectrum : :class:`numpy.ndarray`
               Intensity vector of the EPR signal.

    Returns
    -------
    spectrumc : :class:`numpy.ndarray`
                Intensity vector of the EPR signal convoluted with a Gaussian


    Notes
    -------
    The function convolutes the EPR signal with a Gaussian function of a
    given line-width. The function uses a FFT approach to calculate the
    signal convolution.


    Examples
    --------

    >>> import numpy as np
    >>> import cwEPRsim.Tools as tool
    >>> field, spc = np.loadtxt("EPR_spectrum.txt")
    >>> FWHM = 0.2
    >>> spcc = tool.convolution_G(FWHM, field, spc)
    """
    npoints = len(spectrum)
    mid = npoints//2
    std = (width/(2*np.sqrt(2*np.log(2))))
    G = np.exp(-0.5*(field-field[mid])**2/(std**2))
    spectrum_conv = np.convolve(spectrum, G, mode='same')
    spectrum = spectrum_conv
    return spectrum_conv


def normalize2area(spectrum, Harmonic=1):
    """
    Normalization of EPR spectra to their area

    Parameters
    ----------
    spectrum : :class:`numpy.ndarray`
               Intensity vector of the EPR signal.

    Harmonic : :class:`int`
               0: absorptive EPR signal, 1: first derivative (default)

    Returns
    -------
    spectrumn : :class:`numpy.ndarray`
                Intensity vector of the EPR signal normalized to the area

    Notes
    -------

    Normalizes EPR signals to their area. The function can be applied
    to absorptive EPR signals and EPR signals measured as first derivative.

    Examples
    --------

    Normalization of a EPR spectra, measured as first deriative

    >>> import cwEPRsim.EPRsim as sim
    >>> import cwEPRsim.Tools as tool
    >>> Param = sim.Parameters()
    >>> field, spc = sim.simulate(Param)
    >>> spcn = sim.normalize2area(spc)
    """

    if Harmonic == 1:
        pass
        spectrum = spectrum/np.sum(np.absolute(np.cumsum(spectrum)))
    else:
        spectrum = spectrum/np.sum(np.absolute(spectrum))
    return spectrum


def generalized_Pascal(n, I):
    """
    Generates the values of a row of a generalized Pascal triangle

    Parameters
    ----------
    n :     :class:`int`
            Number of equivalent nuclei.

    I :    :class:`int`
           Nuclear spin quantum number

    Returns
    -------
    P :    :class:`numpy.ndarray`
            n-th row of the generalized Pascal triangle for a spin quantum
            number I.


    Notes
    -------

    Uses a general approach to calulate the relevant row of a generalized
    Pascal triangle, a generalization of a Pascal triangle which is only
    valid for a spin quantum number of I = 1/2. The function returns
    the intesity pattern of n equivalent coupling spins with
    spin quantum number I.

    Note that the rows start with n = 1 (no coupling).


    Examples
    --------

    Intenisty pattern of a spin (I=1/2) coupled to 3 equivalent spins.

    >>> import cwEPRsim.Tools as tool
    >>> n = 3
    >>> I = 0.5
    >>> P = tool.generalized_Pascal(n, I)
    >>> print(P)
    [1. 3. 3. 1.]


    Intenisty pattern of a spin (I=1) coupled to 3 equivalent spins.

    >>> import cwEPRsim.Tools as tool
    >>> n = 3
    >>> I = 1
    >>> P = tool.generalized_Pascal(n, I)
    >>> print(P)
    [1. 3. 6. 7. 6. 3. 1.]

    Intenisty pattern of a spin (I=3/2) coupled to 4 equivalent
    spins.

    >>> import cwEPRsim.Tools as tool
    >>> n = 4
    >>> I = 3/2
    >>> P = tool.generalized_Pascal(n, I)
    >>> print(P)
    [ 1.  4. 10. 20. 31. 40. 44. 40. 31. 20. 10.  4.  1.]
    """
    if n == 0:
        return [1.0]

    s0 = int(2*I*n+1)
    A = np.zeros((n, s0))
    A[0, 0:int(2*I+1)] = 1
    I2 = 2*I
    for i in range(1, n):
        for j in range(0, s0):
            if j+I2 >= s0:
                ub = int(s0-1)
            else:
                ub = int(j+I2)
            if j-I2 < 0:
                lb = 0
            else:
                lb = int(j-I2)
            A[i, j] = np.sum(A[i-1, lb:ub-int(2*I)+1])
    return A[n-1, :]


def gyro2gn(gyro):
    """
    Convertes a nuclear gyromagnetic ratio into a nuclear g-factor

    Parameters
    ----------
    gyro :   :class:`float` or :class:`numpy.ndarray`
             Gyromagnetic ratio of an isotope, given in MHz/T.
             Alternatively, an array of gyromagnetic ratios can be passed.

    Returns
    -------
    gn :    :class:`float` or :class:`numpy.ndarray`
            Nuclear g-factor for the isotope(s) (dimensionless).


    Notes
    -------
    The function convertes the gyromagnetic ratio into a
    nuclear g-factor using:

    .. math:: g_{\\mathrm{n}} =
              \\dfrac{10^{6}\\cdot\\gamma\\dot h}{\\beta_{\\mathrm{n}}} \,

    Where :math:`\\gamma` is the gyromagnetic ratio, :math:`h` is the
    Planck constant, :math:`\\beta_{\\mathrm{n}}` is the nuclear Bohr
    magneton. :math:`g_{\\mathrm{n}}` is the obtained nuclear g-factor.

    Examples
    --------

    Example of a gyromagnetic ratio of a 1H nucleus
    (42.57747876 MHz/T).

    >>> import cwEPRsim.Tools as tool
    >>> gyro = 42.57747876
    >>> gn = tool.gyro2gn(gyro)
    >>> print(gn)
    5.585694680337019
    """
    con = physical_constants()
    gn = (1e6*gyro*con.h)/con.beta_n
    return gn


def add_noise(spc, SNR):
    """
    Adds Gaussian noise to a signal

    Parameters
    ----------
    spc : :class:`numpy.ndarray`
          Intensity vector of the EPR signal.

    SNR : :class:`float`
          Signal-to-noise ratio with respect to the intensity.

    Returns
    -------
    spc :  :class:`numpy.ndarray`
           Intensity vector of the EPR signal with addition of noise.


    Notes
    -------
    Adds Gaussian noise (white noise to a signal). The noise level is
    reference to the maximum intesity of the noise-free signal. SNR decribes
    the relative singal-to-noise ratio where the input variable
    SNR is the standard deviation
    :math:`\\sigma`
    of the probability density :math:`p(x)`, normalized to the
    maximum intesity of the signal:

    .. math:: p(x) = \\dfrac{1}{\\sqrt{2\\pi \\sigma^2}}\\cdot\\mathrm{exp}
              \\left[-\\dfrac{x^2}{2\\sigma^2}\\right] \, .

    The vector :math:`p(x)` is added to the signal.

    Examples
    --------

    Example of adding noise for a final SNR of 10.

    >>> import cwEPRsim.EPRsim as sim
    >>> import cwEPRsim.Tools as tool
    >>> Range = [339, 353]
    >>> mwFreq = 9.7
    >>> g = [2.0086, 2.0064, 2.0021]
    >>> A = [17.4, 20, 105]
    >>> Nucs = "14N"
    >>> lw = [0.1, 0.1]
    >>> tcorr = 7e-10
    >> Param2 = sim.Parameters(mwFreq=mwFreq, A=A, Nucs=Nucs, g=g, Range=Range,
                              tcorr=tcorr, lw=lw)
    >>> field, spc = sim.simulate(Param2)
    >>> tool.plot(field, spc)
    >>> spc = tool.add_noise(spc, 10)
    >>> tool.plot(field, spc)

    .. figure:: ../docfigures/Noisfree.png
        :width: 60 %
        :align: center

    .. figure:: ../docfigures/SNR_10.png
        :width: 60 %
        :align: center

    Example for adding noise directly by passing the information to
    the simulation routine.

    >>> import cwEPRsim.EPRsim as sim
    >>> import cwEPRsim.Tools as tool
    >>> Range = [339, 353]
    >>> mwFreq = 9.7
    >>> g = [2.0086, 2.0064, 2.0021]
    >>> A = [17.4, 20, 105]
    >>> Nucs = "14N"
    >>> lw = [0.1, 0.1]
    >>> tcorr = 7e-10
    >>> SNR = 30
    >> Param = sim.Parameters(mwFreq=mwFreq, A=A, Nucs=Nucs, g=g, Range=Range,
                              tcorr=tcorr, lw=lw, SNR=SNR)
    >>> field, spc = sim.simulate(Param)
    >>> tool.plot(field, spc)

    .. figure:: ../docfigures/SNR_30.png
        :width: 60 %
        :align: center
    """
    maxSNR = np.max(np.abs(spc))
    lenghtspc = len(spc)
    noisvec = np.random.normal(0, maxSNR/SNR, lenghtspc)
    spc += noisvec
    return spc


def Eulermatrix(phi, theta, psi=0):
    """
    input: tensor, phi, theta

    output: rotatedTensor

    Algorithm:
    Euler transformation using y-convention. The euler matrix is set up
    with the given angles. Phi and theta is necessary, psi is optional.
    The Euler matrix is returned
    """
    # Allocations
    cosphi = np.cos(phi)
    sinphi = np.sin(phi)
    costhet = np.cos(theta)
    sinthet = np.sin(theta)
    cospsi = np.cos(psi)
    sinpsi = np.sin(psi)
    eulermatrix = np.zeros((3, 3))
    # Set up the full 3-dimensional Euler matrix
    eulermatrix[0][0] = cosphi*costhet*cospsi - sinphi*sinpsi
    eulermatrix[0][1] = -cosphi*costhet*sinpsi - sinphi*cospsi
    eulermatrix[0][2] = cosphi*sinthet
    eulermatrix[1][0] = sinphi*costhet*cospsi + cosphi*sinpsi
    eulermatrix[1][1] = - sinphi*costhet*sinpsi + cosphi*cospsi
    eulermatrix[1][2] = sinphi*sinthet
    eulermatrix[2][0] = -sinthet*cospsi
    eulermatrix[2][1] = sinthet*sinpsi
    eulermatrix[2][2] = costhet
    return eulermatrix


def tensor_rotation(tensor, eulermatrix):
    """
    input: tensor, eulermatrix

    output: rotatedTensor
    The euler matrix O of the SO(3) Group in y-convention is set up in already
    multiplicated form. Than the orthogonal similarity transformation of the
    tensor T is carried out:
    T' = O^(-1)*T*O    with  O^(-1) = O^T
    The rotated tensor is returned

    (c) Stephan Rein, 31.10.2017
    """
    # Final two sided matrix multiplication (similarity transformation)
    eulermatrix_transpose = eulermatrix.T
    rotatedTensor = eulermatrix_transpose@(tensor@eulermatrix)

    return rotatedTensor


#************************************************************************
# Plot Results
#************************************************************************

def plot(field, spectrum, font=11):
    """
    Plots cw-EPR data

    Parameters
    ----------
    field :  :class:`numpy.ndarray`
             Intensity vector of the EPR signal.

    spectrum : :class:`numpy.ndarray`
               Magnetic field vector of the EPR signal in mT.

    font : :class:`int`
           Font size of the axis labels.

    Notes
    -------
    The function relies on the matplotlib data visualization framework.
    The function creates high-quality figures of cwEPR spectra. Margins
    are set adaptively. The fontsize can be adjusted if necessary.
    The function normalizes the signal to the absolute maximum.


    Examples
    --------

    Standard plot of a simulated EPR spectrum, using standard parameters.

    >>> import cwEPRsim.EPRsim as sim
    >>> import cwEPRsim.Tools as tool
    >>> Param = sim.Parameters()
    >>> field, spc = sim.simulate(Param)
    >>> tool.plot(field, spc)

    .. figure:: ../docfigures/Plot_1_example.png
        :width: 60 %
        :align: center

    Plot of a spimple isotropic nixtroxide EPR spectrum with adjusted
    label-fontsize.

    >>> import cwEPRsim.EPRsim as sim
    >>> import cwEPRsim.Tools as tool
    >>> Param = sim.Parameters()
    >>> Param.Nucs = "14N"
    >>> Param.A  = 50
    >>> field, spc = sim.simulate(Param)
    >>> tool.plot(field, spc, 14)

    .. figure:: ../docfigures/Plot_2_example.png
        :width: 60 %
        :align: center

    """

    try:
        spectrum = spectrum/max(abs(spectrum))
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.margins(x=0, y=0.08)
        ax.set_xlabel(r'$magnetic\ field\, / \,  \rm{mT}$', fontsize=font)
        ax.set_ylabel(r'$Normalised\ Intensity$', fontsize=font)
        ax.plot(field, spectrum)
        plt.show()
    except:
        print("Matplotlib is not available")
        pass
    return
