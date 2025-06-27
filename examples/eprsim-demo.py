""" demo of basic functionality. Originally from https://pypi.org/project/EPRsim/ """


import EPRsim_src as es
# Note this only works on non-packaged version, as it references the *folder* EPRsim.
# In the packaged version, should use this instead:
# from EPRsim import EPRsim, Tools
import numpy as np
import matplotlib.pyplot as plt


sim = es.EPRsim # sim = EPRsim #if packaged
tool = es.Tools # tool = Tools #if packaged

run_all = True

# Simple example for the simulation of an isotropic nitroxide spectrum.
if run_all:
    P = sim.Parameters()
    P.Range = [335, 350]
    P.mwFreq = 9.6
    P.g = 2.002
    P.A = 45.5
    P.Nucs = "N"
    P.lw = [0.2, 0.2]
    P.motion = "fast"
    B0, spc, flag = sim.simulate(P)
    tool.plot(B0, spc)

# Simple example for the simulation of an anisotropic nitroxide spectrum (only 14N) in the fast-motion regime.
if run_all:
    Ra = [335, 350]
    freq = 9.6
    g = [2.0083, 2.0061, 2.0022]
    A = [12, 13, 110]
    Nucs = "14N"
    lw = [0.2, 0.2]
    tcorr = 1e-10
    motion = "fast"
    Param = sim.Parameters(
        Range=Ra, g=g, A=A, Nucs=Nucs, mwFreq=freq, lw=lw, tcorr=tcorr, motion=motion
    )
    B0, spc, flag = sim.simulate(Param)
    tool.plot(B0, spc)

# Simple example for the simulation of an anisotropic nitroxide spectrum (only 14N) in the solid-state regime.
if run_all:
    P = sim.Parameters()
    P.Range = [335, 350]
    P.mwfreq = 9.6
    P.g = [2.0083, 2.0061, 2.0022]
    P.A = [[12, 13, 110], [20, 30, 30]]
    P.Nucs = "14N,H"
    P.lw = [0.5, 0.2]
    P.motion = "solid"
    B0, spc, flag = sim.simulate(P)
    tool.plot(B0, spc)


# Simple example for the simulation of an anisotropic nitroxide spectrum (only 14N) in the solid-state regime, coupled to an additional hydrogen nucleus.
if run_all:
    P = sim.Parameters()
    P.Range = [335, 350]
    P.mwfreq = 9.6
    P.g = [2.0083, 2.0061, 2.0022]
    P.A = [[12, 13, 110], [20, 30, 30]]
    P.Nucs = "14N,H"
    P.lw = [0.5, 0.2]
    P.motion = "solid"
    B0, spc, flag = sim.simulate(P)
    tool.plot(B0, spc)


# Simple example for the simulation of two radical species.
if run_all:
    P = sim.Parameters()
    P.Range = [335, 350]
    P.mwfreq = 9.6
    P.g = [2.0083, 2.0061, 2.0022]
    P.A = [12, 13, 110]
    P.Nucs = "14N"
    P.lw = [0.5, 0.2]
    P.motion = "solid"
    P2 = sim.Parameters()
    P2.Range = [335, 350]
    P2.mwfreq = 9.6
    P2.g = 2.0003
    P2.lw = [0.3, 0.0]
    P2.motion = "solid"
    P2.weight = 0.1
    B0, spc, flag = sim.simulate([P, P2])
    tool.plot(B0, spc)

# spin-polarized triplet spectrum
if run_all:
    P = sim.Parameters()
    P.S = 1
    P.Range = [130, 450]
    P.mwfreq = 9.6
    P.g = 2
    P.lw = [4, 1]
    P.D = [-1400, 20]
    P.Population = [0.2, 0.3, 0.4]
    P.Harmonic = 0
    B0, spc, flag = sim.simulate(P)
    tool.plot(B0, spc)
