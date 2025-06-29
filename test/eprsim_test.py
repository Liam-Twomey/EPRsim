from eprsim import EPRsim#, Tools
from pytest import raises
import numpy as np
#from test.data_helper import *
import matplotlib.pyplot as plt

sim =  EPRsim
#tool = Tools
dataPath = "simfiles/"
writeData = False # used for writing new comparison data for testing.
lsqTol = 0.05 # difference tolerance in simulation

### Helper Functions -- not read by pytest because don't match "test_*"

def sim_io_helper(testName:str,B0:np.ndarray, spc:np.ndarray,fPath:str, write=False):
	"""
	A helper function for pytest tests to enable reading and writing of simulated
	spectra.
	
	Parameters
	----------
	testName: str
		Name of the test running this function
	B0: np.ndarray
		Abscissa vector of the simulation
	spc: np.ndarray
		Signal vector of the simulation
	fPath: str
		subdirectory of cwd to read/write files from
	write: bool = False
		Whether to write to the file. Should remain false unless simulation methods
		are updated and data needs to be overwritten.

	Returns
	-------
	Status: tuple
		(B0 status, spc status). If reading, statuses will be the data, if writing
		each value will be "written".
	"""
	BFile   = fPath+testName+"_B.bin"
	spcFile = fPath+testName+"_spc.bin"
	if write:
		B0.tofile(BFile)
		spc.tofile(spcFile)
		return ("written", "written")
	else:
		Bt   = np.fromfile(BFile)
		spct = np.fromfile(spcFile) 
		return (Bt,spct)

def sim_diff(BSim,SSim,BRef,SRef,Tol=0.001):
	"""
	Checks for a difference between two simulated signals.
	"""
	# Field checks
	assert(np.shape(BRef) == np.shape(BSim))
	assert(max(BRef) == max(BSim))
	assert(min(BRef) == min(BSim))
	# Signal checks
	assert(np.shape(SRef) == np.shape(SSim))
	# Calculate RMSE and normalize to max val of simulation
	diffMetric = np.sqrt(np.mean((SSim-SRef)**2))/max(SSim) 
	print("Quality:",diffMetric)
	assert(diffMetric < Tol)

### Test functions 

def test_iso_nitrox():
	"""Simple example for the simulation of an isotropic nitroxide spectrum."""
	P = sim.Parameters()
	P.Range = [335, 350]
	P.mwFreq = 9.6
	P.g = 2.002
	P.A = 45.5
	P.Nucs = "N"
	P.lw = [0.2, 0.2]
	P.motion = "fast"
	B0, spc, flag = sim.simulate(P)
	#tool.plot(B0, spc, fignum=1)
	BRef, spcRef = sim_io_helper("iso_nitrox",B0,spc,dataPath, writeData)
	sim_diff(B0,spc,BRef,spcRef) if not writeData else None

def test_aniso_N_fm():
	"""Simple example for the simulation of an anisotropic nitroxide spectrum (only 15N) in the fast-motion regime."""
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
	#tool.plot(B0, spc, fignum=2)
	BRef,spcRef = sim_io_helper("aniso_N_fm",B0,spc,dataPath, writeData)
	sim_diff(B0,spc,BRef,spcRef) if not writeData else None

def test_aniso_N_ss():
	"""Simple example for the simulation of an anisotropic nitroxide spectrum (only 14N) in the solid-state regime."""
	P = sim.Parameters()
	P.Range = [335, 350]
	P.mwfreq = 9.6
	P.g = [2.0083, 2.0061, 2.0022]
	P.A = [12, 13, 110]
	P.Nucs = "14N"
	P.lw = [0.5, 0.2]
	P.motion = "solid"
	B0, spc, flag = sim.simulate(P)
	#tool.plot(B0, spc, fignum=3)
	BRef, spcRef = sim_io_helper("aniso_N_ss",B0,spc,dataPath, writeData)
	sim_diff(B0,spc,BRef,spcRef) if not writeData else None

def test_aniso_NH_ss():
	"""Simple example for the simulation of an anisotropic nitroxide spectrum (only 14N) in the solid-state regime, coupled to an additional hydrogen nucleus."""
	P = sim.Parameters()
	P.Range = [335, 350]
	P.mwfreq = 9.6
	P.g = [2.0083, 2.0061, 2.0022]
	P.A = [[12, 13, 110], [20, 30, 30]]
	P.Nucs = "14N,H"
	P.lw = [0.5, 0.2]
	P.motion = "solid"
	B0, spc, flag = sim.simulate(P)
	#tool.plot(B0, spc, fignum=4)
	BRef, spcRef = sim_io_helper("aniso_NH_ss",B0,spc,dataPath, writeData)
	sim_diff(B0,spc,BRef,spcRef) if not writeData else None

def test_multicomponent_ss(): 
	"""Simple example for the simulation of two radical species."""
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
	#tool.plot(B0, spc, fignum=5)
	BRef,spcRef = sim_io_helper("multicomponent_ss",B0,spc,dataPath, writeData)
	sim_diff(B0,spc,BRef,spcRef) if not writeData else None

def test_triplet():
	"""spin-polarized triplet spectrum."""
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
	#tool.plot(B0, spc, fignum=6)
	BRef,spcRef = sim_io_helper("triplet",B0,spc,dataPath, writeData)
	sim_diff(B0,spc,BRef,spcRef) if not writeData else None

def test_falseNegative():
	"""Tests that the checking for the difference between spectra is working
	by feeding it two different spectra.

	"""
	P = sim.Parameters()
	P.Range = [335, 350]
	P.mwFreq = 9.6
	P.g = 2.002
	P.A = 45.5
	P.Nucs = "N"
	P.lw = [0.2, 0.2]
	P.motion = "fast"
	B0, spc, flag = sim.simulate(P)
	
	#load *different* datafile 
	BRef, spcRef = sim_io_helper("aniso_N_ss",B0,spc,dataPath)
	#These are different speectra, so it should not pass the assertion
	with raises(AssertionError):
		sim_diff(B0,spc,BRef,spcRef) if not writeData else None
