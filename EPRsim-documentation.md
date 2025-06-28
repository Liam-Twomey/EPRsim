# EPRsim documentation
> This information was copied from the original PyPI page, rearranged for clarity by LPT.

Open-source simulation package for cw-EPR spectra. EPRsim has been developed in the group of Prof. Dr. Stefan Weber at the University of Freiburg, Freiburg im Breisgau, Germany, during the last couple of years. EPRsim was developed by Stephan Rein. The program uses several concepts described in [^1]. EPRsim is open-source and available free of charge.

Find the full documentation at the link below: https://www.radicals.uni-freiburg.de/de/software
## Installation

Install EPRsim via pip:
```shell
pip install eprsim
```

## Using EPRsim
### Getting started
Import it as package when running Python.
```python
import EPRsim.EPRsim as sim
```

Define parameters and run the simulation by invoking the simulate() function of EPRsim.

```python
Param = sim.Parameters()
B, spc, flag = sim.simulate(Param)
```

### The Parameters object

A Parameters object can be instantiated with `Param = sim.Parameters(parameter1 = value1, ... , parameterN = valueN)`

The Parameter syntax was kept similar to the one used in EasySpin [1], to make it Optional Parameters (with their defaults):

| Parameter | Default 	| Meaning					|
|-----------|-----------|---------------------------
| mwFreq	| 9.6		| microwave frequency (GHz)|
| A			| None		| Hyperfine couplings (MHz)|
| abund_threshold| 0.0001| Threshold for isotope mixtures|
| D			| None		| Zero-field splitting (MHz) |
| g			| 2.0023193 | g-tensor
|Harmonic|1|Harmonic of the spectrum|
|J|None|Exchange coupling|
|tcorr|None|Rotational correlation time in ns|
|logtcorr|None|Decadic logarithm of tcorr|
|lw|[0.1,0.1]|Line-widths (Gaussian, Lorentzian)|
|ModAmp|0|Modulation amplitude|
|motion|‘solid’|Motional regime|
|mwPhase|0|Microwave phase offset|
|n|1|Number of equivalent nuclei|
|nKnots|12|Initial number of theta values|
|Nucs|None|Isotope specification|
|Points|1024|Number of points|
|Range|[330,360]|Magentic field range in mT|
|S|0.5|Electron spin quantum number|
|SNR|None|Signal-to-noise ratio|
|verbosity|True|Print output information|
|weight|1|Weighting (for multiple species)|
|gFrame|None|Euler angles for the g tensor|
|AFrame|None|Euler angles for the A tensors|
|DFrame|None|Euler angles for the D tensor|
|Temperature|300|Experimental temperature|
|Population|None|Zero-field populations|
|LevelSelect|5e-5|Threshold for level selection|

### The simulate() function

`fieldnumpy, spcnumpy, flaglist = sim.simulate(Param)`

Arguments
- **Parametersobject:** Object containing all simulation parameters.

Returns:

- **fieldnumpy** *(ndarray):* Magnetic field vector
- **spcnumpy** *(ndarray)*: Intesity vector of the cw-EPR signal
- **flaglist** *(int)*: Flags with warning codes (below)

	0. Everything is alright
	1. Solid-state is not possible due to too large matrix dimension.
	2. Fast-motion/iso is not possible due to S > 1/2.

## Examples

Simple example for the simulation of an isotropic nitroxide spectrum.
```python
import EPRsim.EPRsim as sim
P = sim.Parameters()
P.Range = [335 ,350]
P.mwFreq = 9.6
P.g = 2.002
P.A = 45.5
P.Nucs = 'N'
P.lw = [0.2, 0.2]
P.motion = 'fast'
B0, spc, flag = sim.simulate(P)
```
Simple example for the simulation of an anisotropic nitroxide spectrum (only 14N) in the fast-motion regime.
```python
import EPRsim.EPRsim as sim
Ra = [335 ,350]
freq = 9.6
g = [2.0083, 2.0061, 2.0022]
A = [12, 13, 110]
Nucs = '14N'
lw = [0.2, 0.2]
tcorr = 1e-10
motion = 'fast'
Param = sim.Parameters(Range=Ra, g=g, A=A, Nucs=Nucs, mwFreq=freq, lw=lw, tcorr=corr,motion=motion)
B0, spc, flag = sim.simulate(Param)
```

Simple example for the simulation of an anisotropic nitroxide spectrum (only 14N) in the solid-state regime.

```python
import EPRsim.EPRsim as sim
import EPRsim.Tools as tool
P = sim.Parameters()
P.Range = [335 ,350]1] in a fixed-point iteration. Anisotropic line-broadening effects in the fast-motion regime are calculated via the Kivelson formula [
P.mwfreq = 9.6
P.g = [2.0083, 2.0061, 2.0022]
P.A = [[12, 13, 110], [20, 30, 30]]
P.Nucs = '14N,H'
P.lw = [0.5, 0.2]
P.motion = 'solid'
B0, spc, flag = sim.simulate(P)
tool.plot(B0, spc)
```

Simple example for the simulation of an anisotropic nitroxide spectrum (only 14N) in the solid-state regime, coupled to an additional hydrogen nucleus.

```python
import EPRsim.EPRsim as sim
import EPRsim.Tools as tool
P = sim.Parameters()
P.Range = [335 ,350]
P.mwfreq = 9.6
P.g = [2.0083, 2.0061, 2.0022]
P.A = [[12, 13, 110], [20, 30, 30]]
P.Nucs = '14N,H'
P.lw = [0.5, 0.2]
P.motion = 'solid'
B0, spc, flag = sim.simulate(P)
tool.plot(B0, spc)
```

Simple example for the simulation of two radical species.

```python
import EPRsim.EPRsim as sim
import EPRsim.Tools as tool
P = sim.Parameters()
P.Range = [335 ,350]
P.mwfreq = 9.6
P.g = [2.0083, 2.0061, 2.0022]
P.A = [12, 13, 110]
P.Nucs = '14N'
P.lw = [0.5, 0.2]
P.motion = 'solid'
P2 = sim.Parameters()
P2.Range = [335 ,350]
P2.mwfreq = 9.6
P2.g = 2.0003
P2.lw = [0.3, 0.0]
P2.motion = 'solid'
P2.weight = 0.1
B0, spc, flag = sim.simulate([P, P2])
tool.plot(B0, spc)
```

Simple example for the simulation of a spin-polarized triplet spectrum.

```python
import EPRsim.EPRsim as sim
import EPRsim.Tools as tool
P = sim.Parameters()
P.S = 1
P.Range = [130 ,450]
P.mwfreq = 9.6
P.g = 2
P.lw = [4, 1]
P.D = [-1400, 20]
P.Population = [0.2, 0.3, 0.4]
P.Harmonic = 0
B0, spc, flag = sim.simulate(P)
tool.plot(B0, spc)
```
### Technical details

**Main function:** for the simulation of cw-EPR in different motional regimes (isotropic, fast-motion and solid state) All spectra are simulated as field sweep spectra.

**Isotropic/fast-motion:** For the fast-motion regime/isotropic limit, the program solves the implicit Breit-Rabi formula [^1] in a fixed-point iteration. Anisotropic line-broadening effects in the fast-motion regime are calculated via the Kivelson formula [^2]. Currently, Euler angles between tensors are ignored by the algorithm! All tensors (only relevant for fast-motion) need to be in their principal axis system and colinear to each other.

**Solid-state:** In the solid-state regime, the program uses a full matrix diagonalization algorithm. Therefore, only spin systems with a Hilbert space dimension of dim(H) < 512 can be calculated. The powder average is partially generated by interpolation of eigenvalues and transition probabilitites (similar to [^1]). The interpolation level is automatically set by the program. The solid state algorithm treats arbitrary spin systems as long as the Hilbert space dimension is within the threshold. Spin-polarization can be defined (withing the electronic sublevels) as zero-field populations. The program constructs (sparse) density matrices out of the zero-field eigenvectors, to efficiently calculate the population transformation from zero field to high field. Per default, the program calculates with thermal equilibrium. Nuclear quadrupolar couplings (for I > 0.5) are currently not implemented.
## Properties

EPRsim provides:

- Simulation for cw-EPR spectra in the solid-state limit and fast-motion regime
- Flexible simualtion options
- Highly-optimized performance of the simulation algorithm
- Various EPR-data processing function
- Open-source
## Feedback

We are eager to hear about your experiences with GloPel. You can email me at stephan.rein@physchem.uni-freiburg.de.

## References

[^1]: S. Stoll, A. Schweiger, J. Magn. Reson., 2006, 178, 42-55

[^2]: N. M. Atherton, Principles of Electron Spin Resonance, 1993
Acknowledgement

A number of people have helped shaping EPRsim and the ideas behind. First and foremost, Prof. Dr. Stefan Weber and Dr. Sylwia Kacprzak (now Bruker Biospin) were for years the driving force behind EPRsim.
