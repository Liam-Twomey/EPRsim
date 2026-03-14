'''
ESFit is responsible for least-squares fitting of a simulation function to data.

(c) Liam Twomey, 2026
'''
import numpy as np
from scipy import optimize as opt # cut these down to only used parts later.

class fitOpt():
	'''
	Arguments	
	----------
	algorithm: str
		Specify the optimization algorithm: ``simplex``, ``levmar``, ``montecarlo``,
		``genetic``, ``grid``, or ``swarm``.
	target: str
		Specify the target function: ``fn``, ``int``, ``dint``, ``diff``, ``fft``
	autoscale: str
		Property to reference for autoscaling: lsq, maxabs, or none. Defaults to
		lsq for EPRsim simulations, none in all other scenarios.
	outputs: tuple of ints
		``(nOut, iOut)``. ``nOut`` being the number of outputs of the simulation
		function, and ``iOut`` being which of these outputs to use for fitting.
	mask: :class:`np.ndarray`
		Boolean array, same dimensions as the data vector, to adjust what regions
		are fitted.
	weight: :class:`np.ndarray`
		Array of weights to use 

	Notes
	-----
	Each argument is set to an attribute of the same name.

	'''
	def __init__(self,algorithm:str='nelder-mead', target:str='fn',
		autoscale:str='lsq', outputs:tuple = (1,1), mask:(np.ndarray|None)=None,
		weight:float=1):
		#raise NotImplementedError("ESFit is still a WIP. Please try again later.")
		self.algorithm	= algorithm.lower()
		self.target		= target.lower()
		self.autoscale	= autoscale.lower()
		self.outputs	= outputs
		# the selected dimension cannot be greater than the number of dims
		assert (self.outputs[0] >= self.outputs[1])
		self.mask		= mask
		self.weight		= weight


class fitter():
	'''
	The ESFit function performs least-squares fitting on EPR and other data.
	This version is reimplemented from Stefan Stoll's EasySpin.

	Arguments
	---------
	pfit: :class:`np.ndarray`

	pnames: list of str

	pfit_full: ??

	argsfit: ??

	pstd: ??

	ci95: 

	cov: :class:`np.ndarray`

	corr: :class:`np.ndarray`

	p_start: 

	fitraw:

	fit:

	scale:

	baseline:

	mask:

	residuals:

	ssr:

	rmsd: 

	bestFitHistory:

	Attributes
	----------
	'''

	class fit():
		'''
		An object representing the result of fitting by esfit.fitter()

		Parameters
		----------
		data: :class:`np.ndarray`
			Array of experimental data to be fit.
		fn: str
			Simulation function.
		par: :class"`dict`
			Experimental parameters with which `data` was collected.
		vary: :class:`dict`
			Allowed variation in each value of `par`.
		FitOpt: object
			A fit options object.

		..	lbnd: :class:`dict` 
			minimum value for each parameter under optimization
			ubnd: :class:`dict` 
			maximum value for each parameter under optimization
		
		Attributes
		----------

		Notes
		-----
		- scipy.optimization `docs`_
		- I think par, var, lb, and ub can be collapsed into a singe object.
		  Maybe a np.ndarray with columns [par, vary, min, max]?
		- That's the cleanest on the back end, but what will be quickest to
		  write scripts for?
		- Also look at `PyO3`_ and `numeris`_ for fast optimization

		.. _docs: https://docs.scipy.org/doc/scipy/tutorial/optimize.html
		.. _PyO3: https://github.com/PyO3/pyo3 
		.. _numeris: https://crates.io/crates/numeris
		'''
		def __init__(self,pfit:np.ndarray, pnames:list, pfit_full, argsfit, pstd,
			ci95, cov, corr, p_start, fitraw,fit,scale, baseline,mask, residuals,
			ssr,rmsd,bestFitHistory):
			#raise NotImplementedError("ESFit is still a WIP. Please try again later.")
			self.pfit			 = pfit 
			self.pnames			 = pnames
			self.pfit_full		 = pfit_full
			self.argsfit		 = argsfit
			self.pstd			 = pstd
			self.ci95			 = ci95 
			self.cov			 = cov
			self.corr			 = corr
			self.p_start		 = p_start 
			self.fitraw			 = fitraw
			self.fit			 = fit
			self.scale			 = scale
			self.baseline		 = baseline
			self.mask			 = mask
			self.residuals		 = residuals
			self.ssr			 = ssr
			self.rmsd			 = rmsd 
			self.bestFitHistory  = bestFitHistory

	def __init__(self, data:np.ndarray, fn:str, par:dict, vary:dict, fopt:fitOpt):#,
		#lbnd, ubnd):
		self.data	= data
		self.fn		= fn
		self.par	= par
		self.vary	= vary
		#self.lbnd	 = lbnd
		#self.ubnd	 = ubnd
		self.fopt	= fopt
		# Quality checks
		if self.fopt.mask is not None:
			# This may actually need to match only one dimension of data
			assert self.fopt.mask.shape == self.data.shape


	def _sim(self):
		'''
		Function to be fit
		'''
		return

	def _lsq(self):
		'''
		Perform least-squares step, checking goodness of fit for optimization.
		'''
		opt.least_squares(sim)

	def start(self):
		'''
		General goal of the fitter is to:

		1. Calculate a simulation within the bounds of self.par +/- self.vary
		2. Determine the deviation (lsq or otherwise) of the simulation
		3. Iterate an optimization over 1 and 2, minimizing the deviation

		Notes
		-----
		
		For faster optimization, look into `pymoo`_, `optimparallel`_,
		and ``multiprocessing.Pool`` (in that order).

		.. _pymoo: https://pymoo.org
		.. _optimparallel: https://pypi.org/project/optimparallel/
		'''

		opt.minimize(self.lsq,method=self.fopt.algorithm)
		return
