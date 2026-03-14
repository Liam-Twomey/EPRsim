'''
ESFit is responsible for least-squares fitting of a simulation function to data.
'''
import numpy as np
from scipy import optimize as opt

class esfit():
    '''
    The ESFit function performs least-squares fitting on EPR and other data.
    This version is reimplemented from Stefan Stoll's EasySpin.
    
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
    lbnd: ...
        ...?
    ubnd: ...
        ...
    FitOpt: object
        A fit options object.
    
    Attributes
    ----------

    Notes
    -----
    - Scipy optimization [docs](https://docs.scipy.org/doc/scipy/tutorial/optimize.html)
    '''

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

        Attributes
        ----------
        Each argument is set to an attribute of the same name.
        '''
        def __init__(self,algorithm:str='simplex', target:str='fn',
            autoscale:str='lsq', outputs:tuple = (1,1), mask:(np.ndarray|None)=None,
            weight:float=1):
            #raise NotImplementedError("ESFit is still a WIP. Please try again later.")
            self.algorithm  = algorithm
            self.target     = target
            self.autoscale  = autoscale
            self.outputs    = outputs
            self.mask       = mask
            self.weight     = weight


    class fit():
        '''
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
        def __init__(self,pfit:np.ndarray, pnames:list, pfit_full, argsfit, pstd, ci95, cov, corr,
               p_start, fitraw,fit,scale, baseline,mask, residuals,ssr,rmsd,bestFitHistory):
            #raise NotImplementedError("ESFit is still a WIP. Please try again later.")
            self.pfit            = pfit 
            self.pnames          = pnames
            self.pfit_full       = pfit_full
            self.argsfit         = argsfit
            self.pstd            = pstd
            self.ci95            = ci95 
            self.cov             = cov
            self.corr            = corr
            self.p_start         = p_start 
            self.fitraw          = fitraw
            self.fit             = fit
            self.scale           = scale
            self.baseline        = baseline
            self.mask            = mask
            self.residuals       = residuals
            self.ssr             = ssr
            self.rmsd            = rmsd 
            self.bestFitHistory  = bestFitHistory

    def __init__(self, data:np.ndarray, fn:str, par:dict, vary:dict, lbnd, ubnd):
        self.data   = data
        self.fn     = fn
        self.par    = par
        self.vary   = vary
        self.lbnd   = lbnd
        self.ubnd   = ubnd
        #raise NotImplementedError("ESFit is still a WIP. Please try again later.")
        self.fopt = self.fitOpt()
        match self.fopt.algorithm.lower():
            case ('nelder-mead' | 'simplex'):
                raise NotImplementedError("ESFit is still a WIP. Please try again later.")
            case 'bfgs':
                raise NotImplementedError("ESFit is still a WIP. Please try again later.")
            case ('ncg' | 'newton-cg'):
                raise NotImplementedError("ESFit is still a WIP. Please try again later.")
            case ('trust-ncg' | 'tncg'):
                raise NotImplementedError("ESFit is still a WIP. Please try again later.")
            case ('trust-kyrlov' | 'krylov'):
                raise NotImplementedError("ESFit is still a WIP. Please try again later.")
            case 'trust-exact':
                raise NotImplementedError("ESFit is still a WIP. Please try again later.")
            case 'trust-constr':
                raise NotImplementedError("ESFit is still a WIP. Please try again later.")
            case 'slsqp':
                raise NotImplementedError("ESFit is still a WIP. Please try again later.")
