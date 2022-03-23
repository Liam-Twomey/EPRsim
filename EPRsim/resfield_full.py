# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 14:41:32 2017

@author: stephan
"""


import numpy as np
from scipy import interpolate
from EPRsim.Hamiltonian_Eig import define_nKnots_pattern
import scipy.sparse as sparse
try:
    from numba import int64, float64, int32, jit, types
    Numba = 1
except ImportError:
    Numba = 0
import EPRsim.Tools as tool

# *****************************************************************************
# Physical constants and unit conversion factors + global default settings
# *****************************************************************************
# Load physical constans
con = tool.physical_constants()


def stick_spectrum_calculation(Par):
    """
    Calculates the resonance field and intensities from the eigenvalues/vectors

    Parameters
    ----------
    Par :     :class:`object`
              Object with all user-defined parameters.


    See Also
    --------
    HF_Eig : Function for field-dependent eigendecomposition of the Hamiltonian

    Notes
    -----
    The functions takes the eigenvalues and eigenvectors to calculate resonance
    fields and the corresponding intensities for each transition. To reduce
    the computing time, a transition selection is carried out. First, all
    transitions are removed for which for all orientations the resoance
    frequency is outside the magentic field range. Furhtermore, all transitions
    are removed with small (over different orientations averaged) transition
    probabilities.
    For the remaining transitions the resonance fiels and the corresponding
    intensities are calculated for all orientations. If spin-polarization
    is present, the zero-field populations are transformed to high-field
    populations for all transitions. Otherwise, a thermal equilibrium
    is assumed and the populations are distributed via a Boltzmann function.
    Finally, a post-selection of the transition is carried out and transitions
    with a small (averaged) intensity are removed. The triangular distribution
    of orientation is subsequently interpolated to a regular grid over theta
    and phi.
    The resonance fields for all kept transitions and all orientations are
    saved in the member variable Par.res, while the corresponding intensities
    are saved in the member variable Par.intensity.
    """
    #  Define knots pattern
    Knots_theta_vec = define_nKnots_pattern(Par)
    # Initialization of dimension
    S_sp_x_y = preparations(Par)
    # Check if its spin-polarized. Transfer it to local bool to spare time
    popu, ispopu, rho_0 = check_spin_polarization(Par)
    # Thow all fully off Par.field range transitions
    arg = (Par.nKnots, Par.mwFreq, Par.trans_dim, Par.Hilbert_dim, Par.phinKnots,
           Par.eigval, Par.n_explicit)
    Delta_t, signum = preselect_off_res(*arg)
    Par.trans_dim = len(Delta_t)
    if Par.trans_dim > 30 and Par.Point_Group != "Dhinfty":
        args = (Par, Delta_t, signum, S_sp_x_y, Knots_theta_vec)
        Delta_t, signum, Par.trans_dim = preselect_to_probability(*args)
    arg = (Par, S_sp_x_y, Knots_theta_vec, Delta_t, ispopu, signum,
           rho_0, popu)
    res, intensity, Warning_counter = resonance_loop(*arg)
    args = (Par, intensity, Warning_counter, res, Knots_theta_vec)
    postprocess_resonances(*args)
    if int(np.sum(Warning_counter)) > Par.Transdim:
        Par.field_warning = True
    return Par.intensity, Par.res, Par


def preparations(Par):
    """
    Defines some dimension (transition dimension, Hamiltonian dimension)
    and creates a sparse Matrix out of the x,y-Pauli matrices (Sx+Sy) if
    the dimension of the Pauli matrices is larger than 96 or spin polarization
    is present
    """
    # Initialization of dimension
    Par.trans_dim = int((len(Par.eigvec)-1)*len(Par.eigvec)//2)
    Par.all_trans_dim = Par.trans_dim
    Par.Hilbert_dim = len(Par.Pauli[0, 0, :])
    # Ceate sparse matrices for transition probablities
    if Par.ispopu or Par.Hilbert_dim > 96:
        S_spx = sparse.csr_matrix(Par.S_tot[0])
        S_spy = sparse.csr_matrix(Par.S_tot[1])
    else:
        S_spx = Par.S_tot[0]
        S_spy = Par.S_tot[1]
    S_sp_x_y = S_spx+S_spy
    return S_sp_x_y


def dec_preselect_off_res():
    """
    Creates a decorator for the function preselect_off_res(). If Numba
    is available, the decorator is a Numba jit (just-in-time compilation),
    otherwise the identity decorator is called (this decorator does nothing).
    The jit decorator is called with defined static types to improve the
    performance. Cashing is enabled to avoid avoid compilation each time
    the Python program is invoked.
    """
    if Numba == 1:
        return (jit(types.Tuple((float64[:, :, :, :], int32[:, :]))
                (int64, float64, int64, int64, int64,
                 float64[:, :, :, :], int64),
                nopython=True, cache=True))
    else:
        return dec_identity


def dec_find_resonances():
    """
    Creates a decorator for the function find_resonances(). If Numba
    is available, the decorator is a Numba jit (just-in-time compilation),
    otherwise the identity decorator is called (this decorator does nothing).
    The jit decorator is called with defined static types to improve the
    performance. Cashing is enabled to avoid avoid compilation each time
    the Python program is invoked.
    """
    if Numba == 1:
        return (jit(types.Tuple((float64, float64, float64))(float64[:],
                float64[:], int32), nopython=True, cache=True))
    else:
        return dec_identity


def dec_thermal_popdiff():
    """
    Creates a decorator for the function thermal_popdiff(). If Numba
    is available, the decorator is a Numba jit (just-in-time compilation),
    otherwise the identity decorator is called (this decorator does nothing).
    The jit decorator is called with defined static types to improve the
    performance. Cashing is enabled to avoid avoid compilation each time
    the Python program is invoked.
    """
    if Numba == 1:
        return (jit(float64(float64[:], float64[:],
                int32, int32, float64, float64),
                nopython=True, cache=True))
    else:
        return dec_identity


def dec_identity(ob):
    """
    Identity decorator. If Numba is not available, this decorator is used. The
    decorator does nothing.
    """
    return ob


def population_trans(rho_0, eig, popu):
    """
    Calculates the population transformation from zero-field populations
    to high field populations. Uses sparse density matrix--vector products
    for high-performant calculations.
    """
    Pop = abs(np.vdot(eig, rho_0.dot(eig)))
    return Pop


def transition_probability(eigvec1, eigvec2, S_x_y):
    TransProb = abs(np.vdot(eigvec1, S_x_y.dot(eigvec2)))**2
    return TransProb


@dec_preselect_off_res()
def preselect_off_res(nKnots, mwFreq, trans_dim, Hilbert_dim, phinKnots,
                      eigval_hf, n_explicit):
    signum = np.zeros((trans_dim, 2), dtype=np.int32)
    Delta_t = np.zeros((trans_dim, nKnots, phinKnots, n_explicit))
    p = 0
    for s in range(0, Hilbert_dim):
        for t in range(0, s):
            k = (eigval_hf[s, :, :, :]-eigval_hf[t, :, :, :])
            if not (np.all((k > mwFreq) | (k == 0))):
                if not (np.all((k < mwFreq) | (k == 0))):
                    Delta_t[p] = k-mwFreq
                    signum[p] = [t, s]
                    p += 1
    Delta_t = np.real(Delta_t[0:p])
    return Delta_t, signum


def preselect_to_probability(Par, Delta_t, signum, S_sp_x_y, Knots_theta_vec):
    """
    The function selects transitions according to there transition
    probabilities.

    """
    index = np.zeros(Par.trans_dim, dtype=np.int32)
    prob_tmp = np.zeros((Par.trans_dim, Par.nKnots, Par.phinKnots))
    Par.field_length = len(Delta_t[0, 0, 0, :])-1
    ind = np.round(np.linspace(0, 1, 3)*Par.nKnots)
    ind[2] = ind[2]-1
    ind = ind.astype(int)
    ind2 = np.round(np.linspace(0, 1, 4)*Par.nKnots)
    ind2[3] = ind2[3]-1
    ind2_st = ind2.astype(int)
    for k in range(0, 3):
        if k == 0:
            for g in range(0, Par.trans_dim):
                if np.mod(k, 2) == 0:
                    index[g] = 0
                else:
                    index[g] = Par.field_length
                t1 = signum[g][0]
                t2 = signum[g][1]
                eigv1 = Par.eigvec[t1][0, 0, index[g]]
                eigv2 = Par.eigvec[t2][0, 0, index[g]]
                prob = transition_probability(eigv1, eigv2, S_sp_x_y)
                prob_tmp[g, 0, 0] = prob
        else:
            for q in range(0, k+2):
                if k+1 == 2:
                    ind2 = ind
                else:
                    ind2 = ind2_st
                for i in range(0, Par.trans_dim):
                    if np.mod(k, 2) == 0:
                        index = 0
                    else:
                        index = Par.field_length
                    t1 = signum[i][0]
                    t2 = signum[i][1]
                    eigv1 = Par.eigvec[t1][ind2[q], ind[k], index]
                    eigv2 = Par.eigvec[t2][ind2[q], ind[k], index]
                    prob = transition_probability(eigv1, eigv2, S_sp_x_y)
                    prob_tmp[i, ind2[q], ind[k]] = prob
    maxtransprob = Par.LevelSelect*(np.max(np.sum(np.sum(np.absolute(prob_tmp),
                                                         axis=2), axis=1)))
    for s in range(Par.trans_dim-1, -1, -1):
        if np.sum(np.absolute(prob_tmp[s, :, :])) < maxtransprob:
            Delta_t = np.delete(Delta_t, s, 0)
            signum = np.delete(signum, s, 0)
    Par.trans_dim = len(Delta_t)
    return Delta_t, signum, Par.trans_dim


def resonance_loop(Par, S_sp_x_y, Knots_theta_vec, Delta_t, ispopu, signum,
                   rho_0, popu):
    """
    Function with the loop for finding resonance fiels and intensities to find
    the required information (resonance fields, intensities) for generating
    a stick spectrum. Linear interpolation is for the resonance fields and
    the corresponding intensities. Linear extrapolation is used if the
    resonance is partially (for some orientations) out of the defined
    magentic field range.
    """
    thermal_energy = (Par.T*con.kb)/con.h
    # Allocation of all vectors which should be filled
    intensity = np.zeros((Par.trans_dim, Par.nKnots, Par.phinKnots), order='C')
    res = np.zeros((Par.trans_dim, Par.nKnots, Par.phinKnots), order='C')
    Warning_counter = np.zeros((Par.trans_dim, 2))
    # Extrapolative Par.field if necessary
    stepsize = (Par.field[len(Par.field)-1]-Par.field[0])/(len(Par.field)-1)
    Par.field_extra = np.append(np.append(Par.field[0]-10*stepsize, Par.field),
                                Par.field[len(Par.field)-1]+10*stepsize)
    t1 = []
    t2 = []
    for i in range(0, Par.trans_dim):
        t1.append(signum[i][0])
        t2.append(signum[i][1])
    for i in range(0, Par.trans_dim):
        t1.append(signum[i][0])
        t2.append(signum[i][1])
    # MAIN LOOP FOR FINDING ALL RESONANCE FIELDS
    for k in range(0, Par.nKnots):
        for q in range(0, Knots_theta_vec[k]):
            for i in range(0, Par.trans_dim):
                index = np.argmax(Delta_t[i, k, q, :] > 0)
                Delta = Delta_t[i, k, q, :]
                if index > 0:
                    res[i][k, q], steep, ediff = find_resonance(Par.field,
                                                                Delta, index)
                    fac = ((con.beta*(Par.field[index]-Par.field[index-1])) /
                           (con.h*ediff*1e3))
                    prob1 = transition_probability(Par.eigvec[t1[i]][k, q, index-1],
                                                   Par.eigvec[t2[i]][k, q, index-1],
                                                   S_sp_x_y)
                    prob2 = transition_probability(Par.eigvec[t1[i]][k, q, index],
                                                   Par.eigvec[t2[i]][k, q, index], S_sp_x_y)

                    prob = (1-steep)*prob1+(steep)*prob2
                    Warning_counter[i][0] = 1
                    if ispopu:
                        eigv_1 = (Par.eigvec[t1[i]][k, q, index-1]*(1-steep)+
                                  Par.eigvec[t1[i]][k, q, index]*(steep))
                        eigv_2 = (Par.eigvec[t2[i]][k, q, index-1]*(1-steep)+
                                  Par.eigvec[t2[i]][k, q, index]*(steep))
                else:
                    if min(Delta) > 0:
                        ind = len(Par.field)-1
                        steep = 1
                    else:
                        ind = 0
                        steep = 0
                    index = ind
                    spl = interpolate.splrep(Par.field, Delta, k=1)
                    extrapol = interpolate.splev(Par.field_extra, spl, ext=0)
                    spl2 = interpolate.splrep(Par.field_extra, extrapol, k=3)
                    kt = interpolate.sproot(spl2, mest=1)
                    if len(kt) == 1:
                        res[i][k, q] = kt
                    else:
                        res[i][k, q] = 0
                    fac = 0.5
                    Warning_counter[i][1] = 1
                    eigv1 = Par.eigvec[t1[i]][k, q, ind]
                    eigv2 = Par.eigvec[t2[i]][k, q, ind]
                    prob = transition_probability(eigv1, eigv2, S_sp_x_y)
                    if ispopu:
                        eigv_1 = eigv1
                        eigv_2 = eigv2
                # Spare some calculations for very small transition prob.
                if prob < 1e-6:
                        Warning_counter[i][1] = 0
                else:
                    if ispopu:
                        eigv_1 = eigv_1/np.linalg.norm(eigv_1)
                        eigv_2 = eigv_2/np.linalg.norm(eigv_2)
                        pop1 = population_trans(rho_0[0, k, q], eigv_1, popu)
                        pop2 = population_trans(rho_0[0, k, q], eigv_2, popu)
                        popdiff = np.asscalar(np.real(pop2-pop1))
                    else:
                        popdiff = thermal_popdiff(Par.eigval[:, k, q, index-1],
                                                  Par.eigval[:, k, q, index],
                                                  t1[i], t2[i], thermal_energy,
                                                  steep)
                    intensity[i, k, q] = fac*prob*popdiff
    return res, intensity, Warning_counter


@dec_find_resonances()
def find_resonance(field, Delta, index):
    """
    Function that finds the resonance fields out of the eigenvalue information.
    Linear interpolation between eigenvalue points is used to find the
    resonance fields.
    """
    fielddiff = field[index]-field[index-1]
    ediff = (Delta[index]-Delta[index-1])
    steep = -Delta[index-1]/(ediff)
    res = field[index-1]+steep*fielddiff
    return res, steep, ediff


@dec_thermal_popdiff()
def thermal_popdiff(eigval, eigval2, t1, t2, thermal_energy, steep):
    """
    Calculates the thermal equilibrium populations. This is necessary for
    system with large zero-field splittings or low temperatures.
    For systems where the resonance frequency is substantially smaller than
    the thermal enery, the high-temperature approximation is used.
    Otherwise the fully Boltzman factor is calculated.
    """
    levels = (1-steep)*eigval+(steep)*eigval2
    ekbt = (levels[t2]-levels[t1])/thermal_energy
    if ekbt < 0.1:
        poppges = np.sum(1-levels/thermal_energy)
        popdiff = (1-ekbt)/poppges
    else:
        poppges = np.sum(np.exp(levels/thermal_energy))
        popdiff = np.exp(ekbt)/poppges
    return popdiff


def postprocess_resonances(Par, intensity, Warning_counter, res,
                           Knots_theta_vec):
    """
    input: Par,intensity,Warning_counter,res,Knots_theta_vec

    output: res_new, intensity_new, Warning_counter

    Algorithm:
    1. post_selection()
    Checks the average absolute intensity of a transitions. If it falls under
    a threshold of 5e-4 of the maximal averaged absolute intensity of all
    transitions its deleted. This spares some time in the reconstruction of
    the final stick spectrum. The Warning counter is also deleted for this
    transtition, what means it does not matter if it was a transition which
    was partially out of bounce.
    2. equiv_grid()
    Makes cubic spline interpolation to get the resonance Par.fields and
    intensities for a square grid over theta and phi. This is done for
    each transition.
    (c) Stephan Rein, 29.11.2017
    """
    # Define transition dimension
    intensity, res = post_selection(Par, intensity, res, Warning_counter)
    # Equivalent grid of theta and phi
    equiv_grid(Par, res, intensity, Knots_theta_vec)
    return


def post_selection(Par, intensity, res, Warning_counter):
    """
    Makes a post-selection of resonances. If the transition probability of
    a transition, summed over all orientations is small, the transition
    is removed from the resonance fields and intensities list.
    """
    # Define transition dimension
    Par.Transdim = len(intensity[:, 0, 0])
    # Kick out all transitions with an average intensity smaller than
    maxv = np.max((np.sum(np.sum(np.absolute(intensity), axis=2), axis=1)))
    maxtransprob = Par.LevelSelect*maxv
    for s in range(Par.Transdim-1, -1, -1):
        if np.amax(np.sum(np.absolute(intensity[s, :, :]))) < maxtransprob:
            intensity = np.delete(intensity, s, 0)
            res = np.delete(res, s, 0)
            Warning_counter = np.delete(Warning_counter, s, 0)
    Par.Transdim = len(intensity[:, 0, 0])
    Par.Warning_counter = Warning_counter
    return intensity, res


def equiv_grid(Par, res, intensity, Knots_theta_vec):
    """
    Transforms the triangular grid into a regular grid over theta and phi.
    The final results for the intensities and resonance fields are saved
    in the Par object.
    """
    # Allocation
    Par.intensity = np.zeros((Par.Transdim, Par.nKnots, Par.phinKnots))
    Par.res = np.zeros((Par.Transdim, Par.nKnots, Par.phinKnots))
    # Fill the first orientation vectorized
    Par.intensity[:, 0, :] = np.reshape(np.repeat(Par.intensity[:, 0, 0],
                                        Par.phinKnots, axis=0),
                                        Par.intensity[:, 0, :].shape)
    Par.res[:, 0, :] = np.reshape(np.repeat(res[:, 0, 0],
                                  Par.phinKnots, axis=0),
                                  Par.res[:, 0, :].shape)
    k2 = np.linspace(0, Par.nOctants*np.pi/2, Par.phinKnots, endpoint=True)
    # Cubic spline interpolation to bring it to the square grid over angles
    for k in range(1, Par.nKnots):
        k1 = np.linspace(0, Par.nOctants*np.pi/2, Knots_theta_vec[k],
                         endpoint=True)
        for i in range(0, Par.Transdim):
            spl = interpolate.splrep(k1, intensity[i, k, 0:Knots_theta_vec[k]],
                                     k=3)
            Par.intensity[i, k, :] = interpolate.splev(k2, spl)
            spl2 = interpolate.splrep(k1, res[i, k, 0:Knots_theta_vec[k]], k=3)
            Par.res[i, k, :] = interpolate.splev(k2, spl2)
    return


def check_spin_polarization(Par):
    """
    The function checks if a spin-polarized system should be calculated
    and sets some default parameters.
    """
    if hasattr(Par, 'Population'):
        popu = Par.Population
        ispopu = True
        rho_0 = Par.rho_0
    else:
        popu = None
        ispopu = False
        rho_0 = None
    return popu, ispopu, rho_0
