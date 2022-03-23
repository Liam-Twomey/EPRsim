# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 12:59:48 2017

@author: stephan
"""
import numpy as np
import EPRsim.Nucdic as nucdic
from EPRsim.Hamiltonian_Point_Group import Symmetry_Group
import EPRsim.Tools as tool

# *****************************************************************************
# Physical constants and unit conversion factors + global default settings
# *****************************************************************************
# Load physical constans
con = tool.physical_constants()


def convert_user_input_and_Set_up_defaults(Par, SimPar):
    """
    Convertes the user input into tensors and arrays used in the simulation

    Parameters
    ----------
    Par :     :class:`object`
              Object with all user-defined parameters.

    SimPar : :class:`object`
             Object with all simulation parameters.

    Returns
    -------
    Par :     :class:`object`
              Object with all user-defined parameters.

    SimPar : :class:`object`
             Object with all simulation parameters.


    See Also
    --------
    Parameters : Class of Parameters (with a full list of parameters)

    Notes
    -----
    The functions defines addtional parameters in the Par object,
    necessary for the simulation. For example from the information of
    Par.A the tensor(s) Par.A_tensors are generated. Essential paramaeters
    like the

    Examples
    --------

    Example for getting back the hyperfine tensors, out of the input
    information as well as the point group of the Hamiltonian

    >>> import EPRsim.EPRsim as sim
    >>> import EPRsim.Presettings as pre
    >>> P = sim.Parameters()
    >>> P.Range = [335 ,350]
    >>> P.mwFreq = 9.6
    >>> P.g = 2.002
    >>> P.A = [[12, 120], [100, 100]]
    >>> P.Nucs = '14N,1H'
    >>> Val = sim.Validate_Parameters(P)
    >>> P, B  pre.convert_user_input_and_Set_up_defaults(P, Val.Sim_objects[0])
    >>> P.A_tensor
    [[[1.2e+07 0.0e+00 0.0e+00]
    [0.0e+00 1.2e+07 0.0e+00]
    [0.0e+00 0.0e+00 1.2e+08]]
    [[1.0e+08 0.0e+00 0.0e+00]
    [0.0e+00 1.0e+08 0.0e+00]
    [0.0e+00 0.0e+00 1.0e+08]]]
    >>> P.Point_Group
    Dhinfty
    """
    Par.warning = 0
    # Convert the microwave frequency to Hz
    Par.mwFreq = Par.mwFreq*con.GHz2Hz
    # Spectrum exceeds field range boolean
    Par.S = Par.S*1.0
    # Some simulation thresholds
    threshold_settings(Par)
    # D-Tensor
    define_D_tensors(Par)
    # Get number of coupled electrons
    get_number_of_coupled_electrons(Par)
    # Define dimension of hyperfines and Hyperfine Tensors
    define_hyperfines(Par, SimPar)
    if Par.warning == 1:
        return  Par, SimPar
    # g-Tensor
    Par.g_tensor = create_g_or_hyperfine_tensors(Par.g, Par.coupled_e_dim)
    zero_field_populations(Par)
    # Do rotations to general Frame
    rotate_coordsystems(Par)
    symmetry(Par)
    allowed_trans(Par)
    return Par, SimPar


def allowed_trans(Par):
    Par.allowed_NMR_trans = Par.e_dimension*(Par.dim_nuc_tot-1)
    Par.allowed_EPR_trans = (Par.e_dimension-1)*Par.dim_nuc_tot
    return


def get_number_of_nuclei(Par):
    if Par.A is not None:
        Par.A = np.asarray(Par.A)
        Par.A_tensor = np.zeros((Par.A.ndim, 3, 3))
        if Par.A.ndim > 1:
            Par.number_of_nuclei = Par.A.shape[0]
        else:
            Par.number_of_nuclei = 1
    else:
        Par.number_of_nuclei = 0
    return


def get_number_of_coupled_electrons(Par):
    if isinstance(Par.S, float):
        Par.coupled_e_dim = 1
        Par.e_dimension = int(2*Par.S+1)
    else:
        Par.coupled_e_dim = len(Par.S)
        Par.e_dimension = 1
        for i in range(0, Par.coupled_e_dim):
            Par.e_dimension *= int(2*Par.S[i]+1)
    return


def get_full_nuclear_dimension(Par):
    if hasattr(Par, 'ENucCoupling') and Par.SepHilbertspace:
        Par.dim_nuc = np.ones((Par.coupled_e_dim, Par.number_of_nuclei),
                              dtype=np.int)
        Par.dim_nuc_tot = 0
        for s in range(0, Par.coupled_e_dim):
            for i in range(0, Par.number_of_nuclei):
                if Par.ENucCoupling[s, i]:
                    nucstring = Par.Nucs[i]  # each nuclei seperately
                    Par.dim_nuc[s, i] = nucdic.nuclear_properties(nucstring)[1]
            Par.dim_nuc_tot_tmp = int(np.prod(Par.dim_nuc[s]))
            if Par.dim_nuc_tot < Par.dim_nuc_tot_tmp:
                Par.dim_nuc_tot = Par.dim_nuc_tot_tmp  # full nuclear dimension
    else:
        Par.dim_nuc = np.zeros(Par.number_of_nuclei, dtype=np.int)
        for i in range(0, Par.number_of_nuclei):
            if Par.number_of_nuclei > 1:
                nucstring = Par.Nucs[i]  # each nuclei seperately
            else:
                nucstring = Par.Nucs
            Par.dim_nuc[i] = nucdic.nuclear_properties(nucstring)[1]
        Par.dim_nuc_tot = int(np.prod(Par.dim_nuc))  # full nuclear dimension
    return


def create_g_or_hyperfine_tensors(input_array, dim):
    def return_diags(vector):
        try:
            arraylen = len(vector)
        except TypeError:
            arraylen = 1
        if arraylen == 1:
            v1, v2, v3 = vector, vector, vector
        elif arraylen == 2:
            v1, v2, v3 = vector[0], vector[0], vector[1]
        else:
            v1, v2, v3 = vector[0], vector[1], vector[2]
        outpout_array = np.diag(np.array([v1, v2, v3]))
        return outpout_array

    outpout_array = np.zeros((dim, 3, 3))
    if dim == 1:
        outpout_array[0] = return_diags(input_array)
    else:
        for i in range(0, dim):
            outpout_array[i] = return_diags(input_array[i])
    return outpout_array


def threshold_settings(Par):
    if Par.nKnots < 8:
        Par.nKnots = 8
    if not hasattr(Par, 'LevelSelect'):
        Par.LevelSelect = 5e-5
    if not hasattr(Par, 'SepHilbertspace'):
        Par.SepHilbertspace = False
    if hasattr(Par, 'Population'):
        Par.ispopu = True
    else:
        Par.ispopu = False

    if hasattr(Par, 'Temperature'):
        Par.T = Par.Temperature
    else:
        Par.T = 300
    return


def define_hyperfines(Par, SimPar):
    if Par.A is not None:
        Par.Nucs = SimPar._Nucsvec
        get_equivalent_nucs(Par)
        get_number_of_nuclei(Par)
        if Par.number_of_nuclei > 1:
            Par.Nucs = Par.Nucs.split(',')
        else:
            Par.Nucs = Par.Nucs
        get_full_nuclear_dimension(Par)
        if Par.dim_nuc_tot*Par.e_dimension > 512:
            Par.warning = 1
            return
        Par.A_tensor = create_g_or_hyperfine_tensors(Par.A,
                                                     Par.number_of_nuclei)
        Par.A_tensor = Par.A_tensor*con.MHz2Hz
    else:
        Par.dim_nuc_tot = 1
        Par.number_of_nuclei = 0
    return


def get_equivalent_nucs(Par):
    try:
        maxn = max(Par.n)
    except TypeError:
        maxn = Par.n
    if maxn > 1:
        if isinstance(Par.A, float) or isinstance(Par.A, int):
            Par.A = [Par.A]
        try:
            ntenss = len(Par.n)
        except TypeError:
            ntenss = 1
            Par.n = [Par.n]
            Par.A = [Par.A]
            if hasattr(Par, 'AFrame'):
                Par.AFrame = [Par.AFrame]
        if hasattr(Par, 'AFrame'):
            AFrame_tmp = []
        Par.A_tmp = []
        for i in range(0, ntenss):
            for k in range(Par.n[i]):
                Par.A_tmp.append(Par.A[i])
                if hasattr(Par, 'AFrame'):
                    AFrame_tmp.append(Par.AFrame[i])
        Par.A = Par.A_tmp
        if hasattr(Par, 'AFrame'):
            Par.AFrame = AFrame_tmp
    return


def define_D_tensors(Par):
    def arrange_D(D):
        if float(abs(D[0])) < 0.01:
            D[0] = 0.01
        DT = np.array([[-1/3*D[0]+D[1], 0, 0], [0, -1/3*D[0]-D[1], 0],
                      [0, 0, 2/3*D[0]]])
        DT = DT*con.MHz2Hz
        return DT

    if Par.D is not None:
        Par.D = np.asarray(Par.D)*1.0
        Par.D_tensor = arrange_D(Par.D)
    # D-Tensor _pair
    if Par.DPair is not None:
        Par.DPair = Par.DPair*1.0
        Par.DPair_tensor = arrange_D(Par.DPair)
    return


def define_J(Par):
    if Par.J is not None:
        Par.J_tensor = Par.J*1e6*np.eye(3)
    return


def symmetry(Par):
    Symmetry_Group(Par)
    if Par.Point_Group == "O3":
        Par.nKnots = 4
    Par.phinKnots = Par.nOctants*Par.nKnots
    return


def rotate_coordsystems(Par):
    if hasattr(Par, 'gFrame'):
        eulermatrix = tool.Eulermatrix(Par.gFrame[2], Par.gFrame[1],
                                       Par.gFrame[0])
        Par.g_tensor = tool.tensor_rotation(Par.g_tensor, eulermatrix)
    if hasattr(Par, 'AFrame'):
        if Par.number_of_nuclei == 1:
                eulermatrix = tool.Eulermatrix(Par.AFrame[2], Par.AFrame[1],
                                               Par.AFrame[0])
                Par.A_tensor = tool.tensor_rotation(Par.A_tensor, eulermatrix)
        else:
            Par.AFrame = np.asarray(Par.AFrame)
            for i in range(Par.number_of_nuclei):
                eulermatrix = tool.Eulermatrix(Par.AFrame[i, 2],
                                               Par.AFrame[i, 1],
                                               Par.AFrame[i, 0])
                Par.A_tensor[i, :] = tool.tensor_rotation(Par.A_tensor[i, :],
                                                          eulermatrix)
    if hasattr(Par, 'DFrame'):
        eulermatrix = tool.Eulermatrix(Par.DFrame[2], Par.DFrame[1],
                                       Par.DFrame[0])
        Par.D_tensor = tool.tensor_rotation(Par.D_tensor, eulermatrix)
    return


def zero_field_populations(Par):
    # Set up spinpolarization for electron sublevels
    def renormalize(Pop):
        Pop = Pop/sum(Pop)
        Pop = 1.0-Pop
        Pop = Pop/sum(Pop)
        return Pop

    def triplet(Par):
        Par.Triplet = True
        Par.Population_t = np.zeros(3)
        Par.Population_t[0:3] = Par.Population[1:4]
        Par.Population_t = renormalize(Par.Population_t)
        return

    if hasattr(Par, 'Population'):
        if isinstance(Par.Population, str):
            if Par.Population == "S" or Par.Population == "s":
                Par.Singlet = True
            elif Par.Population == "T" or Par.Population == "t":
                    Par.Population = ["T", 1/3.0, 1/3.0, 1/3.0]
                    triplet(Par)
        elif isinstance(Par.Population, list):
            if isinstance(Par.Population[0], str):
                if Par.Population[0] == "T" or Par.Population[0] == "t":
                    triplet(Par)
            else:
                Par.Population = np.asarray(Par.Population)
                Par.Singlet = False
                Par.Population = renormalize(Par.Population)
    return
