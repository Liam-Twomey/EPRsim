# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 14:35:34 2017

@author: Stephan Rein, University of Freiburg, 2017
"""

import numpy as np


def check_tensor(tensor, ani, inft, iso):
    k = check_off_diag(tensor)
    if not k:
        ani = True
    else:
        q, p = check_diag(tensor)
        if not q:
            inft = False
        if not p:
            iso = False
    return ani, inft, iso


def mul_and_add_tensors(mul, add, tensor):
    mul = mul@tensor
    add += tensor
    return mul, add


def Symmetry_Group(Par):
    """
   Determin the symmetry group of a (pre-processed) Hamiltonian

    Parameters
    ----------
    Par :     :class:`object`
              Object with all user-defined parameters.


    See Also
    --------
    Parameters : Class of Parameters (with a full list of parameters)


    Notes
    -----
    Works fully on call on reference on the Par object. The function
    categrorized the Hamiltonian in different point groups. The symmetry
    group determines how many octants of a hemisphere are calculated
    and how many points, at maxmimum, of theta and phi are caluclated.

    ========  ================  ============   ====================
    sym_grou    Point group        Octants        max nphi/ntheta
    ========  ================  ============   ====================
     -1             O3                1               4/4
      0           Dhinfty             1          400/4  or  4/400
      1             D2h               1              400/400
      2             C2h               2              800/400
      4             Ci                4             1600/400
    ========  ================  ===========    ====================



    Examples
    --------

    Example for getting back the hyperfine tensors, out of the input
    information as well as the point group of the Hamiltonian

    >>> import numpy as np
    >>> import EPRsim.Presettings as pre
    >>> import EPRsim.Hamiltonian_Point_Group as gr
    >>> P = sim.Parameters()
    >>> P.Range = [335 ,350]
    >>> P.mwFreq = 9.6
    >>> P.g = 2.002
    >>> P.A = [12, 20, 120]
    >>> P.Nucs = '14N'
    >>> Val = sim.Validate_Parameters(P)
    >>> P.A_tensor = [np.diag(np.array([12, 20, 120]))]
    >>> P.g_tensor = [np.diag(np.array([2.002, 2.002, 2.002]))]
    >>> P.number_of_nuclei = 1
    >>> P.coupled_e_dim = 1
    >>> gr.Symmetry_Group(P)
    D2h
    """

    # Important Symmetry Point Groups
    sym_grou = {4: "Ci", 2: "C2h",  1: "D2h", 0: "Dhinfty", -1: "O3"}
    sym_num = {"Ci": 4, "C2h": 2, "D2h": 1, "Dhinfty": 1, "O3": 1}
    n = 0
    # Assume the best symmetry Group (O3) at the beginning.
    isD2h = True
    inft = True
    iso = True
    ani = False
    mul = np.eye(3)
    add = np.zeros((3, 3))
    if Par.D is not None:
        ani, inft, iso = check_tensor(Par.D_tensor, ani, inft, iso)
        mul, add = mul_and_add_tensors(mul, add, Par.D_tensor)
    if Par.DPair is not None:
        ani, inft, iso = check_tensor(Par.DPair_tensor, ani, inft, iso)
        mul, add = mul_and_add_tensors(mul, add, Par.DPair_tensor)
    # Hyperfine tensors
    if Par.A is not None:
        if Par.number_of_nuclei == 1:
            ani, inft, iso = check_tensor(Par.A_tensor[0], ani, inft, iso)
            mul, add = mul_and_add_tensors(mul, add, Par.A_tensor[0])
        else:
            for t in range(0, Par.number_of_nuclei):
                ani, inft, iso = check_tensor(Par.A_tensor[t], ani, inft, iso)
                mul, add = mul_and_add_tensors(mul, add, Par.A_tensor[t])
    # g tensors
    if Par.coupled_e_dim == 1:
        ani, inft, iso = check_tensor(Par.g_tensor[0], ani, inft, iso)
        mul, add = mul_and_add_tensors(mul, add, Par.g_tensor[0])
    else:
        ani, inft, iso = check_tensor(Par.g_tensor[0], ani, inft, iso)
        for t in range(0, Par.coupled_e_dim):
            ani, inft, iso = check_tensor(Par.g_tensor[t], ani, inft, iso)
            mul, add = mul_and_add_tensors(mul, add, Par.g_tensor[t])
    # Multiplication of the interactions tensors
    ani, inft, iso = check_tensor(mul, ani, inft, iso)
    ani, inft, iso = check_tensor(add, ani, inft, iso)
    # Evaluate final symmetry Group
    if isD2h:
        if inft:
            if iso:
                n = -1
            else:
                n = 0
        else:
            n = 1
    else:
        n = 2
    if ani:
        n = 4
    Par.Point_Group = sym_grou[n]
    Par.nOctants = sym_num[Par.Point_Group]
    # Immediatelly lift the degenraty of the D-tensor (break the symmetry)
    # at zero fuield. Lift it for 10 kHz.
    if Par.D is not None:
        if Par.D[1] == 0:
            Par.D_tensor[0, 0] += 10000
            Par.D_tensor[1, 1] -= 10000
    # Make permutation of the Hamiltonian for D2H:
    if Par.Point_Group == "Dhinfty":
        permutate_tensor_elements(Par, add)
    return


def check_off_diag(tensor):
    iszero = True
    for i in range(0, 2):
        for j in range(i+1, 3):
            if abs(tensor[i, j]) > 1e-8:
                iszero = False
    return iszero


def check_diag(tensor):
    issame = False
    allsame = False
    if (abs(tensor[0, 0]-tensor[1, 1]) < 1e-8 or
       abs(tensor[0, 0]-tensor[2, 2]) < 1e-8 or
       abs(tensor[2, 2]-tensor[1, 1]) < 1e-8):
        issame = True
    if (abs(tensor[0, 0]-tensor[1, 1]) < 1e-8 and
       abs(tensor[0, 0]-tensor[2, 2]) < 1e-8 and
       abs(tensor[2, 2]-tensor[1, 1]) < 1e-8):
        allsame = True
    return issame, allsame


def permutate_tensor_elements(Par, add):
    # Make permutation of the Hamiltonian for D2H:
    if add[0, 0] == add[1, 1]:
        transformat = np.eye(3)
    elif add[0, 0] == add[2, 2]:
        transformat = np.zeros((3, 3))
        transformat[0, 1] = 1
        transformat[1, 2] = 1
        transformat[2, 0] = 1
    else:
        transformat = np.zeros((3, 3))
        transformat[0, 2] = 1
        transformat[1, 0] = 1
        transformat[2, 1] = 1
    if Par.coupled_e_dim == 1:
        Par.g_tensor = Par.g_tensor@transformat
    else:
        for t in range(0, Par.coupled_e_dim):
            Par.g_tensor[t] = Par.g_tensor[t]@transformat
    if Par.D is not None:
        Par.D_tensor = Par.D_tensor@transformat
    if Par.DPair is not None:
        Par.DPair_tensor = Par.DPair_tensor@transformat
    # Hyperfine tensors
    if Par.A is not None:
        if Par.number_of_nuclei == 1:
            Par.A_tensor[0] = Par.A_tensor[0]@transformat
        else:
            for t in range(0, Par.number_of_nuclei):
                Par.A_tensor[t] = Par.A_tensor[t]@transformat
    return
