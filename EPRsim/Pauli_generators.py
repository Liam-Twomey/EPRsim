# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 10:05:21 2017

@author: Stephan
"""

import numpy as np


def Pauli_matrices(dim):
    """
    Pauli matrix generator


    Parameters
    ----------
    dim :  :class:`int`
           Dimension of the Pauli matrices (dim x dim).

    Returns
    -------
    j_x :  :class:`numpy.ndarray`
           Pauli matrix :math:`\\sigma_x`.

    j_y :  :class:`numpy.ndarray`
           Pauli matrix :math:`sigma_y`.

    j_z :  :class:`numpy.ndarray`
           Pauli matrix :math:`sigma_z`.

    Notes
    -------
    Calculates the three Pauli matrices (:math:`x, y, z`)
    for an arbitrary Hilbert
    space dimension. Note that the dimension for a specific spin quantum number
    :math:`S` is given as: dim = :math:`2S+1`.

    Examples
    --------
    Example for Pauli matrices of a spin :math:`S = 1` system.

    >>> import cwEPRsim.Useful_Spin_Mathematics as US
    >>> j_x, j_y, j_z = US.Pauli_matrices(3)
    >>> print(j_x)
    [[0.        +0.j 0.70710678+0.j 0.        +0.j]
    [0.70710678+0.j 0.        +0.j 0.70710678+0.j]
    [0.        +0.j 0.70710678+0.j 0.        +0.j]]
    >>> print(j_y)
    [[0.-0.j         0.-0.70710678j 0.-0.j        ]
    [0.+0.70710678j 0.-0.j         0.-0.70710678j]
    [0.-0.j         0.+0.70710678j 0.-0.j        ]]
    >>> print(j_z)
    [[ 1.+0.j  0.+0.j  0.+0.j]
    [ 0.+0.j  0.+0.j  0.+0.j]
    [ 0.+0.j  0.+0.j -1.+0.j]]
    """
    d = np.zeros(int(dim-1))
    j_x = np.zeros((dim, dim))
    for i in range(0, dim-1):
        d[i] = np.sqrt((i+1)*(dim-1-(i)))
    j_x = np.zeros((dim, dim), dtype=np.complex, order='C')
    j_y = np.zeros((dim, dim), dtype=np.complex, order='C')
    j_z = np.zeros((dim, dim), dtype=np.complex, order='C')
    for i in range(0, dim):
        j_z[i, i] = dim/2.0-i-0.5
        if i+1 <= dim-1:
            j_x[i, i+1] = 0.5*d[i]
            j_y[i, i+1] = 0.5*1j*d[i]
    j_x = j_x+np.transpose(np.conjugate(j_x))
    j_y = j_y+np.transpose(np.conjugate(j_y))
    return j_x, j_y, j_z


def create_Pauli_matrices(Par, ZE=False):
    S_p = np.zeros((Par.coupled_e_dim, 3, Par.e_dimension*Par.dim_nuc_tot,
                   Par.e_dimension*Par.dim_nuc_tot), dtype=np.complex64,
                   order='C')
    Par.S_tot = np.zeros((3, Par.e_dimension, Par.e_dimension),
                         dtype=np.complex64, order='C')
    Par.S_tot_without_nuc = 0
    Par.S_pure_ele = np.zeros((Par.coupled_e_dim, 3, Par.e_dimension,
                               Par.e_dimension), dtype=np.complex64, order='C')
    for i in range(0, Par.coupled_e_dim):
        if Par.coupled_e_dim == 1:
            i_x_tmp, i_y_tmp, i_z_tmp = Pauli_matrices(Par.e_dimension)
            Par.S_tot_without_nuc = np.array([i_x_tmp, i_y_tmp, i_z_tmp])
            S_p[0, 0] = np.kron(i_x_tmp, np.eye(Par.dim_nuc_tot))
            S_p[0, 1] = np.kron(i_y_tmp, np.eye(Par.dim_nuc_tot))
            S_p[0, 2] = np.kron(i_z_tmp, np.eye(Par.dim_nuc_tot))
            Par.S_pure_ele[0, 0] = i_x_tmp
            Par.S_pure_ele[0, 1] = i_y_tmp
            Par.S_pure_ele[0, 2] = i_z_tmp
            Par.S_tot = S_p[0, :, :]
        else:
            dim = int(Par.S[i]*2+1)
            i_x_tmp, i_y_tmp, i_z_tmp = Pauli_matrices(dim)
            for k in range(0, Par.coupled_e_dim):
                dim = int(Par.S[k]*2+1)
                if k != i:
                    if k < i:
                        i_x_tmp = np.kron(np.eye(dim), i_x_tmp)
                        i_y_tmp = np.kron(np.eye(dim), i_y_tmp)
                        i_z_tmp = np.kron(np.eye(dim), i_z_tmp)
                    else:
                        i_x_tmp = np.kron(i_x_tmp, np.eye(dim))
                        i_y_tmp = np.kron(i_y_tmp, np.eye(dim))
                        i_z_tmp = np.kron(i_z_tmp, np.eye(dim))
            Par.S_tot_without_nuc += np.array([i_x_tmp, i_y_tmp, i_z_tmp])
            Par.S_pure_ele[i, 0] = i_x_tmp
            Par.S_pure_ele[i, 1] = i_y_tmp
            Par.S_pure_ele[i, 2] = i_z_tmp
            S_p[i, 0] = np.kron(i_x_tmp, np.eye(Par.dim_nuc_tot))
            S_p[i, 1] = np.kron(i_y_tmp, np.eye(Par.dim_nuc_tot))
            S_p[i, 2] = np.kron(i_z_tmp, np.eye(Par.dim_nuc_tot))
            Par.S_tot = np.sum(S_p, axis=0)
    return S_p


# *****************************************************************************
# CURRENTLY NOT USED MODULE FOR SPIN CORRELATED RADICAL PAIRS
# *****************************************************************************
def create_Pauli_matrices_Nuc(Par, ZE=False):
    # reate Pauli matrices
    dim_tot = Par.e_dimension*Par.dim_nuc_tot
    i_p = np.zeros((Par.number_of_nuclei,3,dim_tot,dim_tot),
                   dtype=np.complex64,order='C')

    dim_electron = Par.e_dimension #electron spin dimension
    for i in range(0,Par.number_of_nuclei):
        i_x_tmp,i_y_tmp,i_z_tmp = Pauli_matrices(Par.dim_nuc[i])
        if Par.number_of_nuclei == 1:
            i_p[i,0] = np.kron(np.eye(dim_electron),i_x_tmp)
            i_p[i,1] = np.kron(np.eye(dim_electron),i_y_tmp)
            i_p[i,2] = np.kron(np.eye(dim_electron),i_z_tmp)
        else:
            for k in range(0,Par.number_of_nuclei):
                dim =  Par.dim_nuc[k]
                if k != i:
                    if k < i:
                        i_x_tmp = np.kron(np.eye(dim),i_x_tmp)
                        i_y_tmp = np.kron(np.eye(dim),i_y_tmp)
                        i_z_tmp = np.kron(np.eye(dim),i_z_tmp)
                    else:
                        i_x_tmp = np.kron(i_x_tmp,np.eye(dim))
                        i_y_tmp = np.kron(i_y_tmp,np.eye(dim))
                        i_z_tmp = np.kron(i_z_tmp,np.eye(dim))
            i_p[i,0] = np.kron(np.eye(dim_electron),i_x_tmp)
            i_p[i,1] = np.kron(np.eye(dim_electron),i_y_tmp)
            i_p[i,2] = np.kron(np.eye(dim_electron),i_z_tmp)
    return i_p


def create_seperate_Pauli_matrices_Nuc(Par, ZE = False):
    #Create Pauli matrices
    dim_tot = Par.e_dimension*Par.dim_nuc_tot
    i_p = np.zeros((Par.coupled_e_dim,Par.number_of_nuclei,3,dim_tot,dim_tot),
                   dtype=np.complex64,order='C')

    dim_electron = Par.e_dimension #electron spin dimension
    for s in range(0,Par.coupled_e_dim):
        for i in range(0,Par.number_of_nuclei):
            dim_tmp = Par.dim_nuc[s,i]
            if dim_tmp > 1:
                i_x_tmp,i_y_tmp,i_z_tmp = Pauli_matrices(Par.dim_nuc[s,i])
                if Par.number_of_nuclei == 1:
                    i_p[i,0] = np.kron(np.eye(dim_electron),i_x_tmp)
                    i_p[i,1] = np.kron(np.eye(dim_electron),i_y_tmp)
                    i_p[i,2] = np.kron(np.eye(dim_electron),i_z_tmp)
                else:
                    for k in range(0,Par.number_of_nuclei):
                        dim =  Par.dim_nuc[s,k]
                        if k != i:
                            if k < i:
                                i_x_tmp = np.kron(np.eye(dim),i_x_tmp)
                                i_y_tmp = np.kron(np.eye(dim),i_y_tmp)
                                i_z_tmp = np.kron(np.eye(dim),i_z_tmp)
                            else:
                                i_x_tmp = np.kron(i_x_tmp,np.eye(dim))
                                i_y_tmp = np.kron(i_y_tmp,np.eye(dim))
                                i_z_tmp = np.kron(i_z_tmp,np.eye(dim))
                    if len(i_z_tmp[:,0]) != Par.dim_nuc_tot:
                        dimfac = Par.dim_nuc_tot//len(i_z_tmp[:,0])
                    else:
                        dimfac = 1
                    i_p[s,i,0] = np.kron(np.eye(dim_electron*dimfac),i_x_tmp)
                    i_p[s,i,1] = np.kron(np.eye(dim_electron*dimfac),i_y_tmp)
                    i_p[s,i,2] = np.kron(np.eye(dim_electron*dimfac),i_z_tmp)
    return i_p