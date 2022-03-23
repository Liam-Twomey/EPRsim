# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 14:35:34 2017

(c) Stephan Rein, 29.11.2017
"""

# External libraries
import numpy as np
import math as math
from scipy import linalg as LAS
from EPRsim.Pauli_generators import (create_Pauli_matrices,
                                     create_Pauli_matrices_Nuc,
                                     create_seperate_Pauli_matrices_Nuc)
import EPRsim.Tools as tool
import EPRsim.Nucdic as Nucdic
try:
    from numba import complex64, int32, jit
    Numba = 1
except ImportError:
    Numba = 0


def dec_eigvector_phase():
    """
    Creates a decorator for the function preselect_off_res(). If Numba
    is available, the decorator is a Numba jit (just-in-time compilation),
    otherwise the identity decorator is called (this decorator does nothing).
    The jit decorator is called with defined static types to improve the
    performance. Cashing is enabled to avoid avoid compilation each time
    the Python program is invoked.
    """
    if Numba == 1:
        return (jit(complex64[:, :](complex64[:, :], int32),
                nopython=True, cache=True))
    else:
        return dec_identity


def dec_identity(ob):
    """
    Identity decorator. If Numba is not available, this decorator is used. The
    decorator does nothing.
    """
    return ob


def create_Bilinear_Hamiltonian(Tensor, S, I, eulermatrix):
    """
    input: Tensor,Par.Pauli,Par.I, phi, theta, psi = 0

    output: Ham

    Algorithm:
    Create an arbitrary bilinear interaction Hamiltonian between the spin
    vector Par.Pauli and Par.I and the interaction matrix Tensor.
    Rotates tensor for a given theta, phi and psi angle (default is zero).
    T' = Tensor(phi,theta,psi)
    H = S*T'*Par.I
    Returns the bilinear Zeeman Hamiltonian.

    (c) Stephan Rein, 31.10.2017
    """
    # Tensor rotation
    Tensor_rot = tool.tensor_rotation(Tensor, eulermatrix)
    Ham = 0
    # Set up the full bilinear Hamiltonian
    for i in range(0, 3):
        for j in range(0, 3):
            Ham = np.add(Ham, np.dot(S[i], Tensor_rot[i, j]*I[j]))
    return Ham


def create_linear_Hamiltonian(Tensor, S, eulermatrix, z=True):
    """
    input: Tensor, Par.Pauli, phi, theta, z = True, psi = 0

    output: Ham

    Algorithm:
    Create an linear interaction Hamiltonian (usually the Zeeman interaction)
    between the spin vector Par.Pauli and the interaction matrix Tensor and a
    vector with 2 cartesian orientations set to zero.
    This can be either the x,y components (the z-component remains) or
    the y,z components (the x-component remains). The unecessary matrix
    multiplications are spared compared to a bilinear interaction.
    Rotates tensor for a given theta, phi and psi angle (default is zero).
    T' = Tensor(phi,theta,psi)
    H = S*T'*u_i
    with u_i the i-th component of a canoncial vector (either x or z)
    Returns the linear Zeeman Hamiltonian for a unit vector on one side.

    (c) Stephan Rein, 31.10.2017
    """
    # Tensor rotation
    Tensor_rot = tool.tensor_rotation(Tensor, eulermatrix)
    if z:
        # Set up the Hamiltonian for perpendicular magnetic Par.field
        Ham = S[0]*Tensor_rot[0, 2]+S[1]*Tensor_rot[1, 2]+S[2]*Tensor_rot[2, 2]
    else:
        # Set up the Hamiltonian for parallel magnetic Par.field
        Ham = S[0]*Tensor_rot[0, 0]+S[1]*Tensor_rot[1, 0]+S[2]*Tensor_rot[2, 0]
    return Ham


def create_Isotrop_Hamiltonian(S, sc):
    """
    input: S, zeeman

    output: Ham

    Algorithm:
    Create an linear interaction Hamiltonian with a scalar Par.field sc.
    Of course no rotations are necessary.
    H = S*sc
    (c) Stephan Rein, 31.10.2017
    """
    # Scalar multiplication with the Pauli matrix
    Ham_zee = S*sc
    return Ham_zee

def define_nKnots_pattern(Par):
    """
    input: Par,Par

    output: Par.nKnots_theta_vec

    Algorithm:
    Creates an integer rounded sinusoidal grid for a given number of theta
    steps. This depends additionally from the symmetry factor of the
    Hamiltonian. It can be one octant, two or the full hemisphere. This
    is defined by the Par.phinKnots which can be 1 (one octant), 2 or
    4 (full hemisphere).
    The grid is set up with the number of phi values n for a given theta
    (here on the exapmle of Par.phinKnots = 1)
    n = round(sin(theta)+0.5)
    The +0.5 forces the round function to always round up.
    The mimimal number for phi values for a given theta is 4. This ensures
    a later stable cubic spline where 4 points are necessary.
    The maximal number is the number of theta values.
    Returns an integer vector with the numbers of phi values which should be
    calculated for a given theta value.

    (c) Stephan Rein, 31.10.2017
    """
    sinvector = np.zeros(Par.nKnots, dtype=int)
    sinvector[0] = 2
    for i in range(1, Par.nKnots):
        theta = (math.pi/2)*((i)/(Par.nKnots-1))
        sinvector[i] = np.int(round(math.sin(theta)*Par.phinKnots+0.4999))
        if sinvector[i] < 4:
            sinvector[i] = 4
        if sinvector[i] > Par.phinKnots:
            sinvector[i] = Par.phinKnots
        if Par.Point_Group == "Dhinfty":
            sinvector[i] = 4
    Par.nKnots_theta_vec = sinvector
    return Par.nKnots_theta_vec


def define_I(Par, dimension):
    """
    Defines the Nuclear spin Pauli matrices.
    """
    if Par.A is not None:
        rand = np.linspace(0, 1000, dimension, endpoint=False)
        randmatrix = np.diag(rand)
        if hasattr(Par, 'ENucCoupling') and Par.SepHilbertspace:
            Par.I = create_seperate_Pauli_matrices_Nuc(Par, ZE=False)
        else:
            Par.I = create_Pauli_matrices_Nuc(Par)
    else:
        randmatrix = 0
    return randmatrix


def ZFS_Hamiltonian(Par):
    """  zero_field_eigendecomposition_and_Hamiltonian_setup(Par, Par)

    input: Par, Par

    output: Par.eigvec , Par.Ham_ZFS, Par.Ham_FD

    Algorithm:
    Creates the Hamiltonian for a given set of orientations M
    H(M) = H_EZ_k(M) + H_DD(M) + H_HF_ki(M) + H_NZ_i
    for all k electron zeeman interactions EZ. All dipolar couplings DD,
    all hyperfine interations HF and all nuclear zeeman interactions NZ. The
    nuclear zeeman interaction is treated isotropically. The gyromagnetic
    ratio is taken from the library in Nuc_zeemandict().
    All Par.field dependned interactions are multiplied with a unit vector
    instead of the magnetic Par.field vector.
    The hamiltonian is given back as two independent Hamiltonians:
    H1 = H_ZFS and H2 = H_FD
    Where H_ZFS is Par.field independent and H_FD is Par.field dependent.
    Both Hamiltonians are returned independently.
    Additionlly the zero-Par.field eigenvectors are returned. They are zero
    for the Boltzman case as they are not calculated. They are fully
    calculated if spin-polarization is switched on as they are later used
    for the pojection of the population at zero-Par.field to high Par.field.
    The diagonalization at the zero-fied hamiltonian are only carried out if
    spin-polarization is enabled.

    (c) Stephan Rein, 31.10.2017
    """
    # Conversion constants
    mT2Hzperg = 1e9*tool.mT2GHz(1.0, 1.0)
    # Allocations
    Par.Pauli = create_Pauli_matrices(Par)
    dimension = len(Par.Pauli[0, 0, :])
    # Interval for phi due to symmetry of the Hamiltonian
    phiKnots = Par.phinKnots
    # Define theta/phi grid
    Knots_theta_vec = define_nKnots_pattern(Par)
    # Allocation of the Hamiltonian
    Par.Ham_ZFS = np.zeros((1, Par.nKnots, phiKnots, dimension, dimension),
                           dtype=np.complex64, order='C')
    Ham = np.zeros((dimension, dimension), dtype=np.complex64, order='C')
    Par.Ham_FD = np.zeros((1, Par.nKnots, phiKnots, dimension, dimension),
                          dtype=np.complex64, order='C')
    randmatrix = define_I(Par, dimension)
    # Run the Loop over theta and phi for setting up the Hamiltonian
    for k in range(0, Par.nKnots):
        theta = (math.pi/2)*((k)/(Par.nKnots-1))
        for q in range(0, Knots_theta_vec[k]):
            phi = (Par.nOctants*np.pi/2)*((q)/(Knots_theta_vec[k]*1.0-1))
            Ham = 0
            eulermatrix = tool.Eulermatrix(phi, theta)
            if Par.D is not None:
                for i in range(0, Par.coupled_e_dim):
                    Ham_tmp = create_Bilinear_Hamiltonian(Par.D_tensor,
                                                          Par.Pauli[i],
                                                          Par.Pauli[i],
                                                          eulermatrix)
                    Ham += Ham_tmp
            if Par.DPair is not None:
                Ham_tmp = create_Bilinear_Hamiltonian(Par.DPair_tensor,
                                                      Par.Pauli[0],
                                                      Par.Pauli[1],
                                                      eulermatrix)
                Ham += 2*Ham_tmp
            if Par.J is not None:
                Ham_tmp = create_Bilinear_Hamiltonian(Par.J_tensor,
                                                      Par.Pauli[0],
                                                      Par.Pauli[1],
                                                      eulermatrix)
                Ham -= Ham_tmp
            # Add the Zeeman interaction
            if Par.A is not None:
                for i in range(0, Par.number_of_nuclei):
                    if Par.number_of_nuclei > 1:
                        nucstring = Par.Nucs[i]  # each nuclei seperately
                    else:
                        nucstring = Par.Nucs
                    for s in range(0, Par.coupled_e_dim):
                        # Define coupling to different nuclei for electrons
                        if hasattr(Par, 'ENucCoupling'):
                            if Par.ENucCoupling[s, i]:
                                if Par.SepHilbertspace:
                                    arg = [Par.A_tensor[i], Par.Pauli[s],
                                           Par.I[s, i], eulermatrix]
                                    Ham_tmp = create_Bilinear_Hamiltonian(*arg)
                                    Ham += Ham_tmp
                                else:
                                    arg = [Par.A_tensor[i], Par.Pauli[s],
                                           Par.I[i], eulermatrix]
                                    Ham_tmp = create_Bilinear_Hamiltonian(*arg)
                                    Ham += Ham_tmp
                        else:
                            arg = [Par.A_tensor[i], Par.Pauli[s], Par.I[i],
                                   eulermatrix]
                            Ham_tmp = create_Bilinear_Hamiltonian(*arg)
                            Ham += Ham_tmp
                    zeeman = 1.0*Nucdic.nuclear_properties(nucstring)[0]
                    if Par.coupled_e_dim > 1:
                        pass
                    else:
                        Ham_tmp = create_Isotrop_Hamiltonian(Par.I[i, 2],
                                                             zeeman)
                        Par.Ham_FD[0, k, q, :, :] = (Par.Ham_FD[0, k, q, :, :]+
                                                     Ham_tmp)
            # Add the Zeeman interaction
            for i in range(0, Par.coupled_e_dim):
                args = [Par.g_tensor[i, :, :]*mT2Hzperg, Par.Pauli[i],
                        eulermatrix]
                Ham_tmp = create_linear_Hamiltonian(*args)
                Par.Ham_FD[0, k, q] = Par.Ham_FD[0, k, q]+Ham_tmp
            Ham += randmatrix
            Par.Ham_ZFS[0, k, q] = Ham
    zero_field_diag(Par, Knots_theta_vec, phiKnots)
    return


def zero_field_diag(Par, Knots_theta_vec, phiKnots):
    """
    input: Par, Knots_theta_vec, phiKnots

    The function diagonalizes the zero-field Hamiltonian to get the
    eigenvectors as projection operators for the zero-field density matrix.
    The function stores the zero-field density matrices as a list of
    zero-field matrices
    """
    rho_0 = np.zeros((1, Par.nKnots, phiKnots), dtype=object)
    # Do only calculate the zero Par.field eigenvectors for spin-polarization
    if Par.ispopu:
        rho_init, rho_0_tmp = set_up_density_mat(Par)
        # Run the Loop over theta and phi for setting up the ZF density matrix
        for k in range(0, Par.nKnots):
            theta = (math.pi/2)*((k)/(Par.nKnots-1))
            for q in range(0, Knots_theta_vec[k]):
                phi = (Par.nOctants*np.pi/2)*((q)/(Knots_theta_vec[k]*1.0-1))
                Ham = 0
                eulermatrix = tool.Eulermatrix(phi, theta)
                args = [None, Par.Pauli[0], Par.Pauli[0], eulermatrix]
                if not Par.Singlet and not Par.Triplet:
                    if Par.D is not None:
                        args[0] = Par.D_tensor
                        Ham_tmp = create_Bilinear_Hamiltonian(*args)
                        Ham += Ham_tmp
                    if Par.DPair is not None:
                        args[0] = Par.D_tensor
                        Ham_tmp = create_Bilinear_Hamiltonian(*args)
                        Ham += 2*Ham_tmp
                    if Par.J is not None:
                        args[0] = Par.J_tensor
                        Ham_tmp = create_Bilinear_Hamiltonian(*args)
                        Ham_tmp = create_Bilinear_Hamiltonian(*args)
                        Ham -= Ham_tmp
                    # Add the hyperfine interaction
                    if Par.A is not None:
                        for i in range(0, Par.number_of_nuclei):
                            for s in range(0, Par.coupled_e_dim):
                                # Coupling different nuclei for diff. electrons
                                if hasattr(Par, 'ENucCoupling'):
                                    if Par.ENucCoupling[s, i]:
                                        if Par.SepHilbertspace:
                                            ags = (Par.A_tensor[i],
                                                   Par.Pauli[s],
                                                   Par.I[s, i], eulermatrix)
                                        else:
                                            ags = (Par.A_tensor[i],
                                                   Par.Pauli[s],
                                                   Par.I[i], eulermatrix)
                                else:
                                    ags = (Par.A_tensor[i], Par.Pauli[s],
                                           Par.I[i], eulermatrix)
                                Ham_tmp = create_Bilinear_Hamiltonian(*ags)
                                Ham += Ham_tmp
                    w, v = LAS.eigh(Ham, turbo=True, check_finite=True)
                    w = np.real(w)
                    tmp = np.sum(v.real, axis=0)
                    s = np.where(tmp < 0)[0]
                    v[:, s] = -1*v[:, s]
                    rho_0[0, k, q] = v.dot((rho_0_tmp).dot((v.T).conj()))
                else:
                    rho_0[0, k, q] = rho_init
        Par.rho_0 = rho_0
    return


def set_up_density_mat(Par):
    """
    The function sets up a zero-field density matrix rho_init for
    either a spin-polarized radical pair (Par.Singlet), for spin-polarized
    radical pairs in the triplet state (Par.Triplet), or for arbitrary
    spin-temperature (provided as vector Par.Population).
    This density matrix is later rotated by the zero-field eigenvectors.
    """
    rho_init, rho_0_tmp1 = 0, 0
    if Par.ispopu and Par.Singlet:
        rho_init = np.zeros((4, 4))
        rho_init[1, 1] = 1/2.0
        rho_init[1, 2] = -1/2.0
        rho_init[2, 1] = -1/2.0
        rho_init[2, 2] = 1/2.0
        rho_init = (1.0/Par.dim_nuc_tot)*np.kron(rho_init,
                                                 np.eye(Par.dim_nuc_tot))
    elif Par.ispopu and Par.Triplet:
        rho_init = np.zeros((4, 4))
        rho_init[0, 0] = Par.Population_t[0]
        rho_init[3, 3] = Par.Population_t[2]
        rho_init[1, 2] = Par.Population_t[1]/2.0
        rho_init[2, 1] = Par.Population_t[1]/2.0
        rho_init[1, 1] = Par.Population_t[1]/2.0
        rho_init[2, 2] = Par.Population_t[1]/2.0
        rho_init = (1.0/Par.dim_nuc_tot)*np.kron(rho_init,
                                                 np.eye(Par.dim_nuc_tot))
    elif Par.ispopu and not Par.Singlet:
        rho_0_tmp1 = np.zeros((Par.e_dimension, Par.e_dimension))
        for i in range(0, Par.e_dimension):
            rho_0_tmp1[i, i] = Par.Population[i]
        rho_0_tmp1 = (1.0)*np.kron(rho_0_tmp1, np.eye(Par.dim_nuc_tot))
    return rho_init, rho_0_tmp1


def FD_diagonalization(nKnots, Knots_theta_vec, phiKnots, Ham_ZFS, field,
                       dimension, Ham_FD, ispopu=False):
    """
    Core diagonalization routine for a field-dependent Hamiltonian

    Parameters
    ----------
    nKnots :  :class:`int`
              Number of theta points

    Knots_theta_vec : :class:`np.array`
                      Vector with theta values

    phiKnots : :class:`int`
               Number of phi points

    Ham_ZFS : :class:`object`
              Multidimensional array with all field-inddependent Hamiltionian
              parts

    field : :class:`np.array`
            Array with field points

    dimension: :class:`int`
               Dimension of the eigenvalues per field point and orientation

    Ham_FD: :class:`np.array`
            Multidimensional array with all field-dependent Hamiltionian parts

    Returns
    -------
    eigvec :  :class:`np.array`
              Multidimensional arrray with eigenvectors

    eigval : :class:`object`
              Multidimensional arrray with eigenvalues

    Notes
    -----
    Algorithm:
    Carries out a diagonalization for given number of field points, for
    each orientation (theta, phi). Uses a  divide and conquer algorithm from
    Lapack for the hermitian eigendecomposition.
    The eigenvectors are rephased by a phase of :math:`\\pi` if the sum
    over the real part of the eigenvector is negative.
    Returns the eigenvalues and vectors for all provided field points and
    orientations.
    """
    # Allocations
    n_explicit = len(field)
    eigvec = np.zeros((dimension, nKnots, phiKnots, n_explicit, dimension),
                      dtype=np.complex64, order='C')
    eigval = np.zeros((dimension, nKnots, phiKnots, n_explicit),
                      dtype=np.float64, order='C')
    # Calculate eigendecomposition  (use divide and conquer)
    for k in range(0, nKnots):
        for q in range(0, Knots_theta_vec[k]):
            for m in range(0, n_explicit):
                Hamilonian = Ham_ZFS[k, q]+Ham_FD[k, q]*field[m]
                w, v = LAS.eigh(Hamilonian, overwrite_a=True,
                                turbo=True, check_finite=False)
                if ispopu:
                    eigvec[:, k, q, m, :] = phase_eigenvectors(v, dimension)
                else:
                    eigvec[:, k, q, m, :] = np.transpose(v)
                eigval[:, k, q, m] = w
    return eigvec, eigval


@dec_eigvector_phase()
def phase_eigenvectors(v, dimension):
    # Previous method for phase change in eigenvectors
    for j in range(0, dimension):
        s = np.where(np.absolute((v[:, j])) > 1e-8)[0]
        if np.abs(v[s[0], j].real) > np.abs(v[s[0], j].imag):
            if np.real(v[s[0], j]) < 0:
                v[:, j] = 1.0*v[:, j]
        else:
            if np.imag(v[s[0], j]) < 0:
                v[:, j] = -1.0*v[:, j]
    v = np.transpose(v)
    return v


def HF_Eig(Par):
    """
    Core algorithm for field bisection and eigenvalue determination

    Parameters
    ----------
    Par :     :class:`object`
              Object with all user-defined parameters.

    Notes
    -----
    This function is the core function for high-field diagonalization.
    A field bisection is carried out to use a sufficient number of field
    points but not more than necessary. The more detailed algorithm is:

        1. Carries out a diagonalization at 3 Par.field points using the
        subroutine FD_diagonalization()

        2. Bisectioning of the field if the linearity is not fulfilled
        (in general when the deviation from linearity falls beyond a given
        treshold) Carries out again diagonalizations by calling
        FD_diagonalization() for the new field points.

    The FD_diagonalization() is the core function for the diagonalization
    for one magnetic field position and all orientations.

    The resulting field points, used in the diagonalization are saved in
    the member variable Par.field. The eigenvector and eigenvalues
    are stored in the member variable Par.eigvec and Par.eigval, respectively.
    """

    # Create angle pattern
    Knots_theta_vec = define_nKnots_pattern(Par)
    # Allocations
    dimension = len(Par.Pauli[0, 0, :])
    Par.field = np.linspace(Par.Range[0], Par.Range[1], Par.n_explicit,
                            endpoint=True)
    # Initial calcualtion at three magnetic Par.field points
    Par.eigvec, Par.eigval = FD_diagonalization(Par.nKnots, Knots_theta_vec,
                                                Par.phinKnots, Par.Ham_ZFS[0],
                                                Par.field, dimension,
                                                Par.Ham_FD[0], Par.ispopu)
    eigaverage = (Par.eigval[:, :, :, 2]+Par.eigval[:, :, :, 0])*0.5
    eigdiffmax = np.amax(eigaverage-Par.eigval[:, :, :, 1])
    j = 1
    counter = 0
    threshold = 2e4*(Par.Range[1]-Par.Range[0])
    while eigdiffmax > threshold and counter < 512:
        eigdiffmax = 0
        for i in range(0, j, 2):
            # Make bisection if for one transition and one orientation
            # is larger than the threshold.
            eigaverage = (Par.eigval[:, :, :, i]+Par.eigval[:, :, :, i+2])*0.5
            eigdiff = np.amax(eigaverage-Par.eigval[:, :, :, i+1])
            if eigdiff > threshold:
                field_tmp1 = (Par.field[i]+Par.field[i+1])/2
                field_tmp2 = (Par.field[i+1]+Par.field[i+2])/2
                field_tmp = np.array([field_tmp1, field_tmp2])
                eigvec, eigval = FD_diagonalization(Par.nKnots,
                                                    Knots_theta_vec,
                                                    Par.phinKnots,
                                                    Par.Ham_ZFS[0],
                                                    field_tmp,
                                                    dimension,
                                                    Par.Ham_FD[0], Par.ispopu)
                eigdiff = np.amax(eigaverage-eigval[:, :, :, 0])
                Par.field = np.append(Par.field, field_tmp)
                Par.eigval = np.append(Par.eigval, eigval, axis=3)
                Par.eigvec = np.append(Par.eigvec, eigvec, axis=3)
                idx = np.argsort(Par.field, axis=-1, kind='quicksort')
                Par.field = Par.field[idx]
                Par.eigvec = Par.eigvec[:, :, :, idx]
                Par.eigval = Par.eigval[:, :, :, idx]
                eigaverage = (Par.eigval[:, :, :, i] +
                              Par.eigval[:, :, :, i+2])*0.5
                eigdiff = np.amax(eigaverage-Par.eigval[:, :, :, i+1])
                if eigdiff > eigdiffmax:
                    eigdiffmax = eigdiff
                eigaverage = (Par.eigval[:, :, :, i+2] +
                              Par.eigval[:, :, :, i+4])*0.5
                eigdiff = np.amax(eigaverage-Par.eigval[:, :, :, i+3])
                if eigdiff > eigdiffmax:
                    eigdiffmax = eigdiff
                Par.n_explicit = len(Par.field)
                j += 2
                counter += 1
    return
