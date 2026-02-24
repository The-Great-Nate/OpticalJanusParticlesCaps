#
# Hydrodynamic things
#

import numpy as np
import itertools as it
from scipy.spatial.transform import Rotation
from pyoptics import constants as ct


def func3(a, r):
    return ((a ** 2) / (r ** 2)) + 1

def func4(a, b, r):
    return (a * b) / (r ** 2)

def Djj(particle_radius):  # For Diffusion
    djj = (ct.k_B * ct.temperature) / (6 * np.pi * ct.viscosity * particle_radius)
    D = np.zeros([3, 3])
    np.fill_diagonal(D, djj)
    return D

def Djk(x, y, z, r):
    D = np.zeros([3, 3])
    D[0][0] = func3(x, r)
    D[1][1] = func3(y, r)
    D[2][2] = func3(z, r)
    D[0][1] = func4(x, y, r)
    D[1][0] = D[0][1]
    D[0][2] = func4(x, z, r)
    D[2][0] = D[0][2]
    D[1][2] = func4(y, z, r)
    D[2][1] = D[1][2]
    return ((ct.k_B * ct.temperature) / (8 * np.pi * ct.viscosity * r)) * D

def Djk_rp(rvec, r, particle_radius):
    #
    # Expression from Ermak and McCammon
    #
    factor = 2*(particle_radius/r)**2
    D = np.eye(3)*(1+factor/3) + np.outer(rvec,rvec)*(1-factor)/r**2
    #
    return ((ct.k_B * ct.temperature) / (8 * np.pi * ct.viscosity * r)) * D


def oseen_tensor(array_of_positions,number_of_particles,particle_radius):
    #
    # D_matrix is the Dij referred to in Ermak and McCammon.  The Djk_array and
    # Djj_array respectively are for the off-diagonal and on-diagonal terms.
    #
    D_matrix = np.zeros([number_of_particles, number_of_particles], dtype=object)
    Djj_array = Djj(particle_radius)
    for i in range(number_of_particles):
        D_matrix[i,i] = Djj_array
        for j in range(i+1,number_of_particles):
            displacement = array_of_positions[j]-array_of_positions[i]
            distance = np.linalg.norm(displacement)
            Djk_array = Djk(
                displacement[0],
                displacement[1],
                displacement[2],
                distance)
            D_matrix[i,j] = Djk_array
            D_matrix[j,i] = Djk_array
    return D_matrix


def rp_tensor(array_of_positions,number_of_particles,particle_radius):
    #
    # Rotne-Prager tensor.
    #
    # D_matrix is the Dij referred to in Ermak and McCammon.  The Djk_array and
    # Djj_array respectively are for the off-diagonal and on-diagonal terms.
    #
    D_matrix = np.zeros([number_of_particles, number_of_particles], dtype=object)
    Djj_array = Djj(particle_radius)
    for i in range(number_of_particles):
        D_matrix[i,i] = Djj_array
        for j in range(i+1,number_of_particles):
            displacement = array_of_positions[j]-array_of_positions[i]
            distance = np.linalg.norm(displacement)
            Djk_array = Djk_rp(displacement, distance, particle_radius)
            D_matrix[i,j] = Djk_array
            D_matrix[j,i] = Djk_array
    return D_matrix


def blake_tensor(array_of_positions,number_of_particles,particle_radius, height=0.0):
    #
    # Rotne-Prager-Blake tensor.
    # As given in Yann von Hansen, Michael Hinczewski and Roland R. Netz,
    # J. Chem. Phys. 134, 235102 (2011).
    # The Djk_array and Djj_array respectively are for the off-diagonal
    # and on-diagonal terms.
    #
    # height is the z position of the wall
    #
    prefix = (ct.k_B * ct.temperature) / (4 * np.pi * ct.viscosity)
    a2 = particle_radius**2
    D_matrix = np.zeros([number_of_particles, number_of_particles], dtype=object)
    Djj_array = Djj(particle_radius)
    for i in range(number_of_particles):
        #
        # Self mobility from approx in eqs (8) & (9)
        #
        Djj_temp = np.zeros((3,3))
        z = abs(array_of_positions[i,2]-height) # separation z to the wall;
        az = particle_radius/z  # Particle CANNOT coincide with position of wall;
        az3 = az**3   # will need additional potential to keep particle from wall.
        Djj_temp[0,0] = Djj_array[0,0]*(1-9*az/16+az3/8)
        Djj_temp[1,1] = Djj_temp[0,0]
        Djj_temp[2,2] = Djj_array[2,2]*(1-9*az/8+az3/2)
        D_matrix[i,i] = Djj_temp
        #
        # Off diagonal terms from eqs A2 for RP and A4 to A8 for the Blake terms
        #
        for j in range(number_of_particles):
            if j != i:
                #
                # First the two Rotne-Prager terms
                #
                displacement = array_of_positions[i]-array_of_positions[j]
                jtemp = np.asarray((array_of_positions[j,0],array_of_positions[j,1],2*height-array_of_positions[j,2]))
                disp_image = array_of_positions[i]-jtemp
                distance = np.linalg.norm(displacement)
                dist_image = np.linalg.norm(disp_image)
                Djk_array = Djk_rp(displacement, distance, particle_radius) - Djk_rp(disp_image, dist_image, particle_radius)
                #
                # Next the Blake terms, each one is different
                #
                Rx = disp_image[0]
                Ry = disp_image[1]
                Rz = disp_image[2]
                Rx2 = Rx*Rx
                Ry2 = Ry*Ry
                Rz2 = Rz*Rz
                R = dist_image
                R2 = R*R
                R3 = R2*R
                R5 = R3*R2
                R7 = R5*R2
                zi = abs(array_of_positions[i,2]-height)
                zj = abs(array_of_positions[j,2]-height)
                #
                Djk_array[0,0] += prefix*(-(zi*zj/R3)*(1-3*Rx2/R2) + (a2*Rz2/R5)*(1-5*Rx2/R2))
                Djk_array[1,1] += prefix*(-(zi*zj/R3)*(1-3*Ry2/R2) + (a2*Rz2/R5)*(1-5*Ry2/R2))
                Djk_array[2,2] += prefix*( (zi*zj/R3)*(1-3*Rz2/R2) - (a2*Rz2/R5)*(3-5*Rz2/R2))
                #
                tempxy = prefix*(3*zi*zj*Rx*Ry/R5 - 5*a2*Rx*Ry*Rz2/R7)
                Djk_array[0,1] += tempxy
                Djk_array[1,0] += tempxy
                #
                Djk_array[0,2] += prefix*((zj*Rx/R3)*(1-3*zi*Rz/R2) - (a2*Rx*Rz/R5)*(2-5*Rz2/R2))
                Djk_array[1,2] += prefix*((zj*Ry/R3)*(1-3*zi*Rz/R2) - (a2*Ry*Rz/R5)*(2-5*Rz2/R2))
                #
                Djk_array[2,0] += prefix*((zj*Rx/R3)*(1+3*zi*Rz/R2) - (5*a2*Rx*Rz*Rz2/R7))
                Djk_array[2,1] += prefix*((zj*Ry/R3)*(1+3*zi*Rz/R2) - (5*a2*Ry*Rz*Rz2/R7))
                #
                D_matrix[i,j] = Djk_array
    return D_matrix


def diffusion_matrix(array_of_positions, particle_radius, tensor_choice='OSEEN', wall_height=0.0):
    """
    This version designed to return either the Oseen matrix, the Rotne-Prager matrix or the Rotne-Prager-Blake matrix.  The Rotne-Prager-Blake matrix is specifically for use with a wall perpendicular to the z-axis.
    Inputs:
        array_of_positions (number_of_particles,3)
        particle_radius
    The assumption in all the methods is that all particles have the same radius.
    Outputs:
        Diffusion matrix in two formats:
            D(3*number_of_particles,3*number_of_particles): for use with the force vector to generate the set of displacements.
            D_matrix(number_of_particles,number_of_particles): where each entry is a 3x3 tensor, which is used in the constraints calculation.
    """
    #
    # It's worth noting the 3x3 matrices are symmetric, as is the final diffusion matrix.
    # The concatenation at the end appears to rely on this.
    # NOTE: This will be a problem for the Blake tensor which is not symmetric.
    #
    number_of_particles = len(array_of_positions)
    #
    # This is where we need to choose between methods to fill the diffusion matrix.
    # These functions return an (NxN) array of (3x3) matrices.
    #
    if tensor_choice == 'OSEEN':
        D_matrix = oseen_tensor(array_of_positions,number_of_particles,particle_radius)
    elif tensor_choice == 'RP':
        D_matrix = rp_tensor(array_of_positions,number_of_particles,particle_radius)
    elif tensor_choice == 'BLAKE':
        D_matrix = blake_tensor(array_of_positions,number_of_particles,particle_radius,wall_height)
    else:
        print("Hydrodynamic tensor choice not recognised")
    #
    # Flatten the diffusion matrix into (3Nx3N) for use in random number generator
    # and elsewhere.
    #
    temporary_array = np.zeros(number_of_particles, dtype=object)
    for i in range(number_of_particles):
        temporary_array[i] = np.concatenate(D_matrix[i])
    D = np.hstack(temporary_array)
    return D,D_matrix


def diffusion_matrix_OLD(array_of_positions, particle_radius):
    """
    This version designed to return either the Oseen matrix, the Rotne-Prager matrix or the Rotne-Prager-Blake matrix.  The Rotne-Prager-Blake matrix is specifically for use with a wall perpendicular to the z-axis.
    Inputs:
    array_of_positions (number_of_particles,3)
    particle_radius
    The assumption in all the methods is that all particles have the same radius.
    Outputs:
    Diffusion matrix in two formats:
    D(3*number_of_particles,3*number_of_particles): for use with the force vector to generate the set of displacements.
    D_matrix(number_of_particles,number_of_particles): where each entry is a 3x3 tensor, which is used in the constraints calculation.
    """
    #
    # This is a complicated way of making a list of displacements for use in the
    # mobility tensor calculation.  The method relies on knowing the order the iterator
    # will produce the list so they can be mapped into the diffusion matrix.  If I were
    # writing this I would just use two nested loops.  It's worth noting the 3x3
    # matrices are symmetric, as is the final diffusion matrix.  The concatenation at
    # the end appears to rely on this.
    #
    list_of_displacements = [u - v for u, v in it.combinations(array_of_positions, 2)]
    array_of_displacements = np.zeros(len(list_of_displacements), dtype=object)
    for i in range(len(list_of_displacements)):
        array_of_displacements[i] = list_of_displacements[i]
    array_of_distances = np.array([np.linalg.norm(w) for w in array_of_displacements])
    number_of_particles = len(array_of_positions)
    number_of_displacements = len(array_of_displacements)
    #
    # D_matrix is the Dij referred to in Ermak and McCammon.  The Djk_array and
    # Djj_array respectively are lists for the off-diagonal and on-diagonal terms.
    #
    D_matrix = np.zeros([number_of_particles, number_of_particles], dtype=object)
    Djk_array = np.zeros(
        number_of_displacements, dtype=object
    )  # initialize array to store D_jk matrices
    Djj_array = np.zeros(
        number_of_particles, dtype=object
    )  # initialize array to store D_jj matrices
    iu = np.triu_indices(number_of_particles, 1)
    di = np.diag_indices(number_of_particles)
    #
    # This is where we need to choose between methods to fill the diffusion matrix.
    #
    for i in range(number_of_displacements):
        Djk_array[i] = Djk(
            array_of_displacements[i][0],
            array_of_displacements[i][1],
            array_of_displacements[i][2],
            array_of_distances[i],
        )
    for i in range(number_of_particles):
        Djj_array[i] = Djj(particle_radius)
    #
    # Fill in and symmetrize the resulting diffusion array, and perform the
    # necessary concatenations, and return.
    #
    D_matrix[iu] = Djk_array
    D_matrix.T[iu] = D_matrix[iu]
    D_matrix[di] = Djj_array
    temporary_array = np.zeros(number_of_particles, dtype=object)
    for i in range(number_of_particles):
        temporary_array[i] = np.concatenate(D_matrix[i])
    D = np.hstack(temporary_array)
    return D,D_matrix

def grand_diffusion_matrix(array_of_positions,
                           particle_radius,
                           tensor_choice='RP',
                           wall_height=0.0):
    if tensor_choice.upper() == 'BLAKE':
        raise ValueError("Coupled BD not compatable with walls.")
    number_of_particles = len(array_of_positions)
    N = number_of_particles
    kT = ct.k_B * ct.temperature
    eta = ct.viscosity
    a = particle_radius
    #TT same as non-grnad
    D_TT_flat, D_TT_blocks = diffusion_matrix(
        array_of_positions,
        particle_radius,
        tensor_choice=tensor_choice,
        wall_height=wall_height
    )

    D6 = np.zeros((6*N, 6*N), dtype=float)
    D6[0:3*N, 0:3*N] = D_TT_flat

    for i in range(N):

        D_RR_self = (kT / (8.0 * np.pi * eta * a**3)) * np.eye(3)

        row_rot_i = 3*N + 3*i
        D6[row_rot_i:row_rot_i+3,
           row_rot_i:row_rot_i+3] = D_RR_self

        for j in range(i+1, N):

            rvec = array_of_positions[j] - array_of_positions[i]
            r = np.linalg.norm(rvec)

            #we were getting division by 0 errors that hadnt shown up in the translational hydro, i suspect this is because of the 1/r^3 terms. There may be a better method than this
            r_min = 2.0 * a
            if r < r_min:
                r = r_min

            rhat = rvec / np.linalg.norm(rvec)  

            #RR
            pref_RR = kT / (16.0 * np.pi * eta * r**3)
            D_RR = pref_RR * (np.eye(3) - 3.0 * np.outer(rhat, rhat))

            row_rot_j = 3*N + 3*j

            D6[row_rot_i:row_rot_i+3,
               row_rot_j:row_rot_j+3] = D_RR

            D6[row_rot_j:row_rot_j+3,
               row_rot_i:row_rot_i+3] = D_RR

            #TR
            C = np.array([[0,        -rhat[2],  rhat[1]],
                          [rhat[2],   0,       -rhat[0]],
                          [-rhat[1],  rhat[0],  0       ]])

            pref_TR = kT / (8.0 * np.pi * eta * r**2)

            D_TR = pref_TR * C
            #RT should just be TR transposed
            D_RT = D_TR.T

            row_trans_i = 3*i
            row_trans_j = 3*j

            D6[row_trans_i:row_trans_i+3,
               row_rot_j:row_rot_j+3] = D_TR

            D6[row_trans_j:row_trans_j+3,
               row_rot_i:row_rot_i+3] = -D_TR

            D6[row_rot_j:row_rot_j+3,
               row_trans_i:row_trans_i+3] = D_RT

            D6[row_rot_i:row_rot_i+3,
               row_trans_j:row_trans_j+3] = -D_RT

    D6 = 0.5 * (D6 + D6.T)

    return D6


def coupled_bd_stepper(position_vectors,
                             array_of_rotated_dipoles,
                             total_force_array,
                             total_torques,
                             particle_radius,
                             timestep,
                             tensor_choice='RP',
                             brownian=False):
    
    N = len(position_vectors)
    if tensor_choice.upper() == 'BLAKE':
        raise ValueError("Coupled BD not compatable with walls.")

    D6 = grand_diffusion_matrix(position_vectors, particle_radius,
                                tensor_choice=tensor_choice,
                                wall_height=0.0)   #cant actually use walls but diffusion_matrix expects some value

    

    F_trans = np.hstack(total_force_array).astype(float)   
    T_rot   = np.hstack(total_torques).astype(float)       
    Fg = np.concatenate([F_trans, T_rot])                  

    kT = ct.k_B * ct.temperature
    sumDijFj = (1.0 / kT) * (D6 @ Fg)   
    drift = sumDijFj * timestep

    if brownian:

        cov = 2.0 * D6 * timestep

        #this was my first attempt at stopping the division by 0 at low seperations 
        eps = 1e-18 * np.trace(cov) / cov.shape[0]
        cov += eps * np.eye(6 * N)

        L = np.linalg.cholesky(cov)
        xi = np.random.normal(size=6 * N)

        R = L @ xi

    else:
        R = np.zeros(6 * N)

    delta_y = drift + R  

    delta_trans = delta_y[:3 * N].reshape((N, 3))
    new_positions = position_vectors + delta_trans


    for p in range(N):
        rotvec = delta_y[3 * N + 3 * p : 3 * N + 3 * p + 3]
        rot = Rotation.from_rotvec(rotvec)
        array_of_rotated_dipoles[p] = rot.apply(array_of_rotated_dipoles[p])

    return new_positions, array_of_rotated_dipoles

def trans_bd_hi(position_vectors, radius, total_force_array, number_of_particles, timestep, tensor_choice='OSEEN', wall_height=0.0):
    """
    'trans_bd_hi'
    Performs one time step of translational Brownian dynamics with HI.
    """
    mean = np.zeros(number_of_particles * 3)
    D,_ = diffusion_matrix(position_vectors, radius, tensor_choice, wall_height)
    F = np.hstack(total_force_array)
    #
    cov = 2 * timestep * D
    R = np.random.multivariate_normal(mean, cov)
    SumDijFj = (1 / (ct.k_B * ct.temperature)) * np.dot(D, F)
    positions_stacked = np.hstack(position_vectors)
    new_positions = positions_stacked + SumDijFj * timestep + R
    new_positions_list = np.hsplit(new_positions, number_of_particles)
    new_positions_array = np.zeros((number_of_particles,3), dtype=np.float64)
    for j in range(len(new_positions_list)):
        new_positions_array[j] = new_positions_list[j]
    return new_positions_array

def trans_bd_hi_no_brownian(position_vectors, radius, total_force_array, number_of_particles, timestep, tensor_choice='OSEEN', wall_height=0.0):
    """
    'trans_bd_hi_no_brownian'
    Performs one time step of dynamics with HI.
    """
    mean = np.zeros(number_of_particles * 3)
    D,_ = diffusion_matrix(position_vectors, radius, tensor_choice, wall_height)
    F = np.hstack(total_force_array)
    SumDijFj = (1 / (ct.k_B * ct.temperature)) * np.dot(D, F)
    positions_stacked = np.hstack(position_vectors)
    new_positions = positions_stacked + SumDijFj * timestep
    new_positions_list = np.hsplit(new_positions, number_of_particles)
    new_positions_array = np.zeros((number_of_particles,3), dtype=np.float64)
    for j in range(len(new_positions_list)):
        new_positions_array[j] = new_positions_list[j]
    return new_positions_array


def trans_bd_shake_hi(position_vectors, radius, total_force_array, number_of_particles, timestep, separation_list, particle_separation, constraints, tensor_choice='OSEEN', wall_height=0.0):
    """
    'trans_bd_shake_hi'
    Performs one time step of translational Brownian dynamics with constrained HI.
    """
    mean = np.zeros(number_of_particles * 3)

    #particle_separation = 0.95e-6
    D,D_matrix = diffusion_matrix(position_vectors, radius, tensor_choice, wall_height)
    F = np.hstack(total_force_array)
    #
    cov = 2 * timestep * D
    R = np.random.multivariate_normal(mean, cov)
    SumDijFj = (1 / (ct.k_B * ct.temperature)) * np.dot(D, F)
    current_positions = np.copy(position_vectors)
    r_0_list = np.copy(current_positions)
        
    for k in range(number_of_particles):
        index1 = k*3
        index2 = index1 + 3

        #Equation A1

        r_0_k = current_positions[k]
        R_k = R[index1:index2]
        SumDkjFk = SumDijFj[index1:index2]

        #print('init change',time_step * SumDkjFk + R_k)
        r_k_prime = r_0_k + timestep * SumDkjFk + R_k
        current_positions[k] = r_k_prime

    r_prime_list = np.copy(current_positions)
    #print("Constraints",constraints)
    for j in constraints:
        #print("J",j)
        r_prime_change = np.zeros([number_of_particles,3])

        for p in j:
            #print("p",p)
            m = p[0]
            n = p[1]
            k_a = p[0]
        
            if n != m:
                squared_dif = particle_separation**2 - (np.linalg.norm(r_prime_list[m] - r_prime_list[n]))**2
                r_m_n_0_c = r_0_list[m] - r_0_list[n]
                separation_list.append(np.linalg.norm(r_m_n_0_c))
                D_i_n = D_matrix[k_a][n]
                D_m_m = D_matrix[m][m]
                D_i_m = D_matrix[k_a][m]
                D_m_n = D_matrix[m][n]

                current_positions[k_a] += (squared_dif * np.dot((D_i_m-D_i_n),r_m_n_0_c)/(4*np.dot(np.dot(r_m_n_0_c,(D_m_m-D_m_n)),r_m_n_0_c)))

                r_prime_change[k_a] += (squared_dif * np.dot((D_i_m-D_i_n),r_m_n_0_c)/(4*np.dot(np.dot(r_m_n_0_c,(D_m_m-D_m_n)),r_m_n_0_c)))

        r_prime_list += r_prime_change
            
    return current_positions, separation_list