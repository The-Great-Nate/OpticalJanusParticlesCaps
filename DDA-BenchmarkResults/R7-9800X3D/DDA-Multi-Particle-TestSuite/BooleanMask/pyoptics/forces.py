#
# Force calculations of various sorts.
#
import numpy as np
import itertools as it
from scipy.spatial import KDTree


def displacement_matrix(array_of_positions):
    number_of_particles = len(array_of_positions)
    #print("x",number_of_dipoles)
    list_of_displacements = [u - v for u, v in it.combinations(array_of_positions, 2)]
    array_of_displacements = np.zeros(len(list_of_displacements), dtype=object)
    for i in range(len(list_of_displacements)):
        array_of_displacements[i] = list_of_displacements[i]
    displacement_matrix = np.zeros([number_of_particles, number_of_particles], dtype=object)
    iu = np.triu_indices(number_of_particles, 1)
    displacement_matrix[iu] = array_of_displacements
    displacement_matrix.T[iu] = -array_of_displacements
    return displacement_matrix

def non_bonded(constant1, constant2, r):
#    r_max = 2.2 * particle_radius
    r_abs = np.linalg.norm(r)
    if r_abs > constant2:
        print("Eeek!! Too Big! r_abs = ", r_abs)
        r_abs = constant2  # capping the force
    force = - constant1 * (1.0-r_abs/constant2)*r/r_abs
    return force

def non_bonded_array(number_of_particles, displacements_matrix, positions, exceptions):
    #
    # This version using a DPD-like non-bonded force expression.
    #
    displacements_matrix_T = np.transpose(displacements_matrix)
    #
    ConstantA = 1e-12 # 1pN reasonable max force
    ConstantB = 0.6e-6 # guess at nearest separation
#    nonbonded_force_matrix = np.zeros(
#        [number_of_particles, number_of_particles], dtype=object
#    )
    nonbonded_force_array = np.zeros((number_of_particles,3),dtype=np.float64)
    #
    # The force expression is not defined beyond r = ConstantB, therefore
    # start by finding list of beads within that distance
    #
    mytree = KDTree(positions)
    numbers = mytree.query_ball_point(positions,ConstantB,return_length=False)
    #print(numbers)
    for i in range(number_of_particles):
        for j in range(len(numbers[i])):
            #print(i,j,numbers[i][j])
            k = numbers[i][j]
            if exceptions is not None:
                exception = exceptions[i,k]
            else:
                exception = False
            #if i<92 and k<92:
            #    exception = True
            if i != k and exception==False:
                #print("Nonbond:",i,k)
                nonbonded_force_array[i] += non_bonded(ConstantA, ConstantB,
                    displacements_matrix_T[i][k])
#    temp = np.sum(nonbonded_force_matrix, axis=1)
#    for i in range(number_of_particles):
#        nonbonded_force_array[i] = temp[i]
    return nonbonded_force_array

def non_bonded_OLD(constant1, constant2, r, particle_radius):
    r_max = 2.2 * particle_radius
    r_abs = np.linalg.norm(r)
    if r_abs < r_max:
        print("Eeek!! r_abs = ", r_abs)
        r_abs = r_max  # capping the force
    force = np.array(
        [- constant1 * constant2 * np.exp(-constant2 * r_abs)
            * (r[i] / r_abs)
            for i in range(3)])
    return force

def non_bonded_array_OLD(number_of_particles, displacements_matrix, particle_radius):
    #
    displacements_matrix_T = np.transpose(displacements_matrix)
    #
    ConstantA = 1.0e-18
    ConstantB = 5e7
    nonbonded_force_matrix = np.zeros(
        [number_of_particles, number_of_particles], dtype=object
    )
    for i in range(number_of_particles):
        for j in range(number_of_particles):
            if i == j:
                nonbonded_force_matrix[i][j] = np.asarray([0, 0, 0],dtype=np.float64)
            else:
                nonbonded_force_matrix[i][j] = non_bonded(
                    ConstantA,
                    ConstantB,
                    displacements_matrix_T[i][j],
                    particle_radius,
                )
    nonbonded_force_array = np.zeros((number_of_particles,3),dtype=np.float64)
    temp = np.sum(nonbonded_force_matrix, axis=1)
    for i in range(number_of_particles):
        nonbonded_force_array[i] = temp[i]
    return nonbonded_force_array



def buckingham_force(Hamaker, constant1, constant2, r, particle_radius):
    r_max = 2.2 * particle_radius
    r_abs = np.linalg.norm(r)
    if r_abs < r_max:
        print("Eeek!! r_abs = ", r_abs)
        r_abs = r_max  # capping the force

    force = np.array(
        [
            -(
                constant1 * constant2 * np.exp(-constant2 * r_abs)
                - (
                    (32 * Hamaker * (particle_radius ** 6))
                    / (3 * (r_abs ** 3) * (r_abs ** 2 - 4 * (particle_radius ** 2)) ** 2)
                )
            )
            * (r[i] / r_abs)
            for i in range(3)
        ]
    )
    # force = np.zeros(3)
    # force[0] = -(constant1*constant2*np.exp(-constant2*r_abs) - ((32*Hamaker*(dipole_radius**6))/(3*(r_abs**3)*(r_abs**2 - 4*(dipole_radius**2))**2)))*(r[0]/r_abs)
    # force[1] = -(constant1*constant2*np.exp(-constant2*r_abs) - ((32*Hamaker*(dipole_radius**6))/(3*(r_abs**3)*(r_abs**2 - 4*(dipole_radius**2))**2)))*(r[1]/r_abs)
    # force[2] = -(constant1*constant2*np.exp(-constant2*r_abs) - ((32*Hamaker*(dipole_radius**6))/(3*(r_abs**3)*(r_abs**2 - 4*(dipole_radius**2))**2)))*(r[2]/r_abs)

    return force


def spring_force_OLD(constant1, r, particle_radius):
    #print("Dipole Radius:",dipole_radius)
    factor = 2
    r_abs = np.linalg.norm(r)
    #force = [constant1 * (r_abs - 2 * particle_radius) * (r[i] / r_abs) for i in range(3)]
    force = np.zeros(3)
    force[0] = constant1*(r_abs-factor*particle_radius)*(r[0]/r_abs)
    force[1] = constant1*(r_abs-factor*particle_radius)*(r[1]/r_abs)
    force[2] = constant1*(r_abs-factor*particle_radius)*(r[2]/r_abs)

    return force

def spring_force(stiffness, r, separation):
    # this version using supplied stiffness and separation
    r_abs = np.linalg.norm(r)
    force = np.zeros(3)
    force = stiffness*(r_abs-separation)*r/r_abs
    return force


def driving_force(constant1, r, particle_radius):
    #print("Dipole Radius:",dipole_radius)
    r_abs = np.linalg.norm(r)
    force = [constant1 * (r_abs - 2 * particle_radius) * (r[i] / r_abs) for i in range(3)]
    # force = np.zeros(3)
    # force[0] = constant1*(r_abs-2*dipole_radius)*(r[0]/r_abs)
    # force[1] = constant1*(r_abs-2*dipole_radius)*(r[1]/r_abs)
    # force[2] = constant1*(r_abs-2*dipole_radius)*(r[2]/r_abs)

    return force


def bending_force(bond_stiffness, ri, rj, rk):
    rij = rj - ri
    rik = rk - ri
    rij_abs = np.linalg.norm(rij)
    rik_abs = np.linalg.norm(rik)
    rijrik = rij_abs * rik_abs
    rij2 = rij_abs * rij_abs
    rik2 = rik_abs * rik_abs
    costhetajik = np.dot(rij, rik) / rijrik
    force = np.zeros([3, 3])
    i = 1
    force[i] = bond_stiffness * (
        (rik + rij) / rijrik - costhetajik * (rij / rij2 + rik / rik2)
    )
    force[i - 1] = bond_stiffness * (costhetajik * rij / rij2 - rik / rijrik)
    force[i + 1] = bond_stiffness * (costhetajik * rik / rik2 - rij / rijrik)
    #    print(force)
    return force



def buckingham_force_array(number_of_particles, displacements_matrix, particle_radius):
    #
    displacements_matrix_T = np.transpose(displacements_matrix)
    #
    Hamaker = 0
    ConstantA = 1.0e23
    ConstantB = 2.0e8  # 4.8e8
    buckingham_force_matrix = np.zeros(
        [number_of_particles, number_of_particles], dtype=object
    )
    for i in range(number_of_particles):
        for j in range(number_of_particles):
            if i == j:
                buckingham_force_matrix[i][j] = [0, 0, 0]
            else:
                buckingham_force_matrix[i][j] = buckingham_force(
                    Hamaker,
                    ConstantA,
                    ConstantB,
                    displacements_matrix_T[i][j],
                    particle_radius,
                )
    buckingham_force_array = np.zeros((number_of_particles,3),dtype=np.float64)
    temp = np.sum(buckingham_force_matrix, axis=1)
    for i in range(number_of_particles):
        buckingham_force_array[i] = temp[i]
    return buckingham_force_array


def spring_force_array_OLD(number_of_particles, displacements_matrix, particle_radius, stiffness):
    #number_of_dipoles = len(array_of_positions)
    #displacements_matrix = displacement_matrix(array_of_positions)
    displacements_matrix_T = np.transpose(displacements_matrix)
    #stiffness = 1.0e-4
    spring_force_matrix = np.zeros([number_of_particles, number_of_particles], dtype=object)
    for i in range(number_of_particles):
        for j in range(number_of_particles):
            spring_force_matrix[i][j] = np.zeros(3)
    # this code for pairs
#    for i in range(0,number_of_dipoles,2):
#        j = i + 1
#        spring_force_matrix[i][j] = spring_force(stiffness, displacements_matrix_T[i][j], dipole_radius)
#    for i in range(1,number_of_dipoles,2):
#        j = i - 1
#        spring_force_matrix[i][j] = spring_force(stiffness, displacements_matrix_T[i][j], dipole_radius)
    # this code for all in a line
    for i in range(0,number_of_particles-1):
        j = i + 1
        spring_force_matrix[i][j] = spring_force(stiffness, displacements_matrix_T[i][j], particle_radius)
    for i in range(1,number_of_particles):
        j = i - 1
        spring_force_matrix[i][j] = spring_force(stiffness, displacements_matrix_T[i][j], particle_radius)
    # this code for all in two lines
    """
    for i in range(0,number_of_dipoles//2-1):
        j = i + 1
        spring_force_matrix[i][j] = spring_force(stiffness, displacements_matrix_T[i][j], dipole_radius)
    for i in range(1,number_of_dipoles//2):
        j = i - 1
        spring_force_matrix[i][j] = spring_force(stiffness, displacements_matrix_T[i][j], dipole_radius)
    for i in range(number_of_dipoles//2,number_of_dipoles-1):
        j = i + 1
        spring_force_matrix[i][j] = spring_force(stiffness, displacements_matrix_T[i][j], dipole_radius)
    for i in range(number_of_dipoles//2+1,number_of_dipoles):
        j = i - 1
        spring_force_matrix[i][j] = spring_force(stiffness, displacements_matrix_T[i][j], dipole_radius)
    """
    #print("Spring force array shape before:",spring_force_matrix.shape)
    #print(spring_force_matrix)
    spring_force_array = np.zeros((number_of_particles,3),dtype=np.float64)
    #print("Springs",spring_force_array)
    #print("Spring force array shape after:",spring_force_array.shape)
    temp = np.sum(spring_force_matrix, axis=1)
    for i in range(number_of_particles):
            spring_force_array[i] = temp[i]
    return spring_force_array


def spring_force_array2_OLD(number_of_particles, displacements_matrix, particle_radius, stiffness, spring_list):
    """
    This version takes a list of springs, denoted by a pair of particle indices.
    """
    displacements_matrix_T = np.transpose(displacements_matrix)
    spring_force_matrix = np.zeros([number_of_particles, number_of_particles], dtype=object)
    for i in range(number_of_particles):
        for j in range(number_of_particles):
            spring_force_matrix[i][j] = np.zeros(3)
    # this code for list of pairs
    for ij in range(len(spring_list)):
        i = spring_list[ij,0]
        j = spring_list[ij,1]
        spring_force_matrix[i][j] = spring_force(stiffness, displacements_matrix_T[i][j], particle_radius)
        i = spring_list[ij,1]
        j = spring_list[ij,0]
        spring_force_matrix[i][j] = spring_force(stiffness, displacements_matrix_T[i][j], particle_radius)
        #print("ij, i, j: ",ij,i,j)
    # this code for all in a line
    #for i in range(0,number_of_particles-1):
    #    j = i + 1
    #    spring_force_matrix[i][j] = spring_force(stiffness, displacements_matrix_T[i][j], particle_radius)
    #for i in range(1,number_of_particles):
    #    j = i - 1
    #    spring_force_matrix[i][j] = spring_force(stiffness, displacements_matrix_T[i][j], particle_radius)
    # this code for all in two lines
    """
    for i in range(0,number_of_dipoles//2-1):
        j = i + 1
        spring_force_matrix[i][j] = spring_force(stiffness, displacements_matrix_T[i][j], dipole_radius)
    for i in range(1,number_of_dipoles//2):
        j = i - 1
        spring_force_matrix[i][j] = spring_force(stiffness, displacements_matrix_T[i][j], dipole_radius)
    for i in range(number_of_dipoles//2,number_of_dipoles-1):
        j = i + 1
        spring_force_matrix[i][j] = spring_force(stiffness, displacements_matrix_T[i][j], dipole_radius)
    for i in range(number_of_dipoles//2+1,number_of_dipoles):
        j = i - 1
        spring_force_matrix[i][j] = spring_force(stiffness, displacements_matrix_T[i][j], dipole_radius)
    """
    #print("Spring force array shape before:",spring_force_matrix.shape)
    #print(spring_force_matrix)
    spring_force_array = np.zeros((number_of_particles,3),dtype=np.float64)
    #print("Springs",spring_force_array)
    #print("Spring force array shape after:",spring_force_array.shape)
    temp = np.sum(spring_force_matrix, axis=1)
    for i in range(number_of_particles):
            spring_force_array[i] = temp[i]
    return spring_force_array

def spring_force_array(number_of_particles, displacements_matrix, separations, stiffnesses, spring_list):
    """
    This version takes a list of springs, denoted by a pair of particle indices.
    """
    displacements_matrix_T = np.transpose(displacements_matrix)
    spring_force_matrix = np.zeros([number_of_particles, number_of_particles], dtype=object)
    for i in range(number_of_particles):
        for j in range(number_of_particles):
            spring_force_matrix[i][j] = np.zeros(3)
    # this code for list of pairs
    for ij in range(len(spring_list)):
        i = spring_list[ij,0]
        j = spring_list[ij,1]
        spring_force_matrix[i][j] = spring_force(stiffnesses[ij], displacements_matrix_T[i][j], separations[ij])
        spring_force_matrix[j][i] = -spring_force_matrix[i][j]
        #i = spring_list[ij,1]
        #j = spring_list[ij,0]
        #spring_force_matrix[i][j] = spring_force(stiffnesses[ij], displacements_matrix_T[i][j], separations[ij])

    #print("Spring force array shape before:",spring_force_matrix.shape)
    #print(spring_force_matrix)
    spring_force_array = np.zeros((number_of_particles,3),dtype=np.float64)
    #print("Springs",spring_force_array)
    #print("Spring force array shape after:",spring_force_array.shape)
    temp = np.sum(spring_force_matrix, axis=1)
    for i in range(number_of_particles):
            spring_force_array[i] = temp[i]
    return spring_force_array



def driving_force_array(number_of_particles, displacements_matrix, particle_radius):
    #number_of_dipoles = len(array_of_positions)
    #displacements_matrix = displacement_matrix(array_of_positions)
    displacements_matrix_T = np.transpose(displacements_matrix)
    driver = 3.0e-7#6
    driving_force_array = np.zeros(number_of_particles, dtype=object)
    for i in range(0,number_of_particles,2):
        j = i+1
        driving_force_array[i] = driving_force(driver, displacements_matrix_T[i][j], particle_radius)
        driving_force_array[j] = driving_force_array[i]
    return driving_force_array


def bending_force_array_OLD(number_of_particles, array_of_positions, particle_radius, bond_stiffness):
    #number_of_dipoles = len(array_of_positions)
    #bond_stiffness = BENDING
    bending_force_matrix = np.zeros([number_of_particles], dtype=object)
    bending_force_temp = np.zeros([3], dtype=object)
    for i in range(1, number_of_particles - 1):
        bending_force_temp = bending_force(
            bond_stiffness,
            array_of_positions[i],
            array_of_positions[i - 1],
            array_of_positions[i + 1],
        )
        #        print(bending_force_temp)
        bending_force_matrix[i - 1] += bending_force_temp[0]
        bending_force_matrix[i] += bending_force_temp[1]
        bending_force_matrix[i + 1] += bending_force_temp[2]
    """
    for i in range(number_of_dipoles//2 + 1, number_of_dipoles - 1):
        bending_force_temp = bending_force(
            bond_stiffness,
            array_of_positions[i],
            array_of_positions[i - 1],
            array_of_positions[i + 1],
        )
        #        print(bending_force_temp)
        bending_force_matrix[i - 1] += bending_force_temp[0]
        bending_force_matrix[i] += bending_force_temp[1]
        bending_force_matrix[i + 1] += bending_force_temp[2]
    """
    #    print("Springs",spring_force_array)
    bending_force_array = np.zeros((number_of_particles,3),dtype=np.float64)
    for i in range(number_of_particles):
            bending_force_array[i] = bending_force_matrix[i]
    return bending_force_array
    

def bending_force_array2_OLD(number_of_particles, array_of_positions, particle_radius, bond_stiffness, bending_list):
    #number_of_dipoles = len(array_of_positions)
    #bond_stiffness = BENDING
    bending_force_matrix = np.zeros([number_of_particles], dtype=object)
    bending_force_temp = np.zeros([3], dtype=object)
    for ijk in range(len(bending_list)):
        j = bending_list[ijk,0]
        i = bending_list[ijk,1]
        k = bending_list[ijk,2]
    #for i in range(1, number_of_particles - 1):
        bending_force_temp = bending_force(
            bond_stiffness,
            array_of_positions[i],
            array_of_positions[j],
            array_of_positions[k],
        )
        #        print(bending_force_temp)
        bending_force_matrix[j] += bending_force_temp[0]
        bending_force_matrix[i] += bending_force_temp[1]
        bending_force_matrix[k] += bending_force_temp[2]
    """
    for i in range(number_of_dipoles//2 + 1, number_of_dipoles - 1):
        bending_force_temp = bending_force(
            bond_stiffness,
            array_of_positions[i],
            array_of_positions[i - 1],
            array_of_positions[i + 1],
        )
        #        print(bending_force_temp)
        bending_force_matrix[i - 1] += bending_force_temp[0]
        bending_force_matrix[i] += bending_force_temp[1]
        bending_force_matrix[i + 1] += bending_force_temp[2]
    """
    #    print("Springs",spring_force_array)
    bending_force_array = np.zeros((number_of_particles,3),dtype=np.float64)
    for i in range(number_of_particles):
            bending_force_array[i] = bending_force_matrix[i]
    return bending_force_array
    

def bending_force_array(number_of_particles, array_of_positions, bond_angles, bending_stiffness, bending_list):
    #
    # bond_stiffness is now an array
    # bond angle not yet used
    #
    bending_force_matrix = np.zeros([number_of_particles], dtype=object)
    bending_force_temp = np.zeros([3], dtype=object)
    for ijk in range(len(bending_list)):
        j = bending_list[ijk,0]
        i = bending_list[ijk,1] # central bead
        k = bending_list[ijk,2]
    #
        bending_force_temp = bending_force(
            bending_stiffness[ijk],
            array_of_positions[i],
            array_of_positions[j],
            array_of_positions[k],
        )
        #        print(bending_force_temp)
        bending_force_matrix[j] += bending_force_temp[0]
        bending_force_matrix[i] += bending_force_temp[1]
        bending_force_matrix[k] += bending_force_temp[2]

    #    print("Springs",spring_force_array)
    bending_force_array = np.zeros((number_of_particles,3),dtype=np.float64)
    for i in range(number_of_particles):
            bending_force_array[i] = bending_force_matrix[i]
    return bending_force_array
    




def gravity_force_array(number_of_particles, array_of_positions, particle_radius):
    #number_of_beads = len(array_of_positions)
    #bond_stiffness = BENDING
    gravity_force_matrix = np.zeros([number_of_particles], dtype=object)
    gravity_force_temp = np.zeros([3], dtype=object)
    gravity_force_temp[2] = -1e-12#-9.81*mass
    for i in range(number_of_particles):
        gravity_force_matrix[i] = gravity_force_temp
    return gravity_force_matrix


def wall_force(z, wall_position,force_max,order,zrange):
    """
    Mainly for use with the Blake tensor hydrodynamics near a wall.
    Using a simple polynomial centred on the wall, repulsive in both directions.
    Parameters for the potential are:
        wall_position: position along the z axis
        force_max: maximum force actually at the wall - force is scaled to this value
        order: positive integer - the actual force will go as z^(2*order-1)
        zrange: effective range of the force in either direction from the wall i.e. force
            will be force_max at the wall and zero at wall_position+/-zrange.
    """
    
    force =-force_max*((z-wall_position-zrange)/zrange)**(2*order-1)

    return force


def wall_force_array(num_particles, positions, wall_position, force_max, order, zrange):
    """
    Mainly for use with the Blake tensor hydrodynamics near a wall.
    Using a simple polynomial centred on the wall, repulsive in both directions.
    Parameters for the potential are:
        wall_position: position along the z axis
        force_max: maximum force actually at the wall - force is scaled to this value
        order: positive integer - the actual force will go as z^(2*order-1)
        zrange: effective range of the force in either direction from the wall i.e. force
            will be force_max at the wall and zero at wall_position+/-zrange.
        zthick: the central thickness of the wall that particles cannot enter.  This is normally
            equal to the radius of the particle.  In this region the force is capped.
    """
    zthick = 0.45e-6 # temporary fix to radius of particle
    wall_forces = np.zeros((num_particles,3))
    for i in range(num_particles):
        zi = positions[i,2]
        if zi < wall_position:
        # lower side of wall
            if zi > wall_position-zthick:
                wall_forces[i,2] = -force_max
            elif zi > wall_position-zrange-zthick:
                wall_forces[i,2] = wall_force(zi+zthick,wall_position,-force_max,order,-zrange)
        else:
        # upper side of wall
            if zi < wall_position+zthick:
                wall_forces[i,2] = force_max
            elif zi < wall_position+zrange+zthick:
                wall_forces[i,2] = wall_force(zi-zthick,wall_position,force_max,order,zrange)
    return wall_forces
