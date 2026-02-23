# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 10:59:31 2018
Original author: Chaoyi Zhang
This version March 2025
Extensively reworked by SH
"""

import sys
import numpy as np
import time
import h5py
import pyoptics.forces as pyf
import matplotlib.pyplot as plt
from pyoptics import dipoles
from pyoptics import readyaml
from pyoptics import hydro
from pyoptics import constants as ct
from scipy.spatial.transform import Rotation

def perform_simulation(number_of_particles, positions, sphere_radius, dipole_radius):

    position_vectors = positions
    #print(positions)
    separation_list = []
    #positions_list = []

    #=========================================================
    # Set up some simulation constants and parameters
    #=========================================================
    number_of_timesteps = frames
    #spring_stiffness = parameters.spring_stiffness
    bending_stiffness = parameters.bending_stiffness
    wall_position = parameters.wall_position
    wall_force_max = parameters.wall_force_max
    wall_order = parameters.wall_order
    wall_zrange = parameters.wall_zrange
    #=========================================================
    #include_dynamics = simulation.include_dynamics
    dynamics_method = simulation.dynamics_method
    hi_method = simulation.hi_method
    include_nonbonded = simulation.include_nonbonded
    include_springs = simulation.include_springs
    include_bending = simulation.include_bending
    include_driving = simulation.include_driving
    include_wall = simulation.include_wall
    include_gravity = simulation.include_gravity
    #=========================================================
    topology = options.topology
    pair_connections = topology.get_bond_connections()
    pair_separations = topology.get_separations()
    pair_stiffnesses = topology.get_bond_stiffness()
    pair_constraints = topology.get_constraints()
    constrained_separation = topology.particle_separation
    angle_connections = topology.get_angle_connections()
    angle_values = topology.get_angles()
    bending_stiffnesses = topology.get_angle_stiffness()
    #print("Spring info:")
    #print(pair_connections)
    #print(pair_separations)
    #print(pair_stiffnesses)
    #print(pair_constraints)
    #print(angle_connections)
    #print(angle_values)
    #print(bending_stiffnesses)
    #=========================================================
    # If springs are in use, we probably want to omit non-bonded
    # interactions for the connected pairs.  This is just to
    # produce the list of exceptions.
    #=========================================================
    #
    nonbonded_exceptions = None
    if pair_connections is not None:
        number_of_exceptions = len(pair_connections)
        nonbonded_exceptions = np.zeros((number_of_particles,number_of_particles),dtype=np.bool_) # all False
        for ij in range(number_of_exceptions):
            i = pair_connections[ij,0]
            j = pair_connections[ij,1]
            nonbonded_exceptions[i,j] = True
            nonbonded_exceptions[j,i] = True
    #print(nonbonded_exceptions)
    #=========================================================
    #
    vectors_list = []
    vectors_array = np.zeros(number_of_timesteps, dtype=object)
    temp_array1 = np.zeros(number_of_timesteps, dtype=object)
    #
    # Generate a list of dipoles for one sphere
    #
    dipole_primitive, _ = dipoles.sphere_positions(sphere_radius, dipole_radius)
    print("dipole primitive shape")
    print(dipole_primitive.shape)
    print("dipole primitive type")
    print(type(dipole_primitive))
    
    '''
    Now calculate new dipole rotations and store it
    So store dipole_primitive and tipole_primitive t-1 separate for each particle
    
    '''

    
    #
    #=========================================================
    # Output options
    #=========================================================

    #if excel_output==True:
    optpos = np.zeros((frames,n_particles,3))
    dipole_positions_all = np.zeros((frames,dipole_primitive.shape[0]*n_particles,3))
    if include_force==True:
        optforce = np.zeros((frames,n_particles,3))
    else:
        optforce = None
    if include_couple==True:
        optcouple = np.zeros((frames,n_particles,3))
    else:
        optcouple = None

    #=========================================================
    # Main simulation loop
    #=========================================================
    for i in range(number_of_timesteps):
        #
        # Pass in list of dipole positions to generate total dipole array;
        # All changes inside optical_force_array().
        #
        #optical,couples = optical_force_array(position_vectors, E0, dipole_radius, dipole_primitive)
        if i == 0:
            rotated_dipoles = dipole_primitive
            dipole_above_zero = []
            dipole_below_zero = []
            for di_ind, dipole in enumerate(rotated_dipoles):
                if dipole[0] > 0:
                    dipole_above_zero.append(di_ind)
                else:
                    dipole_below_zero.append(di_ind)
                    
            
        #test_rot = Rotation.from_rotvec([np.pi/3000, 0, 0])
        #test_rot.as_matrix()
        #rotated_dipoles = test_rot.apply(rotated_dipoles)
        #print("rotated_dipoles type")
        #print(type(rotated_dipoles))
        #print("rotated_dipoles shape")
        #print(rotated_dipoles.shape)
        #print(f"BEFORE FRED {dipole_above_zero}")
        optical, torques, couples, dipole_positions = dipoles.py_optical_force_torque_array(position_vectors, dipole_radius, rotated_dipoles, inverse_polarizability, beam_collection, inverse_polarizability_2, dipole_above_zero, dipole_below_zero)
        

        
        #couples = None
        #include_couple==False
        #if excel_output==True:
        for j in range(n_particles):
            for k in range(3):
                optpos[i,j,k] = position_vectors[j][k]
        if include_force==True:
            for j in range(n_particles):
                for k in range(3):
                    optforce[i,j,k] = optical[j][k]
        if include_couple==True:
            for j in range(n_particles):
                for k in range(3):
                    optcouple[i,j,k] = couples[j][k] + torques[j][k]

        
        #    print(i,optical[1])
        #
        # Generate displacements matrix
        #
        displacements_matrix = pyf.displacement_matrix(position_vectors)
        #
        # Compute relevant force components
        #
        total_force_array = optical
        if include_nonbonded == True:
            #buckingham = pyf.buckingham_force_array(number_of_particles, displacements_matrix, radius)
            #total_force_array += buckingham
            nonbonded = pyf.non_bonded_array(number_of_particles, displacements_matrix, position_vectors, nonbonded_exceptions)
            total_force_array += nonbonded
        if include_springs == True:
            if pair_connections is not None:
                spring = pyf.spring_force_array(number_of_particles, displacements_matrix, pair_separations, pair_stiffnesses, pair_connections)
                total_force_array += spring
            else:
                sys.exit("Simulation error: spring forces selected but no connection information given.")

        if include_driving == True:
            driver = pyf.driving_force_array(number_of_particles, displacements_matrix, radius)
            total_force_array += driver
        
        if include_wall == True:
            wall_force = pyf.wall_force_array(number_of_particles, position_vectors, wall_position, wall_force_max, wall_order, wall_zrange)
            total_force_array += wall_force
        
        if include_bending == True:
            if angle_connections is not None:
                bending = pyf.bending_force_array(number_of_particles, position_vectors, angle_values, bending_stiffnesses, angle_connections)
                total_force_array += bending
            else:
                sys.exit("Simulation error: bending forces selected but no connection information given.")
        
        if include_gravity == True:
            gravity = pyf.gravity_force_array(number_of_particles, position_vectors, radius)
            total_force_array += gravity
        #
        # Brownian dynamics with hydrodynamics (translational) and choice of Oseen, Rotne-Prager, and Rotne-Prager-Blake tensor
        #
        if dynamics_method=='BD_TRANS_HI' or dynamics_method=='OSEEN': # OSEEN Deprecated
            new_positions = hydro.trans_bd_hi(position_vectors, radius, total_force_array, number_of_particles, timestep, tensor_choice=hi_method, wall_height=wall_position)

        #
        # Brownian dynamics with constraints (SHAKE_HI) hydrodynamics (translational) and choice of Oseen, Rotne-Prager, and Rotne-Prager-Blake tensor.  The Blake tensor also requires a wall to be specified on the z-axis.
        #
        elif dynamics_method=='BD_TRANS_SHAKE_HI' or dynamics_method=='SHAKE_HI': # SHAKE_HI Deprecated
            new_positions, final_separation_list = hydro.trans_bd_shake_hi(position_vectors, radius, total_force_array, number_of_particles, timestep, separation_list, constrained_separation, pair_constraints, tensor_choice=hi_method, wall_height=wall_position)

        #
        # General translational HI code here +/- constraints +/- wall NO LONGER NEEDED
        #
        #elif dynamics_method=='TRANS_HI':
        #    new_positions, separation_list = hydro.translational_method(position_vectors, radius, total_force_array, number_of_particles, timestep, separation_list)

        #
        # Force scan dimer_z_theta
        #
        elif dynamics_method=='DIMER_Z_THETA':
            print("Rotational scan")
            iz = i//(number_of_timesteps+1) #101
            itheta = i%(number_of_timesteps+1) #101
            xpos = position_vectors[0,0]#3.286e-6 #6e-6
            #dtheta = 26.565*np.pi/180#2*np.pi/number_of_timesteps
            dtheta = 2*np.pi/number_of_timesteps
            z_com = (position_vectors[0,2]+position_vectors[1,2])/2
            theta = (itheta+1)*dtheta # start at 0 degrees
            radial_separation = 0.0#1.5e-8
            myradius = radius + radial_separation
            new_positions = np.zeros((number_of_particles,3),dtype=np.float64)
            new_positions[1,0] = xpos
            new_positions[0,0] = xpos
            new_positions[1,1] = myradius*np.cos(theta)
            new_positions[0,1] = myradius*np.cos(theta+np.pi)
            new_positions[1,2] = z_com+myradius*np.sin(theta)
            new_positions[0,2] = z_com+myradius*np.sin(theta+np.pi)
            #print(new_positions)
            # rotate the primitive
            rotation_matrix = np.eye(3)
            rotation_matrix[1,1] = np.cos(dtheta)
            rotation_matrix[2,2] = np.cos(dtheta)
            rotation_matrix[1,2] = -np.sin(dtheta)
            rotation_matrix[2,1] = np.sin(dtheta)
            dipole_primitive_new = np.dot(rotation_matrix,dipole_primitive.T).T
            dipole_primitive = dipole_primitive_new.copy(order='C')
            #print(dipole_primitive)

        #
        # Force scan dimer_z_x
        #
        elif dynamics_method=='DIMER_Z_X':
            nmp = 10
            iz = i//nmp
            ix = i%nmp
            offz = 3.1e-6 #5.75e-6
            xpos = 2.3e-6 #6e-6
            dz = 10e-8
            dx = dz
            z_com = (iz-nmp/2)*dz+offz
            x_com = (ix-nmp/2)*dx+xpos
            theta = np.pi#0.0#(itheta+1)*dtheta-np.pi/2 # start at -90 degrees
            radial_separation = 0.0#1.5e-8
            myradius = radius + radial_separation
            new_positions = np.zeros((number_of_particles,3),dtype=np.float64)
            new_positions[0,0] = x_com
            new_positions[1,0] = x_com
            new_positions[0,1] = myradius*np.cos(theta)
            new_positions[1,1] = myradius*np.cos(theta+np.pi)
            new_positions[0,2] = z_com+myradius*np.sin(theta)
            new_positions[1,2] = z_com+myradius*np.sin(theta+np.pi)

        elif dynamics_method=='X_SCAN':
            dx = 10e-9
            dy = 0.0
            dz = 0.0
            new_positions = np.zeros((number_of_particles,3),dtype=np.float64)
            for iscan in range(number_of_particles):
                new_positions[iscan,0] = position_vectors[iscan,0] + dx
                new_positions[iscan,1] = position_vectors[iscan,1] + dy
                new_positions[iscan,2] = position_vectors[iscan,2] + dz

        elif dynamics_method=='XZ_BINDING':
            dz = 100e-9
            dy = 0.0
            dx = 0.0
            new_positions = np.zeros((number_of_particles,3),dtype=np.float64)
            for iscan in range(number_of_particles):
                if iscan==0:
                    new_positions[iscan,0] = position_vectors[iscan,0]
                    new_positions[iscan,1] = position_vectors[iscan,1]
                    new_positions[iscan,2] = position_vectors[iscan,2]
                else:
                    new_positions[iscan,0] = position_vectors[iscan,0] + dx
                    new_positions[iscan,1] = position_vectors[iscan,1] + dy
                    new_positions[iscan,2] = position_vectors[iscan,2] + dz


        elif dynamics_method=='DIMER_Z':
            nmp = 100
            iz = i#//nmp
            ix = nmp/2
            offz = 2.5e-6 #5.75e-6
            xpos = 2.3e-6 #6e-6
            dz = 5e-8
            dx = dz
            z_com = (iz-nmp/2)*dz+offz
            x_com = (ix-nmp/2)*dx+xpos
            theta = 0.0#np.pi#0.0#(itheta+1)*dtheta-np.pi/2 # start at -90 degrees
            radial_separation = 0.0#1.5e-8
            myradius = radius + radial_separation
            new_positions = np.zeros((number_of_particles,3),dtype=np.float64)
            new_positions[0,0] = x_com
            new_positions[1,0] = x_com
            new_positions[0,1] = myradius*np.cos(theta)
            new_positions[1,1] = myradius*np.cos(theta+np.pi)
            new_positions[0,2] = z_com+myradius*np.sin(theta)
            new_positions[1,2] = z_com+myradius*np.sin(theta+np.pi)


        elif dynamics_method == "ONLY_OPTICS":
            new_positions = hydro.trans_bd_hi_no_brownian(position_vectors, radius, total_force_array, number_of_particles, timestep, tensor_choice=hi_method, wall_height=wall_position)
            pass

        elif dynamics_method == "NOTHING":
            new_positions = position_vectors
            pass
        '''
        It's not complicated to store the amount of rotation in each timestep
        and rot commutitivititytyty should not matter bc timestep is small
        
        Can this work if we never store actual rotation of lab frame of each particle
        at a given timestep? But only store change in rot per timestep
        
        1 Set of dipole pos for each particle
        Every timestep - apply neccessary rotation to that object
        When forces are calculated, make ref to rhat object.
        '''
        position_vectors = new_positions
        dipole_positions = rotated_dipoles[None, :, :] + position_vectors[:, None, :] # there'll be a shape mismatch. But transform the dipole pos. to lab frame
        dipole_positions_all[i] = dipole_positions.reshape(-1, 3) # Makes sure there's always 3 rows
        #dipole_positions_all
        vectors_list.append(position_vectors)  # returns list of position vector arrays of all particles
        # print("Positions:",vectors_list)
        if i%10 == 0:
            print("Step ",i)
            print(i,optical[0])
            #print(position_vectors)
    for k in range(number_of_timesteps):
        vectors_array[k] = vectors_list[k]
        temp_array1[k] = np.hstack(vectors_array[k])

    xyz_list1 = np.vsplit(np.vstack(temp_array1).T, number_of_particles)
    
    if dynamics_method=='BD_TRANS_SHAKE_HI' or dynamics_method=='SHAKE_HI':
        print("Mean separation: ",np.mean(final_separation_list))
        #print(separation_list)
    return xyz_list1,optpos,optforce,optcouple,dipole_positions_all,dipole_above_zero,dipole_below_zero



###################################################################################
# Start of program
###################################################################################

if int(len(sys.argv)) != 2:
    sys.exit("Usage: python {} <FILESTEM>".format(sys.argv[0]))

filestem = sys.argv[1]
filename_vtf = filestem+".vtf"
filename_xl = filestem+".xlsx"
filename_hdf = filestem+".h5"
filename_yaml = filestem+".yml"
#===========================================================================
# Read the yaml file into a system parameter dictionary
#===========================================================================
options = readyaml.Options(filestem)
print(options)
#===========================================================================
# Read simulation parameters
#===========================================================================
parameters = options.parameters
wavelength = parameters.wavelength
dipole_radius = parameters.dipole_radius
timestep = parameters.time_step
#===========================================================================
# Read simulation options
#===========================================================================
simulation = options.simulation
frames = options.simulation.frames
print("Frames:",frames)
#===========================================================================
# Read output options
#===========================================================================
output = options.output
vmd_output = output.vmd_output
excel_output = output.excel_output
hdf_output = output.hdf_output
include_force = output.include_force
include_couple = output.include_couple
#===========================================================================
# Read display options
#===========================================================================
display = options.display
#===========================================================================
# Read beam options and create beam collection
#===========================================================================
beam_collection = options.beam_collection
#===========================================================================
# Read particle options and create particle collection
#===========================================================================
particle_collection = options.particle_collection
print(particle_collection.num_particles)
n_particles = particle_collection.num_particles

ref_ind = particle_collection.get_refractive_indices()
particle_types = particle_collection.get_particle_types()
colors = particle_collection.get_particle_colours()
vtfcolors = particle_collection.get_particle_vtfcolours()
radii = particle_collection.get_particle_radii()
radius = radii[0] # because we cannot handle variable radii yet.
density = particle_collection.get_particle_density()
rho = density[0] # not yet implemented.
positions = particle_collection.get_particle_positions()

for i in range(n_particles):
    print(i,particle_types[i],ref_ind[i],colors[i],radii[i],density[i],positions[i])
#===========================================================================
# Set up particle polarisabilities and other spurious options
#===========================================================================`
ep1 = ref_ind**2
#mass = (4/3)*rho*np.pi*radius**3
#gravity = np.zeros(3,dtype=np.float64)
#gravity[1] = -9.81*mass
print("dipole radius is:",dipole_radius,type(dipole_radius))

nm = 1.333 # water
epm = nm*nm
k = nm*2 * np.pi / wavelength
#k = 2 * np.pi / wavelength

a0 = (4 * 6 * ct.eps0) * (dipole_radius ** 3) * ((ep1 - epm) / (ep1 + 2*epm)) # Corrected formula
a = a0 / (1 - (2 / 3) * 1j * k ** 3 * a0/(4*np.pi*ct.eps0))

polarizability = a
print("Polarisability:")
print(polarizability)
inverse_polarizabilitys = (1.0+0j)/polarizability # added this for the C++ wrapper (Chaumet's alpha bar)
print("Inverse Polarisability:")
print(inverse_polarizabilitys)

inverse_polarizability = np.ascontiguousarray(inverse_polarizabilitys[:,0])
inverse_polarizability_2 = np.ascontiguousarray(inverse_polarizabilitys[:,1])
print(inverse_polarizability)
print(inverse_polarizability_2)
'''
ep2 = ref_ind[1]**2
a0_2 = (4 * 6 * ct.eps0) * (dipole_radius ** 3) * ((ep2 - epm) / (ep2 + 2*epm)) # Corrected formula
a_2 = a0_2 / (1 - (2 / 3) * 1j * k ** 3 * a0_2/(4*np.pi*ct.eps0))
polarizability_2 = a_2
print("Polarisability_2:")
print(polarizability_2)
inverse_polarizability_2 = np.array([(1.0+0j)/polarizability_2])
print("Inverse Polarisability_2:")
print(inverse_polarizability_2)
print(type(inverse_polarizability_2))
'''
z_offset = wavelength / 4.0 # needed for odd order Bessel beams
z_offset = 0.0 # for most other situations

#===========================================================================
# Perform the simulation
#===========================================================================

initialT = time.time()
particles,optpos, optforces,optcouples,dipole_positions,dipole_above_zero,dipole_below_zero = perform_simulation(n_particles, positions, radius, dipole_radius)
finalT = time.time()
print("Elapsed time: {:8.6f} s".format(finalT-initialT))
_, number_of_dipoles = dipoles.sphere_positions(radius, dipole_radius)
file = open("R7-9800X3D-DDA-BooleanMask.txt", "a")
file.write(f"{number_of_dipoles}\t{finalT-initialT}\n")
file.close()

#===========================================================================
# This code for matplotlib animation
#===========================================================================
dynamics_method = simulation.dynamics_method
if display.show_output==True:

    fig,ax = display.plot_intensity(beam_collection)

    parpole_ani = display.plot_parpoles(fig,ax,particles,dipole_positions,radius,colors,dipole_above_zero,dipole_below_zero)
    #particle_ani = display.animate_particles(fig,ax,particles,radius,colors)
    #dipole_ani = display.animate_dipoles(fig,ax,dipole_positions,radius,colors)
    #particle_ani = display.animate_particles(fig,ax,particles,radius,colors)
    #print(particles)
    # ax.set_xlim(-7E-6,7E-6)
    # ax.set_ylim(-7E-6,7E-6)
    parpole_ani.save(f"dipole_radius_{dipole_radius}-dynamics_method_{dynamics_method}.mp4", dpi = 300, fps=60)
    print("==========================")
    plt.show()
    print("--------------------------")
    #print(dipole_positions)
#===========================================================================
# Write out data to files
#===========================================================================

if vmd_output==True:
    output.make_vmd_file(filename_vtf,n_particles,frames,timestep,particles,optpos,beam_collection,finalT-initialT,radius,dipole_radius,z_offset,particle_types,vtfcolors)

if excel_output==True:
    output.make_excel_file(filename_xl,n_particles,frames,timestep,particles,optpos,include_force,optforces,include_couple,optcouples)

if hdf_output==True:
    h5file = h5py.File(filename_hdf,'w')
    h5file.create_dataset('positions',data=optpos)
    #h5file.create_dataset('particles',data=particles)
    if include_force==True:
        h5file.create_dataset('forces',data=optforces)
    if include_couple==True:
        h5file.create_dataset('couples',data=optcouples)
    h5file.close()

