#
# Python wrapper for Dipoles.cpp
#
""" Python wrapper for the C++ shared library Dipoles.  This is for
anything related to the dipole generation and optical force calculation."""
import sys, platform
import ctypes, ctypes.util
import numpy.ctypeslib
import numpy as np
from pyoptics import beams
import pyoptics

###################################################################################

# Find the library and load it
#dipoles_path = ctypes.util.find_library(pyoptics.__path__[0]+"/dipoleslib")
#if not dipoles_path:
#    print("Unable to find the specified Dipoles library.")
#    sys.exit()
#
#try:
#    Dipoles = ctypes.CDLL(dipoles_path)
#except OSError:
#    print("Unable to load the Dipoles C++ library")
#    sys.exit()

try:
    Dipoles = ctypes.cdll.LoadLibrary("./pyoptics/dipoleslib.so")
except OSError:
    print("Unable to load the Dipoles C++ library")
    sys.exit()

#
#
# Python helper functions
#
def sphere_positions(sphere_radius, dipole_radius):
    """
    Dipoles displaced 0.5 radius from origin - must be even number,
    e.g. works best if you have say 16 dipoles fitting within sphere diameter.
    This model consistent with ADDA and MOSH-DDA.
    """
    dipole_diameter = 2*dipole_radius
    dd2 = dipole_diameter**2
    sr2 = sphere_radius**2
    print(sphere_radius,dipole_radius)
#    num = int(sphere_radius//dipole_diameter)
    num = int(sphere_radius/dipole_diameter+0.5)
    print("Dipoles per diameter:",2*num)
    number_of_dipoles = 0
    points = []
    for i in range(-num,num):
        x = (i+0.5)*dipole_diameter
        i2 = x*x
        for j in range(-num,num):
            y = (j+0.5)*dipole_diameter
            j2 = y*y
            for k in range(-num,num):
                z = (k+0.5)*dipole_diameter
                k2 = z*z
                rad2 = (i2+j2+k2)
                if rad2 < sr2:
                    number_of_dipoles += 1
                    points.append([x,y,z])
    pts = np.asarray(points)
#    number_of_dipoles = 0
#    for i in range(-num,num+1):
#        i2 = i*i
#        x = i*dipole_diameter
#        for j in range(-num,num+1):
#            j2 = j*j
#            y = j*dipole_diameter
#            for k in range(-num,num+1):
#                k2 = k*k
#                z = k*dipole_diameter
#                rad2 = (i2+j2+k2)*dd2
#                if rad2 < sr2:
#                    pts[number_of_dipoles][0] = x+1e-20
#                    pts[number_of_dipoles][1] = y+1e-20
#                    pts[number_of_dipoles][2] = z
#                    number_of_dipoles += 1
    print(number_of_dipoles," dipoles generated")
    #print(pts)
    return pts
    
def sphere_positions_indexed(sphere_radius, dipole_radius):
    """
    Dipoles displaced 0.5 radius from origin - must be even number,
    e.g. works best if you have say 16 dipoles fitting within sphere diameter.
    This model consistent with ADDA and MOSH-DDA.
    This version returns an indexing array as well.
    """
    dipole_diameter = 2*dipole_radius
    dd2 = dipole_diameter**2
    sr2 = sphere_radius**2
    print(sphere_radius,dipole_radius)
#    num = int(sphere_radius//dipole_diameter)
    num = int(sphere_radius/dipole_diameter+0.5)
    print("Dipoles per diameter:",2*num)
    number_of_dipoles = 0
    points = []
    ijkarray = -np.ones((2*num,2*num,2*num),dtype=int)
    for i in range(-num,num):
        x = (i+0.5)*dipole_diameter
        i2 = x*x
        for j in range(-num,num):
            y = (j+0.5)*dipole_diameter
            j2 = y*y
            for k in range(-num,num):
                z = (k+0.5)*dipole_diameter
                k2 = z*z
                rad2 = (i2+j2+k2)
                if rad2 < sr2:
                    ijkarray[i+num,j+num,k+num] = number_of_dipoles
                    number_of_dipoles += 1
                    points.append([x,y,z])
    pts = np.asarray(points)
#                    number_of_dipoles += 1
    print(number_of_dipoles," dipoles generated")
    return pts,ijkarray
    
def sphere_positions_OLD(sphere_radius, dipole_radius):
    dipole_diameter = 2*dipole_radius
    dd2 = dipole_diameter**2
    sr2 = sphere_radius**2
    print(sphere_radius,dipole_radius)
    num = int(sphere_radius//dipole_diameter)
    number_of_dipoles = 0
    for i in range(-num,num+1):
        i2 = i*i
        for j in range(-num,num+1):
            j2 = j*j
            for k in range(-num,num+1):
                k2 = k*k
                rad2 = (i2+j2+k2)*dd2
                if rad2 < sr2:
                    number_of_dipoles += 1
    pts = np.zeros((number_of_dipoles, 3))
    number_of_dipoles = 0
    for i in range(-num,num+1):
        i2 = i*i
        x = i*dipole_diameter
        for j in range(-num,num+1):
            j2 = j*j
            y = j*dipole_diameter
            for k in range(-num,num+1):
                k2 = k*k
                z = k*dipole_diameter
                rad2 = (i2+j2+k2)*dd2
                if rad2 < sr2:
                    pts[number_of_dipoles][0] = x+1e-20
                    pts[number_of_dipoles][1] = y+1e-20
                    pts[number_of_dipoles][2] = z
                    number_of_dipoles += 1
    print(number_of_dipoles," dipoles generated")
    return pts
    
#
# Python wrapper functions
#
def py_grad_E_cc(position, polarisation, kvec):
    """
    position: x, y, z coordinates of point (double precision);
    polarisation: complex vector;
    kvec: scalar wave vetor (should be corrected for medium)
    gradEE: a complex array to receive the gradients.
    """
    dgradEE = np.zeros(18,dtype=np.float64)
    polarisation_unwrap = polarisation.view(dtype=np.float64).reshape((6,1)).flatten()
    Dipoles.grad_E_cc(position, polarisation_unwrap, kvec, dgradEE)
    gradEE = dgradEE.view(dtype=np.complex128).reshape((3,3))

#    print("position:",position)
#    print("polarisation:",polarisation)
#    print("kvec:",kvec)
#    print("gradEE:",gradEE)

    return gradEE

#==============================================================================
# Wrapper for new optical force array code
#==============================================================================

def py_optical_force_array(array_of_particles,dipole_radius,dipole_primitive,inverse_polarisation, beam_collection):
    """
    Wrapper function for the new optical force code.
    """
    num_particles = len(array_of_particles)
    num_dipoles = len(dipole_primitive)
#    print(array_of_particles.shape)
#    print(dipole_primitive.shape)
    forces = np.zeros((num_particles,3),dtype=np.float64)
    inv_polar_unwrap = inverse_polarisation.view(dtype=np.float64).reshape((num_particles*2,1)).flatten()
    #print(inv_polar_unwrap)
    Dipoles.optical_force_array(array_of_particles, num_particles, dipole_radius, dipole_primitive, num_dipoles, inv_polar_unwrap, beam_collection, forces)
    return forces
#==============================================================================
# Wrapper for new optical force array code
#==============================================================================

def py_optical_force_array_precomp(array_of_particles,dipole_radius,dipole_primitive,dpl_moments, beam_collection):
    """
    Wrapper function for the new optical force code.
    """
    num_particles = len(array_of_particles)
    num_dipoles = len(dipole_primitive)
    num_moments = len(dpl_moments)
#    print(array_of_particles.shape)
#    print(dipole_primitive.shape)
    forces = np.zeros((num_particles,3),dtype=np.float64)
    dpl_moment_unwrap = dpl_moments.view(dtype=np.float64).reshape((num_moments,6))#.flatten()
    print(dpl_moment_unwrap)
    print("Number of moments:",num_moments)
    Dipoles.optical_force_array_precomp(array_of_particles, num_particles, dipole_radius, dipole_primitive, num_dipoles, dpl_moment_unwrap, num_moments, beam_collection, forces)
    return forces
#==============================================================================
# Wrapper for new optical force and torque array code
#==============================================================================

def py_optical_force_torque_array(array_of_particles,dipole_radius,dipole_primitive,inverse_polarisation, beam_collection):
    """
    Wrapper function for the new optical force code, including torques (= r X F), and couples (= p X E).
    """
    #print(f"====================\n{array_of_particles}\n====================")
    #print(len(array_of_particles))
    num_particles = len(array_of_particles)
    num_dipoles = len(dipole_primitive)
    num_parpoles = num_particles*num_dipoles
#    print(array_of_particles.shape)
#    print(dipole_primitive.shape)
    forces = np.zeros((num_particles,3),dtype=np.float64)
    torques = np.zeros((num_particles,3),dtype=np.float64)
    couples = np.zeros((num_particles,3),dtype=np.float64)
    dipole_positions = np.zeros((num_parpoles,3),dtype=np.float64)
    inv_polar_unwrap = inverse_polarisation.view(dtype=np.float64).reshape((num_particles*2,1)).flatten() # Split real & imaginary bit into two separate entries
    #print(f"Inv_polar_unwrap: {inv_polar_unwrap}")
    #print(f"dipoles before\n{dipole_positions}")
    Dipoles.optical_force_torque_array(array_of_particles, num_particles, dipole_radius, dipole_primitive, num_dipoles, inv_polar_unwrap, beam_collection, forces, torques, couples, dipole_positions)
    #print("=========================")
    #print(f"dipoles after\n{dipole_positions}")
    return forces, torques, couples, dipole_positions
#
# General dipole function interfaces
#
grad_E_cc = Dipoles.grad_E_cc
grad_E_cc.argtypes = [numpy.ctypeslib.ndpointer(dtype=np.float64, ndim=1, shape=(3), flags='C_CONTIGUOUS'), numpy.ctypeslib.ndpointer(dtype=np.float64, ndim=1, shape=(6), flags='C_CONTIGUOUS'), ctypes.c_double, numpy.ctypeslib.ndpointer(dtype=np.float64, ndim=1, shape=(18), flags='C_CONTIGUOUS'),]

ND_POINTER_1 = np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags="C")
ND_POINTER_2 = np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags="C")

optical_force_array = Dipoles.optical_force_array

optical_force_array.argtypes = [ND_POINTER_2,ctypes.c_int,ctypes.c_double,ND_POINTER_2,ctypes.c_int, ND_POINTER_1,ctypes.POINTER(beams.BEAM_COLLECTION),ND_POINTER_2]

optical_force_array_precomp = Dipoles.optical_force_array_precomp

optical_force_array_precomp.argtypes = [ND_POINTER_2,ctypes.c_int,ctypes.c_double,ND_POINTER_2,ctypes.c_int, ND_POINTER_2,ctypes.c_int,ctypes.POINTER(beams.BEAM_COLLECTION),ND_POINTER_2]

#optical_force_array.restype = ctypes.POINTER(ctypes.c_double)
optical_force_torque_array = Dipoles.optical_force_torque_array

optical_force_torque_array.argtypes = [ND_POINTER_2,ctypes.c_int,ctypes.c_double,ND_POINTER_2,ctypes.c_int, ND_POINTER_1,ctypes.POINTER(beams.BEAM_COLLECTION),ND_POINTER_2,ND_POINTER_2,ND_POINTER_2,ND_POINTER_2]
#
#
###################################################################################
