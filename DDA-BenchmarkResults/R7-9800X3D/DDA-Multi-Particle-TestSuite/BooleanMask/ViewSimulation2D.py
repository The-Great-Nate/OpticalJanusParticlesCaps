# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 12:52:31 2022
Animated trajectories
"""

import sys
import h5py
from operator import pos
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
from cycler import cycler
import cmath
import itertools as it
from numpy.core import numerictypes
from scipy.special import j0, j1, jvp, jv
from numpy import sin, cos, pi, arctan2
import time
from pyoptics import beams
#import Dipoles
from pyoptics import readyaml
import ctypes
import datetime
import matplotlib.animation as animation
#import pyvista as pv
#from pyvistaqt import BackgroundPlotter
import pandas as pd

def init():
    for trajectory in trajectories:
        trajectory.set_data([], [])
    return trajectories


# Function for updating the positions of the meshes and light in the update loop as a callback
def update_scene():
    global frame
    #frame1 = frame
    #frame2 = frame1+frame_interval
    frame_previous = frame
    frame+=frame_interval
    if frame>=frame_max:
        frame=frame-(frame_max-frame_min)
    # We call inplace=True to update the changes to the mesh directly
    for i in range(n_particles):
        spheres[i].translate(particles[frame][i]-particles[frame_previous][i],inplace=True)
    p.add_text(f"Iteration: {frame}", name='time-label')#,color="#FFFFFF")
    # We update the whole plotter
    p.update()



###################################################################################
# Start of program
###################################################################################

if int(len(sys.argv)) != 2:
    sys.exit("Usage: python {} <FILESTEM>".format(sys.argv[0]))

filestem = sys.argv[1]
filename_vtf = filestem+".vtf"
filename_xl = filestem+".xlsx"
filename_hdf = filestem+".h5"
#filename_yaml = filestem+".yml"
#===========================================================================
# Read the yaml file into a system parameter dictionary
#===========================================================================
project = readyaml.Options(filestem)

#sys_params = ReadYAML.load_yaml(filename_yaml)
#print(sys_params)
#===========================================================================
# Parse the sys_params yaml file
#===========================================================================
#beaminfo = ReadYAML.read_section(sys_params,'beams')
#displayinfo = ReadYAML.read_section(sys_params,'display')
#particleinfo = ReadYAML.read_section(sys_params,'particles')
#===========================================================================
# Read simulation parameters (this should be done externally)
#===========================================================================
wavelength = project.parameters.wavelength
dipole_radius = project.parameters.dipole_radius
time_step = project.parameters.time_step
frames = project.simulation.frames
vmd_output = project.output.vmd_output
excel_output = project.output.excel_output
hdf_output = project.output.hdf_output
include_force = project.output.include_force
include_couple = project.output.include_couple
#===========================================================================
# Read beam options and create beam collection
#===========================================================================
beam_collection = project.beam_collection
#===========================================================================
# Read particle options and create particle collection
#===========================================================================
particle_collection = project.particle_collection
print(particle_collection.num_particles)
n_particles = particle_collection.num_particles
#c = 3e8
#n1 = 3.9
#n1a = 1.5
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

#===============================================================
# Read the excel file with spins
#===============================================================
if hdf_output==True:
    h5_file = h5py.File(filename_hdf, 'r')
    print(h5_file.keys())
    particles = h5_file['positions'][()]
    #particles = h5_file['particles'][()]
    h5_file.close()
    print(particles.shape)
    #print(positions.shape)

elif excel_output==True:
    xl=pd.ExcelFile(filename_xl)
    dfl = pd.read_excel(filename_xl,sheet_name=xl.sheet_names[0])
    npdata = dfl.to_numpy()
    print(npdata.shape)
    particles = np.zeros((frames,n_particles,3),dtype=np.float64)
    for i in range(n_particles):
        particles[:,i,0:3] = npdata[:,1+i*3:1+(i+1)*3]

else:
    sys.exit("Unable to find input file:", filename_xl,"or",filename_hdf)

###################################################################################
# Plot the field intensity for the beam configuration
###################################################################################
#fig = plt.figure()
#
frame_interval = project.display.frame_interval
frame_min = project.display.frame_min
frame_max = project.display.frame_max
max_size = project.display.max_size
print("FRAME INTERVAL",frame_interval)
frame = frame_min # initial value
#
lower = -max_size
upper = -lower
nsize = 100
spacing = (upper-lower)/(nsize-1)
ndim = np.array([nsize,nsize,nsize])
origins = np.array([lower,lower,lower])
if beam_collection.BEAM_ARRAY[0].beamtype == beams.BEAMTYPE_BESSEL:
    ndim[2] = nsize // 10
    origins[2] = 0.0
    print(ndim)
if beam_collection.BEAM_ARRAY[0].beamtype == beams.BEAMTYPE_GAUSS_CSP:
    ndim[0] = nsize // 2
    ndim[1] = nsize // 2
    origins[0] = lower / 2
    origins[1] = lower / 2

#grid = plot_intensity_xyz_volume(ndim, origins, spacing, beam_collection)

###################################################################################
# Reverse the order of the array
###################################################################################

positions = np.zeros((n_particles,3,frames))
for i in range(frames):
    for j in range(n_particles):
        for k in range(3):
            positions[j][k][i] = particles[i,j,k]

###################################################################################
# Do the animation
###################################################################################


fig,ax = project.display.plot_intensity(beam_collection)

project.display.animate_particles(fig,ax,positions,radius,colors)
