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
import pyvista as pv
from pyvistaqt import BackgroundPlotter
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


def plot_intensity_xyz_volume(ndims, origins, spacing, beam_collection):
#    iscale = 5
#    nx = ndim
#    ny = ndim
#    nz = ndim//iscale
    Ex = np.zeros((ndims[2],ndims[1],ndims[0]), dtype=complex)
    Ey = np.zeros((ndims[2],ndims[1],ndims[0]), dtype=complex)
    Ez = np.zeros((ndims[2],ndims[1],ndims[0]), dtype=complex)
#    I = []
    E = np.zeros(3,dtype=np.complex128)

    x = np.linspace(origins[0], origins[0]+spacing*ndims[0], ndims[0])
    y = np.linspace(origins[1], origins[1]+spacing*ndims[1], ndims[1])
    z = np.linspace(origins[2], origins[2]+spacing*ndims[2], ndims[2])
#    sp = (upper-lower)/(ndim-1)
    grid = pv.ImageData(
        dimensions=ndims,
        spacing=(spacing,spacing,spacing),
        origin=origins
    )
    
    for ij in range(ndims[2]):
        for j in range(ndims[1]):
            for i in range(ndims[0]):
                beams.all_incident_fields((x[i], y[j], z[ij]), beam_collection, E)
                Ex[ij][j][i] = E[0]
                Ey[ij][j][i] = E[1]
                Ez[ij][j][i] = E[2]

    I = np.square(np.abs(Ex)) + np.square(np.abs(Ey)) + np.square(np.abs(Ez))

    I0 = np.max(I)
    I = I / I0
    grid['intensity'] = I.flatten()

    return grid


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
fig = plt.figure()
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

grid = plot_intensity_xyz_volume(ndim, origins, spacing, beam_collection)

###################################################################################
# Read the atom types and set up the colours
###################################################################################
"""
"""

print("# Particle types: ",particle_types)
print("# Refractive indices: ", ref_ind)
marker_size = 1.25*10.0*(radius/200e-9)*(5e-6/upper) # Size 10 for 200nm radius and upper limit 5 microns.

###################################################################################
# Do the animation
###################################################################################


if __name__ == '__main__':
    # We create a plotter for non-blocking visualization, using pyvistaqt. If we don't need interactivity or animations we can just call pyvista.Plotter
    p = BackgroundPlotter()#window_size=(600,400))
    
    # We set the camera position to be in the xy plane
    p.camera_position = 'xy'
    # We set the front and back clipping planes
    p.camera.clipping_range = (0, 10)

    # We create the beads
    spheres=[]
    for i in range(n_particles):
        spheres.append(pv.Sphere(radius = radius*1.2, center=particles[frame][i]))
        p.add_mesh(spheres[i], color=colors[i], roughness=0.2, diffuse=1, opacity=1.0)#0.99)

    ## We add the laser volume
    ##p.add_volume(grid,name='laser',
    ##cmap="reds",#gray",#inferno",
    ##opacity=[0,0.25])#"sigmoid_3")#'linear')
    
    opaque = 0.75
    single_slice_x = grid.slice(normal=[0,1,0],origin=[0,0,0])
    p.add_mesh(single_slice_x,opacity=opaque)#, cmap=cmap)
    #single_slice_z1 = grid.slice(normal=[0,0,1],origin=[0,0,-upper/2])
    #p.add_mesh(single_slice_z1,opacity=opaque)#, cmap=cmap)
    single_slice_z2 = grid.slice(normal=[0,0,1],origin=[0,0,0])
    p.add_mesh(single_slice_z2,opacity=opaque)#, cmap=cmap)
    #single_slice_z3 = grid.slice(normal=[0,0,1],origin=[0,0,upper/2])
    #p.add_mesh(single_slice_z3,opacity=opaque)#, cmap=cmap)
#slices = grid.slice_along_axis(n=5, axis="z")
#p.add_mesh(slices,opacity=0.75)

#p = pv.Plotter()
    p.add_mesh(grid.outline(), color="k")
# Add a beam direction arrow
    arrows = pv.Arrow(start=[0,0,lower*1.25],direction=[0,0,1],scale=2*upper*1.25,shaft_radius=0.005,tip_radius=0.02,tip_length=0.125)
    #p.add_mesh(arrows,color="FFFFFF")
    xm = 7e-6
    #wall = pv.Cube(bounds=(-xm,xm,-xm,xm,2e-6,3e-6))
    #p.add_mesh(wall,color="FF4444",opacity=0.5)

    
    # We create a callback for the plotter and we set the function that will be run in the update loop. We also set the interval of the update
    # It is important to set the interval explicitly.
    p.add_callback(update_scene, interval=50)

    # We show the plotter and call p.app.exec_(), so the plotter stays open. This is important when running the pyvistaqt plotter.
    # If your visualization automatically closes add the last line
    p.remove_scalar_bar()
    #p.background_color = 'k'
    p.show()
    p.show_axes()
    p.camera.parallel_projection = True
    p.view_isometric()# forces the view
    p.camera.zoom(1.5)

    p.app.exec_()

