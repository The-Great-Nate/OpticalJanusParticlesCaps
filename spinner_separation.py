import numpy as np
seperations = np.linspace(1,100,12)

for separator in seperations:
    filename = f"OPTICAL_spinner_separation{separator}.yml"
    
    content = f"""simulation:
  frames: 20000
  include_dynamics: True
  dynamics_method: COUPLED_ROTATION
  include_springs: False
  include_bending: False
  include_nonbonded: True
  include_optics: True
  include_rotation: True
  optics_method: SH_DDA
parameters:
  wavelength: 0.785e-6 # meters in vacuum
  dipole_radius: 45e-9 # meters # meters
  time_step: 1e-5 # seconds
output:
  vmd_output: False
  excel_output: True
  hdf_output: True
  include_force: True
  include_couple: True
  include_orientation: True
display:
  show_output: True
  frame_interval: 10
  max_size: 20e-6 # range will be 2 times this
  resolution: 201 # number of points in each direction of plot
  frame_min: 0 # starting frame for animation
  frame_max: 20000 # will default to number of frames
  z_offset: 0.0e-6
beams:
  beam_1:
    beamtype: BEAMTYPE_RICHARDS_WOLF_ANALYTIC
    E0: 2.5e7
    kt_by_kz: 0.3
    order: 0
    gouy: 0
    w0: 1.0
    sigma: 1.0
    zernike: 0 0 0
    nm: 1.333
    translation: -{separator}e-6 0 0
    rotation: 0 0
    jones: POLARISATION_LCP
    NA: 1.25
    numkpoints: 100
  beam_2:
    beamtype: BEAMTYPE_RICHARDS_WOLF_ANALYTIC
    E0: 2.5e7
    kt_by_kz: 0.3
    order: 0
    gouy: 0
    w0: 1.0
    sigma: 1.0
    zernike: 0 0 0
    nm: 1.333
    translation: {separator}e-6 0 0
    rotation: 0 0
    jones: POLARISATION_RCP
    NA: 1.25
    numkpoints: 100
particles:
  default_radius: 2e-07
  default_material: FusedSilica
  particle_list:
    particle_0:
      material1: Air
      material2: Air
      coords: 0.0 -4.0e-6 0.0
    particle_1:
      material1: Silicon
      material2: Air
      coords: -{separator} 0.0 0.0
    particle_2:
      material1: Silicon
      material2: Air
      coords: {separator} 0.0 0.0
"""

    with open(filename, "w") as f:
        f.write(content)

print("YAML files generated successfully.")