"""
Service functions for reading the parameters in
the yaml file for a given simulation.
"""
import os
import yaml
import sys

import pyoptics
from pyoptics import simulation
from pyoptics import output
from pyoptics import parameters
from pyoptics import display
from pyoptics import beams
from pyoptics import particles
from pyoptics import topology

def load_yaml(filename_yaml):
    yamlpath = './'+filename_yaml
    check_file = os.path.isfile(yamlpath)
    if not check_file:
        print("Unable to find configuration file: ",yamlpath)
        sys.exit()
    with open(yamlpath, 'r') as yamlfile:
        sys_params = yaml.safe_load(yamlfile)
    return sys_params
    
def load_defaults(filename_yaml):
    yamlpath = pyoptics.__path__[0]+'/'+filename_yaml
    check_file = os.path.isfile(yamlpath)
    if not check_file:
        print("Unable to find defaults file: ",yamlpath)
        sys.exit()
    with open(yamlpath, 'r') as yamlfile:
        sys_defaults = yaml.safe_load(yamlfile)
    return sys_defaults
    
def load_particle_types(filename_yaml):
    yamlpath = pyoptics.__path__[0]+'/'+filename_yaml
    check_file = os.path.isfile(yamlpath)
    if not check_file:
        print("Unable to find defaults file: ",yamlpath)
        sys.exit()
    with open(yamlpath, 'r') as yamlfile:
        sys_defaults = yaml.safe_load(yamlfile)
    return sys_defaults
    

def read_section(sys_params, section):
    sectioninfo = sys_params.get(section,None)
    return sectioninfo


class Options(object):
    num_simulations = 0

    def __init__(self, filestem):
        filename_default = "defaults.yml"
        filename_particle_types = "particle_types.yml"
        filename_yaml = filestem+".yml"
        self.name = filename_yaml
        Options.num_simulations += 1
        sys_params = load_yaml(filename_yaml)
        sys_defaults = load_defaults(filename_default)
        self.beaminfo = read_section(sys_params,'beams')
        self.beamdefaults = read_section(sys_defaults,'beams')
        self.paraminfo = read_section(sys_params,'parameters')
        self.paramdefaults = read_section(sys_defaults,'parameters')
        self.simulationinfo = read_section(sys_params,'simulation')
        self.simulationdefaults = read_section(sys_defaults,'simulation')
        self.displayinfo = read_section(sys_params,'display')
        self.displaydefaults = read_section(sys_defaults,'display')
        self.outputinfo = read_section(sys_params,'output')
        self.outputdefaults = read_section(sys_defaults,'output')
        self.particleinfo = read_section(sys_params,'particles')
        self.particledefaults = read_section(sys_defaults,'particles')
        self.particletypes = load_particle_types(filename_particle_types)
        self.topologyinfo = read_section(sys_params,'topology')
        self.topologydefaults = read_section(sys_defaults,'topology')

        self.simulation = simulation.SimulationObject(self.simulationinfo,self.simulationdefaults)
        self.parameters = parameters.ParameterObject(self.paraminfo,self.paramdefaults)
        self.output = output.OutputObject(self.outputinfo,self.outputdefaults)
        self.display = display.DisplayObject(self.displayinfo,self.displaydefaults,self.simulation.frames)
        self.beam_collection = beams.create_beam_collection(self.beaminfo,self.parameters.wavelength)
        self.particle_collection = particles.ParticleCollection(self.particleinfo,self.particledefaults,self.particletypes)
        #print(self.particletypes)
        self.topology = topology.TopologyObject(self.topologyinfo,self.topologydefaults, self.particle_collection)


    def __del__(self):
        Options.num_simulations -= 1
    
    
