"""
Particle creation routines
"""
import numpy as np

class ParticleCollection (object):
    num_particles = 0
    
    def __init__(self,particleinfo,particledefaults,particlespecs):
        if particleinfo==None:
            # Set defaults
            self.default_material = particledefaults['default_material']
            self.default_radius = float(particledefaults['default_radius'])
            self.default_density = float(particledefaults['default_density'])
            # Create default particle array (one particle)
            ParticleCollection.num_particles = 1
            self.particle_type = np.asarray([self.default_material])
            self.particle_radius = np.asarray([self.default_radius])
            self.particle_indices = np.asarray([self.get_particle_spec(particlespecs,self.default_material,'refractive_index')],dtype=complex)
            self.particle_colour = np.asarray([self.get_particle_spec(particlespecs,self.default_material,'display_colour')])
            self.particle_vtfcolour = np.asarray([self.get_particle_spec(particlespecs,self.default_material,'vtf_colour_name')])
            self.particle_density = np.asarray([self.default_density])
            self.particle_positions = np.zeros((1,3),dtype=float)
        else:
            # Read from file
            self.default_material = particleinfo.get('default_material',particledefaults['default_material'])
            self.default_radius = float(particleinfo.get('default_radius',particledefaults['default_radius']))
            self.default_density = float(particleinfo.get('default_density',particledefaults['default_density']))
            self.particle_list = particleinfo.get('particle_list',None)
            if self.particle_list==None or self.particle_list==False:
                # Set defaults
                ParticleCollection.num_particles = 1
                self.particle_type = np.asarray([self.default_material])
                self.particle_radius = np.asarray([self.default_radius])
                self.particle_density = np.asarray([self.default_density])
                self.particle_indices = np.asarray([self.get_particle_spec(particlespecs,self.default_material,'refractive_index')],dtype=complex)
                self.particle_colour = np.asarray([self.get_particle_spec(particlespecs,self.default_material,'display_colour')])
                self.particle_vtfcolour = np.asarray([self.get_particle_spec(particlespecs,self.default_material,'vtf_colour_name')])
                self.particle_positions = np.zeros((1,3),dtype=float)
            else:
                # Read individual particles
                i=0
                self.particle_type = []
                self.particle_radius = []
                self.particle_colour = []
                self.particle_vtfcolour = []
                self.particle_positions = []
                self.particle_indices = []
                self.particle_density = []
                for newparticle in self.particle_list:
                    particle = self.particle_list[newparticle]
                    print("Loading particle",particle)
                    if particle != None:
                        self.particle_type.append(particle.get('material',self.default_material))
                        self.particle_radius.append(float(particle.get('radius',self.default_radius)))
                        self.particle_density.append(float(particle.get('density',self.default_density)))
                        self.altcolour = bool(particle.get('altcolour',False))
                        if self.altcolour==False:
                            self.particle_colour.append(self.get_particle_spec(particlespecs,self.particle_type[i],'display_colour'))
                        else:
                            self.particle_colour.append(self.get_particle_spec(particlespecs,self.particle_type[i],'alt_display_colour'))
                        self.particle_vtfcolour.append(self.get_particle_spec(particlespecs,self.particle_type[i],'vtf_colour_name'))
                        self.particle_indices.append(self.get_particle_spec(particlespecs,self.particle_type[i],'refractive_index'))
                        self.coords = particle.get('coords',"0.0 0.0 0.0")
                        self.fields = self.coords.split(" ")
                        if self.fields[0]=="None":
                            self.particle_positions.append(np.array((0.0,0.0,0.0),dtype=np.float64))
                        else:
                            self.particle_positions.append(np.array((0.0,0.0,0.0),dtype=np.float64))
                            for j in range(min(len(self.fields),3)):
                                self.particle_positions[i][j] = float(self.fields[j])
                    else:
                        # not sure this part is ever used
                        self.particle_type.append(self.default_material)
                        self.particle_radius.append(self.default_radius)
                        self.particle_density.append(self.default_density)
                        self.particle_indices.append(self.get_particle_spec(particlespecs,self.default_material,'refractive_index'))
                        self.particle_colour.append(self.get_particle_spec(particlespecs,self.default_material,'display_colour'))
                        self.particle_vtfcolour.append(self.get_particle_spec(particlespecs,self.default_material,'vtf_colour_name'))
                        self.particle_positions.append(np.array((0.0,0.0,0.0),dtype=np.float64))
                    i+=1
                ParticleCollection.num_particles = i

    def get_particle_spec(self,particle_specs,name,property):
        """
        Returns the values from the new particle spec dictionary.
        """
        value = particle_specs[name][property]
        if property == 'refractive_index':
            return_value = complex(value)
        elif property == 'density':
            return_value = float(value)
        else:
            return_value = value
        return return_value

    def get_refractive_indices(self):
        return np.asarray(self.particle_indices,dtype=complex)
        
    def get_particle_types(self):
        return np.asarray(self.particle_type)
        
    def get_particle_colours(self):
        return np.asarray(self.particle_colour)
        
    def get_particle_vtfcolours(self):
        return np.asarray(self.particle_vtfcolour)
        
    def get_particle_radii(self):
        return np.asarray(self.particle_radius,dtype=float)

    def get_particle_density(self):
        return np.asarray(self.particle_density,dtype=float)

    def get_particle_positions(self):
        return np.asarray(self.particle_positions,dtype=float).reshape((ParticleCollection.num_particles,3))
    
