#======================================================================
# Topology Object Class
# Created 5-Apr-25 SH
# Generates the connections for bonds and bends within a simulated
# system.  Bonds connect pairs of beads, either through springs or
# constraints.  Bends connect linear triplets of connected beads with
# bending springs.
# Topology can be stored in separate file and read in as needed.
# Also produces lists of faces for flexible membrane-like systems.
#======================================================================

import sys
import numpy as np

class TopologyObject (object):
    """
    Loads the topology section of the YAML file if it exists.  Loads 
    individual values and assigns defaults if no value given.
    """

    def __init__(self,topologyinfo,topologydefaults,particle_collection):
        #print(topologydefaults)
        #print(topologyinfo)
        #
        # Type of topology to use, chosen from list
        # {DIMERS, TRIMERS, TETRAHEDRA, LINEAR, GENERAL, FROM_FILE}
        #
        self.topotype = str(self.set_topology_value(topologyinfo,topologydefaults,'topotype'))
        #
        # Default spring stiffness to use for any connected particles
        #
        self.default_spring_stiffness = float(self.set_topology_value(topologyinfo,topologydefaults,'default_spring_stiffness'))
        #
        # Default bending stiffness to use for any connected particles
        #
        self.default_bending_stiffness = float(self.set_topology_value(topologyinfo,topologydefaults,'default_bending_stiffness'))
        #
        # Particle separation in constraints method - all particles same separation,
        # given, not determined by initial positions.
        #
        self.particle_separation = float(self.set_topology_value(topologyinfo,topologydefaults,'particle_separation'))
        #
        numparts = particle_collection.num_particles
        print(numparts,"particles defined in collection")
        
        if self.topotype == 'None' or self.topotype == 'NONE':
            # Nothing defined
            self.bond_constraints = None
            self.bond_connections = None
            self.bond_stiffness = None
            self.bond_lengths = None
            self.angle_connections = None
            self.angle_stiffness = None
            self.angle_values = None
            self.faces = None
            print("No topology defined")

        elif self.topotype == 'DIMERS':
            # bond constraints
            # bond connections
            # bond stiffnesses
            # bond lengths
            self.angle_connections = None
            self.angle_stiffness = None
            self.angle_values = None
            self.faces = None
            #
            print("Defining DIMER connections")
            self.bond_constraints = []
            if numparts//2 == numparts/2:
                self.bond_connections = np.zeros((numparts//2,2),dtype=np.int64)
                self.bond_lengths = np.zeros(numparts//2,dtype=np.float64)
                self.bond_stiffness = np.zeros(numparts//2,dtype=np.float64)
                j = 0
                for i in range(0,numparts,2):
                    self.bond_connections[j,0] = i
                    self.bond_connections[j,1] = i + 1
                    self.bond_stiffness[j] = self.default_spring_stiffness
                    self.bond_lengths[j] = np.linalg.norm(particle_collection.particle_positions[i]-particle_collection.particle_positions[i+1])
                    self.bond_constraints.append([[i,i+1],[i+1,i]])
                    j += 1
            else:
                sys.exit("Error defining topology: \'Number of particles must be even for DIMERS topology\'.")

        elif self.topotype == 'TRIMERS':
            # bond constraints
            # bond connections
            # bond stiffnesses
            # bond lengths
            # faces
            self.angle_connections = None
            self.angle_stiffness = None
            self.angle_values = None
            #
            self.bond_constraints = []
            print("Defining TRIMER connections")
            if numparts//3 == numparts/3:
                self.bond_connections = np.zeros((numparts,2),dtype=np.int64)
                self.bond_lengths = np.zeros(numparts,dtype=np.float64)
                self.bond_stiffness = np.zeros(numparts,dtype=np.float64)
                self.faces = np.zeros((numparts//3,3),dtype=np.int64)
                #
                for i in range(0,numparts,3):
                    self.bond_constraints.append([[i,i+1],[i+1,i]])
                    self.bond_constraints.append([[i,i+2],[i+2,i]])
                    self.bond_constraints.append([[i+1,i+2],[i+2,i+1]])
                    #self.faces[i//3,0] = i
                    #self.faces[i//3+1,0] = i+1
                    #self.faces[i//3+2,0] = i+2
                    for j in range(3):
                        self.faces[i//3+j,j] = i + j
                        self.bond_connections[i+j,0] = i + j
                        self.bond_connections[i+j,1] = i + (j + 1)%3
                        self.bond_stiffness[i+j] = self.default_spring_stiffness
                        self.bond_lengths[i+j] = np.linalg.norm(particle_collection.particle_positions[i+j]-particle_collection.particle_positions[i+(j+1)%3])
            else:
                sys.exit("Error defining topology: \'Number of particles must be multiple of 3 for TRIMERS topology\'.")

        elif self.topotype == 'TETRAHEDRA':
            # bond constraints
            # bond connections
            # bond stiffnesses
            # bond lengths
            # faces
            self.angle_connections = None
            self.angle_stiffness = None
            self.angle_values = None
            #
            self.bond_constraints = []
            print("Defining TETRAHEDRAL connections")
            if numparts//4 == numparts/4:
                self.bond_connections = np.zeros(((3*numparts)//2,2),dtype=np.int64)
                self.bond_lengths = np.zeros((3*numparts)//2,dtype=np.float64)
                self.bond_stiffness = np.zeros((3*numparts)//2,dtype=np.float64)
                self.faces = np.zeros((numparts,3),dtype=np.int64)
                #
                l = 0
                for i in range(0,numparts,4):
                    for j in range(3):
                        for k in range(j,3):
                            #print(i,j,k,i+j,i+k+1)
                            self.bond_constraints.append([[i+j,i+k+1],[i+k+1,i+j]])
                            self.bond_connections[l,0] = i + j
                            self.bond_connections[l,1] = i + k + 1
                            self.bond_lengths[l] = np.linalg.norm(particle_collection.particle_positions[i+j] - particle_collection.particle_positions[i+k+1])
                            self.bond_stiffness[l] = self.default_spring_stiffness
                            l += 1
                    #
                    for j in range(4):
                        self.faces[i//4+j,0] = i+j
                        self.faces[i//4+j,1] = i+(j+1)%4
                        self.faces[i//4+j,2] = i+(j+2)%4
            else:
                sys.exit("Error defining topology: \'Number of particles must be multiple of 4 for TETRAHEDRA topology\'.")

        elif self.topotype == 'LINEAR':
            # bond constraints
            # bond connections
            # bond stiffnesses
            # bond lengths
            # angle connections
            # angle stiffnesses
            # angle values
            self.faces = None
            #
            self.bond_constraints = []
            print("Defining LINEAR connections")

            if numparts > 2:
                self.bond_connections = np.zeros((numparts-1,2),dtype=np.int64)
                self.bond_lengths = np.zeros(numparts-1,dtype=np.float64)
                self.bond_stiffness = np.zeros(numparts-1,dtype=np.float64)
                self.angle_connections = np.zeros((numparts-2,3),dtype=np.int64)
                self.angle_values = np.zeros(numparts-2,dtype=np.float64)
                self.angle_stiffness = np.zeros(numparts-2,dtype=np.float64)
                #
                for i in range(numparts-1):
                    self.bond_connections[i,0] = i
                    self.bond_connections[i,1] = i + 1
                    self.bond_stiffness[i] = self.default_spring_stiffness
                    self.bond_lengths[i] = np.linalg.norm(particle_collection.particle_positions[i]-particle_collection.particle_positions[i+1])
                    self.bond_constraints.append([[i,i+1],[i+1,i]])

                for i in range(numparts-2):
                    self.angle_connections[i,0] = i
                    self.angle_connections[i,1] = i + 1
                    self.angle_connections[i,2] = i + 2
                    self.angle_stiffness[i] = self.default_bending_stiffness
                    vec1 = particle_collection.particle_positions[i] - particle_collection.particle_positions[i+1]
                    vec2 = particle_collection.particle_positions[i+2] - particle_collection.particle_positions[i+1]
                    nvec1 = np.linalg.norm(vec1)
                    nvec2 = np.linalg.norm(vec2)
                    costheta = np.dot(vec1,vec2)/(nvec1*nvec2)
                    self.angle_values[i] = np.arccos(costheta)
                    #
            else:
                sys.exit("Error defining topology: \'Number of particles must be > 2 for LINEAR topology\'.")

        elif self.topotype == 'GENERAL' or self.topotype == 'FROM_FILE':
            # Optionally:
            #   bond connections
            #   bond stiffnesses
            #   bond lengths
            # Optionally:
            #   angle connections
            #   angle stiffnesses
            #   angle values
            # Optionally:
            #   faces
            #
            self.bond_constraints = None
            print("Defining GENERAL connections")
            
            if self.topotype == 'FROM_FILE':
                #
                # Load the topology here
                #
                print("check")
            
            #print(topologyinfo['connections'])
            #
            # Check we have some connections
            #
            query_bonds = topologyinfo.get('connections',None)
            if query_bonds != None:
                numconnects = len(query_bonds)
                #numconnects = len(topologyinfo['connections'])
                self.bond_connections = np.zeros((numconnects,2),dtype=np.int64)
                self.bond_lengths = np.zeros(numconnects,dtype=np.float64)
                self.bond_stiffness = np.zeros(numconnects,dtype=np.float64)
                #
                i = 0
                connects = topologyinfo['connections']
                for connect in connects:
                    # split into fields and read two connections, a stiffness and a length
                    self.fields = connects[connect].split(" ")
                    if len(self.fields) >= 2:
                        self.bond_connections[i,0] = int(self.fields[0])
                        self.bond_connections[i,1] = int(self.fields[1])
                    else:
                        sys.exit("Error defining topology: \'Connection info incorrect:\'  {:s}.".format(connect,connects[connect]))
                    if len(self.fields) >= 3:
                        self.bond_stiffness[i] = float(self.fields[2])
                    if len(self.fields) == 4:
                        self.bond_lengths[i] = float(self.fields[3])
                    i+=1
                #print(self.bond_connections)
            else:
                self.bond_connections = None
                self.bond_stiffness = None
                self.bond_lengths = None
                print("No bond information found.")
            #
            #
            # Check we have some angles
            #
            query_angles = topologyinfo.get('angles',None)
            if query_angles != None:
                numangles = len(query_angles)
                #numangles = len(topologyinfo['angles'])
                self.angle_connections = np.zeros((numangles,3),dtype=np.int64)
                self.angle_values = np.zeros(numangles,dtype=np.float64)
                self.angle_stiffness = np.zeros(numangles,dtype=np.float64)
                #
                i = 0
                angles = topologyinfo['angles']
                for angle in angles:
                    # split into fields and read three connections, a stiffness and an angle
                    self.fields = angles[angle].split(" ")
                    if len(self.fields) >= 3:
                        self.angle_connections[i,0] = int(self.fields[0])
                        self.angle_connections[i,1] = int(self.fields[1])
                        self.angle_connections[i,2] = int(self.fields[2])
                    else:
                        sys.exit("Error defining topology: \'Angle info incorrect:\'  {:s}.".format(angle,angles[angle]))
                    if len(self.fields) >= 4:
                        self.angle_stiffness[i] = float(self.fields[3])
                    if len(self.fields) == 5:
                        self.angle_values[i] = float(self.fields[4])
                    i+=1
                print(self.angle_connections)
            else:
                self.angle_connections = None
                self.angle_stiffness = None
                self.angle_values = None
                print("No angle information found.")
            #
            #
            # Check we have some faces
            #
            query_faces = topologyinfo.get('faces',None)
            if query_faces != None:
                numfaces = len(query_faces)
                self.faces = np.zeros((numfaces,3),dtype=np.int64)
                #
                i = 0
                faceinfo = topologyinfo['faces']
                for face in faceinfo:
                    # split into fields and read three connections
                    self.fields = faceinfo[face].split(" ")
                    if len(self.fields) == 3:
                        self.faces[i,0] = int(self.fields[0])
                        self.faces[i,1] = int(self.fields[1])
                        self.faces[i,2] = int(self.fields[2])
                    else:
                        sys.exit("Error defining topology: \'Face info incorrect:\'  {:s}.".format(face,faceinfo[face]))
                    i+=1
                #print(self.faces)
            else:
                self.faces = None
                print("No face information found.")

        elif self.topotype == 'FROM_FILE':
            # bond connectivity
            # angle constraints
            # faces
            # bond constraints
            a=1
        else:
            # topology not recognised
            #print("all OK")
            sys.exit("Error reading topology: \'Topology not recognised\':\n  check you are using one of {DIMERS, TRIMERS, TETRAHEDRA, LINEAR, GENERAL, FROM_FILE}.")

    def set_topology_value(self,topologyinfo,topologydefaults,name):
        if topologyinfo==None:
            topology_value = topologydefaults[name]
        else:
            topology_value = topologyinfo.get(name,topologydefaults[name])
        return topology_value

    # methods to return:
    # (a) list of bonds: <#0> <#1> <length> [<stiffness>]
    # (b) list of angles: <#0> <#1> <#2> <angle> [<stiffness>] ! angle not used yet
    # (c) list of faces: <#0> <#1> <#2>
    # (d) list of constraints: ((<#0> <#1>),(<#1> <#0>)) <length>

    def get_bond_connections(self):
        return self.bond_connections

    def get_bond_stiffness(self):
        return self.bond_stiffness

    def get_separations(self):
        return self.bond_lengths
        
    def get_constraints(self):
        return self.bond_constraints

    def get_angle_connections(self):
        return self.angle_connections

    def get_angle_stiffness(self):
        return self.angle_stiffness

    def get_angles(self):
        return self.angle_values

    def get_faces(self):
        return self.faces
