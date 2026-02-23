

class ParameterObject (object):

    def __init__(self,paraminfo,paramdefaults):
        self.wavelength = float(self.set_param_value(paraminfo,paramdefaults,'wavelength'))
        self.dipole_radius = float(self.set_param_value(paraminfo,paramdefaults,'dipole_radius'))
        self.time_step = float(self.set_param_value(paraminfo,paramdefaults,'time_step'))
        #self.spring_stiffness = float(self.set_param_value(paraminfo,paramdefaults,'spring_stiffness'))
        self.bending_stiffness = float(self.set_param_value(paraminfo,paramdefaults,'bending_stiffness'))
        self.wall_position = float(self.set_param_value(paraminfo,paramdefaults,'wall_position'))
        self.wall_force_max = float(self.set_param_value(paraminfo,paramdefaults,'wall_force_max'))
        self.wall_order = float(self.set_param_value(paraminfo,paramdefaults,'wall_order'))
        self.wall_zrange = float(self.set_param_value(paraminfo,paramdefaults,'wall_zrange'))


    def set_param_value(self,paraminfo,paramdefaults,name):
        if paraminfo==None:
            param_value = paramdefaults[name]
        else:
            param_value = paraminfo.get(name,paramdefaults[name])
        return param_value
