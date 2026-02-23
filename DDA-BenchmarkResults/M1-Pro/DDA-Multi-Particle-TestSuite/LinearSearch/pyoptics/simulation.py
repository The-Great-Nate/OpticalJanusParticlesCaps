

class SimulationObject (object):

    def __init__(self,simulationinfo,simulationdefaults):
        self.frames = int(self.set_simulation_value(simulationinfo,simulationdefaults,'frames'))
        self.include_dynamics = bool(self.set_simulation_value(simulationinfo,simulationdefaults,'include_dynamics'))
        self.dynamics_method = str(self.set_simulation_value(simulationinfo,simulationdefaults,'dynamics_method'))
        self.hi_method = str(self.set_simulation_value(simulationinfo,simulationdefaults,'hi_method'))
        self.include_springs = bool(self.set_simulation_value(simulationinfo,simulationdefaults,'include_springs'))
        self.include_bending = bool(self.set_simulation_value(simulationinfo,simulationdefaults,'include_bending'))
        self.include_driving = bool(self.set_simulation_value(simulationinfo,simulationdefaults,'include_driving'))
        self.include_wall = bool(self.set_simulation_value(simulationinfo,simulationdefaults,'include_wall'))
        self.include_gravity = bool(self.set_simulation_value(simulationinfo,simulationdefaults,'include_gravity'))
        self.include_nonbonded = bool(self.set_simulation_value(simulationinfo,simulationdefaults,'include_nonbonded'))
        self.include_optics = bool(self.set_simulation_value(simulationinfo,simulationdefaults,'include_optics'))
        self.optics_method = str(self.set_simulation_value(simulationinfo,simulationdefaults,'optics_method'))

    def set_simulation_value(self,simulationinfo,simulationdefaults,name):
        if simulationinfo==None:
            simulation_value = simulationdefaults[name]
        else:
            simulation_value = simulationinfo.get(name,simulationdefaults[name])
        return simulation_value
