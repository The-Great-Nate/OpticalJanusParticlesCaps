"""
Make a display object for animations etc
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation

from pyoptics import beams

class DisplayObject (object):

    def __init__(self,displayinfo,displaydefaults,frames):
        self.show_output = bool(self.set_display_value(displayinfo,displaydefaults,'show_output'))
        self.max_size = float(self.set_display_value(displayinfo,displaydefaults,'max_size')) # range will be 2 times this
        self.resolution = int(self.set_display_value(displayinfo,displaydefaults,'resolution')) # number of points in each direction of plot
        self.z_offset = float(self.set_display_value(displayinfo,displaydefaults,'z_offset')) # this is the z value for the intensity plot
        self.frame_min = abs(int(self.set_display_value(displayinfo,displaydefaults,'frame_min'))) # starting frame for animation
        self.frame_max = min(frames,int(self.set_display_value(displayinfo,displaydefaults,'frame_max')))
        self.frame_interval = int(self.set_display_value(displayinfo,displaydefaults,'frame_interval'))
    def rotate_vector_by_quaternion(self, v, q):
        w = q[0]
        u = q[1:4]

        uv = np.cross(u, v)
        uuv = np.cross(u, uv)

        return v + 2 * (w * uv + uuv)

    def set_display_value(self,displayinfo,displaydefaults,name):
        if displayinfo==None:
            display_value = displaydefaults[name]
        else:
            display_value = displayinfo.get(name,displaydefaults[name])
        return display_value

    def plot_intensity(self, beam_collection):
        nx = self.resolution
        ny=nx
        Ex = np.zeros((nx, ny), dtype=complex)
        Ey = np.zeros((nx, ny), dtype=complex)
        Ez = np.zeros((nx, ny), dtype=complex)
        z = self.z_offset
        #I = []
        E = np.zeros(3,dtype=np.complex128)
        #fig, ax = plt.subplots(1, num_plots)
        upper = self.max_size
        lower = -upper
        x = np.linspace(lower, upper, nx)
        y = np.linspace(lower, upper, ny)
        X, Y = np.meshgrid(x, y)
        for j in range(ny):
            for i in range(nx):
                beams.all_incident_fields((x[i], y[j], self.z_offset), beam_collection, E)
                Ex[j][i] = E[0]
                Ey[j][i] = E[1]
                Ez[j][i] = E[2]

        I = np.square(np.abs(Ex)) + np.square(np.abs(Ey)) + np.square(np.abs(Ez))
        print(np.shape(I))
#        for j in range(ny):
#            print(j,I[0][100][j])
        I0 = np.max(I)
#            ax.axis('equal')

        fig = plt.figure()
        ax = plt.axes(xlim=(lower, upper), ylim=(lower, upper))

        ax.set_aspect('equal','box')
        #cs=ax.contourf(X, Y, I, cmap=cm.viridis, levels=30)
        #cs=ax.imshow(I[k],cmap=cm.summer)
        extents = (lower,upper,lower,upper)
        cs=ax.imshow(I,cmap=cm.Blues,vmin=0.0,vmax=I0,origin="lower",extent=extents)
        ax.set_aspect('equal','box')
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        cbar = fig.colorbar(cs)
        # ax.set_title("z = {:.1e}".format(z[k]))
        return fig,ax


    def plot_parpoles(self,fig,ax,particle_positions,dipole_positions,radius,colors, dipole_above_zero,dipole_below_zero, allqs):
        n_particles = len(colors)
        marker_size_particles = (1.25*10.0*(radius/200e-9)*(5e-6/self.max_size))
        self.particles_positions = particle_positions

        self.particles_trajectories = []
        for i in range(n_particles):
            marker = ax.plot([], [], marker="o", markersize=marker_size_particles,
                            c=colors[i], mec='white', mew=0.75, alpha=1,
                            animated=True)[0]
            self.particles_trajectories.append(marker)

        n_frames, n_dipoles, _ = dipole_positions.shape
        n_particles_dipoles = len(colors)
        self.dipoles_positions = dipole_positions
        self.dipoles_per_particle = n_dipoles // n_particles_dipoles
        marker_size_dipoles = ((1.25*10.0*(radius/200e-9)*(5e-6/self.max_size))/4)

        self.dipoles_trajectories = []
        for i in range(n_particles_dipoles):
            for die in range(self.dipoles_per_particle):
                colors = colors
                if die in dipole_above_zero:
                    col = colors[i]
                elif die in dipole_below_zero:
                    col = colors[i]
                marker = ax.plot([], [], marker="s", markersize=marker_size_dipoles,
                                c=col, alpha=1, animated=True, linestyle="None")[0]
                self.dipoles_trajectories.append(marker)
        axis_length = radius * 2
        axis_colors = ['r', 'g', 'b']

        self.body_axes = []
        self.title_text = ax.text(0.5, 1.02, '', transform=ax.transAxes, ha='center', va='bottom', fontsize=14, animated=True)


        for p in range(n_particles):
            axes_for_particle = []
            for c in axis_colors:
                line = ax.plot([], [], lw=1.5, c=c, animated=True)[0]
                axes_for_particle.append(line)
            self.body_axes.append(axes_for_particle)
        
        def init_anim():
            for trajectory in self.particles_trajectories + self.dipoles_trajectories:
                trajectory.set_data([], [])
            self.title_text.set_text('')
            return [self.title_text] + self.particles_trajectories + self.dipoles_trajectories

        def animate(fff):
            p = 0
            frames = self.frame_min + fff * self.frame_interval
            start = max(frames - 2, 0)
            end = frames
            real_n_particles = int(n_particles/2) #because n_particles is actually the number of colours now
            
            self.title_text.set_text(f'Timestep = {frames}')
            for p in range(real_n_particles):
                q = allqs[frames, p]  

                
                ex = np.array([1.,0.,0.])
                ey = np.array([0.,1.,0.])
                ez = np.array([0.,0.,1.])

                ex_r = self.rotate_vector_by_quaternion(ex, q)
                ey_r = self.rotate_vector_by_quaternion(ey, q)
                ez_r = self.rotate_vector_by_quaternion(ez, q)

                pos = particle_positions[p][:, frames]

                rotated_axes = [ex_r, ey_r, ez_r]

                for i in range(3):
                    arrow_start = pos
                    arrow_end = pos + axis_length * rotated_axes[i]
                    self.body_axes[p][i].set_data(
                        [arrow_start[0], arrow_end[0]],
                        [arrow_start[1], arrow_end[1]]
                    )
            
            for trajectory, particle in zip(self.particles_trajectories, self.particles_positions):
                trajectory.set_data(
                    particle[0, frames - 2 : frames], particle[1, frames - 2 : frames]
                )
            
            for i, trajectory in enumerate(self.dipoles_trajectories):
                particle_idx = i // self.dipoles_per_particle
                sub_idx = i % self.dipoles_per_particle
                points = self.dipoles_positions[start:end,
                                                particle_idx * self.dipoles_per_particle + sub_idx,
                                                :]
                x = points[:, 0].reshape(-1)
                y = points[:, 1].reshape(-1)
                trajectory.set_data(x, y)

            return (self.particles_trajectories + self.dipoles_trajectories + [axis for particle_axes in self.body_axes for axis in particle_axes])

        ani = animation.FuncAnimation(fig, animate, init_func=init_anim,
                                    frames=(self.frame_max - self.frame_min)//self.frame_interval,
                                    interval=25, blit=True)

        return ani


    def animate_dipoles(self,fig,ax,positions,radius,colors):
        
        n_frames, self.dipoles, _ = positions.shape
        self.positions = positions
        particles = len(colors)
        self.dipoles_per_particle = self.dipoles // particles 
        
        self.trajectories = []
        marker_size = ((1.25*10.0*(radius/200e-9)*(5e-6/self.max_size))/4)*8 # Size 10 for 200nm radius and upper limit 5 microns,
        for i in range(particles):
            for die in range(self.dipoles):
                marker = ax.plot([], [], marker="s", markersize=marker_size, c=colors[i], alpha=0.1, animated=True, linestyle="None",)[0]
                self.trajectories.append(marker)
        
        
        #self.particles = positions
        #self.trajectories = [ax.plot([], [], markersize=marker_size, marker="o", c=colors[i], mec='white', mew=0.75, alpha=1, animated=True)[0] for i in np.arange(n_particles)]

        ani = animation.FuncAnimation(fig, self.animate_dip, init_func=self.init_anim, frames=(self.frame_max-self.frame_min) // self.frame_interval, interval=25, blit=True)
        return ani
        #plt.show()
# writer = animation.PillowWriter(fps=30)

# ani.save("bessel-ang-mom-test.gif", writer=writer)
# =====================================

    def init_anim(self):
        for trajectory in self.trajectories:
            trajectory.set_data([], [])
        return self.trajectories
    
    def animate_dip(self, fff):

        frame = self.frame_min + fff * self.frame_interval
        start = max(frame - 2, 0)
        end = frame

        for i, trajectory in enumerate(self.trajectories):
            start_idx = i * self.dipoles_per_particle
            end_idx   = start_idx + self.dipoles_per_particle
            points = self.positions[start:end, start_idx:end_idx, :] 
            x = points[:,:, 0]
            y = points[:,:, 1]

            trajectory.set_data(x, y)

        return self.trajectories

    '''
    def animate_dip(self,fff):
        frames = self.frame_min + fff * self.frame_interval
        #slice_from = max(frames-2,0)
        #slice_to = frames
        for trajectory, particle in zip(self.trajectories, self.particles):
            trajectory.set_data(
                particle[0, frames - 2 : frames], particle[1, frames - 2 : frames]
            )
        return self.trajectories
    '''

    def animate_particles(self,fig,ax,positions,radius,colors):
            n_particles = len(colors)
            marker_size = (1.25*10.0*(radius/200e-9)*(5e-6/self.max_size)) # Size 10 for 200nm radius and upper limit 5 microns.
            self.particles = positions
            self.trajectories = [ax.plot([], [], markersize=marker_size, marker="o", c=colors[i], mec='white', mew=0.75, alpha=0.5, animated=True)[0] for i in np.arange(n_particles)]

            ani = animation.FuncAnimation(fig, self.animate, init_func=self.init_anim, frames=(self.frame_max-self.frame_min) // self.frame_interval, interval=25, blit=True)
            return ani
            #plt.show()
    # writer = animation.PillowWriter(fps=30)

    # ani.save("bessel-ang-mom-test.gif", writer=writer)
    # =====================================


    def init_anim(self):
        for trajectory in self.trajectories:
            trajectory.set_data([], [])
        return self.trajectories


    def animate(self,fff):
        frames = self.frame_min + fff * self.frame_interval
        for trajectory, particle in zip(self.trajectories, self.particles):
            trajectory.set_data(
                particle[0, frames - 2 : frames], particle[1, frames - 2 : frames]
            )
        return self.trajectories





