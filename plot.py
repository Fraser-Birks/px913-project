import netCDF4 as NC
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
from matplotlib.ticker import (StrMethodFormatter, AutoLocator, AutoMinorLocator)
from matplotlib.animation import FuncAnimation

# read in data from NetCDF file
dat = NC.Dataset("particle_simulation_data","r", format="NETCDF4")

# read variable values for the various plots from dat
part_x = dat.variables['part_x'][:]
part_y = dat.variables['part_y'][:]
Vx = dat.variables['Vx'][:]
Vy = dat.variables['Vy'][:]
acc_x = dat.variables['ax'][:]
acc_y = dat.variables['ay'][:]

cell_x = dat.variables['cell_x'][:]
cell_y = dat.variables['cell_y'][:]
Ex = dat.variables['Ex'][:]
rho = dat.variables['charge_density'][:]

# array for major ticks to ensure they range from -1 to 1
ticks = np.linspace(-1,1,5)


# plot consisting of 4 subplots
fig,axs = plt.subplots(2,2)
fig.suptitle('Various outputs from particle moving in electric field simulation')

# subplot1 - particle positions
# get the positions, velcoties and accelerations for 20 evenly spaced points
x = part_x[:,0][0:-1:250]
y = part_y[:,0][0:-1:250]
vel_x = Vx[:,0][0:-1:250]
vel_y = Vy[:,0][0:-1:250]
acc_u = acc_x[:,0][0:-1:250]
acc_v = acc_y[:,0][0:-1:250]
# plot particle path
axs[0,0].plot(part_x[:,0],part_y[:,0], label=r'position', alpha=0.8)
# plot velocitiy and accerlation for 20 points as arrows
axs[0,0].quiver(x,y,vel_x,vel_y, width=0.004, headwidth=4, headlength=6, headaxislength=5, color='green', label=r'velocity')
axs[0,0].quiver(x,y,acc_u,acc_v, width=0.004, headwidth=4, headlength=6, headaxislength=5, color='purple', label=r'acceleration')
axs[0,0].scatter(x,y, s=6, c='black')
# customise graph 
axs[0,0].set_title('Particle path', size=10)
axs[0,0].set_xlabel(r'$x$ [arbitrary units]', size=9)
# axs[0,0].set_xlim(cell_x[0], cell_x[-1])
axs[0,0].xaxis.set_major_locator(AutoLocator())
axs[0,0].xaxis.set_minor_locator(AutoMinorLocator())
axs[0,0].tick_params(axis='x', labelsize=7)
axs[0,0].set_ylabel(r'$y$ [arbitrary units]', size=9)
# axs[0,0].set_ylim(cell_y[0], cell_y[-1])
axs[0,0].yaxis.set_major_locator(AutoLocator())
axs[0,0].yaxis.set_minor_locator(AutoMinorLocator())
axs[0,0].tick_params(axis='y', labelsize=7)
axs[0,0].legend(fontsize=6)


# subplot 2 - animation of particle position
#Note that this animation code draws heavily on the 
#code at the URLs below. It does not copy any of them exactly, but is instead a
#mix of all of them.
#https://stackoverflow.com/questions/29832055/animated-subplots-using-matplotlib
#https://brushingupscience.com/2016/06/21/matplotlib-animations-the-easy-way/
#https://stackoverflow.com/questions/17212722/matplotlib-imshow-how-to-animate
#https://matplotlib.org/stable/api/animation_api.html#

intervaltime = 10 #ms, how long in ms between each animation frame

#Plot initial graphs:
#First graph is a 2D line of the particle displacement from the beginning of the simulation to the current timestep
#Second graph is a scatter point of the particle position as it moves through the simulation.
time_axis = np.linspace(0, len(part_x[:,0])-1, len(part_x[:,0])-1).astype(int)
path, = axs[0,1].plot(part_x[:,0][0], part_y[:,0][0], label='particle path')
scatpoint, = axs[0,1].plot(part_x[:,0][0], part_y[:,0][0], marker=".", label='particle')

#Set axes labels and legend etc.
axs[0,1].set_title('Animation of particle path', size=10)
axs[0,1].set_xlim(cell_x[0], cell_x[-1])
axs[0,1].set_xlabel(r'$x$ [arbitrary units]', size=9)
axs[0,1].set_xticks(ticks)
axs[0,1].xaxis.set_minor_locator(AutoMinorLocator())
axs[0,1].tick_params(axis='x', labelsize=7)
axs[0,1].set_ylim(cell_x[0], cell_x[-1])
axs[0,1].set_ylabel(r'$y$ [arbitrary units]', size=9)
axs[0,1].set_yticks(ticks)
axs[0,1].yaxis.set_minor_locator(AutoMinorLocator())
axs[0,1].tick_params(axis='y', labelsize=7)
axs[0,1].legend(fontsize=6)

#Create a list which holds both the things we wish to animate, both the
#2D grid above and the scattered point.
graphlines = [path,scatpoint]


def update(t):
    #Define an update function for the animation. This is what is called each frame
    #to update the graph. As blitting is set to on (it has to be to use reasonable computing power)
    #This function needs to return a single array containing the objects to animate with new data set.
    #Function takes as an input argument the timestep it is being called at.
    #Set the data of each element of graphlines to the corresponding value of t passed to the function
    graphlines[0].set_data(part_x[:,0][0:t], part_y[:,0][0:t])
    graphlines[1].set_data(part_x[:,0][t], part_y[:,0][t]) 
    return graphlines


#Create an animation. Arguments to pass are:
#figure that holds the axes
#the update function to be called each timestep
#interval - the time in milliseconds each between each frame
#frames - an array to iterate over, with each value in the array being passed as the
#sole argument of the update function the each animation frame as the animation progresses.
#In this case, this array is the time axis, which holds values between 1 and the number of timesteps.
#Blit - this is set to true here to turn on blitting - something which makes animations much
#less compuationally expensive. Read more about blitting: https://matplotlib.org/stable/api/animation_api.html#
ani = FuncAnimation(fig, update, interval=intervaltime, frames=time_axis, blit=True)


# subplot 3 - Ex
# pcolor plot with cell_x and cell_y as x and y axes
pc = axs[1,1].pcolor(cell_x, cell_y, Ex)
# customise plot 
axs[1,1].set_title('Ex component of the electric field', size=10)
axs[1,1].set_xlabel(r'$x$ [arbitrary units]', size=9)
axs[1,1].set_xlim(cell_x[0], cell_x[-1])
axs[1,1].set_xticks(ticks)
axs[1,1].xaxis.set_minor_locator(AutoMinorLocator())
axs[1,1].tick_params(axis='x', labelsize=7)
axs[1,1].set_ylabel(r'$y$ [arbitrary units]', size=9)
axs[1,1].set_ylim(cell_y[0], cell_y[-1])
axs[1,1].set_yticks(ticks)
axs[1,1].yaxis.set_minor_locator(AutoMinorLocator())
axs[1,1].tick_params(axis='y', labelsize=7)
# customise colorbar
ticklist = np.linspace(Ex.min(), Ex.max(), 8)
cmap = mpl.cm.viridis
cbar = axs[1,1].figure.colorbar(pc, ax=axs[1,1], cmap=cmap, ticks=ticklist, label='Strength of electric field')
cbar.ax.tick_params(labelsize=7)
cbar.ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.3f}"))


# subplot 4 - charge density rho
# pcolor plot with cell_x and cell_y as x and y axes
pc = axs[1,0].pcolor(cell_x, cell_y, rho)
# customise plot
axs[1,0].set_title('Charge density of the field', size=10)
axs[1,0].set_xlabel(r'$x$ [arbitrary units]', size=9)
axs[1,0].set_xlim(cell_x[0], cell_x[-1])
axs[1,0].set_xticks(ticks)
axs[1,0].xaxis.set_minor_locator(AutoMinorLocator())
axs[1,0].tick_params(axis='x', labelsize=7)
axs[1,0].set_ylabel(r'$y$ [arbitrary units]',size=9)
axs[1,0].set_ylim(cell_y[0], cell_y[-1])
axs[1,0].set_yticks(ticks)
axs[1,0].yaxis.set_minor_locator(AutoMinorLocator())
axs[1,0].tick_params(axis='y', labelsize=7)
# customise colorbar
ticklist = np.linspace(rho.min(), rho.max(), 8)
cmap = mpl.cm.viridis
cbar = axs[1,0].figure.colorbar(pc, ax=axs[1,0], cmap=cmap, ticks=ticklist, label='Charge density')
cbar.ax.tick_params(labelsize=7)
cbar.ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.3f}"))


# plot the graph with tight layout, show and save as .png
plt.tight_layout()
plt.show()
# plt.savefig('Mini_project_visualisation.png')
# plt.close()
