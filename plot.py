import netCDF4 as NC
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
from matplotlib.ticker import (StrMethodFormatter, AutoMinorLocator)

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
x = part_x[:,0][0:-1:500]
y = part_y[:,0][0:-1:500]
vel_x = Vx[:,0][0:-1:500]
vel_y = Vy[:,0][0:-1:500]
acc_u = acc_x[:,0][0:-1:500]
acc_v = acc_y[:,0][0:-1:500]
# plot particle path
axs[0,0].plot(part_x[:,0],part_y[:,0], label=r'position', alpha=0.8)
# plot velocitiy and accerlation for 20 points as arrows
axs[0,0].quiver(x,y,vel_x,vel_y, width=0.004, headwidth=4, headlength=6, headaxislength=5, color='green', label=r'velocity')
axs[0,0].quiver(x,y,acc_u,acc_v, width=0.004, headwidth=4, headlength=6, headaxislength=5, color='purple', label=r'acceleration')
axs[0,0].scatter(x,y, s=6, c='black')
# customise graph 
axs[0,0].set_title('Particle path', size=10)
axs[0,0].set_xlabel(r'$x$ [arbitrary units]', size=9)
axs[0,0].set_xlim(cell_x[0], cell_x[-1])
axs[0,0].set_xticks(ticks)
axs[0,0].xaxis.set_minor_locator(AutoMinorLocator())
axs[0,0].tick_params(axis='x', labelsize=7)
axs[0,0].set_ylabel(r'$y$ [arbitrary units]', size=9)
axs[0,0].set_ylim(cell_y[0], cell_y[-1])
axs[0,0].set_yticks(ticks)
axs[0,0].yaxis.set_minor_locator(AutoMinorLocator())
axs[0,0].tick_params(axis='y', labelsize=7)
axs[0,0].legend(fontsize=6)


# subplot 2 - Ex
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
cbar = axs[1,1].figure.colorbar(pc, ax=axs[1,1], cmap=cmap, ticks=ticklist)
cbar.ax.tick_params(labelsize=7)
cbar.ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.3f}"))


# subplot 3 - charge density rho)
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
cbar = axs[1,0].figure.colorbar(pc, ax=axs[1,0], cmap=cmap, ticks=ticklist)
cbar.ax.tick_params(labelsize=7)
cbar.ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.3f}"))


# plot the graph with tight layout, show and save as .png
plt.tight_layout()
plt.show()
plt.savefig('Mini_project_visualisation.png', format='png')
plt.close()
