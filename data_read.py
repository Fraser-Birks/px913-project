import netCDF4 as NC
import matplotlib.pyplot as plt
dat = NC.Dataset("particle_simulation_data","r", format="NETCDF4")


part_x = dat.variables['part_x'][:]
part_y = dat.variables['part_y'][:]
Vx = dat.variables['Vx'][:]
Vy = dat.variables['Vy'][:]
ax = dat.variables['ax'][:]
ay = dat.variables['ay'][:]

time_axis = dat.variables['time'][:] #grab time axes

print(part_x, part_y)
plt.plot(part_x[:,0],part_y[:,0])
#plt.plot(part_x[:,1],part_y[:,1])
plt.show()