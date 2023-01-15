# px913-project
Fraser and Chantal px913 project. <br> <br>
Distribution of work is: <br> 

Chantal: Generation of potential grid<br>
Fraser: Particle motion solver and data writing<br>
Chantal: Visualisation<br>

## Quick run:
To run the basic simulation with the 'single' charge density format and visualisation, please just run the build script ./build_particle_solver_netcdf.
It is possible to change runtime parameters in this buildscript to load the 'double' or 'null' formats.
Additionally, if you navigate to the folder 'multiple_particles' the run_multi_particle build script in there allows for the simulation and visualisation of a grid of particles moving in the 'double' regime.

## Helper files provided:
Note that the modules for writing to command line (command_line.f90) and the module for generating the grid axis which is used in the file writing (create_axis.f90) are not original work, and have been provided by Heather and Chris.

## Note on collaboration:
Whilst the broad distribution of work is that layed out above, both of the authors collaborated together on much of the code.