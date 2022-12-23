# px913-project
Fraser and Chantal px913 project. <br> <br>
Distribution of work is: <br> 

Chantal: Generation of potential grid<br>
Fraser: Particle motion solver and data writing<br>
Chantal: Visualisation<br>

####Work completed so far<br>
##Essential work:<br>
Field solver: A field solver has been written which takes a provided potential grid and turns it into an electric field.<br>
Particle solver: A particle solver has been implemented which uses velocity verlet to propagate one (or more) particles a single timestep. This can then be called from the main program.<br>
File writing: A file writer has been created that has functions to create and write to a netcdf file which holds all of the necessary information for the entire simulation (about the charge density, potential, field and particle positions, velocities and accelerations). An interface was created for the module which allows for the same function call to be able to write information about a single particle or an array of particles.<br>
Data structures: A module was created which holds all of the necessary global data used in all the simulation modules. This file is USED everywhere in the code.<br>
<br>
##Non-essential work.
Test potential generator: A file was written that could generate 2 simple numerical potentials for testing; one which was constant everywhere and one which constantly increases from left to write. Particles were observed to move in parabolas in this potential, which is the expected result.<br>
<br>
<br>
##Helper files provided:
Note that the modules for writing to command line (command_line.f90) and the module for generating the grid axis which is used in the file writing (create_axis.f90) are not original work, and have been provided by Heather and Chris.
