MODULE particle_write_netcdf
  !Note that the syntax of this module especially the report and exit statements
  !are taken from the provided example, which itself is based on  https://people.sc.fsu.edu/~jburkardt/f_src/netcdf/netcdf.html
  !which is accessible using https://web.archive.org/web/20190623025346/http://people.sc.fsu.edu/~jburkardt/f_src/netcdf/netcdf.html (accessed in Nov 2021)
  USE data_structures
  USE domain_tools
  USE netcdf
  IMPLICIT NONE


  INTERFACE write_to_file
    MODULE PROCEDURE write_single_particle_to_file

    MODULE PROCEDURE write_particle_array_to_file
  END INTERFACE

  CONTAINS
  SUBROUTINE initialise_file(filename, ierr)
    !Subroutine for file initialisation. Creates a file, defines the dimensions and axes
    !and defines the variables that the data will be written to and which dimensions they
    !correspond to.

    !Function takes as input arguments the name of the file it is creating,
    !ierr as an integer to store error values 

    INTEGER,INTENT(OUT) :: ierr
    CHARACTER(LEN=*),INTENT(IN) :: filename

    !The file needs to have 9 dimensions, x, y, part_x, part_y, vx, vy, ax, ay and time, so ndims = 9
    INTEGER, PARAMETER :: ndims = 4
    
    !griddims holds the length of each dimension for easy initialisation
    INTEGER, DIMENSION(ndims)  :: griddims


    !dim_ids holds the ids of dimensions once they have been initialised.
    INTEGER, DIMENSION(ndims) :: dim_ids
    
    !The names of the dimensions used in the file
    CHARACTER(LEN=1), DIMENSION(ndims) :: dims=(/"x", "y", "p","t" /)
    

    INTEGER, DIMENSION(N_particles) :: particle_axis


    INTEGER :: file_id, i, var_id, var_id_x, var_id_y, var_id_part

    !The length of each dimension. In this case, x is of length Nx
    !y is of length Ny and t is the record dimension (has unlimited length).
    !the particles axis has length equal to the number of particles and
    !allows values of position, velocity, acceleration etc to be written for multiple particles.

    !The t dimension can be extended as long as the program is run,
    !giving additional flexibility.
    griddims = (/Nx,Ny,N_particles,NF90_UNLIMITED/)

    !In the following bits of code, if any errors are generated, they are caught
    !and the subroutine ends with the error printed to console.

    !Create the file
    ierr = nf90_create(filename, NF90_CLOBBER, file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN 
    END IF

    !Define each grid size dimension, x, y and t, and their corresponding sizes.
    DO i = 1,ndims
      ierr = nf90_def_dim(file_id, dims(i), griddims(i), dim_ids(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO

    !Define axis variables - this is cell_x,cell_y,t and particles here.
    !keep hold of x and y axis variable ids for writing, others do not matter.
    !x
    ierr = nf90_def_var(file_id, 'cell_x', NF90_DOUBLE, dim_ids(1), var_id_x)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !y
    ierr = nf90_def_var(file_id, 'cell_y', NF90_DOUBLE, dim_ids(2), var_id_y)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !particles
    ierr = nf90_def_var(file_id, 'particles', NF90_INT, dim_ids(3), var_id_part)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF


    !Define 1D array to hold time axes, it grows each step
    ierr = nf90_def_var(file_id, 'time', NF90_INT, dim_ids(4), var_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF


    !Define variables to hold the required in the grid.

    !define variable to hold charge density along x and y 
    ierr = nf90_def_var(file_id, 'charge_density', NF90_DOUBLE, dim_ids(1:2), var_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !define variable to hold potential along x and y 
    ierr = nf90_def_var(file_id, 'potential', NF90_DOUBLE, dim_ids(1:2), var_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !define variable to hold Ex along x and y
    ierr = nf90_def_var(file_id, 'Ex', NF90_DOUBLE, dim_ids(1:2), var_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !define variable to hold Ey along x and y
    ierr = nf90_def_var(file_id, 'Ey', NF90_DOUBLE, dim_ids(1:2), var_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !define variables to hold the x, y positions, Vx, Vy velocities and ax, ay accelerations for each particle, 
    !using particles and time dimensions

    !part_x
    ierr = nf90_def_var(file_id, 'part_x', NF90_DOUBLE, dim_ids(3:4), var_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF


    !part_y
    ierr = nf90_def_var(file_id, 'part_y', NF90_DOUBLE, dim_ids(3:4), var_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF


    !Vx
    ierr = nf90_def_var(file_id, 'Vx', NF90_DOUBLE, dim_ids(3:4), var_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !Vy
    ierr = nf90_def_var(file_id, 'Vy', NF90_DOUBLE, dim_ids(3:4), var_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !ax
    ierr = nf90_def_var(file_id, 'ax', NF90_DOUBLE, dim_ids(3:4), var_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !ay
    ierr = nf90_def_var(file_id, 'ay', NF90_DOUBLE, dim_ids(3:4), var_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !Put file into write mode
    ierr = nf90_enddef(file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !Create the x and y axes along corresponding grid dimensions
    !Create full x axis array:
    !(note that we need to set n_ghosts to 0 
    !as we are only saving axis values in the bulk of the array.)
    CALL create_axis(x_axis,Nx,x_axis_range)
    !Create full y axis array
    CALL create_axis(y_axis,Ny,y_axis_range)
    !Create particle axis
    DO i = 1,N_particles
      particle_axis(i) = i
    END DO

    !write axes to file
    !write x
    ierr = nf90_put_var(file_id, var_id_x, x_axis)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !write y
    ierr = nf90_put_var(file_id, var_id_y, y_axis)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !write particle axis
    ierr = nf90_put_var(file_id, var_id_part, particle_axis)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ! Close the file
    ierr = nf90_close(file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
  END SUBROUTINE

  SUBROUTINE write_single_particle_to_file(filename,curr_timestep,single_part,ierr)

  !This function opens the file, writes the particle position, velocity, acceleration
    !and the current corresponding timestep to the file.

    !This function takes the following arguments:
    !The current timestep - which is the value to write to the time axis
    !ierr - an integer to hold error codes.
    !part - an array of particle types
    !filename - the filename to write to.
    
    INTEGER, INTENT(IN) :: curr_timestep
    INTEGER, INTENT(OUT) :: ierr
    TYPE(particle),INTENT(IN) :: single_part
    CHARACTER(LEN=*),INTENT(IN) :: filename


    INTEGER :: file_id
    INTEGER,DIMENSION(7) :: var_ids
    !The different start and count variables are to hold 
    !the positions to start writing in each dimension (start)
    !and how many datapoints to write in each dimension (count).
    INTEGER,DIMENSION(2) :: start_write_part
    INTEGER,DIMENSION(1) :: startT,countT,timestepArr
   
    !As above, if any errors are returned by netcdf, they are caught and printed.
    !open file to write
    ierr = nf90_open(filename,NF90_WRITE,file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

        
    !Fetch variable ids of all variables we need to write to in the file,
    !such that we can write to them in next step.

    !Fetch ID of the variable holding the particle x position
    ierr = nf90_inq_varid(file_id,'part_x',var_ids(1))
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !Fetch ID of the variable holding the particle y position
    ierr = nf90_inq_varid(file_id,'part_y',var_ids(2))
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !Fetch ID of the variable holding the particle x velocity
    ierr = nf90_inq_varid(file_id,'Vx',var_ids(3))
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !Fetch ID of the variable holding the particle y velocity
    ierr = nf90_inq_varid(file_id,'Vy',var_ids(4))
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !Fetch ID of the variable holding the particle x acceleration
    ierr = nf90_inq_varid(file_id,'ax',var_ids(5))
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !Fetch ID of the variable holding the particle y acceleration
    ierr = nf90_inq_varid(file_id,'ay',var_ids(6))
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !fetch the time axis variable
    ierr = nf90_inq_varid(file_id,'time',var_ids(7))
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !Define the point to start writing in the file (first item on each array for the corresponding slice timestep)
    start_write_part = (/1,curr_timestep/)

    !write particle x
    ierr = nf90_put_var(file_id, var_ids(1), single_part%position(1), start_write_part)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !write particle y
    ierr = nf90_put_var(file_id, var_ids(2), single_part%position(2), start_write_part)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !write particle x velocity
    ierr = nf90_put_var(file_id, var_ids(3), single_part%velocity(1), start_write_part)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF


    !write particle y velocity
    ierr = nf90_put_var(file_id, var_ids(4), single_part%velocity(2), start_write_part)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

      !write particle x acceleration
    ierr = nf90_put_var(file_id, var_ids(5), single_part%acceleration(1), start_write_part)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !write particle y acceleration
    ierr = nf90_put_var(file_id, var_ids(6), single_part%acceleration(2), start_write_part)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF


    !Extend the time axis by the value corresponding to the current timestep.
    !To write to the file, the timestep has to be in array form,
    !turn it into a single element array
    timestepArr = (/curr_timestep/)

    !Define where to start and how many counts to write in the file, once again these
    !need to be arrays. Start writing at the point in the list corresponding to
    !timestep, and write a single value.

    startT = (/curr_timestep/)
    countT = (/1/)
    ierr = nf90_put_var(file_id, var_ids(7),timestepArr,startT,countT)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !close the file
    ierr = nf90_close(file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

  END SUBROUTINE

  SUBROUTINE write_particle_array_to_file(filename,curr_timestep,part_arr,ierr)
    !This function opens the file, writes the particle position, velocity, acceleration
    !and the current corresponding timestep to the file.

    !This function takes the following arguments:
    !The current timestep - which is the value to write to the time axis
    !ierr - an integer to hold error codes.
    !part - an array of particle types
    !filename - the filename to write to.

    INTEGER, INTENT(IN) :: curr_timestep
    INTEGER, INTENT(OUT) :: ierr
    TYPE(particle),INTENT(IN),DIMENSION(:) :: part_arr
    CHARACTER(LEN=*),INTENT(IN) :: filename

    INTEGER :: file_id, loop
    INTEGER,DIMENSION(7) :: var_ids
    !The different start and count variables are to hold 
    !the positions to start writing in each dimension (start)
    !and how many datapoints to write in each dimension (count).
    INTEGER,DIMENSION(2) :: start_write_part
    INTEGER,DIMENSION(1) :: startT,countT,timestepArr
   
    !As above, if any errors are returned by netcdf, they are caught and printed.
    !open file to write
    ierr = nf90_open(filename,NF90_WRITE,file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

        
    !Fetch variable ids of all variables we need to write to in the file,
    !such that we can write to them in next step.

    !Fetch ID of the variable holding the particle x position
    ierr = nf90_inq_varid(file_id,'part_x',var_ids(1))
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !Fetch ID of the variable holding the particle y position
    ierr = nf90_inq_varid(file_id,'part_y',var_ids(2))
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !Fetch ID of the variable holding the particle x velocity
    ierr = nf90_inq_varid(file_id,'Vx',var_ids(3))
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !Fetch ID of the variable holding the particle y velocity
    ierr = nf90_inq_varid(file_id,'Vy',var_ids(4))
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !Fetch ID of the variable holding the particle x acceleration
    ierr = nf90_inq_varid(file_id,'ax',var_ids(5))
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !Fetch ID of the variable holding the particle y acceleration
    ierr = nf90_inq_varid(file_id,'ay',var_ids(6))
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !fetch the time axis variable
    ierr = nf90_inq_varid(file_id,'time',var_ids(7))
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !Define the point to start writing in the file (first item on each array for the corresponding slice timestep)
    !Define the count (how many items to write in each dimension, in this case the entire slice for a single timestep)
    DO loop = 1,N_particles
      !write for all particles at the timestep in question.
      start_write_part = (/loop,curr_timestep/)

      !write particle x
      ierr = nf90_put_var(file_id, var_ids(1), part_arr(loop)%position(1), start_write_part)
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF

      !write particle y
      ierr = nf90_put_var(file_id, var_ids(2), part_arr(loop)%position(2), start_write_part)
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF

      !write particle x velocity
      ierr = nf90_put_var(file_id, var_ids(3), part_arr(loop)%velocity(1), start_write_part)
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF


      !write particle y velocity
      ierr = nf90_put_var(file_id, var_ids(4), part_arr(loop)%velocity(2), start_write_part)
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF

      !write particle x acceleration
      ierr = nf90_put_var(file_id, var_ids(5), part_arr(loop)%acceleration(1), start_write_part)
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF

      !write particle y acceleration
      ierr = nf90_put_var(file_id, var_ids(6), part_arr(loop)%acceleration(2), start_write_part)
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO

    !Extend the time axis by the value corresponding to the current timestep.
    !To write to the file, the timestep has to be in array form,
    !turn it into a single element array
    timestepArr = (/curr_timestep/)

    !Define where to start and how many counts to write in the file, once again these
    !need to be arrays. Start writing at the point in the list corresponding to
    !timestep, and write a single value.

    startT = (/curr_timestep/)
    countT = (/1/)
    ierr = nf90_put_var(file_id, var_ids(7),timestepArr,startT,countT)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !close the file
    ierr = nf90_close(file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

  END SUBROUTINE
  
  SUBROUTINE write_metadata(filename,rundata,ierr)
    !This subroutine exists purely to write in the file metadata as global file attributes
    !It takes the type rundata as an argument, which should contain within it
    !all of the metadata the user wishes to write to the file.

    !filename is the name of the file to write to
    !ierr is a value to hold the error codes
    !rundata is of the type program_metadata, and holds metadata the user defines.
    CHARACTER(LEN=*),INTENT(IN) :: filename
    TYPE(program_metadata), INTENT(IN) :: rundata
    INTEGER, INTENT(OUT):: ierr
    INTEGER :: file_id
    !open file to write
    ierr = nf90_open(filename,NF90_WRITE,file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !Place file into define mode, such that each attribute can be added globally.
    ierr = nf90_redef(file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !Write each rundata metadata component to the file.
    !----ADD ADDITIONAL METADATA WRITING STEPS HERE IF rundata has more attributes added!-----
    !use NF90_GLOBAL to write these as global attributes.

    !Write generating file name and version to file
    ierr = nf90_put_att(file_id, NF90_GLOBAL,'name_and_version',TRIM(rundata%namestring))
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !add full string of command line arguments used to file.
    ierr = nf90_put_att(file_id, NF90_GLOBAL,'command_line_args',TRIM(rundata%commandlineargs))
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ! Close the file
    ierr = nf90_close(file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
  END SUBROUTINE

END MODULE