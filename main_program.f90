PROGRAM MAIN
    USE data_structures
    USE testpotentialgenerator
    USE particlevelocitysolver
    USE domain_tools
    USE particle_write_netcdf
    USE chargedensitygenerator
    USE fieldsolver
    USE command_line
    IMPLICIT NONE
    !TYPE(particle),DIMENSION(:),ALLOCATABLE :: part
    TYPE(particle) :: part
    INTEGER :: timesteps_elapsed
    INTEGER :: ierr
    CHARACTER(LEN=*), PARAMETER :: data_filename = 'particle_simulation_data'
    CHARACTER(LEN = 7) :: problem
    LOGICAL :: success
    !----------------------READ FROM COMMAND LINE----------------------!
    ! Calls and storts all the arguments passed in via the command line
    ! using the module provided by H Ratcliffe
    CALL parse_args
    ! Checks each argument provided is of the right type and ensures that both
    ! lower and upper case variable names are recognised
    ! Returns 'Could not read in <argument> if any of the arguments could not
    ! be read in successfully 
    success = get_arg('nx', Nx) .OR. get_arg('Nx', Nx)
    IF (success) THEN 
        PRINT*, 'Nx :', Nx
    ELSE
        PRINT*, 'Could not read in nx'
        ERROR STOP
    END IF 
    success = get_arg('ny', Ny) .OR. get_arg('Ny', Nx)
    IF (success) THEN 
        PRINT*, 'Ny :', Ny
    ELSE
        PRINT*, 'Could not read in ny, '
        ERROR STOP
    END IF 
    success = get_arg('problem', problem)
    IF (success) THEN 
        PRINT*, 'problem :', TRIM(problem)
    ELSE
        PRINT*, 'Could not read in problem'
        ERROR STOP
    END IF 
    !-----------------END READ-----------------------!
    
    !-------------SAVE META DATA FOR WRITING LATER-----!
    meta_data%namestring = 'build_particle_solver_netcdf'
    WRITE(meta_data%commandlineargs,*) 'Nx=',Nx,' Ny=',Ny,' problem=',&
                                    problem
    !---------------------------------------------------!
    !create axis
    x_axis_range = (/-1.0_REAL64,1.0_REAL64/)
    y_axis_range = (/-1.0_REAL64,1.0_REAL64/)
    nghosts = 0

    !set number of particles
    N_particles = 1
    !IF (N_particles > 1) THEN
        !ALLOCATE(part(1:N_particles),stat=ierr)
        !IF (ierr/=0) stop 'Error allocating particles, check that they have been intialised correctly'
    !END IF

    !set initial positions and velocities of particles
    !part%position = (/0.1_REAL64,-0.1_REAL64/)
    !part%velocity = (/0.0_REAL64,0.01_REAL64/)
    !part(2)%position = (/0.1_REAL64,0.1_REAL64/)
    !part(2)%velocity = (/0.0_REAL64,-0.01_REAL64/)

    !find dx,dy
    dx = (x_axis_range(2)-x_axis_range(1))/(REAL(Nx,kind=REAL64))
    dy = (y_axis_range(2)-y_axis_range(1))/(REAL(Ny,kind=REAL64))

    !how long to run simulation for (timesteps)
    total_time = 1000

    !how long is each timestep
    dt = 0.01
    !Allocate the charge density grid
    ALLOCATE(rho(1:Nx,1:Ny))
    !generate charge density
    IF (problem == 'null') THEN
        CALL generate_null_charge_density
        part%position = (/0.0_REAL64,0.0_REAL64/)
        part%velocity = (/0.1_REAL64,0.1_REAL64/)
    ELSE IF (problem == 'single') THEN
        CALL generate_single_charge_density
        part%position = (/0.1_REAL64,0.0_REAL64/)
        part%velocity = (/0.0_REAL64,0.0_REAL64/)
    ELSE IF (problem == 'double') THEN
        CALL generate_double_charge_density
        part%position = (/0.0_REAL64,0.5_REAL64/)
        part%velocity = (/0.0_REAL64,0.0_REAL64/)
    ELSE 
        PRINT*, problem, 'is not a valid input'
        ERROR STOP
    END IF
    !solve for the potential
    CALL solve_gauss_seidel()
    !initialise the file to write in
    CALL initialise_file(data_filename,ierr)
    !generate constant gradient potential (left to right)
    !CALL generate_const_grad_potential(-0.1_REAL64)
    !get electric field
    CALL get_field()
    !Write charge density, potential and field to file
    CALL write_grid_vars(data_filename,ierr)

    !initialise particle(s) (takes array of particles or single particle)
    CALL initialise_particles(part)



    timesteps_elapsed = 0

    !Speeding up program by not putting if statement into loop.
    DO WHILE((.NOT.(particle_out_of_bounds)).AND.(timesteps_elapsed<=total_time))
        timesteps_elapsed = timesteps_elapsed + 1
        !Call the propagate function to advance all particles one timestep
        !takes either an array of particles or a single particle
        CALL propagate(part)
        !Write the particle data to file. part can be either a single particle
        !or an array of particles (an interface is used to allow either.)
        CALL write_to_file(data_filename,timesteps_elapsed,part,ierr)
    END DO

    !print final velocity,position and acceleration
    PRINT*,'                   x component                 y component'
    PRINT*,'final velocity',part%prev_velocity
    PRINT*,'final position',part%prev_position
    PRINT*,'final acceleration',part%prev_acceleration
    PRINT*,'final timesteps',timesteps_elapsed

    !Write metadata to file
    CALL write_metadata(data_filename,ierr)
END PROGRAM