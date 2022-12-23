PROGRAM MAIN
    USE data_structures
    USE testpotentialgenerator
    USE particlevelocitysolver
    USE domain_tools
    USE particle_write_netcdf
    IMPLICIT NONE
    !TYPE(particle),DIMENSION(:),ALLOCATABLE :: part
    TYPE(particle) :: part
    INTEGER :: timesteps_elapsed
    INTEGER :: ierr
    CHARACTER(LEN=*), PARAMETER :: data_filename = 'particle_simulation_data'

    !set grid dimensions
    Nx = 20
    Ny = 20
    
    !create axis
    x_axis_range = (/-1.0_REAL64,1.0_REAL64/)
    y_axis_range = (/-1.0_REAL64,1.0_REAL64/)
    nghosts = 0
    !set number of particles
    N_particles = 1
    !ALLOCATE(part(N_particles))
    !find dx,dy
    dx = (x_axis_range(2)-x_axis_range(1))/(REAL(Nx,kind=REAL64))
    dy = (y_axis_range(2)-y_axis_range(1))/(REAL(Ny,kind=REAL64))

    !how long to run simulation for (timesteps)
    total_time = 100

    !how long is each timestep
    dt = 0.01
    
    CALL initialise_file(data_filename,ierr)
    !generate constant gradient potential (left to right)
    CALL generate_const_grad_potential(0.01_REAL64)
    !get electric field
    CALL get_field()
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
END PROGRAM