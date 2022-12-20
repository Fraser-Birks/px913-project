PROGRAM MAIN
    USE data_structures
    USE testpotentialgenerator
    USE particlevelocitysolver
    IMPLICIT NONE
    TYPE(particle) :: part
    INTEGER :: timesteps_elapsed
    !set grid dimensions
    Nx = 100
    Ny = 100
    
    !set number of particles
    N_particles = 1

    !find dx,dy assuming grid coords range between -1 and 1
    dx = 2.0_REAL64/(REAL(Nx,kind=REAL64))
    dy = 2.0_REAL64/(REAL(Ny,kind=REAL64))

    !how long to run simulation for (timesteps)
    total_time = 100

    !how long is each timestep
    dt = 0.01

    !generate constant gradient potential (left to right)
    CALL generate_const_grad_potential(0.01_REAL64)

    !get electric field
    CALL get_field()

    !initialise particle(s) (takes array of particles or single particle)
    CALL initialise_particles(part)

    timesteps_elapsed = 0

    !Speeding up program by not putting if statement into loop.
    DO WHILE((.NOT.(particle_out_of_bounds)).AND.(timesteps_elapsed<=total_time))
        CALL propagate(part)
        !writing to netCDF file goes here
        timesteps_elapsed = timesteps_elapsed + 1
        PRINT*,timesteps_elapsed
    END DO

    !print final velocity,position and acceleration
    PRINT*,'                   x component                 y component'
    PRINT*,'final velocity',part%prev_velocity
    PRINT*,'final position',part%prev_position
    PRINT*,'final acceleration',part%prev_acceleration
END PROGRAM