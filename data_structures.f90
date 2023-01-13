!Use this module as a place to define all of the relevant data structures in each program
!this will make all the structures here accessible inside every function, which means we don't have to 
!do as much passing things around (as fortran is pass-by-reference, not pass-by-value)

MODULE data_structures
    USE iso_fortran_env

    !physical constants
    REAL(REAL64), PARAMETER :: epsilon0 = 1
    REAL(REAL64), PARAMETER :: m = 1
    REAL(REAL64), PARAMETER :: q = -1


    TYPE :: particle
        REAL(REAL64), DIMENSION(2) :: position
        REAL(REAL64), DIMENSION(2) :: velocity
        REAL(REAL64), DIMENSION(2) :: acceleration
        REAL(REAL64), DIMENSION(2) :: prev_position
        REAL(REAL64), DIMENSION(2) :: prev_velocity
        REAL(REAL64), DIMENSION(2) :: prev_acceleration
    END TYPE

    !Define the type to hold the metadata about the program
    !This type must be passed to the function write_metadata in order
    !to write the variables included within it to the file.
    TYPE :: program_metadata
        !Namestring includes the name of the code used to run the file,
        !as well as the version of the code.
        CHARACTER(LEN=30) :: namestring
        
        !Command line args is a single string which includes
        !all of the args that need to be specified on the command line
        !to reproduce the result shown in the data
        CHARACTER(LEN=2000) :: commandlineargs
    END TYPE

    REAL(REAL64), DIMENSION(:,:,:), ALLOCATABLE :: field !the placeholder to store the field once it's generated.
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: potential !this is placeholder array created to store a dummy potential
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: rho !allocate rho (the charge density)
    REAL(REAL64) :: dy,dx,dt !variables which say how much of a step along the x axis and y axis each grid unit is and how large the timestep is.
    REAL(REAL64) :: total_time !variable which says how long the total time in the program is
    INTEGER :: Nx, Ny !How large the grid size will be
    INTEGER :: N_particles !How many particles we're simulating.
    LOGICAL :: particle_out_of_bounds

    REAL(REAL64), DIMENSION(:), ALLOCATABLE :: x_axis, y_axis
    INTEGER :: nghosts
    REAL(REAL64), DIMENSION(2) :: x_axis_range, y_axis_range
END MODULE
