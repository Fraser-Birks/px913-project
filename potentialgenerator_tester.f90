PROGRAM potentialgenerator_tester
    USE data_structures
    USE chargedensitygenerator
    USE fieldsolver
    USE iso_fortran_env

    IMPLICIT NONE
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: rho

    !set grid dimensions
    Nx = 20
    Ny = 20

    ALLOCATE(rho(Nx,Ny))

    !create axis
    x_axis_range = (/-1.0_REAL64,1.0_REAL64/)
    y_axis_range = (/-1.0_REAL64,1.0_REAL64/)
    !set number of particles

    !find dx,dy
    dx = (x_axis_range(2)-x_axis_range(1))/(REAL(Nx,kind=REAL64))
    dy = (y_axis_range(2)-y_axis_range(1))/(REAL(Ny,kind=REAL64))

    rho = generate_single_charge_density()
    CALL solve_gauss_seidel(rho)

END PROGRAM 