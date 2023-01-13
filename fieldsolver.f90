MODULE fieldsolver
    USE data_structures
    IMPLICIT NONE

    CONTAINS

    SUBROUTINE solve_gauss_seidel(charge_density)
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: charge_density 
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: rho, init_potential, err, num_dev
        REAL(REAL64) :: tot_err, drms, norm_err, tol
        INTEGER :: i, j !loop variables
        
        ! intialise rho with the initial charge density taken from the chargedensitygenerator
        rho = charge_density
        ! allocate memory to the various arrays
        ! the 2D arrays range from 0 to Nx+1 and 0 to Ny+1
        ! the 1D arrays range to the total number of values within th potential array
        ALLOCATE(init_potential(0:Nx+1, 0:Ny+1))
        ALLOCATE(potential(0:Nx+1, 0:Ny+1))
        ALLOCATE(err(1:Nx, 1:Ny))
        ALLOCATE(num_dev(1:Nx, 1:Ny))

        ! initialise the initial potential, tolerance and noralised error 
        init_potential = 0.5_REAL64
        tol = 0.00001_REAL64
        norm_err = 1.0_REAL64
    
        ! Convergence loop with convergence critetium set by the value of the normalised error being 
        ! below the tolerance
        ! First loops through the i and j to fill the potential array, which is being updated after each iteration
        ! The new potential array is then used to determine the error and the numberical derivative via a seperate loop 
        ! through i and j
        ! These are used to determine the norm_error
        DO WHILE (norm_err >= tol)
            DO j = 1, Ny
                DO i = 1, Nx
                    potential(i,j) = -(rho(i,j)-((init_potential(i+1,j)+init_potential(i-1,j))/dx**2)- &
                    ((init_potential(i,j+1)+init_potential(i,j-1))/dy**2))/((2.0_REAL64/dx**2)+(2.0_REAL64/dy**2))
                    init_potential(i,j) = potential(i,j)
                END DO
            END DO
            DO j = 1,Ny
                DO i = 1,Nx
                    err(i,j) = abs(((potential(i-1,j)-2.0_REAL64*potential(i,j)+potential(i+1,j))/dx**2)+ &
                    (((potential(i,j-1)-2.0_REAL64*potential(i,j)+potential(i,j+1))/dy**2))-rho(i,j))
                    num_dev(i,j) = (((potential(i-1,j)-2.0_REAL64*potential(i,j)+potential(i+1,j))/dx**2)+ &
                    ((potential(i,j-1)-2.0_REAL64*potential(i,j)+potential(i,j+1))/dy**2))**2
                END DO
            END DO
            tot_err = sum(err)
            drms = sqrt((1.0_REAL64/SIZE(err))*sum(num_dev))
            IF (drms == 0.0_REAL64) THEN
                error stop
            END IF
            norm_err = tot_err/drms
            init_potential = potential
            PRINT*, norm_err
        END DO

    END SUBROUTINE
END MODULE