MODULE fieldsolver
    USE data_structures
    IMPLICIT NONE

    CONTAINS

    SUBROUTINE solve_gauss_seidel(charge_density)
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: charge_density 
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: rho, init_potential
        REAL(REAL64), DIMENSION(:), ALLOCATABLE :: err, num_dev
        REAL(REAL64) :: tot_err, drms, norm_err, tol
        INTEGER :: i, j !loop variables
       
        rho = charge_density
        ALLOCATE(init_potential(0:Nx+1, 0:Ny+1))
        ALLOCATE(potential(0:Nx+1, 0:Ny+1))
        ALLOCATE(err(SIZE(potential)))
        ALLOCATE(num_dev(SIZE(potential)))
        init_potential = 0.5_REAL64 ! use 0 as initial guess 
        tol = 0.00001_REAL64
        norm_err = 1.0_REAL64
        PRINT*, norm_err
    
        DO WHILE (norm_err >= tol)
            DO j = 1, Ny
                DO i = 1, Nx
                    potential(i,j) = -(rho(i,j)-((init_potential(i+1,j)+init_potential(i-1,j))/dx**2)- &
                    ((init_potential(i,j+1)+init_potential(i,j-1))/dy**2))/((2.0_REAL64/dx**2)+(2.0_REAL64/dy**2))
                    init_potential(i,j) = potential(i,j)
                    err = abs(((potential(i-1,j)-2.0_REAL64*potential(i,j)+potential(i+1,j))/dx**2.0_REAL64)+ &
                    ((potential(i,j-1)-2.0_REAL64*potential(i,j)+potential(i,j+1)/dy**2.0_REAL64))-rho(i,j))
                    num_dev = abs((potential(i-1,j)-2.0_REAL64*potential(i,j)+potential(i+1,j))/dx**2.0_REAL64)+ &
                    ((potential(i,j-1)-2.0_REAL64*potential(i,j)+potential(i,j+1)/dy**2.0_REAL64))
                END DO
            END DO
            tot_err = sum(err)
            drms = sqrt((1.0_REAL64/SIZE(potential))*sum(num_dev))
            norm_err = tot_err/drms
            init_potential = potential
            PRINT*, norm_err
        END DO

    END SUBROUTINE
END MODULE