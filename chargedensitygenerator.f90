MODULE chargedensitygenerator
    USE data_structures
    IMPLICIT NONE
    SAVE
    CONTAINS

    ! Function to generate a chanrge density field in an
    ! which is 0 everywhere  
    ! Returns rho, a Nx*Ny array of REAL_64 
    FUNCTION generate_null_charge_density() RESULT(rho)
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: rho
        ! Allocate memory
        ALLOCATE(rho(Nx,Ny))

        rho = 0.0_REAL64
    END FUNCTION generate_null_charge_density

    ! Function to return a charge density field with a single
    ! Gaussian peak at the origin.
    ! Reurns rho, a Nx*Ny array of REAL_64
    FUNCTION generate_single_charge_density() RESULT(rho)
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: rho
        REAl(REAL64), DIMENSION(:), ALLOCATABLE :: x, y
        INTEGER :: i, j !loop integers
        ! Allocate memory for the respective arrays
        ALLOCATE(rho(Nx,Ny))
        ALLOCATE(x(Nx))
        ALLOCATE(y(Ny))

        ! set x-axis to go from -x_range to +x-range
        x(1) = x_axis_range(1) + dx/2.0_REAL64
        DO i = 2, Nx
            x(i) = x(i-1) + dx
        END DO
        
        ! set y-axis to go from -y-range to + y-range
        y(1) = y_axis_range(1) + dy/2.0_REAL64
        DO j = 2, Ny
            y(j) = y(j-1) + dy
        END DO

        ! calculate single gaussian
        DO j = 1, Ny
            DO i = 1, Nx
                rho(i,j) = exp(-((x(i)/0.1_REAL64)**2.0_REAL64) - ((y(j))/0.1_REAL64)**2.0_REAL64) 
            END DO
        END DO
    END FUNCTION generate_single_charge_density

    ! Function to return a charge density field with a double
    ! Gaussian peak at the origin.
    ! Reurns rho, a Nx*Ny array of REAL_64
    FUNCTION generate_double_charge_density() RESULT(rho)
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: rho
        REAl(REAL64), DIMENSION(:), ALLOCATABLE :: x, y
        INTEGER :: i, j !loop integers

        ! Allocate memory for the respective arrays
        ALLOCATE(rho(Nx,Ny))
        ALLOCATE(x(Nx))
        ALLOCATE(y(Ny))

        ! set x-axis to go from -x-range to + x-range
        x(1) = x_axis_range(1) + dx/2.0_REAL64
        DO i = 2, Nx
            x(i) = x(i-1) + dx
        END DO
       
        ! set y-axis to go from -y-range to + y-range
        y(1) = y_axis_range(1) + dy/2.0_REAL64
        DO j = 2, Ny
            y(j) = y(j-1) + dy
        END DO

        ! calculate double gaussian
        DO j = 1, Ny
            DO i = 1, Nx
                rho(i,j) = exp(-((x(i)+0.25)/0.1)**2.0_REAL64-((y(j)+0.25)/0.1)**2.0_REAL64)+ &
                exp(-((x(i)+0.75)/0.2)**2.0_REAL64-(((y(j)-0.75)/0.2)**2.0_REAL64)) 
            END DO
        END DO
    END FUNCTION generate_double_charge_density

END MODULE
