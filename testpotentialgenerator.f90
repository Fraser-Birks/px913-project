MODULE testpotentialgenerator
    USE data_structures
    IMPLICIT NONE
    SAVE
    CONTAINS
    SUBROUTINE generate_constant_potential(val)
        REAL(REAL64), INTENT(IN) :: val
        ALLOCATE(potential(Nx,Ny))
        !set all non-ghost cells to the value
        potential(2:Nx+1,2:Ny+1) = val
        !set all ghost cells to 0 in this test case (along edges)
        potential(1,1:Ny+2) = 0
        potential(Nx+2,1:Ny+2) = 0
        potential(1:Nx+2,1) = 0
        potential(1:Nx+2,Ny+2) = 0
    END SUBROUTINE

    SUBROUTINE generate_const_grad_potential(val)
        REAL(REAL64), INTENT(IN) :: val
        INTEGER :: loop
        ALLOCATE(potential(Nx+2,Ny+2))
        !set all non-ghost cells to the value such that it increases left to right
        DO loop = 2,Nx+1
            potential(loop,2:Ny+1) = val*(loop-2/Nx-1)
        END DO 
        !set all ghost cells to 0 in this test case (along edges)
        potential(1,1:Ny+2) = 0
        potential(Nx+2,1:Ny+2) = 0
        potential(1:Nx+2,1) = 0
        potential(1:Nx+2,Ny+2) = 0
    END SUBROUTINE
END MODULE

        