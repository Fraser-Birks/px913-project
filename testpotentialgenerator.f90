MODULE testpotentialgenerator
    USE data_structures
    IMPLICIT NONE
    SAVE
    CONTAINS
    SUBROUTINE generate_constant_potential(val)
        REAL(REAL64), INTENT(IN) :: val
        ALLOCATE(potential(0:Nx+1,0:Ny+1))
        !set all non-ghost cells to the value
        potential(1:Nx,1:Ny) = val
        !set all ghost cells to 0 in this test case (along edges)
        potential(0,1:Ny+1) = 0.0_REAL64
        potential(Nx+1,0:Ny+1) = 0.0_REAL64
        potential(0:Nx+1,0) = 0.0_REAL64
        potential(0:Nx+1,Ny+1) = 0.0_REAL64
    END SUBROUTINE

    SUBROUTINE generate_const_grad_potential(val)
        REAL(REAL64), INTENT(IN) :: val
        INTEGER :: loop
        ALLOCATE(potential(0:Nx+1,0:Ny+1))
        !set all non-ghost cells to the value such that it increases left to right
        DO loop = 1,Nx
            potential(loop,1:Ny) = val*((REAL(loop,kind=REAL64)-1)/(REAL(Nx,kind=real64)))
        END DO 
        !set all ghost cells to 0 in this test case (along edges)
        potential(0,1:Ny+1) = 0.0_REAL64
        potential(Nx+1,0:Ny+1) = 0.0_REAL64
        potential(0:Nx+1,0) = 0.0_REAL64
        potential(0:Nx+1,Ny+1) = 0.0_REAL64
    END SUBROUTINE
END MODULE

        