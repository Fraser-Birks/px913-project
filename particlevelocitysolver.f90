MODULE particlevelocitysolver
    USE data_structures
    IMPLICIT NONE
    CONTAINS

    SUBROUTINE get_field()
        !Subroutine which implements a finite difference solver to get the field across the grid
        !from the defined potential.
        !Note that the variable we are returning for the field here has two slices, the first slice being Ex and the second Ey.
        !just taking a numerical derivative of potential, no need to get fancy with it
        INTEGER :: i, j !loop variables

        ALLOCATE(field(Nx,Ny,2)) !field to return
        DO j = 1,Ny !loop over 2nd - 2nd to last positions in potential grid in both dimension (ignoring ghost cells)
            DO i = 1,Nx
                !get Ex
                field(i,j,1) = (potential(i+1,j)-potential(i-1,j))/(2*dx) !assign electric field to the 1st - Nth positions in the field arrays
                !get Ey
                field(i,j,2) = (potential(i,j+1)-potential(i,j-1))/(2*dy)
            END DO
        END DO
    END SUBROUTINE
    
    IMPURE ELEMENTAL SUBROUTINE initialise_particles(part)
        !Function to get the initial accelerations for the particles
        !Takes following input:
        !part:either an array of particle objects or single particle
        TYPE(particle), INTENT(INOUT) :: part
        INTEGER, DIMENSION(2) :: cell_coords
        REAL(REAL64) :: Ex,Ey
        
        !get initial cell coords for particle
        cell_coords(1) = FLOOR((part%position(1)+1.0_REAL64)/dx)+1
        cell_coords(2) = FLOOR((part%position(2)+1.0_REAL64)/dy)+1
        
        !get the x and y field corresponding to cell_coords
        Ex = field(cell_coords(1),cell_coords(2),1)
        Ey = field(cell_coords(1),cell_coords(2),2)

        !set initial acceleration for particle based on field
        part%acceleration = update_acceleration(Ex,Ey)

    END SUBROUTINE

    FUNCTION update_acceleration(field_valX, field_valY)
        !Function that calculates the acceleration of a particle
        !based on it's position in the electric field
        !Returns the x and y components of the acceleration of a particle.
        REAL(REAL64), INTENT(IN) :: field_valX, field_valY
        REAL(REAL64), DIMENSION(2) :: update_acceleration
        !update acceleration
        update_acceleration(1) = q*(field_valX)/m
        update_acceleration(2) = q*(field_valY)/m
    END FUNCTION


    IMPURE ELEMENTAL SUBROUTINE propagate(part)
        !Subroutine which takes in a given particle (or array of particles)
        !and propagates it (them) one timestep. Using the velocity verlet algorithm. 
        !To read more about the velocity verlet algorithm, see: https://en.wikipedia.org/wiki/Verlet_integration
        TYPE(particle), INTENT(INOUT) :: part
        INTEGER, DIMENSION(2) :: cell_coords
        REAL(REAL64) :: Ex, Ey
        !first, set all 'previous' particle parameters to be the same as the 'current' ones
        !from the previous timestep.
        part%prev_position = part%position
        part%prev_velocity = part%velocity
        part%prev_acceleration = part%acceleration

        !find the new position using velocity verlet (prev_position is nth, new is n+1th)
        part%position = part%prev_position + (part%prev_velocity*dt) + (0.5_REAL64)*(part%prev_acceleration)*(dt**2)

        !find the cell the particle now sits in in the grid and assign the x and y of the cell to cell_coords
        cell_coords(1) = FLOOR((part%position(1)+1.0_REAL64)/dx)+1
        cell_coords(2) = FLOOR((part%position(2)+1.0_REAL64)/dy)+1
            
        !get the x and y field corresponding to cell_coords
        Ex = field(cell_coords(1),cell_coords(2),1)
        Ey = field(cell_coords(1),cell_coords(2),2)

        !get the new accelaration of the particle
        part%acceleration = update_acceleration(Ex,Ey)

        !get the new velocity of the particle
        part%velocity = part%prev_velocity + dt*((part%acceleration + part%prev_acceleration)*0.5_REAL64)

        !Run a check to see if any of the particles are now out of bounds.
        !If they are, then break out of the loop in main_program
        IF (.NOT.(((cell_coords(1) < Nx+1).AND.(cell_coords(1)>0).AND.(cell_coords(2)>0).AND.(cell_coords(2)<Ny+1)))) THEN
            particle_out_of_bounds=.TRUE.
        END IF
    END SUBROUTINE 

END MODULE