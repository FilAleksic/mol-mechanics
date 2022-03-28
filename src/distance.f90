module DistanceModule
    use TypesModule
    implicit none

    PUBLIC :: Rel_Distance, Abs_Distance
    
contains
    ! Calculate relative distance between two atoms.
    ! Relative distance because we don't take the square root
    ! Used for optimization where Abs distance not needed
    real*8 function Rel_Distance(Atom1, Atom2)
        type(Atom), INTENT(IN) :: Atom1, Atom2
        Rel_Distance = (Atom2%x - Atom1%x)**2 + &
                       (Atom2%y - Atom1%y)**2 + &
                       (Atom2%z - Atom1%z)**2
    end function
    
    ! Calculate the absolute distance
    real*8 function Abs_Distance(Atom1, Atom2)
        type(Atom), INTENT(IN) :: Atom1, Atom2
        Abs_Distance = SQRT((Atom2%x - Atom1%x)**2 + &
                       (Atom2%y - Atom1%y)**2 + &
                       (Atom2%z - Atom1%z)**2)
    end function

end module