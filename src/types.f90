module TypesModule
    implicit none

    PUBLIC :: Atom
    PRIVATE
    
    type pnt
        type(Atom), POINTER :: p
    end type

    type Atom
        real*8  :: x, y, z
        character(len=1) :: symbol
        type(pnt), DIMENSION(4) :: CBond
        type(pnt) :: HBond
    end type

    
end module