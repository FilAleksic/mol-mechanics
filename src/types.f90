module TypesModule
    implicit none

    PUBLIC :: Atom, ABCD, pnt
    PRIVATE

    type pnt
        type(Atom), POINTER :: p
    end type

    type Atom
        real*8  :: x, y, z
        character(len=1) :: symbol
        type(pnt), ALLOCATABLE :: Bond(:)
    end type

    type ABCD 
        type(pnt), DIMENSION(4) :: Atom
    end type

    
end module