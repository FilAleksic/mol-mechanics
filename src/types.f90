module TypesModule
    implicit none

    PUBLIC :: Atom, Bonding, pnt
    PRIVATE

    type pnt
        type(Atom), POINTER :: p
    end type

    type Atom
        real*8  :: x, y, z
        character(len=1) :: symbol
        type(pnt), ALLOCATABLE :: Bond(:)
    end type

    type Bonding 
        type(pnt), DIMENSION(4) :: Chain_F
        type(pnt), DIMENSION(3) :: Chain_T
        logical                 :: Del = .false.
    end type

    
end module