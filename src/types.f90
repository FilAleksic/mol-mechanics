module TypesModule
    implicit none

    PUBLIC :: Atom, Bonding, pnt
    PRIVATE

    type pnt
        type(Atom), POINTER :: p
        logical :: Del = .false.
    end type

    type Atom
        real*8  :: x, y, z
        character(len=1) :: symbol
        type(pnt), ALLOCATABLE :: Bond(:)
    end type

    type Bonding
        type(pnt), ALLOCATABLE :: Chain_F(:,:)
        type(pnt), ALLOCATABLE :: Chain_T(:,:)
    end type

    
end module