module BondGroupModule
    use TypesModule
    implicit none
    
contains

    ! Return groups of 4 atoms 
    function Group_ABCD(Atom_Array) result(Out)
        type(Atom), INTENT(IN)  :: Atom_Array
        type(ABCD), ALLOCATABLE :: Out(:), Temp(:)
        
        allocate(Out(1))

    end function

    
end module