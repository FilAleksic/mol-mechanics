module InitModule
    use TypesModule
    use DistanceModule
    implicit none
    
    PUBLIC :: Initialize_Atoms
    PRIVATE

contains

    ! Initialize the atom array by reading the file and assigning the pointers of bonds
    subroutine Initialize_Atoms(Atom_Array)
        type(Atom), TARGET :: Atom_Array(:)
        integer                         :: i, j
        real                            :: Bond_Tresh = 2.5

        do i=1,SIZE(Atom_Array)
            do j=i+1,SIZE(Atom_Array)
                    ! If < bond threshold then we have a bond (Note Rel Distance: Angstrom^2)
                    if(Rel_Distance(Atom_Array(i), Atom_Array(j)) < Bond_Tresh) then
                        Atom_Array(i) = Set_Bond_Pointer(Atom_Array(i), Atom_Array(j))
                        Atom_Array(j) = Set_Bond_Pointer(Atom_Array(j), Atom_Array(i))
                    end if 
            end do
        end do
    end subroutine

    type(Atom) function Set_Bond_Pointer(Atom1, Atom2)
        type(Atom), TARGET :: Atom1, Atom2
        integer :: k
        ! If pointer is not yet pointing to something, assign that pointer to j
        ! For Carbons terminate the loop over the array of pointers once you fill the first pointer
        if(Atom1%symbol == 'C') then
            do k=1,4
                if (.not. associated(Atom1%CBond(k)%p)) then
                    Atom1%CBond(k)%p => Atom2
                    exit
                else
                end if
            end do
        
        else if (Atom1%symbol == 'H') then
            if (.not. associated(Atom1%HBond%p)) then
                Atom1%HBond%p => Atom2
            end if
        end if

        Set_Bond_Pointer = Atom1
    end function

end module