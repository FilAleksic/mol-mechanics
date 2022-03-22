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

        do i=1,SIZE(Atom_Array)-1
            do j=i+1,SIZE(Atom_Array)
                    ! If < bond threshold then we have a bond (Note Rel Distance: Angstrom^2)
                    if(Rel_Distance(Atom_Array(i), Atom_Array(j)) < Bond_Tresh) then
                        call Set_Bond_Pointer(Atom_Array(i), Atom_Array(j))
                        call Set_Bond_Pointer(Atom_Array(j), Atom_Array(i))
                    end if 
            end do
        end do
    end subroutine

    subroutine Set_Bond_Pointer(Atom1, Atom2)
        type(Atom), TARGET :: Atom1, Atom2
        call Push_Back_Atom(Atom1, Atom2)
    end subroutine

    subroutine Push_Back_Atom(Atom1, Bonded)
        type(pnt), ALLOCATABLE :: temp(:)
        type(Atom), TARGET     :: Bonded
        type(Atom), TARGET     :: Atom1
        integer                :: sz
        
        if (allocated(Atom1%Bond)) then
            sz = size(Atom1%Bond)
            allocate(temp(sz))
            temp = Atom1%Bond
            deallocate(Atom1%Bond)
            allocate(Atom1%Bond(sz+1))
            Atom1%Bond(1:sz) = temp
            Atom1%Bond(sz+1)%p => Bonded
        else
            allocate(Atom1%Bond(1))
            Atom1%Bond(1)%p => Bonded
        end if
    end subroutine

end module