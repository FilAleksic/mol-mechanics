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
                        print*, i,j, Atom_Array(i)%symbol, Atom_Array(i)%x
                        Atom_Array(i) = Set_Bond_Pointer(Atom_Array(i), Atom_Array(j))
                        print*, 'Starting J, i'
                        Atom_Array(j) = Set_Bond_Pointer(Atom_Array(j), Atom_Array(i))
                    end if 
            end do
        end do
    end subroutine

    type(Atom) function Set_Bond_Pointer(Atom1, Atom2)
        type(Atom), TARGET :: Atom1, Atom2
        integer            :: k
        type(Atom)         :: temp

        !! NEW CODE THATS CAUSING ISSUES COMPILES BUT DOESNT RUN
        ! Use allocatable array of pointers in type, append a new pointer to it when needed
        ! Allocate it to 1 if it is not yet allocated
        print*, Atom1%symbol, Atom1%x
        print*, allocated(Atom1%Bond)
        if(.not. allocated(Atom1%Bond)) then
            allocate(Atom1%Bond(1))
        else
            temp = Atom1
            print*, 'Allocated'
            deallocate(Atom1%Bond)
            print*, 'Deallocated'
            allocate(Atom1%Bond(SIZE(temp%Bond)+1))
            print*, 'Allocated +1', SIZE(temp%Bond)+1
            Atom1%Bond = temp%Bond
            print*, 'Set to temp'
        end if
        Atom1%Bond(SIZE(Atom1%Bond))%p => Atom2
        print*, 'Set Pointer'
        !! END OF NEW CODE THAT CAUSES ISSUES
        

        !! OLD CODE THAT COMPILES AND RUNS 
        ! If pointer is not yet pointing to something, assign that pointer to j
        ! For Carbons terminate the loop over the array of pointers once you fill the first pointer
        ! if(Atom1%symbol == 'C') then
        !     do k=1,4
        !         if (.not. associated(Atom1%CBond(k)%p)) then
        !             Atom1%CBond(k)%p => Atom2
        !             exit
        !         else
        !         end if
        !     end do
        
        ! else if (Atom1%symbol == 'H') then
        !     if (.not. associated(Atom1%HBond%p)) then
        !         Atom1%HBond%p => Atom2
        !     end if
        ! end if
        !! END OF OLD CODE THAT COMPILES AND RUNS

        Set_Bond_Pointer = Atom1
        print*, 'Return value set'
    end function

end module