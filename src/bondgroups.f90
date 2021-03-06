module BondGroupModule
    use TypesModule
    implicit none

    PUBLIC :: Group_ABCD
    PRIVATE

contains

    subroutine Group_ABCD(Atom_Array, BondGroups)
        type(Atom), TARGET, INTENT(IN) :: Atom_Array(:)
        type(Bonding)                  :: BondGroups, New_Array
        integer                        :: i
        do i=1,SIZE(Atom_Array)
            call Find_Bonding(Atom_Array(i), BondGroups)
        end do

        if (allocated(BondGroups%Chain_F)) then
            call Remove_Duplicates(BondGroups%Chain_F, New_Array%Chain_F)
        end if
        if (allocated(BondGroups%Chain_T)) then
            call Remove_Duplicates(BondGroups%Chain_T, New_Array%Chain_T)
        end if
        BondGroups = New_Array
    end subroutine

    ! Find every 3 and 4 chain starting from a certain atom by iterating over the bonds, and then the bonds of the bonds
    subroutine Find_Bonding(Atom_Array, BondGroups)
        type(Atom), TARGET, INTENT(IN) :: Atom_Array
        type(Bonding)                  :: BondGroups
        type(Bonding)                  :: temp
        type(Atom), TARGET             :: tempAtom, tempAtom2
        type(Atom), POINTER            :: tempPointer
        integer                        :: i, j, k

        if (.not. allocated(temp%Chain_T)) then
            allocate(temp%Chain_T(1,3))
            allocate(temp%Chain_F(1,4))
        end if

        temp%Chain_T(1,1)%p => Atom_Array
        temp%Chain_F(1,1)%p => Atom_Array

        do i=1,SIZE(Atom_Array%Bond)
            temp%Chain_T(1,2)%p => Atom_Array%Bond(i)%p
            temp%Chain_F(1,2)%p => Atom_Array%Bond(i)%p
            tempAtom = Atom_Array%Bond(i)%p
            tempPointer => Atom_Array
            do j=1,SIZE(tempAtom%Bond)
                ! If pointing 1 and 3 are not pointing at the same atom
                if (.not. associated(tempPointer, tempAtom%Bond(j)%p)) then
                    temp%Chain_T(1,3)%p => tempAtom%Bond(j)%p
                    temp%Chain_F(1,3)%p => tempAtom%Bond(j)%p
                    ! If first iteration use optional variable with size of row to allocate in push back bonding subroutine
                    if (i==1) then
                        call Push_Back_Bonding(BondGroups%chain_T, temp%Chain_T, SIZE(temp%Chain_T(1,:)))
                    else
                        call Push_Back_Bonding(BondGroups%chain_T, temp%Chain_T)
                    end if
                    tempAtom2 = tempAtom%Bond(j)%p

                    do k=1,SIZE(tempAtom2%Bond)
                        ! If pointer 2 and 4 are not pointing at the same atom
                        if (.not. associated(Atom_Array%Bond(i)%p, tempAtom2%Bond(k)%p)) then
                            temp%Chain_F(1,4)%p => tempAtom2%Bond(k)%p
                            ! If first iteration use optional variable with size of row to allocate in push back bonding subroutine
                            if (i==1) then
                                call Push_Back_Bonding(BondGroups%chain_F, temp%Chain_F, SIZE(temp%Chain_F(1,:)))
                            else
                                call Push_Back_Bonding(BondGroups%chain_F, temp%Chain_F)
                            end if
                        end if
                    end do

                end if
            end do
        end do 
    end subroutine  

    ! Remove duplicate bondings by checking with the Check_Association function
    subroutine Remove_Duplicates(PointArray, NewArray)
        type(pnt), ALLOCATABLE :: NewArray(:,:), Tester(:), temp(:,:)
        type(pnt)              :: PointArray(:,:)
        integer                :: i, j, count=0, szx, szy
        
        szx = SIZE(PointArray(1,:))
        szy = SIZE(PointArray(:,1))
        allocate(Tester(szx))
        allocate(temp(1,szx))

        do i=1,szy
            if (.not. PointArray(i,1)%Del) then
                temp(1,1:szx) = PointArray(i,1:szx)
                ! If first iteration use optional variable with size of row to allocate in push back bonding subroutine
                if (i==1) then
                    call Push_Back_Bonding(NewArray, temp, szx)
                else 
                    call Push_Back_Bonding(NewArray, temp)
                end if
                do j=i+1,szy
                    Tester(:) = PointArray(j,:)
                    count = Check_Association(PointArray(i,:), Tester)
                    ! If every pointer was matching 1234 = 4321 meaning same link of 4 atoms and can be discarded
                    ! Same for chain of 3 atoms, 123 = 321 so count = 3 when matching
                    if (count == szx) then
                        PointArray(j,:)%Del = .True.
                    end if
                    count = 0
                end do
            end if
        end do
    end subroutine

    ! Appending Algorithm
    subroutine Push_Back_Bonding(BondGroups, NewValue, sz)
        type(pnt), ALLOCATABLE :: BondGroups(:,:), temp(:,:)
        type(pnt)              :: NewValue(:,:)
        integer, OPTIONAL      :: sz
        integer                :: szx, szy

        if (allocated(BondGroups)) then
            szx = SIZE(BondGroups(1,:))
            szy = SIZE(BondGroups(:,1))
            allocate(temp(szy, szx))
            temp(:,:) = BondGroups(:,:)
            deallocate(BondGroups)
            allocate(BondGroups(szy+1, szx))
            BondGroups(1:szy,:) = temp(:,:)
            BondGroups(szy+1,:) = NewValue(1,:)
        else
            allocate(BondGroups(1,sz))
            BondGroups(1,:) = NewValue(1,:)
        end if
    end subroutine

    ! Add 1 to Check_Association if matching pointer
    integer function Check_Association(a, b)
        type(pnt) :: a(:), b(:)
        integer   :: szx, i, j
        szx = SIZE(a)
        j = szx
        Check_Association = 0
        do i=1,szx
            if (associated(a(i)%p, b(j)%p)) then
                Check_Association = Check_Association + 1
            end if
            j = j - 1
        end do

    end function
    
end module