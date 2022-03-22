module BondGroupModule
    use TypesModule
    implicit none

contains

    ! Return groups of 3 and 4 atoms
    subroutine Group_ABCD(Atom_Array, Bonding_Array)
        type(Atom), TARGET, INTENT(IN)  :: Atom_Array(:)
        type(Bonding), ALLOCATABLE :: Bonding_Array(:), New_Array(:)
        integer                    :: i
        do i=1,SIZE(Atom_Array)
            call Find_Bonding(Atom_Array(i), Bonding_Array)
        end do
        call Remove_Duplicate_Pnt(Bonding_Array, New_Array)
        print*, SIZE(New_Array), 'FINAL'
    end subroutine

    subroutine Find_Bonding(Atom_Array, Bonding_Array)
        type(Atom), TARGET, INTENT(IN) :: Atom_Array
        type(Bonding), ALLOCATABLE :: Bonding_Array(:)
        type(Bonding) :: temp
        type(Atom), TARGET :: tempAtom, tempAtom2
        integer :: i, j, k
        type(Atom), POINTER :: x

        temp%Chain_T(1)%p => Atom_Array
        temp%Chain_F(1)%p => Atom_Array

        do i=1,SIZE(Atom_Array%Bond)
            temp%Chain_T(2)%p => Atom_Array%Bond(i)%p
            temp%Chain_F(2)%p => Atom_Array%Bond(i)%p
            tempAtom = Atom_Array%Bond(i)%p
            x => Atom_Array
            do j=1,SIZE(tempAtom%Bond)
                ! If pointing 1 and 3 are not pointing at the same atom
                if (.not. associated(x, tempAtom%Bond(j)%p)) then
                    temp%Chain_T(3)%p => tempAtom%Bond(j)%p
                    temp%Chain_F(3)%p => tempAtom%Bond(j)%p

                    tempAtom2 = tempAtom%Bond(j)%p

                    do k=1,SIZE(tempAtom2%Bond)
                        ! If pointer 2 and 4 are not pointing at the same atom
                        if (.not. associated(Atom_Array%Bond(i)%p, tempAtom2%Bond(k)%p)) then
                            temp%Chain_F(4)%p => tempAtom2%Bond(k)%p
                            call Push_Back_Bonding(Bonding_Array, temp)
                        end if
                    end do

                end if
            end do
        end do
    end subroutine

    ! Appending algorithm
    subroutine Push_Back_Bonding(Bonding_Array, NewValue)
        type(Bonding), ALLOCATABLE :: Bonding_Array(:), temp_Array(:)
        type(Bonding) :: NewValue
        integer :: sz

        if (allocated(Bonding_Array)) then
            sz = SIZE(Bonding_Array)
            allocate(temp_Array(sz))
            temp_Array = Bonding_Array
            deallocate(Bonding_Array)
            allocate(Bonding_Array(sz+1))
            Bonding_Array(1:sz) = temp_Array
            Bonding_Array(sz+1)%Chain_F = NewValue%Chain_F
        else
            allocate(Bonding_Array(1))
            Bonding_Array(1)%Chain_F = NewValue%Chain_F
        end if
    end subroutine

    subroutine Remove_Duplicate_Pnt(Bonding_Array, New_Array)
        type(Bonding) :: Bonding_Array(:)
        type(Bonding), ALLOCATABLE :: New_Array(:)
        type(pnt), ALLOCATABLE :: Tester(:)
        integer, ALLOCATABLE :: perm(:,:)
        integer :: i, j, k, h, m, count = 0

        ! Generate all permutations for an array of numbers between 1 and 4
        ! (hardcoded as this removes duplicates of Chain_F)
        call Permutations(perm, 1,4)
        allocate(Tester(SIZE(perm(1,:))))
        do i=1,SIZE(Bonding_Array)-1
            if (.not. Bonding_Array(i)%Del) then
                ! Append i to the new array
                call Push_Back_Bonding(New_Array, Bonding_Array(i))
                do j=i+1,SIZE(Bonding_Array)
                    ! Y axis of permutations
                    do k=1,SIZE(perm,1)
                        ! X axis of permutations
                        do h=1,SIZE(perm(k,:))
                            Tester(h) = Bonding_Array(j)%Chain_F(perm(k,h))
                        end do
                        ! Iterate over the pointers and check if they are the same
                        do m=1,SIZE(Tester)
                            if (associated(Bonding_Array(i)%Chain_F(m)%p, Tester(m)%p)) then
                                count = count + 1
                            end if
                        end do
                        ! If all the pointers are the same then we are looking at the same grouping of 4 bonds
                        ! Mark j to with Del so it doesn't get iterated, reset count and go to the next iteration of i
                        if (count == SIZE(Tester)) then
                            Bonding_Array(j)%Del = .True.
                        end if
                        count = 0
                    end do
                end do
            end if
        end do
    end subroutine

    ! Permutation algorithm written by Jos Bergevoe
    ! http://computer-programming-forum.com/49-fortran/52bec66e9c6a7145.htm
    subroutine Permutations(Out, minV, maxV)
        integer, INTENT(IN) :: minV, maxV
        integer, ALLOCATABLE :: out(:,:)
        integer :: array(maxV-(minV-1)), i

        do i=minV,maxV-(minV-1)
            if (i==1) then
                array(i) = minV
            else
                array(i)=minv+(i-1)
            end if
        end do

        allocate(out(product(array), size(array)))
        call permutate(array, out)
    end subroutine

    recursive subroutine permutate(E, P)
        integer, intent(in)  :: E(:)       ! array of objects
        integer, intent(out) :: P(:,:)     ! permutations of E
        integer  :: N, Nfac, i, k, S(size(P,1)/size(E), size(E)-1)
        N = size(E); Nfac = size(P,1);
        do i=1,N                           ! cases with E(i) in front
        if( N>1 ) call permutate((/E(:i-1), E(i+1:)/), S)
        forall(k=1:Nfac/N) P((i-1)*Nfac/N+k,:) = (/E(i), S(k,:)/)
        end do
    end subroutine
    
end module