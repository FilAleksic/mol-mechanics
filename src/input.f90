module InputModule
    use TypesModule
    implicit none
    
    PUBLIC
contains
    
    ! Read file and store information in out array
    subroutine Read_File(File, Out_Array)
        character(len=*), INTENT(IN) :: File
        type(Atom), ALLOCATABLE      :: Out_Array(:)
        real*8                       :: x, y, z
        character(len=1)             :: symbol 
        integer                      :: N, i
        open(1, file=File)
        read(1,*) N
        allocate(Out_Array(N))
        do i=1,N
            read(1,*)symbol, x, y, z
            Out_Array(i)%symbol = symbol
            Out_Array(i)%x = x
            Out_Array(i)%y = y
            Out_Array(i)%z = z
        end do
        close(1)
    end subroutine

    subroutine Write_File(File, Atom_Array)
        character(len=*) :: File
        type(Atom)       :: Atom_Array(:)
        integer          :: i

        1 format(i2)
        2 format(a1,999f20.12)
        
        open(1, file=File, status='Replace')
        write(1,1) SIZE(Atom_Array)
        write(1,*) ''
        do i=1,SIZE(Atom_Array)
            write(1,2) Atom_Array(i)%symbol, Atom_Array(i)%x, Atom_Array(i)%y, Atom_Array(i)%z
        end do

    end subroutine
end module