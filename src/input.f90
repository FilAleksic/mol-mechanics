module InputModule
    use TypesModule
    implicit none
    
contains
    
    ! Read file and store information in out array
    subroutine Read_File(File, Out_Array)
        character(len=*), INTENT(IN) :: File
        type(Atom), ALLOCATABLE      :: Out_Array(:)
        real*8                       :: x, y, z
        character(len=1)             :: symbol 
        integer                      :: N, i, j
        open(1, file=File)
        read(1,*) N
        allocate(Out_Array(N))
        do i=1,N
            read(1,*)symbol, x, y, z
            Out_Array(i)%symbol = symbol
            Out_Array(i)%x = x
            Out_Array(i)%y = y
            Out_Array(i)%z = z

            ! Important, set pointers to null, otherwise you get pointers pointing to weird artifacts
            Out_Array(i)%HBond%p => NULL()
            do j=1,4
               Out_Array(i)%CBond(j)%p => NULL()
            end do
        end do
    end subroutine
end module