module InputModule
    use TypesModule
    implicit none
    
contains
    
    subroutine Read_File(File, Out_Array)
        character(len=*), INTENT(IN) :: File
        type(Atom), ALLOCATABLE      :: Out_Array(:)
        real*8                       :: x, y, z
        character(len=1)             :: Symbol 
        integer                      :: N, i
        open(1, file=File)
        read(1,*) N
        allocate(Out_Array(N))
        do i=1,N
            read(1,*)Symbol, x, y, z
            if(Symbol == 'C') then
                Out_Array(i)%C = .true.
            else
                Out_Array(i)%C = .false.
            end if
            Out_Array(i)%x = x
            Out_Array(i)%y = y
            Out_Array(i)%z = z
        end do
    end subroutine
end module