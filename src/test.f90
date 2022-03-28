module TestModule
    use TypesModule
    implicit none
    
    PUBLIC :: Test_Bond
contains

   ! Test if all carbons and hydrogens are bonded
   subroutine Test_Bond(Atom_Array)
        type(Atom), INTENT(IN) :: Atom_Array(:)
        integer                :: i
        do i=1,SIZE(Atom_Array)
            if(Atom_Array(i)%symbol == 'C') then
               if(.not. SIZE(Atom_Array(i)%Bond) == 4) then
                  print*, '!!WARNING!! Carbon with', SIZE(Atom_Array(i)%Bond), 'bonds'
               end if
            else
               if (.not. SIZE(Atom_Array(i)%Bond) == 1) then
                  print*, '!!WARNING!! Hydrogen with', SIZE(Atom_Array(i)%Bond), 'bonds'
               end if
            end if
         end do
   end subroutine

   ! Print the the pointers of the bonds
   subroutine Test_Array_Pointers(Atom_Array)
      type(Atom), INTENT(IN) :: Atom_Array(:)
      integer                :: i, j

      do i=1,SIZE(Atom_Array)
         do j=1,SIZE(Atom_Array(i)%Bond)
            print*, i, Atom_Array(i)%symbol, Atom_Array(i)%Bond(j)%p%symbol, Atom_Array(i)%Bond(j)%p%x
         end do
      end do

   end subroutine
end module