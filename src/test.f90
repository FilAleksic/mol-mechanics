module TestModule
    use TypesModule
    implicit none
    
    PUBLIC :: Test_Bond
contains

   ! Test if all carbons and hydrogens are bonded
   subroutine Test_Bond(Atom_Array)
        type(Atom), INTENT(IN) :: Atom_Array(:)
        integer                :: i, j
        do i=1,SIZE(Atom_Array)
            if(Atom_Array(i)%symbol == 'C') then
               do j=1,4
                  if (.not. associated(Atom_Array(i)%CBond(j)%p)) then
                     print*, '!!WARNING!! Carbon With < 4 Bonds!!'
                  end if
               end do
            else
               if (.not. associated(Atom_Array(i)%HBond%p)) then
                  print*, '!!WARNING!! Non-Bonded Hydrogen!!'
               end if
            end if
         end do
   end subroutine

   ! Print the the pointers of the bonds
   subroutine Test_Array_Pointers(Atom_Array)
      type(Atom), INTENT(IN) :: Atom_Array(:)
      integer                :: i, j

      do i=1,SIZE(Atom_Array)
         if (Atom_Array(i)%symbol == 'C') then
            do j=1,4
               print*, i, Atom_Array(i)%symbol, Atom_Array(i)%CBond(j)%p%symbol, Atom_Array(i)%CBond(j)%p%x &
                     , Atom_Array(i)%CBond(j)%p%y, Atom_Array(i)%CBond(j)%p%z
            end do
         else
            print*, i, Atom_Array(i)%symbol, Atom_Array(i)%HBond%p%symbol, Atom_Array(i)%HBond%p%x &
                  , Atom_Array(i)%HBond%p%y, Atom_Array(i)%HBond%p%z
         end if
      end do

   end subroutine
    
end module