program prog1
   use TypesModule
   use InputModule
   use DistanceModule
   implicit none

   type(Atom), ALLOCATABLE :: Atoms(:)
   integer                 :: i, j
   call Read_File('c4h10.xyz', Atoms)

   do i=1,SIZE(Atoms)-1
      do j=i+1,SIZE(Atoms)
         if(Rel_Distance(Atoms(i), Atoms(j)) < 2.5) then
            print*, Atoms(i)%c, Atoms(j)%C, Rel_Distance(Atoms(i), Atoms(j))
         end if 
      end do
   end do

end program