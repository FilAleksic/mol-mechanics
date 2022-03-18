program prog1
   use TypesModule
   ! use InputModule
   ! use DistanceModule
   use TestModule
   use InitModule
   use InputModule
   implicit none

   type(Atom), TARGET, ALLOCATABLE :: Atoms(:)
   Character(len=*), PARAMETER     :: Input_File = 'c4h10.xyz'

   call Read_File(Input_File, Atoms)

   call Initialize_Atoms(Atoms)
   
   call Test_Bond(Atoms)

   call Test_Array_Pointers(Atoms)
   

end program