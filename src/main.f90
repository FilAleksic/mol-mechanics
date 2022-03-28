program prog1
   use TypesModule
   use InputModule
   use DistanceModule
   use TestModule
   use InitModule
   use BondGroupModule
   use EnergyModule
   implicit none

   type(Atom), TARGET, ALLOCATABLE :: Atoms(:)
   Character(len=*), PARAMETER     :: Input_File = 'c4h10.xyz'
   type(Bonding) :: Bonding_Info
   real*8 :: Energy_Out

   call Read_File(Input_File, Atoms)
   call Initialize_Atoms(Atoms)

   call Test_Bond(Atoms)
   !call Test_Array_Pointers(Atoms)
   call Group_ABCD(Atoms, Bonding_Info)   

   Energy_Out = Energy(Atoms, Bonding_Info)
   print*, Energy_Out
end program