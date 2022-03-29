program prog1
   use TypesModule
   use InputModule
   use TestModule
   use InitModule
   use BondGroupModule
   use MinimizingModule
   implicit none

   type(Atom), TARGET, ALLOCATABLE :: Atoms(:)
   Character(len=*), PARAMETER     :: Input_File = 'ch4.xyz'
   type(Bonding)                   :: Bonding_Info
   real                            :: start, finish

   call cpu_time(start)
   ! Read the input file
   call Read_File(Input_File, Atoms)

   ! Initialize the atoms by figuring out which atoms are bonded
   call Initialize_Atoms(Atoms)

   ! Test if all C have 4 bonds, and H has 1
   call Test_Bond(Atoms)

   ! Create groups of 3 and 4 bonded atoms and store them in Bonding info
   call Group_ABCD(Atoms, Bonding_Info)

   ! Perform metropolis minimizing algorithm on the molecule 
   call Metropolis(Atoms, Bonding_Info)

   ! Write optimized atom coordinates to Output.txt
   call Write_File('Output.txt', Atoms)

   call cpu_time(finish)
   print '("Algorithm time = ",f6.3," seconds.")',finish-start
end program