module MinimizingModule
    use TypesModule
    use EnergyModule
    use TestModule
    implicit none
    
    PUBLIC :: Metropolis
    PRIVATE

contains

    ! Metropolis function
    ! Minimum energy is stored in Energy_Min, and coordinates in Atoms_Min
    ! If no value lower than the minimum is found for 10000 iterations, assume convergence has completed
    subroutine Metropolis(Atom_Array, Bonding_Array)
        type(Atom)                :: Atom_Array(:)
        type(Bonding), INTENT(IN) :: Bonding_Array
        type(Atom), ALLOCATABLE   :: Atoms_New(:), Atoms_Min(:)
        real*8                    :: Energy_Out, Energy_OutNew, EnergyMin=1000
        integer                   :: lastChange=0, i=0

        Energy_Out = Energy(Atom_Array, Bonding_Array)
        print*, Energy_Out, 'Starting Energy'
        do while (lastChange < 5000)
            i = i + 1
            call NewPositions(Atom_Array, Atoms_New)
            call Test_Bond(Atoms_New)
            Energy_OutNew = Energy(Atoms_New, Bonding_Array)
            if (Energy_OutNew - Energy_Out < 0) then 
                if (Energy_OutNew < EnergyMin) then 
                    EnergyMin = Energy_OutNew
                    Atoms_Min = Atoms_New
                    lastChange = 0
                end if
                Atom_Array = Atoms_New
                Energy_Out = Energy_OutNew
            else if (Accept(Energy_OutNew - Energy_Out)) then
                Atom_Array = Atoms_New
                Energy_Out = Energy_OutNew
            end if
            lastChange = lastChange + 1
        end do
        print*, i, 'Iterations'
        print*, EnergyMin, 'Final Energy'
        Atom_Array = Atoms_Min
    end subroutine

    ! Pick a random atom, move it by a random amount
    subroutine NewPositions(Old, New)
        type(Atom)              :: Old(:)
        type(Atom), ALLOCATABLE :: New(:)
        integer                 :: sz, j
        real*8                  ::  rng
        type(Constant)          :: c
        sz = SIZE(Old)
        if (allocated(New)) deallocate(New)
        allocate(New(sz))
        
        call random_number(rng)
        j = 1 + FLOOR(rng * sz)
        New = Old
        New(j)%x = New(j)%x + (c%R * random())
        New(j)%y = New(j)%y + (c%R * random())
        New(j)%z = New(j)%z + (c%R * random())
    end subroutine

    ! Accept values when dEnergy is <0 or with probability based on temperature and Boltzman constant
    logical function Accept(dEnergy)
        real*8, INTENT(IN) :: dEnergy
        real*8             :: q, pa, dEnJ
        type(Constant)     :: c
        if (dEnergy < 0) then
            Accept = .true.
        else
            dEnJ = dEnergy * 4184
            !print*, dEnJ
            pa = exp(-dEnJ /(c%Kb * c%T))
            call random_number(q)
            if (pa > q) then
                Accept = .true.
            else
                Accept = .false.
            end if
        end if
    end function

    ! Return random number between -1 and 1
    real*8 function random()
        real :: x

        call random_number(x)
            x = x * 2
            if (x > 1) then
                call random_number(random) 
                random = -random
            else
                call random_number(random)
            end if
    end function
    
end module