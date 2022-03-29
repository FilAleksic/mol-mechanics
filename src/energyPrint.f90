module EnergyPrintModule
    use TypesModule
    use DistanceModule
    implicit none
    
    PUBLIC :: EnergyP
    PRIVATE

    ! Parameters in kcal, distances in A, angles in degrees
    ! Stretch force constants and equilibrium distances
    real*8, PARAMETER :: KCH = 340, KCC = 310, R0CH = 1.090, R0CC = 1.526
    ! Bending force constants and equilibrium angles 
    real*8, PARAMETER :: KCCC = 40, KCCH = 50, KHCH = 35
    real*8, PARAMETER :: B_Angle = 109.50
    real*8, PARAMETER :: pi = 3.141592653
    ! Torsional parameters
    ! (Note, in the paper a value of V=1.4 is given,
    ! however their value was already divided by 2. For readability 
    ! of the formula I do that in the code so my V here is 2.8)
    integer, PARAMETER :: n = 3
    real*8, PARAMETER  :: V = 2.8, gamma = 0
    ! Electrostatic parameters
    real*8, PARAMETER :: qC = -0.464, qH = 0.116
    real*8, PARAMETER :: eC = 0.1094, eH = 0.0157
    ! Van der Waals parameters 
    real*8, PARAMETER :: RstarC = 1.9080, RstarH = 0.6000


contains

    real*8 function EnergyP(Atom_Array, Bonding_Info)
        type(Atom) :: Atom_Array(:)
        type(Bonding) :: Bonding_Info

        EnergyP = 0
        print*, 'Str'
        EnergyP = EnergyP + Stretch(Atom_Array) + Bending(Bonding_Info) + Torsional(Bonding_Info) &
                          + nbElectrostatic(Atom_Array) + nbVanDerWaals(Atom_Array)
    end function

    real*8 function Stretch(Atom_Array)
        type(Atom) :: Atom_Array(:)
        type(Atom) :: A
        integer :: i, j
        real*8 :: R
        Stretch = 0.0
        do i=1,SIZE(Atom_Array)
            do j=1,SIZE(Atom_Array(i)%Bond)
                A = Atom_Array(i) 
                R = Abs_Distance(A, A%Bond(j)%p)
                ! Note, factor of 0.5 added as all atom pairs are calculated twice (AB and BA)
                ! By simply adding this factor we eliminate this issue
                if (A%symbol // A%Bond(j)%p%symbol == 'CC') then
                    Stretch = Stretch + (0.5 * KCC * (R-R0CC)**2)
                else ! CH Bonds
                    Stretch = Stretch + (0.5 * KCH * (R-R0CH)**2)
                end if
            end do
        end do
        print*, Stretch, 'STRETCH ENERGY'
    end function

    real*8 function Bending(Bonding_Info)
        type(Bonding) :: Bonding_Info
        type(pnt), DIMENSION(3) :: x
        integer :: i
        character(len=3) :: bond
        Bending = 0
        
        do i=1,SIZE(Bonding_Info%Chain_T(:,1))
            x = Bonding_Info%Chain_T(i,:)
            bond = x(1)%p%symbol // x(2)%p%symbol // x(3)%p%symbol
            if (bond == 'CCC') then
                Bending = Bending + KCCC * (Angle(x) - ((B_Angle*pi)/180))
            else if (bond == 'CCH' .or. bond == 'HCC') then
                Bending = Bending + KCCH * (Angle(x) - ((B_Angle*pi)/180))
            else ! HCH Bonds
                Bending = Bending + KHCH * (Angle(x) - ((B_Angle*pi)/180))
            end if
        end do
        print*, Bending, 'BENDING ENERGY'
    end function

    ! Torsional angle ABCD calculated with vector 1 -> BA and vector 2 -> CD
    real*8 function Torsional(Bonding_Info)
        type(Bonding) :: Bonding_Info
        type(pnt), DIMENSION(4) :: x
        integer :: i
        Torsional = 0
        do i=1,SIZE(Bonding_Info%Chain_F(:,1))
            x = Bonding_Info%Chain_F(i,:)
            Torsional = Torsional + 0.5 * V * (1 + cos(n * Angle(x) - gamma))
        end do
        print*, Torsional, 'TORSIONAL ENERGY'
    end function

    real*8 function nbElectrostatic(Atom_Array)
        type(Atom) :: Atom_Array(:)
        integer :: i, j
        real*8 :: rij, eij

        nbElectrostatic = 0
        do i=1,SIZE(Atom_Array)-1
            do j=i+1,SIZE(Atom_Array)
                rij = Abs_Distance(Atom_Array(i), Atom_Array(j))
                if (Atom_Array(i)%symbol // Atom_Array(j)%symbol == 'CC') then
                    eij = sqrt(eC * eC)
                    nbElectrostatic = nbElectrostatic + (qC * qC) / (4 * pi * eij * rij)
                else if (Atom_Array(i)%symbol // Atom_Array(j)%symbol == 'CH' &
                    .or. Atom_Array(i)%symbol // Atom_Array(j)%symbol == 'HC') then
                    eij = sqrt(eC * eH)
                    nbElectrostatic = nbElectrostatic + (qC * qH) / (4 * pi * eij * rij)
                else if (Atom_Array(i)%symbol // Atom_Array(j)%symbol == 'HH') then
                    eij = sqrt(eH * eH)
                    nbElectrostatic = nbElectrostatic + (qH * qH) / (4 * pi * eij * rij)
                end if
            end do
        end do
        print*, nbElectrostatic, 'ELECTROSTATIC ENERGY'
    end function

    real*8 function nbVanDerWaals(Atom_Array)
        type(Atom) :: Atom_Array(:)
        integer :: i, j
        real*8 :: rij, eij, Aij, Bij!, Rstarij

        nbVanDerWaals = 0
        do i=1,SIZE(Atom_Array)-1
            do j=i+1,SIZE(Atom_Array)
                rij = Abs_Distance(Atom_Array(i), Atom_Array(j))
                if (Atom_Array(i)%symbol // Atom_Array(j)%symbol == 'CC') then
                    eij = sqrt(eC * eC)
                    !Rstarij = RstarC + RstarC
                    Aij = eij! * (Rstarij**12)
                    Bij = 2 * eij! * (Rstarij**6)
                    nbVanDerWaals = nbVanDerWaals + (Aij / (rij**12)) - (Bij / (rij**6))
                else if (Atom_Array(i)%symbol // Atom_Array(j)%symbol == 'CH' &
                    .or. Atom_Array(i)%symbol // Atom_Array(j)%symbol == 'HC') then
                    eij = sqrt(eC * eH)
                    !Rstarij = RstarC + RstarH
                    Aij = eij! * (Rstarij**12)
                    Bij = 2 * eij! * (Rstarij**6)
                    nbVanDerWaals = nbVanDerWaals + (Aij / (rij**12)) - (Bij / (rij**6))
                else if (Atom_Array(i)%symbol // Atom_Array(j)%symbol == 'HH') then
                    eij = sqrt(eH * eH)
                    !Rstarij = RstarH + RstarH
                    Aij = eij! * (Rstarij**12)
                    Bij = 2 * eij! * (Rstarij**6)
                    nbVanDerWaals = nbVanDerWaals + (Aij / (rij**12)) - (Bij / (rij**6))
                end if
            end do
        end do
        print*, nbVanDerWaals, 'VDW ENERGY'
    end function

    ! General function to calculate angle between pairs of 3 or 4 atoms
    ! Angle between ABC calculated by taking vector 1 ->BA and Vector 2->BC 
    ! Torsional angle ABCD calculated with vector 1 -> BA and vector 2 -> CD
    ! Dot Product V1 . V2 / Magnitude(V1) * Magnitude(V2) = Cos(Angle)
    real*8 function Angle(Array)
        type(pnt) :: Array(:)
        real*8, ALLOCATABLE :: Vector1(:), Vector2(:)
        integer :: sz
        sz = SIZE(Array)
        allocate(Vector1(sz), Vector2(sz))
        Vector1 = (/Array(1)%p%x - Array(2)%p%x, Array(1)%p%y - Array(2)%p%y, Array(1)%p%z - Array(2)%p%z/)
        Vector2 = (/Array(sz)%p%x - Array(sz-1)%p%x, Array(sz)%p%y - Array(sz-1)%p%y, Array(sz)%p%z - Array(sz-1)%p%z/)

        Angle = acos(DOT_PRODUCT(Vector1, Vector2) / (NORM2(Vector1) * NORM2(Vector2)))
    end function
    
end module