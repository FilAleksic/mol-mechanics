module EnergyModule
    use TypesModule
    use DistanceModule
    implicit none
    
    PUBLIC :: Energy
    PRIVATE
    
contains

    real*8 function Energy(Atom_Array, Bonding_Info)
        type(Atom)    :: Atom_Array(:)
        type(Bonding) :: Bonding_Info

        Energy = 0
        Energy = Energy + Stretch(Atom_Array) + nbElectrostatic(Atom_Array) + nbVanDerWaals(Atom_Array)
        if (allocated(Bonding_Info%Chain_T)) then
            Energy = Energy + Bending(Bonding_Info)
        end if
        if (allocated(Bonding_Info%Chain_F)) then
            Energy = Energy + Torsional(Bonding_Info)
        end if
    end function

    real*8 function Stretch(Atom_Array)
        type(Atom)     :: Atom_Array(:)
        type(Atom)     :: A
        integer        :: i, j
        real*8         :: R
        type(Constant) :: c
        Stretch = 0.0
        do i=1,SIZE(Atom_Array)
            do j=1,SIZE(Atom_Array(i)%Bond)
                A = Atom_Array(i) 
                R = Abs_Distance(A, A%Bond(j)%p)
                ! Note, factor of 0.5 added as all atom pairs are calculated twice (AB and BA)
                ! By simply adding this factor we eliminate this issue
                if (A%symbol // A%Bond(j)%p%symbol == 'CC') then
                    Stretch = Stretch + (0.5 * c%KCC * (R-c%R0CC)**2)
                else ! CH Bonds
                    Stretch = Stretch + (0.5 * c%KCH * (R-c%R0CH)**2)
                end if
            end do
        end do
    end function

    real*8 function Bending(Bonding_Info)
        type(Bonding)           :: Bonding_Info
        type(pnt), DIMENSION(3) :: x
        integer                 :: i
        type(Constant)          :: c
        character(len=3)        :: bond
        Bending = 0
        
        do i=1,SIZE(Bonding_Info%Chain_T(:,1))
            x = Bonding_Info%Chain_T(i,:)
            bond = x(1)%p%symbol // x(2)%p%symbol // x(3)%p%symbol
            if (bond == 'CCC') then
                Bending = Bending + c%KCCC * (Angle(x) - ((c%B_Angle*c%pi)/180))
            else if (bond == 'CCH' .or. bond == 'HCC') then
                Bending = Bending + c%KCCH * (Angle(x) - ((c%B_Angle*c%pi)/180))
            else ! HCH Bonds
                Bending = Bending + c%KHCH * (Angle(x) - ((c%B_Angle*c%pi)/180))
            end if
        end do
    end function

    ! Torsional angle ABCD calculated with vector 1 -> BA and vector 2 -> CD
    real*8 function Torsional(Bonding_Info)
        type(Bonding)           :: Bonding_Info
        type(pnt), DIMENSION(4) :: x
        integer                 :: i
        type(Constant)          :: c
        Torsional = 0
        do i=1,SIZE(Bonding_Info%Chain_F(:,1))
            x = Bonding_Info%Chain_F(i,:)
            Torsional = Torsional + 0.5 * c%V * (1 + cos(c%n * Angle(x) - c%gamma))
        end do
    end function

    real*8 function nbElectrostatic(Atom_Array)
        type(Atom)     :: Atom_Array(:)
        integer        :: i, j
        real*8         :: rij
        type(Constant) :: c

        nbElectrostatic = 0
        do i=1,SIZE(Atom_Array)-1
            do j=i+1,SIZE(Atom_Array)
                rij = Abs_Distance(Atom_Array(i), Atom_Array(j))
                if (Atom_Array(i)%symbol // Atom_Array(j)%symbol == 'CC') then
                    nbElectrostatic = nbElectrostatic + (c%qC * c%qC * c%Ke) / rij
                else if (Atom_Array(i)%symbol // Atom_Array(j)%symbol == 'CH' &
                    .or. Atom_Array(i)%symbol // Atom_Array(j)%symbol == 'HC') then
                    nbElectrostatic = nbElectrostatic + (c%qC * c%qH * c%Ke) / rij
                else if (Atom_Array(i)%symbol // Atom_Array(j)%symbol == 'HH') then
                    nbElectrostatic = nbElectrostatic + (c%qH * c%qH * c%Ke) / rij
                end if
            end do
        end do
    end function

    real*8 function nbVanDerWaals(Atom_Array)
        type(Atom)     :: Atom_Array(:)
        integer        :: i, j
        real*8         :: rij, eij, WDepth, Aij, Bij
        type(Constant) :: c

        nbVanDerWaals = 0
        do i=1,SIZE(Atom_Array)
            do j=i,SIZE(Atom_Array)
                if (i /= j) then
                    rij = Abs_Distance(Atom_Array(i), Atom_Array(j))
                    if (Atom_Array(i)%symbol // Atom_Array(j)%symbol == 'CC') then
                        eij = sqrt(c%eC * c%eC)
                        WDepth = c%WDepthC + c%WDepthC
                    else if (Atom_Array(i)%symbol // Atom_Array(j)%symbol == 'CH' &
                        .or. Atom_Array(i)%symbol // Atom_Array(j)%symbol == 'HC') then
                        eij = sqrt(c%eC * c%eH)
                        WDepth = c%WDepthC + c%WDepthH
                    else if (Atom_Array(i)%symbol // Atom_Array(j)%symbol == 'HH') then
                        eij = sqrt(c%eH * c%eH)
                        WDepth = c%WDepthH + c%WDepthH
                    end if
                    Aij = eij * Wdepth
                    Bij = 2 * eij * WDepth
                    nbVanDerWaals = nbVanDerWaals + (Aij / (rij**12)) - (Bij / (rij**6))
                end if
            end do
        end do
    end function

    ! General function to calculate angle between pairs of 3 or 4 atoms
    ! Angle between ABC calculated by taking vector 1 ->BA and Vector 2->BC 
    ! Torsional angle ABCD calculated with vector 1 -> BA and vector 2 -> CD
    ! Dot Product V1 . V2 / Magnitude(V1) * Magnitude(V2) = Cos(Angle)
    real*8 function Angle(Array)
        type(pnt)           :: Array(:)
        real*8, ALLOCATABLE :: Vector1(:), Vector2(:)
        integer             :: sz
        sz = SIZE(Array)
        allocate(Vector1(sz), Vector2(sz))
        Vector1 = (/Array(1)%p%x - Array(2)%p%x, Array(1)%p%y - Array(2)%p%y, Array(1)%p%z - Array(2)%p%z/)
        Vector2 = (/Array(sz)%p%x - Array(sz-1)%p%x, Array(sz)%p%y - Array(sz-1)%p%y, Array(sz)%p%z - Array(sz-1)%p%z/)

        Angle = acos(DOT_PRODUCT(Vector1, Vector2) / (NORM2(Vector1) * NORM2(Vector2)))
    end function
    
end module