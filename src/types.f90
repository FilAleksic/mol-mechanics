module TypesModule
    implicit none

    PUBLIC

    type pnt
        type(Atom), POINTER :: p
        logical             :: Del = .false.
    end type

    type Atom
        real*8                 :: x, y, z
        character(len=1)       :: symbol
        type(pnt), ALLOCATABLE :: Bond(:)
    end type

    type Bonding
        type(pnt), ALLOCATABLE :: Chain_F(:,:)
        type(pnt), ALLOCATABLE :: Chain_T(:,:)
    end type

    type Constant
        ! Parameters in kcal, distances in A, angles in degrees
        ! Stretch force constants and equilibrium distances
        real*8 :: KCH = 340, KCC = 310, R0CH = 1.090, R0CC = 1.526
        ! Bending force constants and equilibrium angles 
        real*8 :: KCCC = 40, KCCH = 50, KHCH = 35
        real*8 :: B_Angle = 109.50
        real*8 :: pi = 3.141592653
        ! Torsional parameters
        ! (Note, in the paper a value of V=1.4 is given,
        ! however their value was already divided by 2. For readability of the formula I do that in the code so my V here is 2.8)
        integer :: n = 3
        real*8  :: V = 2.8, gamma = 0
        ! Electrostatic parameters. Ke in eV * A * e^-2
        real*8 :: qC = -0.464, qH = 0.116
        real*8 :: Ke = 14.3996
        ! VDW parameters
        real*8 :: WDepthC = 1.9080, WDepthH = 1.4870
        real*8 :: eC = 0.1094, eH = 0.0157
        real*8 :: Avogardro = 6.022140857E23

        ! Minimizing parameters
        real*8 :: T = 297, R = 0.00001, Kb = 1.380649E-23
    end type

    
end module