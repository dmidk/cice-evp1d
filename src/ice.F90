!===============================================================================
module ice_kinds_mod
  implicit none
  public
  save
  integer, parameter :: char_len  = 80, &
                        char_len_long  = 256, &
                        log_kind  = kind(.true.), &
                        int_kind  = selected_int_kind(6), &
                        real_kind = selected_real_kind(6), &
                        dbl_kind  = selected_real_kind(13), &
                        r16_kind  = selected_real_kind(26)
end module ice_kinds_mod
!===============================================================================
module ice_constants
  use ice_kinds_mod
  implicit none
  save
  private
#ifdef _OPENMP_TARGET
  !$omp declare target(c0, c1, p027, p055, p111, p166, c1p5, p2, p222, p25,    &
  !$omp                p333, p5, puny, arlx1i, denom1, Ktens, revp,      &
  !$omp                capping, deltaminEVP, e_factor, epp2i, brlx)
#endif
  real (kind=dbl_kind), parameter, public ::                                   &
    c1p5 = 1.5_dbl_kind, &
    c0   = 0.0_dbl_kind, &
    c1   = 1.0_dbl_kind, &
    c2   = 2.0_dbl_kind, &
    c3   = 3.0_dbl_kind, &
    c4   = 4.0_dbl_kind, &
    c5   = 5.0_dbl_kind, &
    c6   = 6.0_dbl_kind, &
    c8   = 8.0_dbl_kind, &
    c9   = 9.0_dbl_kind, &
    c10  = 10.0_dbl_kind, &
    c12  = 12.0_dbl_kind, &
    c15  = 15.0_dbl_kind, &
    c16  = 16.0_dbl_kind, &
    c20  = 20.0_dbl_kind, &
    c25  = 25.0_dbl_kind, &
    c30  = 30.0_dbl_kind, &
    c100 = 100.0_dbl_kind, &
    c180 = 180.0_dbl_kind, &
    c360 = 360.0_dbl_kind, &
    c365 = 365.0_dbl_kind, &
    c400 = 400.0_dbl_kind, &
    c3600= 3600.0_dbl_kind, &
    c1000= 1000.0_dbl_kind, &
    p001 = 0.001_dbl_kind, &
    p01  = 0.01_dbl_kind, &
    p025 = 0.025_dbl_kind, &
    p1   = 0.1_dbl_kind, &
    p2   = 0.2_dbl_kind, &
    p4   = 0.4_dbl_kind, &
    p5   = 0.5_dbl_kind, &
    p6   = 0.6_dbl_kind, &
    p05  = 0.05_dbl_kind, &
    p15  = 0.15_dbl_kind, &
    p25  = 0.25_dbl_kind, &
    p75  = 0.75_dbl_kind, &
    p166 = c1/c6, &
    p333 = c1/c3, &
    p666 = c2/c3, &
    p111 = c1/c9, &
    p055 = p111*p5, &
    p027 = p055*p5, &
    p222 = c2/c9
!    puny   = 1.0e-11_dbl_kind
  real (kind=dbl_kind), parameter, public :: &
    rhow   = 1026.0_dbl_kind         ! density of seawater (kg/m^3)
end module ice_constants
  module ice_dyn_shared
   use ice_kinds_mod
   use ice_constants
  implicit none
  public set_evp_parameters
  save
! Hardcoded from except for u0, cosw and sinw
  real(kind=dbl_kind), parameter, public            :: &
      u0           = 5e-5_dbl_kind           , & ! hardcoded as in ice_dyn_shared
      cosw         = c1                      , & ! hardcoded as in ice_dyn_shared
      sinw         = c0                      , & ! hardcoded as in ice_dyn_shared
! Default variables from namelist
      elasticDamp  = 0.36_dbl_kind           , &
! variables from namelist
      Ktens        = c0                      , & ! assumes Ktens 0
      revp         = c0                      , & ! assumes revised_evp false
      e_plasticpot = c2                      , &
      e_yieldcurve = c2

! from namelist 
      logical (kind=log_kind), parameter, public :: &
         revised_evp = .false.   ! if true, use revised evp procedure

      integer (kind=int_kind), public        :: &
         ndte
!  Calculated here 
      real(kind=dbl_kind), public            :: &
         e_factor, &       ! assume e_yieldcurve=2 and e_plastpot = 2
         capping,  &
         epp2i    , &
         brlx    , &
         arlx    , &
         arlx1i  , &
         deltaminEVP, & !Read in. Normally set to 10^-9
         denom1
contains
      subroutine set_evp_parameters (dt)
! reduced version. Only contain needed parameters
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

!      character(len=*), parameter :: subname = '(set_evp_parameters)'

      ! elastic time step
!      dtei = real(ndte,kind=dbl_kind)/dt

      ! variables for elliptical yield curve and plastic potential
      epp2i = c1/e_plasticpot**2
      capping = c1  ! This is the default value
      e_factor = e_yieldcurve**2 / e_plasticpot**4
      arlx   = c2 * elasticDamp * real(ndte,kind=dbl_kind)
      arlx1i   = c1/arlx
      brlx   = real(ndte,kind=dbl_kind)
      denom1 = c1/(c1+arlx1i)
      deltaminEVP  = 10e-9_dbl_kind          ! default from namelist
      end subroutine set_evp_parameters
  end module ice_dyn_shared
