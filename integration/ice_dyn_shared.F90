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
      deltaminEVP  = 10e-9_dbl_kind          , & ! default from namelist
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
         capping , & ! Default to 1. In order to mach optimization this is not a 
         epp2i   , &
         brlx    , &
         arlx    , &
         arlx1i  , &
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
      e_factor = e_yieldcurve**2 / e_plasticpot**4
      capping = c1
      arlx   = c2 * elasticDamp * real(ndte,kind=dbl_kind)
      arlx1i   = c1/arlx
      brlx   = real(ndte,kind=dbl_kind)
      denom1 = c1/(c1+arlx1i)

      end subroutine set_evp_parameters
end module ice_dyn_shared

