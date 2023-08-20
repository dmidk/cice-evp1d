  module ice_constants
  use ice_kinds_mod
  implicit none
  public
  save
  real(kind=dbl_kind), parameter :: &
      c0   = 0.0_dbl_kind , &
      c1   = 1.0_dbl_kind , &
      c1p5 = 1.5_dbl_kind , &
      c2   = 2.0_dbl_kind , &
      c3   = 3.0_dbl_kind , &
      c4   = 4.0_dbl_kind , &
      c6   = 6.0_dbl_kind , &
      c9   = 9.0_dbl_kind , &
      p2   = 0.2_dbl_kind , &
      p25  = 0.25_dbl_kind, &
      p5   = 0.5_dbl_kind , &
      p111 = c1/c9        , &
      p166 = c1/c6        , &
      p222 = c2/c9        , &
      p333 = c1/c3        , &
      p055 = p111*p5      , &
      p027 = p055*p5

! This is not in ice_constants. Only a quick fix

      real(kind=dbl_kind), parameter :: puny = 1.0e-11_dbl_kind
   end module ice_constants

