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
