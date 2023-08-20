!===============================================================================
!!#ifdef standalone
module ice_gather_scatter
! This is a very light version of gather_scatter routines
! as it only copies from 3d to 2D in stand alone 
  use ice_kinds_mod
  implicit none
  public
  save
  integer(kind=int_kind), parameter :: max_block=1
  interface gather_global_ext
     module procedure gather_global_ext_dbl,  &
                      gather_global_ext_int,  &
                      gather_global_ext_log
  end interface

 contains

  subroutine gather_global_ext_dbl(ARRAY_G, ARRAY, dst_task, src_dist, spc_val)
   use ice_kinds_mod
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task, including ghost cells.

   integer (kind=int_kind), intent(in) :: &
     dst_task, src_dist   ! task to which array should be gathered

   real(kind=dbl_kind), dimension(:,:,:), intent(in) :: &
     ARRAY      ! array containing horizontal slab of distributed field

   real (kind=dbl_kind), optional :: &
     spc_val

   real (kind=dbl_kind), dimension(:,:),  intent(inout) :: &
     ARRAY_G    ! array containing global horizontal field on dst_task
! not sure if it is standard
    ARRAY_G = ARRAY(:,:,max_block)
   end subroutine gather_global_ext_dbl

   subroutine gather_global_ext_int(ARRAY_G, ARRAY, dst_task, src_dist, spc_val)

!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task, including ghost cells.

   integer (kind=int_kind), intent(in) :: &
     dst_task, src_dist  ! task to which array should be gathered

   integer (kind=int_kind), dimension(:,:,:), intent(in) :: &
     ARRAY      ! array containing horizontal slab of distributed field

   real (kind=int_kind), optional :: &
     spc_val

   integer (kind=int_kind), dimension(:,:), intent(inout) :: &
     ARRAY_G    ! array containing global horizontal field on dst_task

    ARRAY_G = ARRAY(:,:,max_block)
   end subroutine gather_global_ext_int

   subroutine gather_global_ext_log(ARRAY_G, ARRAY, dst_task, src_dist, spc_val)

!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task, including ghost cells.

   integer (int_kind), intent(in) :: &
     dst_task, src_dist   ! task to which array should be gathered

   logical (log_kind), dimension(:,:,:), intent(in) :: &
     ARRAY      ! array containing horizontal slab of distributed field

   logical (log_kind), optional :: &
     spc_val

   logical (log_kind), dimension(:,:), intent(inout) :: &
     ARRAY_G    ! array containing global horizontal field on dst_task

    ARRAY_G = ARRAY(:,:,max_block)
   end subroutine gather_global_ext_log

   end module ice_gather_scatter
