!===============================================================================
! Copyright (C) 2023, Intel Corporation
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!===============================================================================
module my_timer
  use ice_kinds_mod
  use ice_constants, only : c0, c1

  !- directives ----------------------------------------------------------------
  implicit none
  private

  !- module variables ----------------------------------------------------------
  integer(int_kind),   private, parameter :: idlen=10
  character(2), private, parameter        :: clen='10'       ! idlen as a char
  integer(int_kind),   private, parameter :: ntimers_s=80    ! upper bound, serial 
  integer(int_kind),   private, parameter :: ntimers_p=80    ! upper bound, parallel
  integer(int_kind),   private, parameter :: rc = 1          ! error return code
  integer(int_kind),   private, save      :: nthreads
  logical,             private, save      :: timer_active
  integer(int_kind),   private, save      :: ntim_ser=0      ! global counter

! FIXEM: HARDCODED
#if defined (_OPENMP)
  integer(int_kind),   private, save      :: timer_mode=2    ! 2=omp_wtime
  integer(int_kind),   private, save      :: ntim_par=0      ! threadlocal counter
#else
  integer(int_kind),   private, save      :: timer_mode=1    ! 1=date_and_time
#endif
  character(len=idlen), allocatable, save, private :: ccs(:), ccu(:)
  character(len=idlen), allocatable, save, private :: ccp(:,:) 
  integer(int_kind),   private, save      :: nclass
  real(kind=dbl_kind), allocatable, save, private :: t1p(:,:), t2p(:,:)
  real(kind=dbl_kind), allocatable, save, private :: t1s(:), t2s(:)
  real(kind=dbl_kind), parameter, private :: tmin0 = 1000000.0_dbl_kind
  real(kind=dbl_kind), parameter, private :: days = 24.0_dbl_kind*3600.0_dbl_kind
  real(kind=dbl_kind), parameter, private :: hrs  = 3600.0_dbl_kind
  real(kind=dbl_kind), parameter, private :: mins = 60.0_dbl_kind
  real(kind=dbl_kind), parameter, private :: secs = c1
  real(kind=dbl_kind), parameter, private :: msec = 0.001_dbl_kind
!$OMP THREADPRIVATE(ntim_par)

  !- public interface ----------------------------------------------------------
  public :: timer, timer_init, timer_destroy, timer_print

  !- private interface ---------------------------------------------------------
  private :: timer_reset
#if defined (_OPENMP)
  private :: timer_ulist
#endif
contains

!-------------------------------------------------------------------------------

  subroutine timer_init (lflag,nt,ncl)

#if defined (_OPENMP)
    use omp_lib, only : omp_in_parallel
#endif
    implicit none
    logical, intent(in)                     :: lflag
    integer(int_kind), intent(in), optional :: nt
    integer(int_kind), intent(in), optional :: ncl
    integer(int_kind)                       :: s

#if defined (_OPENMP)
    if (omp_in_parallel()) then
      stop 'Calling timer_init() in parrallel region is not allowed'
    endif
#endif

    timer_active = lflag
    if (lflag) then
      if (.not.present(nt)) then
        nthreads=1
      else
        nthreads=nt
      endif
      if (.not.present(ncl)) then
        nclass=1
      else
        nclass=ncl
      endif
      allocate (t1s(ntimers_s), t2s(ntimers_s),                                &
                t1p(ntimers_p,nthreads), t2p(ntimers_p,nthreads),              &
                ccu(ntimers_p),ccs(ntimers_s),ccp(ntimers_p,nthreads),STAT=s)
      if (s /= 0) then
        stop 'Allocation error'
      endif
      write(*,*) 'Setting timer with ', nthreads, ' threads'
    endif
    call timer_reset()

  end subroutine timer_init

!-------------------------------------------------------------------------------

  subroutine timer_reset (nt)

#if defined (_OPENMP)
    use omp_lib, only : omp_in_parallel
#endif

    implicit none
    integer(int_kind), intent(in), optional :: nt
    integer(int_kind)                       :: i,j

#if defined (_OPENMP)
    if (omp_in_parallel()) then
      stop 'Calling timer_reset() in parrallel region is not allowed'
    endif
#endif

    ! initialize the timing
    do j=1,nthreads
      do i=1,ntimers_p
        t1p(i,j) = c0
        t2p(i,j) = c0
        ccp(i,j) = ''
      enddo
    enddo
    do i=1,ntimers_s
      t1s(i) = c0
      t2s(i) = c0
      ccs(i) = ''
    enddo
    if (present(nt)) then
      write(*,*) 'Resetting timer to handle ', nthreads, ' threads'
      nthreads=nt
    endif
    do i=1,ntimers_p
      ccu(i) = ''
    enddo
  end subroutine timer_reset

!-------------------------------------------------------------------------------

  subroutine timer_destroy ()
    implicit none
    integer(int_kind) :: s
   
    if (allocated (t1s)) then   ! admitted, the test is sloppy - so sue me
      deallocate (t1s, t2s, t1p,t2p, ccs, ccp, STAT=s)
      if (s /= 0) then
        stop 'Deallocation error'
      endif
    endif
  end subroutine timer_destroy

!-------------------------------------------------------------------------------
#if defined (_OPENMP)
  subroutine timer_ulist (nu)
    use omp_lib, only : omp_in_parallel
    implicit none
    integer(int_kind),intent(out) :: nu
    integer(int_kind) :: i,j,k
    if (omp_in_parallel()) then
      stop 'Calling timer_ulist() in parrallel region is not allowed'
    endif
    nu=1
    do j=1,nthreads
      do i=1,ntimers_p
        if (ccp(i,j) /= '') then
          k=1
          whileloop: do while (k<nu) 
             if (ccp(i,j) == ccu(k)) exit whileloop
             k=k+1
          enddo whileloop
          if (k==nu) then
            ! unique id
            ccu(nu)=ccp(i,j)
            nu=nu+1
!            write (*,*), 'Unique parallel id:', ccu(nu), nu
          endif
        endif
      enddo
    enddo
  end subroutine timer_ulist 
#endif

  subroutine timer_print (string,id,nnclass)

    implicit none
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: string
    integer(int_kind),       intent(in), optional :: nnclass

    real(kind=dbl_kind), parameter :: tfmtmax = 100000.0_dbl_kind
    integer(int_kind) :: i
#if defined (_OPENMP)
    integer(int_kind) :: j, nu, k
    ! threadlocal stats
    real(kind=dbl_kind) :: t_max(ntimers_p), t_min(ntimers_p),t_ave(ntimers_p)
    integer(int_kind) :: t_hit(ntimers_p), t_max_loc(ntimers_p), t_min_loc(ntimers_p)
#endif

    if (present(nnclass)) then
      if (nnclass > nclass) return

#if defined (_OPENMP)
!$OMP BARRIER
!$OMP MASTER
      call timer_ulist(nu)
      do i=1,nu
        t_max(i) = c0
        t_min(i) = tmin0
        t_ave(i) = c0
        t_hit(i) = 0
        t_max_loc(i)=-1
        t_min_loc(i)=-1
      enddo
      if (present(id)) then
        do j=1,nthreads
          do i=1,ntimers_p
            if (ccp(i,j) == id) then
              kloop: do k=1,nu
                if (ccu(k)==id) then
                  if (t2p(i,j) > t_max(k)) then
                    t_max(k)=t2p(i,j) 
                    t_max_loc(k)=j
                  endif
                  if (t2p(i,j) < t_min(k)) then
                    t_min(k)=t2p(i,j) 
                    t_min_loc(k)=j
                  endif
                  t_ave(k) = t_ave(k) + t2p(i,j)
                  t_hit(k) = t_hit(k)+1
                  exit kloop
                endif
              enddo kloop
            endif
          enddo
        enddo
        do k=1,nu
          if (ccu(k)==id) then
            t_ave(k) = t_ave(k)/nthreads
            if (present(string)) then
              ! Both string and id present for identification:
              write (*,'(a,a1,i4,3f13.3)') string,' ',t_hit(k),t_min(k),    &
                                              t_max(k), t_ave(k)
            else
              write (*,'(a,a,a2,i4,3f13.3)')                                &
                         'Thread timing information (min, max, average) for ', &
                          id, ' :', t_hit(k), t_min(k), t_max(k), t_ave(k)
            endif
            write (*,'(a,a,a2,2i4)')                                        &
                         'Thread timing loc (min, max) for ',                  &
                          id, ' :', t_min_loc(k), t_max_loc(k)
          endif
        enddo
      endif
!$OMP END MASTER
#endif
      ! same round but now for serial timers, only
      if (present(string) .and. present(id)) then
        ! Both string and id present for identification:
        do i = 1, ntim_ser
          if (ccs(i) == id) then
            if (t2s(i) < tfmtmax) then
              write (*,'(a,a1,f10.2,a8)') string, ' ', t2s(i), ' seconds'
            else
              write (*,*) string, ' ', t2s(i), ' seconds'
            endif
            t1s(i)=c0
            t2s(i)=c0
          endif
        enddo
      elseif (present(id)) then
        do i = 1, ntim_ser
          if (ccs(i) == id) then
            write(*,*) 'Timing information for ',id,' : ',ccs(i),t2s(i)
          endif
        enddo
      else
        do i = 1, ntim_ser
          write(*,*) 'Timing information for ', ccs(i), ' : ', t2s(i)
        enddo
      endif
    endif
  end subroutine timer_print

!-------------------------------------------------------------------------------
  subroutine timer (mode,id,lrestart)
    implicit none
    integer(int_kind),       intent(in)           :: mode
    character(len=*), intent(in)           :: id
    logical,          intent(in), optional :: lrestart
    logical    :: llrestart
    if (.not. timer_active) return
    if (len(id) > idlen) then
      stop 'size of id exceeds the allowed size '
    endif
    if (.not.present(lrestart)) then
      llrestart=.true. ! default is to restart the timer
    else
      llrestart=lrestart 
    endif
    select case (timer_mode)
#if defined (_OPENMP)
    case (1)
     call timer_date_and_time(mode,id,llrestart)
    case (2)
     call timer_omp_wtime(mode,id,llrestart)
#else
    case (1)
     call timer_date_and_time(mode,id,llrestart)
    case (2)
     call timer_date_and_time(mode,id,llrestart)
#endif
! ....
    end select
  end subroutine timer

#if defined (_OPENMP)
  subroutine timer_omp_wtime(mode,id,llrestart)
    use omp_lib, only : omp_in_parallel, omp_get_wtime, omp_get_thread_num
    implicit none
    ! input/output parameters
    integer(int_kind),       intent(in) :: mode
    character(len=*), intent(in) :: id
    logical,          intent(in) :: llrestart
    ! local parameters
    integer(int_kind) :: i,tid
    real(kind=dbl_kind)    :: mtime

    if (omp_in_parallel()) then
      tid=omp_get_thread_num()+1
      select case (mode)
      case (1)
        ! start timing of 'id'
        mtime = omp_get_wtime()
        do i = 1, ntim_par
          if (ccp(i,tid) == id) then
            if (llrestart) then
              ! restarting timer
              t1p(i,tid) = mtime
            endif
            return
          endif
        enddo
        if (ntim_par+1 > ntimers_p) then
          stop 'Trying to add more timers than we have reserved space for'
        endif
        ntim_par     = ntim_par + 1
        t1p(ntim_par,tid) = mtime
        t2p(ntim_par,tid) = c0
        ccp(ntim_par,tid) = id
      case (2)
        ! stop timing  of 'id'
        mtime = omp_get_wtime()
        do i = 1, ntim_par
          if (ccp(i,tid) == id) then
            t2p(i,tid) = t2p(i,tid) + mtime - t1p(i,tid)
            return
          endif
        enddo
        if (i == ntim_par+1) then 
         stop 'Trying to stop timer which has not been started'
        endif
      end select
    else
      select case (mode)
      case (1)
        ! start timing of 'id'
        mtime = omp_get_wtime()
        do i = 1, ntim_ser
          if (ccs(i) == id) then
            if (llrestart) then
              ! restarting timer
              t1s(i) = mtime
            endif
            return
          endif
        enddo
        if (ntim_ser+1 > ntimers_s) then
          stop 'Trying to add more timers than we have reserved space for'
        endif
        ntim_ser     = ntim_ser + 1
        t1s(ntim_ser) = mtime
        t2s(ntim_ser) = c0
        ccs(ntim_ser) = id
      case (2)
        ! stop timing  of 'id'
        mtime = omp_get_wtime()
        do i = 1, ntim_ser
          if (ccs(i) == id) then
            t2s(i) = t2s(i) + mtime - t1s(i)
            return
          endif
        enddo
        if (i == ntim_ser+1) then 
         stop 'Trying to stop timer which has not been started'
        endif
      end select
    endif
  end subroutine timer_omp_wtime
#endif

  subroutine timer_date_and_time(mode,id,llrestart)
#if defined (_OPENMP)
    use omp_lib, only : omp_in_parallel, omp_get_thread_num
#endif
    implicit none

    ! input/output parameters
    integer(int_kind),       intent(in) :: mode
    character(len=*), intent(in) :: id
    logical,          intent(in) :: llrestart

    ! local parameters
    integer(int_kind) :: ta(8)
    integer(int_kind) :: i
    real(kind=dbl_kind)    :: mtime

#if defined (_OPENMP)
    integer(int_kind) :: tid
    tid=1
    if (omp_in_parallel()) then
      tid=omp_get_thread_num()+1
      select case (mode)
      case (1)
        ! start timing of 'id'
        call date_and_time(values=ta)
        mtime = ta(3)*(days) + & ! number of days
                ta(5)*(hrs)  + & ! number of hours
                ta(6)*(mins) + & ! number of minutes
                ta(7)*(secs) + & ! number of seconds
                ta(8)*(msec)     ! number of milliseconds
        do i = 1, ntim_par
          if (ccp(i,tid) == id) then
            if (llrestart) then
              ! restarting timer
              t1p(i,tid) = mtime
            endif
            return
          endif
        enddo
        if (ntim_par+1 > ntimers_p) then
          stop 'Trying to add more timers than we have reserved space for'
        endif
        ntim_par     = ntim_par + 1
        t1p(ntim_par,tid) = mtime
        t2p(ntim_par,tid) = c0
        ccp(ntim_par,tid) = id
      case (2)
        ! stop timing  of 'id'
        call date_and_time(values=ta)
        mtime = ta(3)*(days) + & ! number of days
                ta(5)*(hrs)  + & ! number of hours
                ta(6)*(mins) + & ! number of minutes
                ta(7)*(secs) + & ! number of seconds
                ta(8)*(msec)     ! number of milliseconds
        do i = 1, ntim_par
          if (ccp(i,tid) == id) then
            t2p(i,tid) = t2p(i,tid) + mtime - t1p(i,tid)
            return
          endif
        enddo
        if (i == ntim_par+1) then 
         stop 'Trying to stop timer which has not been started'
        endif
      end select
    else
      select case (mode)
      case (1)
        ! start timing of 'id'
        call date_and_time(values=ta)
        mtime = ta(3)*(days) + & ! number of days
                ta(5)*(hrs)  + & ! number of hours
                ta(6)*(mins) + & ! number of minutes
                ta(7)*(secs) + & ! number of seconds
                ta(8)*(msec)     ! number of milliseconds
        do i = 1, ntim_ser
          if (ccs(i) == id) then
            if (llrestart) then
              ! restarting timer
              t1s(i) = mtime
            endif
            return
          endif
        enddo
        if (ntim_ser+1 > ntimers_s) then
          stop 'Trying to add more timers than we have reserved space for'
        endif
        ntim_ser     = ntim_ser + 1
        t1s(ntim_ser) = mtime
        t2s(ntim_ser) = c0
        ccs(ntim_ser) = id
      case (2)
        ! stop timing  of 'id'
        call date_and_time(values=ta)
        mtime = ta(3)*(days) + & ! number of days
                ta(5)*(hrs)  + & ! number of hours
                ta(6)*(mins) + & ! number of minutes
                ta(7)*(secs) + & ! number of seconds
                ta(8)*(msec)     ! number of milliseconds
        do i = 1, ntim_ser
          if (ccs(i) == id) then
            t2s(i) = t2s(i) + mtime - t1s(i)
            return
          endif
        enddo
        if (i == ntim_ser+1) then 
         stop 'Trying to stop timer which has not been started'
        endif
      end select
    endif
#else
    select case (mode)
    case (1)
      ! start timing of 'id'
      call date_and_time(values=ta)
      mtime = ta(3)*(days) + & ! number of days
              ta(5)*(hrs)  + & ! number of hours
              ta(6)*(mins) + & ! number of minutes
              ta(7)*(secs) + & ! number of seconds
              ta(8)*(msec)     ! number of milliseconds
      do i = 1, ntim_ser
        if (ccs(i) == id) then
          if (llrestart) then
            ! restarting timer
            t1s(i) = mtime
          endif
          return
        endif
      enddo
      if (ntim_ser+1 > ntimers_s) then
        stop 'Trying to add more timers than we have reserved space for'
      endif
      ntim_ser     = ntim_ser + 1
      t1s(ntim_ser) = mtime
      t2s(ntim_ser) = c0
      ccs(ntim_ser) = id
    case (2)
      ! stop timing  of 'id'
      call date_and_time(values=ta)
      mtime = ta(3)*(days) + & ! number of days
              ta(5)*(hrs)  + & ! number of hours
              ta(6)*(mins) + & ! number of minutes
              ta(7)*(secs) + & ! number of seconds
              ta(8)*(msec)     ! number of milliseconds
      do i = 1, ntim_ser
        if (ccs(i) == id) then
          t2s(i) = t2s(i) + mtime - t1s(i)
          return
        endif
      enddo
      if (i == ntim_ser+1) then 
         stop 'Trying to stop timer which has not been started'
      endif
    end select
#endif
  end subroutine timer_date_and_time
end module my_timer
