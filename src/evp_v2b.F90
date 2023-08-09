!===============================================================================
! Copyright (C) 2023, Intel Corporation
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!===============================================================================
#ifdef sde
module ssc_fortran
use, intrinsic :: iso_c_binding
interface
   subroutine fortran_sde_start() &
      bind(c, name='fortran_ssc_start')
   end subroutine fortran_sde_start
   subroutine fortran_sde_stop() &
      bind(c, name='fortran_ssc_stop')
   end subroutine fortran_sde_stop
end interface
end module
#endif
!===============================================================================
program ciceevp
#ifdef itt
  use ittnotify
#endif
#ifdef sde
  use ssc_fortran,   only : fortran_sde_start, fortran_sde_stop
#endif
  use ice_kinds_mod, only : int_kind, dbl_kind
  use ice_dyn_shared,only : ndte, set_evp_parameters
  use create_nml,    only : nx_block, ny_block, inpfname, read_nml,            &
                            testscale, na, navel, lindividual, binoutput
  use my_timer,      only : timer, timer_init, timer_print
  use myomp,         only : domp_init
  use numa,          only : numainit_all, numareinit
  use bench,         only : stress, stepu
  use vars,          only : read_nml, alloc_1d_v2, stat_1d, readin_1d_v2,      &
                            stat_out_1d, writeout_1d, replicate,               &
                            uvel, vvel, dxT, dyT, strength,                    &
                            stressp_1, stressp_2, stressp_3, stressp_4,        &
                            stressm_1, stressm_2, stressm_3, stressm_4,        &
                            stress12_1, stress12_2, stress12_3, stress12_4,    &
                            cdn_ocn, aiu, uocn, vocn, waterxU, wateryU,        &
                            forcexU, forceyU, umassdti, fmU, uarear, strintxU, &
                            strintyU, uvel_init, vvel_init, Tbu, Cb,           &
                            ee,ne,se,nw,sw,sse, skipTcell1d,  skipUcell1d,     &
                            str1, str2, str3, str4, str5, str6, str7, str8,    &
                            HTE1d,HTN1d, HTE1dm1,HTN1dm1
  implicit none
  integer (kind=int_kind) :: i, nthreads, myscale, scalefactor
  real    (kind=dbl_kind), parameter :: dt=300._dbl_kind
#ifdef itt
  call itt_pause()
#endif
  !- Initialize openmp ---------------------------------------------------------
  call domp_init(nthreads)
  !- Initialize timers ---------------------------------------------------------
  call timer_init(.true.,nthreads,ncl=10)
  call timer(1,'benchp')
  !--- allocate and fill content into arrays -----------------------------------
  call read_nml()
  call set_evp_parameters(dt)
  call alloc_1d_v2(scalefactor)
!$OMP PARALLEL DEFAULT(shared)
  call numainit_all(1,na*scalefactor,navel*scalefactor) ! numainit_all is for default testscale=1
!$OMP END PARALLEL
  call readin_1d_v2()
  call stat_1d()
  !--- allow simple scaling of the testcase: 2xtestscale -----------------------
  myscale = testscale
  do
    if (myscale <= 1) exit
    write(*,*) 'Bumping the size of the testcase, bumping stat:'
    call replicate()
    call stat_1d()
    myscale = myscale-1
  enddo
  write(*,*) 'BEFORE CALLING --'
  call stat_1d()
  write(*,*) 'BEFORE CALLING --'
  call timer(2,'benchp')
  !--- run kernel --------------------------------------------------------------
  call timer(1,'bench')
#ifdef sde
  call fortran_sde_start()
#endif
#ifdef itt
  call itt_resume()
#endif
!$OMP PARALLEL PRIVATE(i)
  do i = 1, ndte
!DIR$ FORCEINLINE
    call stress (ee, ne, se, 1, na,                                            &
                 uvel, vvel, dxT, dyT, skipTcell1d, strength, HTE1d, HTN1d,    &
                 HTE1dm1,    HTN1dm1,                                          &
                 stressp_1,  stressp_2,  stressp_3,  stressp_4,                &
                 stressm_1,  stressm_2,  stressm_3,  stressm_4,                &
                 stress12_1, stress12_2, stress12_3, stress12_4,               &
                 str1, str2, str3, str4, str5, str6, str7, str8)
!$OMP BARRIER
!DIR$ FORCEINLINE
    call stepu (1, na, cdn_ocn, aiu, uocn, vocn, waterxU, wateryU,             &
                forcexU, forceyU, umassdti, fmU, uarear, strintxU, strintyU,   &
                uvel_init, vvel_init, uvel, vvel,                              &
                str1, str2, str3, str4, str5, str6, str7, str8,                &
                nw, sw, sse, skipUcell1d, Tbu, Cb)
!$OMP BARRIER
  enddo
!$OMP END PARALLEL
#ifdef itt
  call itt_pause()
#endif
#ifdef sde
  call fortran_sde_stop()
#endif
  call timer(2,'bench')
  write(*,*) 'AFTER CALLING --'
  call stat_out_1d()
  write(*,*) 'AFTER CALLING --'
  if ((binoutput) .and. (testscale == 1)) then
    call writeout_1d('b')
  endif
  write(*,*) 'Kernel dim (v2b - SPMD openMP): ', na, navel, ndte
  call timer_print('Benchmark took preparation: ','benchp',1)
  call timer_print('Benchmark took: ','bench',1)

  if (lindividual) then
    call timer(1,'stress')
!$OMP PARALLEL PRIVATE(i)
    do i = 1, ndte
      call stress (ee, ne, se, 1, na,                                          &
                   uvel, vvel, dxT, dyT, skipTcell1d, strength, HTE1d, HTN1d,  &
                   HTE1dm1,    HTN1dm1,                                        &
                   stressp_1,  stressp_2,  stressp_3,  stressp_4,              &
                   stressm_1,  stressm_2,  stressm_3,  stressm_4,              &
                   stress12_1, stress12_2, stress12_3, stress12_4,             &
                   str1, str2, str3, str4, str5, str6, str7, str8)
    enddo
!$OMP END PARALLEL
    call timer(2,'stress')
    call timer_print('Benchmark took stress: ','stress',1)
    call timer(1,'stepu')
!$OMP PARALLEL PRIVATE(i)
    do i = 1, ndte
      call stepu (1, na, cdn_ocn, aiu, uocn, vocn, waterxU, wateryU,           &
                  forcexU, forceyU, umassdti, fmU, uarear, strintxU, strintyU, &
                  uvel_init, vvel_init, uvel, vvel,                            &
                  str1, str2, str3, str4, str5, str6, str7, str8,              &
                  nw, sw, sse, skipUcell1d, Tbu, Cb)
    enddo
!$OMP END PARALLEL
    call timer(2,'stepu')
    call timer_print('Benchmark took stepu: ','stepu',1)
  endif
end program ciceevp
