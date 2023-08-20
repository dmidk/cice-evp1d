!===============================================================================
! Copyright (C) 2023, Intel Corporation
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!===============================================================================
program ciceevp
#ifdef itt
  use ittnotify
#endif
  use ice_kinds_mod, only : int_kind, dbl_kind
  use ice_dyn_shared,only : ndte, set_evp_parameters
#ifdef _OPENMP_TARGET
  use ice_constants, only : arlx1i, denom1, capping, deltaminEVP,e_factor,     &
                            epp2i, brlx, c0
#else
  use ice_constants, only : c0
#endif
  use create_nml,    only : nx_block, ny_block, inpfname, read_nml,            &
                            testscale, na, navel, lindividual, binoutput
  use my_timer,      only : timer, timer_init, timer_print
  use myomp,         only : domp_init
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
  integer (kind=int_kind) :: i, nthreads, myscale, scalefactor, iw
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
  do iw = 1,na*scalefactor
    skipTcell1d(iw)=.false.
    skipUcell1d(iw)=.false.
    ee(iw)=0
    ne(iw)=0
    nw(iw)=0
    se(iw)=0
    sw(iw)=0
    sse(iw)=0
    aiu(iw)=c0
    Cb(iw)=c0
    cdn_ocn(iw)=c0
    dxt(iw)=c0
    dyt(iw)=c0
    fmU(iw)=c0
    forcexU(iw)=c0
    forceyU(iw)=c0
    HTE1d(iw)=c0
    HTE1dm1(iw)=c0
    HTN1d(iw)=c0
    HTN1dm1(iw)=c0
    str1(iw)=c0
    str2(iw)=c0
    str3(iw)=c0
    str4(iw)=c0
    str5(iw)=c0
    str6(iw)=c0
    str7(iw)=c0
    str8(iw)=c0
    strength(iw)= c0
    stress12_1(iw)=c0
    stress12_2(iw)=c0
    stress12_3(iw)=c0
    stress12_4(iw)=c0
    stressm_1(iw)=c0
    stressm_2(iw)=c0
    stressm_3(iw)=c0
    stressm_4(iw)=c0
    stressp_1(iw)=c0
    stressp_2(iw)=c0
    stressp_3(iw)=c0
    stressp_4(iw)=c0
    strintxU(iw)= c0
    strintyU(iw)= c0
    Tbu(iw)=c0
    uarear(iw)=c0
    umassdti(iw)=c0
    uocn(iw)=c0
    uvel_init(iw)=c0
    uvel(iw)=c0
    vocn(iw)=c0
    vvel_init(iw)=c0
    vvel(iw)=c0
    waterxU(iw)=c0
    wateryU(iw)=c0
  enddo
  do iw = 1,navel*scalefactor
    uvel(iw)=c0
    vvel(iw)=c0
    str1(iw)=c0
    str2(iw)=c0
    str3(iw)=c0
    str4(iw)=c0
    str5(iw)=c0
    str6(iw)=c0
    str7(iw)=c0
    str8(iw)=c0
  enddo
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
#ifdef itt
  call itt_resume()
#endif
  do i = 1, ndte
    call stress (ee, ne, se, 1, na,                                            &
                 uvel, vvel, dxT, dyT, skipTcell1d, strength, HTE1d, HTN1d,    &
                 HTE1dm1,    HTN1dm1,                                          &
                 stressp_1,  stressp_2,  stressp_3,  stressp_4,                &
                 stressm_1,  stressm_2,  stressm_3,  stressm_4,                &
                 stress12_1, stress12_2, stress12_3, stress12_4,               &
                 str1, str2, str3, str4, str5, str6, str7, str8)
    call stepu (1, na, cdn_ocn, aiu, uocn, vocn, waterxU, wateryU,             &
                forcexU, forceyU, umassdti, fmU, uarear, strintxU, strintyU,   &
                uvel_init, vvel_init, uvel, vvel,                              &
                str1, str2, str3, str4, str5, str6, str7, str8,                &
                nw, sw, sse, skipUcell1d, Tbu, Cb)
  enddo
#ifdef itt
  call itt_pause()
#endif
  call timer(2,'bench')
  write(*,*) 'AFTER CALLING --'
  call stat_out_1d()
  write(*,*) 'AFTER CALLING --'
  if ((binoutput) .and. (testscale == 1)) then
    call writeout_1d('e')
  endif
  write(*,*) 'Kernel dim (v2e - Fortran concurrent loop): ', na, navel, ndte
  call timer_print('Benchmark took preparation: ','benchp',1)
  call timer_print('Benchmark took: ','bench',1)

  if (lindividual) then
    call timer(1,'stress')
    do i = 1, ndte
      call stress (ee, ne, se, 1, na,                                          &
                   uvel, vvel, dxT, dyT, skipTcell1d, strength, HTE1d, HTN1d,  &
                   HTE1dm1,    HTN1dm1,                                        &
                   stressp_1,  stressp_2,  stressp_3,  stressp_4,              &
                   stressm_1,  stressm_2,  stressm_3,  stressm_4,              &
                   stress12_1, stress12_2, stress12_3, stress12_4,             &
                   str1, str2, str3, str4, str5, str6, str7, str8)
    enddo
    call timer(2,'stress')
    call timer_print('Benchmark took stress: ','stress',1)
    call timer(1,'stepu')
    do i = 1, ndte
      call stepu (1, na, cdn_ocn, aiu, uocn, vocn, waterxU, wateryU,           &
                  forcexU, forceyU, umassdti, fmU, uarear, strintxU, strintyU, &
                  uvel_init, vvel_init, uvel, vvel,                            &
                  str1, str2, str3, str4, str5, str6, str7, str8,              &
                  nw, sw, sse, skipUcell1d, Tbu, Cb)
    enddo
    call timer(2,'stepu')
    call timer_print('Benchmark took stepu: ','stepu',1)
  endif
end program ciceevp
