!===============================================================================
! Copyright (C) 2023, Intel Corporation
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!===============================================================================
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
program cice_evp_kernel
#ifdef itt
  use ittnotify
#endif
#ifdef sde
  use ssc_fortran,   only : fortran_sde_start, fortran_sde_stop
#endif
  use ice_kinds_mod, only : int_kind, dbl_kind
  use ice_boundary,  only : readin_halo, ice_haloUpdate
  use ice_dyn_shared,only : ndte, set_evp_parameters
  use create_nml,    only : nx_block, ny_block, inpfname, read_nml,            & 
                            max_blocks, testscale, lindividual, binoutput
  use my_timer,      only : timer, timer_init, timer_print
  use bench,         only : stress, stepu
  use vars,          only : read_nml, alloc_2d_v0, stat_2d, readin_2d,         &
                            writeout_2d, stat_out_2d, stat_2d_resized,         &
                            resize_2d,                                         &
                            icellt, icellu, indxti, indxtj, indxui, indxuj,    &
                            uvel, vvel, dxT, dyT, dxhy, dyhx, cxp, cyp, cxm,   &
                            cym, DminTarea, strength, stressp_1, stressp_2,    & 
                            stressp_3, stressp_4, stressm_1, stressm_2,        & 
                            stressm_3, stressm_4, stress12_1, stress12_2,      & 
                            stress12_3, stress12_4, strtmp, cdn_ocn,           &
                            aiX, uocn, vocn, waterx, watery,                   & 
                            forcex, forcey, umassdti, fm, uarear, strintx,     &
                            strinty, taubx, tauby, uvel_init, vvel_init, Tbu
  implicit none
  integer (int_kind)     :: iblk, nblocks, ksub, myscale
    real(kind=dbl_kind), parameter :: dt =300._dbl_kind !!!! ADDED TAR MAY NEED TO DO SOMETHING WITH IT!!!! 
#ifdef itt
  call itt_pause()
#endif
  !- Initialize timers ---------------------------------------------------------
  call timer_init(.true.,1,ncl=10)
  call timer(1,'benchp')
  !--- allocate and fill content into arrays -----------------------------------
  call read_nml()
  call set_evp_parameters(dt)
  nblocks=max_blocks
  call alloc_2d_v0()
  call readin_2d()
  if (max_blocks > 1) then
    write(*,*) 'WARNING Running with max_blocks > 1 is not well-tested'
    write(*,*) 'WARNING All binary input files must be dumped from a CICE run with the same number of max_blocks'
    call readin_halo()
  endif
  !--- allow simple scaling of the testcase: 2x2xtestscale ---------------------
  myscale = testscale
  do while (myscale > 1)
    write(*,*) 'Bumping the size of the testcase'
    call resize_2d()
    call stat_2d()
    myscale = myscale-1
  enddo
  write(*,*) 'BEFORE CALLING --'
  call stat_2d()
  if (testscale > 1) then
    call stat_2d_resized()
  endif
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
  do ksub = 1,ndte
!$OMP PARALLEL DO PRIVATE(iblk)
    do iblk = 1, nblocks
      call stress (nx_block*testscale  , ny_block*testscale,                   &
                   icellt(iblk),                                               &
                   indxti      (:,iblk), indxtj      (:,iblk),                 &
                   uvel      (:,:,iblk), vvel      (:,:,iblk),                 &
                   dxT       (:,:,iblk), dyT       (:,:,iblk),                 &
                   dxhy      (:,:,iblk), dyhx      (:,:,iblk),                 &
                   cxp       (:,:,iblk), cyp       (:,:,iblk),                 &
                   cxm       (:,:,iblk), cym       (:,:,iblk),                 &
                                         DminTarea (:,:,iblk),                 &
                   strength  (:,:,iblk),                                       &
                   stressp_1 (:,:,iblk), stressp_2 (:,:,iblk),                 &
                   stressp_3 (:,:,iblk), stressp_4 (:,:,iblk),                 &
                   stressm_1 (:,:,iblk), stressm_2 (:,:,iblk),                 &
                   stressm_3 (:,:,iblk), stressm_4 (:,:,iblk),                 &
                   stress12_1(:,:,iblk), stress12_2(:,:,iblk),                 &
                   stress12_3(:,:,iblk), stress12_4(:,:,iblk),                 &
                   strtmp    (:,:,:) )

      call stepu  (nx_block*testscale,  ny_block*testscale,                    &
                   icellu       (iblk), cdn_ocn (:,:,iblk),                    &
                   indxui     (:,iblk), indxuj    (:,iblk),                    &
                   aiX      (:,:,iblk), strtmp  (:,:,:),                       &
                   uocn     (:,:,iblk), vocn    (:,:,iblk),                    &
                   waterx   (:,:,iblk), watery  (:,:,iblk),                    &
                   forcex   (:,:,iblk), forcey  (:,:,iblk),                    &
                   umassdti (:,:,iblk), fm      (:,:,iblk),                    &
                   uarear   (:,:,iblk),                                        &
                   strintx  (:,:,iblk), strinty (:,:,iblk),                    &
                   taubx    (:,:,iblk), tauby   (:,:,iblk),                    &
                   uvel_init(:,:,iblk), vvel_init(:,:,iblk),                   &
                   uvel     (:,:,iblk), vvel    (:,:,iblk),                    &
                   Tbu      (:,:,iblk))
    enddo
!$OMP END PARALLEL DO
    ! U fields at NE corner - calls ice_haloUpdate, controls bundles and masks
    ! FIXME skip the measurement of THE EXPENSIVE halo update at scale for now
    ! call dyn_haloUpdate (halo_info,halo_info_mask, field_loc_NEcorner,       &
    !                      field_type_vector, uvel, vvel)
    call ice_haloUpdate(uvel)
    call ice_haloUpdate(vvel)
  enddo
#ifdef itt
  call itt_pause()
#endif
#ifdef sde
  call fortran_sde_stop()
#endif
  call timer(2,'bench')
  write(*,*) 'AFTER CALLING --'
  call stat_out_2d()
  write(*,*) 'AFTER CALLING --'
  if ((binoutput) .and. (testscale == 1)) then
    call writeout_2d()
  endif
  if (max_blocks > 1) then
    write(*,*) 'Kernel dim (v0-omp): ', sum(icellt(:)), sum(icellu(:)), ndte, max_blocks
  else
    write(*,*) 'Kernel dim (v0): ', sum(icellt(:)), sum(icellu(:)), ndte, max_blocks
  endif
  call timer_print('Benchmark took preparation: ','benchp',1)
  call timer_print('Benchmark took: ','bench',1)
  if (lindividual) then
    call timer(1,'stress')
    do ksub = 1,ndte
!$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
        call stress (nx_block*testscale  , ny_block*testscale,                 &
                     icellt(iblk),                                             &
                     indxti      (:,iblk), indxtj      (:,iblk),               &
                     uvel      (:,:,iblk), vvel      (:,:,iblk),               &
                     dxT       (:,:,iblk), dyT       (:,:,iblk),               &
                     dxhy      (:,:,iblk), dyhx      (:,:,iblk),               &
                     cxp       (:,:,iblk), cyp       (:,:,iblk),               &
                     cxm       (:,:,iblk), cym       (:,:,iblk),               &
                                           DminTarea (:,:,iblk),               &
                     strength  (:,:,iblk),                                     &
                     stressp_1 (:,:,iblk), stressp_2 (:,:,iblk),               &
                     stressp_3 (:,:,iblk), stressp_4 (:,:,iblk),               &
                     stressm_1 (:,:,iblk), stressm_2 (:,:,iblk),               &
                     stressm_3 (:,:,iblk), stressm_4 (:,:,iblk),               &
                     stress12_1(:,:,iblk), stress12_2(:,:,iblk),               &
                     stress12_3(:,:,iblk), stress12_4(:,:,iblk),               &
                     strtmp    (:,:,:) )
      enddo
!$OMP END PARALLEL DO
    enddo
    call timer(2,'stress')
    call timer_print('Benchmark took stress: ','stress',1)
    call timer(1,'stepu')
    do ksub = 1,ndte
!$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
        call stepu  (nx_block*testscale,  ny_block*testscale,                  &
                     icellu       (iblk), cdn_ocn (:,:,iblk),                  &
                     indxui     (:,iblk), indxuj    (:,iblk),                  &
                     aiX      (:,:,iblk), strtmp  (:,:,:),                     &
                     uocn     (:,:,iblk), vocn    (:,:,iblk),                  &
                     waterx   (:,:,iblk), watery  (:,:,iblk),                  &
                     forcex   (:,:,iblk), forcey  (:,:,iblk),                  &
                     umassdti (:,:,iblk), fm      (:,:,iblk),                  &
                     uarear   (:,:,iblk),                                      &
                     strintx  (:,:,iblk), strinty (:,:,iblk),                  &
                     taubx    (:,:,iblk), tauby   (:,:,iblk),                  &
                     uvel_init(:,:,iblk), vvel_init(:,:,iblk),                 &
                     uvel     (:,:,iblk), vvel    (:,:,iblk),                  &
                     Tbu      (:,:,iblk))
      enddo
!$OMP END PARALLEL DO
    enddo
    call timer(2,'stepu')
    call timer_print('Benchmark took stepu: ','stepu',1)
  endif

end program cice_evp_kernel

