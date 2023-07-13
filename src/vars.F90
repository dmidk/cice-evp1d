!===============================================================================
! Copyright (C) 2023, Intel Corporation
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!===============================================================================
! WARNING: This file is one large piece of survival code to get the kernels setup, ugly as .... but focus was elsewhere
module create_nml
  use ice_kinds_mod
  implicit none
  public
  integer (4) :: nx_block, ny_block, max_blocks, testscale, na, navel, nb, nghost
  logical     :: lconvert21d, lindividual, binoutput
  character (200) :: inpfname
  public :: read_nml
contains
  subroutine read_nml
  use ice_constants, only: ndte
  implicit none
  integer (4) :: nml_err
  integer (4), parameter :: funi=503
  namelist /kernel_nml/ inpfname, ndte, nx_block, ny_block, testscale, na,     &
                        navel, nb, lconvert21d, lindividual, binoutput, max_blocks, nghost
  !default values
  ndte=1
  inpfname='pathtoinputfile'
  nb=-1
  na=-1
  navel=-1
  nx_block=-1
  ny_block=-1
  testscale=1
#ifdef v0
  lconvert21d=.true.
#else
  lconvert21d=.false.
#endif
  lindividual=.false.
  binoutput=.true.
  max_blocks=1 ! ==1, well-tested >1 not so much
  nghost=1 ! Size of ghost zone
  open (funi, file='kernel.nml', status='old',iostat=nml_err)
  if (nml_err .ne. 0) then
    write(*,*) ' Read namelist from kernel.nml'
  endif
  nml_err=9
  do while (nml_err > 0)
    read(funi,nml=kernel_nml,iostat=nml_err)
    if (nml_err > 0) then
        write(*,*) 'ERROR: Can not read namelist: kernel_nml from file kernel.nml'
        stop
    endif
  enddo
  write (*,*)'Namelist kernel_nml from kernel.nml:'
  write (*,*)'Number of evp subcylces    = ',ndte
  write (*,*)'Block size x dim           = ',nx_block
  write (*,*)'Block size y dim           = ',ny_block
  write (*,*)'scale_of_testcase          = ',testscale
  write (*,*)'max_blocks                 = ',max_blocks
  write (*,*)'na                         = ',na
  write (*,*)'navel                      = ',navel
  write (*,*)'nb                         = ',nb
  write (*,*)'lcovert21d                 = ',lconvert21d
  write (*,*)'lindividual                = ',lindividual
  write (*,*)'binoutput                  = ',binoutput
  write (*,*)'nghost                     = ',nghost
  end subroutine read_nml
end module create_nml
!===============================================================================
module ice_dyn_shared
  ! hardcoded hack to allow running with multiple block and OpenMP
  implicit none
  public
  integer(4) :: numLocalBlocks, numLocalCopies
  integer (4), dimension(:), allocatable :: ilo, ihi, jlo, jhi
  integer (4), dimension(:,:), allocatable :: srcLocalAddr, dstLocalAddr
  contains
  subroutine readin_halo()
    implicit none
    integer(4) :: ios, lun 
    character (100)        :: binfile
    binfile = 'base_output_integer_blockindices_ref_v0_2d.bin'
    lun=10
    open(lun,file=binfile, form='unformatted', access='stream', action='read', iostat=ios)
    read(lun,iostat=ios)  numLocalBlocks
    allocate(ilo(numLocalBlocks), ihi(numLocalBlocks), jlo(numLocalBlocks), jhi(numLocalBlocks))
    read(lun,iostat=ios) ilo,jlo,ihi,jhi
    close(lun)
    write(6,*) numLocalBlocks, minval(ilo),maxval(ilo),minval(ihi),maxval(ihi),minval(jlo),maxval(jlo),minval(jhi),maxval(jhi)
    binfile = 'base_output_integer_halo_ref_v0_2d.bin'
    open(lun,file=binfile, form='unformatted', access='stream', action='read', iostat=ios)
    read(lun,iostat=ios)  numLocalCopies
    allocate(srcLocalAddr(3,numLocalCopies), dstLocalAddr(3,numLocalCopies))
    read(lun,iostat=ios)  srcLocalAddr, dstLocalAddr
    close(lun)
    write(6,*) numLocalCopies
  end subroutine readin_halo
  subroutine ice_haloUpdate(array)
    use create_nml,    only : nx_block, ny_block, max_blocks
    real (8), dimension(:,:,:), intent(inout) :: array ! array containing field for which halo needs to be updated
    integer (4) ::           &
      i,j,nmsg,              &! dummy loop indices
      iblk,                  &! block sizes for fill
      iSrc,jSrc,             &! source addresses for message
      iDst,jDst,             &! dest   addresses for message
      srcBlock,              &! local block number for source
      dstBlock               ! local block number for destination
    real (8) :: fill         ! value to use for unknown points

    if (max_blocks == 1) return

    fill = 0.0_8
    do iblk = 1, numLocalBlocks
      do j = 1,1
        array(1:nx_block, jlo(iblk)-j,iblk) = fill
        array(1:nx_block, jhi(iblk)+j,iblk) = fill
      enddo
      do i = 1,1
        array(ilo(iblk)-i, 1:ny_block,iblk) = fill
        array(ihi(iblk)+i, 1:ny_block,iblk) = fill
      enddo
    enddo
    do nmsg=1,numLocalCopies
      iSrc     = srcLocalAddr(1,nmsg)
      jSrc     = srcLocalAddr(2,nmsg)
      srcBlock = srcLocalAddr(3,nmsg)
      iDst     = dstLocalAddr(1,nmsg)
      jDst     = dstLocalAddr(2,nmsg)
      dstBlock = dstLocalAddr(3,nmsg)
      if (srcBlock > 0) then
        if (dstBlock > 0) then
          array(iDst,jDst,dstBlock) = &
          array(iSrc,jSrc,srcBlock)
        else if (dstBlock < 0) then ! tripole copy into buffer
           stop 'should not end here dstBlock<0'
        endif
      else if (srcBlock == 0) then
        array(iDst,jDst,dstBlock) = fill
      endif
    enddo
  end subroutine ice_HaloUpdate 
end module ice_dyn_shared
!===============================================================================
module vars
  use ice_kinds_mod
  use ice_constants, only: ndte,  brlx, arlx1i, denom1, capping, e_factor, epp2i, deltaminEVP
  use create_nml, only : nx_block, ny_block, inpfname, read_nml, max_blocks, na, navel, testscale
  implicit none
  public
  integer (kind=int_kind), public :: nb
#ifdef v0
  ! 1D arrays
  real (kind=dbl_kind), dimension(:), allocatable ::                           &
     Vstressp_1, Vstressp_2, Vstressp_3, Vstressp_4 ,                          &
     Vstressm_1, Vstressm_2, Vstressm_3, Vstressm_4 ,                          &
     Vstress12_1,Vstress12_2,Vstress12_3,Vstress12_4,                          &
     str1, str2, str3, str4, str5, str6, str7, str8
  ! 2D arrays
  integer (kind=int_kind), dimension (:,:), allocatable ::                     &
     indxti, indxtj, indxui, indxuj 
  real (kind=dbl_kind), dimension(:,:,:), allocatable ::                       &
     cdn_ocn,aiX,uocn,vocn,waterx,watery,forcex,forcey,umassdti,fm,uarear,     &
     strintx,strinty,uvel_init,vvel_init,                                      &
     strength, uvel, vvel, dxt, dyt, dxhy, dyhx, cyp, cxp, cym, cxm, DminTarea,&
     stressp_1, stressp_2, stressp_3, stressp_4,                               &
     stressm_1, stressm_2, stressm_3, stressm_4,                               &
     stress12_1,stress12_2,stress12_3,stress12_4,                              &
     strtmp, taubx, tauby, Tbu
#endif
#ifdef v1
  logical(kind=log_kind), allocatable, dimension(:) :: skipTcell1d,skipUcell1d
  integer(kind=int_kind), allocatable, dimension(:) :: ee,ne,se,nw,sw,sse
  integer (kind=int_kind), dimension (:), allocatable ::                     &
     indxti, indxtj, indxui, indxuj 
  real (kind=dbl_kind), dimension(:), allocatable ::                           &
     cdn_ocn,aiX,uocn,vocn,waterx,watery,forcex,forcey,umassdti,fm,uarear,     &
     strintx,strinty,uvel_init,vvel_init,                                      &
     strength, uvel, vvel, dxt, dyt, dxhy, dyhx, cyp, cxp, cym, cxm, DminTarea,&
     stressp_1, stressp_2, stressp_3, stressp_4, stressm_1, stressm_2,         &
     stressm_3, stressm_4, stress12_1, stress12_2, stress12_3, stress12_4,     &
     str1, str2, str3, str4, str5, str6, str7, str8, taubx, tauby, Tbu
#endif
#ifdef v2
  real (kind=dbl_kind),   dimension(:), allocatable ::                         &
    HTE1d,HTN1d, HTE1dm1,HTN1dm1
  logical(kind=log_kind), allocatable, dimension(:) :: skipTcell1d,skipUcell1d
  integer(kind=int_kind), allocatable, dimension(:) :: ee,ne,se,nw,sw,sse
  integer (kind=int_kind), dimension (:), allocatable ::                       &
     indxti, indxtj, indxui, indxuj 
  real (kind=dbl_kind), dimension(:), allocatable ::                           &
     cdn_ocn,aiu,uocn,vocn,waterxU,wateryU,forcexU,forceyU,umassdti,fmU,uarear,&
     strintxU,strintyU,uvel_init,vvel_init,                                    &
     strength, uvel, vvel, dxt, dyt,                                           &
     stressp_1, stressp_2, stressp_3, stressp_4, stressm_1, stressm_2,         &
     stressm_3, stressm_4, stress12_1, stress12_2, stress12_3, stress12_4,     &
     str1, str2, str3, str4, str5, str6, str7, str8, Tbu, Cb
#endif
  integer (kind=int_kind) :: lun
  integer (int_kind), dimension(:), allocatable :: icellt, icellu
  real (kind=dbl_kind), dimension(:), allocatable ::                           &
     Vstrintx,Vstrinty,Vuvel,Vvvel
  contains
  integer(kind=int_kind) function io_new_unit ()
    implicit none
    integer(kind=int_kind) :: i
    logical    :: file_is_open
    integer(kind=int_kind), parameter :: io_min=7, io_max=130
    file_is_open = .true.
    i = io_min-1
    do while (file_is_open)
      i = i+1
      if (i > io_max) then
        stop 'ERROR: Could not find new I/O unit number'
      endif
      inquire(unit=i,opened=file_is_open)
    enddo
    io_new_unit = i
  end function io_new_unit
#ifdef v0
  subroutine alloc_2d_v0()
    implicit none
    integer(kind=int_kind) :: ierr
    allocate( icellt(max_blocks) )
    allocate( icellu(max_blocks) )
    allocate( indxti      (nx_block*ny_block,max_blocks), &
              indxtj      (nx_block*ny_block,max_blocks), &
              indxui      (nx_block*ny_block,max_blocks), &
              indxuj      (nx_block*ny_block,max_blocks), &
              strtmp      (nx_block,ny_block,8),          &
              uvel        (nx_block,ny_block,max_blocks), &
              vvel        (nx_block,ny_block,max_blocks), &
              dxt         (nx_block,ny_block,max_blocks), &
              dyt         (nx_block,ny_block,max_blocks), &
              dxhy        (nx_block,ny_block,max_blocks), &
              dyhx        (nx_block,ny_block,max_blocks), &
              cxp         (nx_block,ny_block,max_blocks), &
              cyp         (nx_block,ny_block,max_blocks), &
              cxm         (nx_block,ny_block,max_blocks), &
              cym         (nx_block,ny_block,max_blocks), &
              DminTarea   (nx_block,ny_block,max_blocks), &
              strength    (nx_block,ny_block,max_blocks), &
              stressp_1   (nx_block,ny_block,max_blocks), &
              stressp_2   (nx_block,ny_block,max_blocks), &
              stressp_3   (nx_block,ny_block,max_blocks), &
              stressp_4   (nx_block,ny_block,max_blocks), &
              stressm_1   (nx_block,ny_block,max_blocks), &
              stressm_2   (nx_block,ny_block,max_blocks), &
              stressm_3   (nx_block,ny_block,max_blocks), &
              stressm_4   (nx_block,ny_block,max_blocks), &
              stress12_1  (nx_block,ny_block,max_blocks), &
              stress12_2  (nx_block,ny_block,max_blocks), &
              stress12_3  (nx_block,ny_block,max_blocks), &
              stress12_4  (nx_block,ny_block,max_blocks), &
              stat=ierr )
    if (ierr/=0) stop 'Error allocating stress'
    allocate( cdn_ocn     (nx_block,ny_block,max_blocks), &
              aiX         (nx_block,ny_block,max_blocks), &
              uocn        (nx_block,ny_block,max_blocks), &
              vocn        (nx_block,ny_block,max_blocks), &
              waterx      (nx_block,ny_block,max_blocks), &
              watery      (nx_block,ny_block,max_blocks), &
              forcex      (nx_block,ny_block,max_blocks), &
              forcey      (nx_block,ny_block,max_blocks), &
              umassdti    (nx_block,ny_block,max_blocks), &
              fm          (nx_block,ny_block,max_blocks), &
              uarear      (nx_block,ny_block,max_blocks), &
              strintx     (nx_block,ny_block,max_blocks), &
              strinty     (nx_block,ny_block,max_blocks), &
              taubx       (nx_block,ny_block,max_blocks), &
              tauby       (nx_block,ny_block,max_blocks), &
              Tbu         (nx_block,ny_block,max_blocks), &
              uvel_init   (nx_block,ny_block,max_blocks), & 
              vvel_init   (nx_block,ny_block,max_blocks), &
              stat=ierr  )
    if (ierr/=0) stop 'Error allocating stepu'
  end subroutine alloc_2d_v0
  subroutine alloc_1d_v0()
    use create_nml, only : na, nb, navel
    implicit none
    integer(kind=int_kind)                              :: ierr
    na=icellt(1)
    nb=icellu(1)

    allocate(                                                                  &
       Vstressp_1(1:na), Vstressp_2(1:na), Vstressp_3(1:na), Vstressp_4(1:na), &
       Vstressm_1(1:na), Vstressm_2(1:na), Vstressm_3(1:na), Vstressm_4(1:na), &
       Vstress12_1(1:na),Vstress12_2(1:na),Vstress12_3(1:na),Vstress12_4(1:na),&
       Vstrintx(1:nb),Vstrinty(1:nb),Vuvel(1:nb), Vvvel(1:nb),                 &
       str1(1:navel), str2(1:navel), str3(1:navel), str4(1:navel)             ,&
       str5(1:navel), str6(1:navel), str7(1:navel), str8(1:navel), stat=ierr   )
    if (ierr/=0) stop 'Error allocating 1D all'
  end subroutine alloc_1d_v0
  subroutine readin_2d()
    implicit none
    integer(kind=int_kind) :: ios
    character (100)        :: binfile
    binfile = 'input_double_2d.bin'
    open(lun,file=binfile, form='unformatted', access='stream', action='read', iostat=ios)
    read(lun,iostat=ios) strength, uvel, vvel, dxt, dyt,                                 &
                         dxhy, dyhx, cxp, cyp, cxm, cym, DminTarea,                      &
                         uarear,cdn_ocn,aiX, uocn, vocn, waterx, watery, forcex, forcey, &
                         umassdti, fm, strintx, strinty, Tbu,                            &         
                         stressp_1, stressp_2, stressp_3, stressp_4,                     &
                         stressm_1, stressm_2, stressm_3, stressm_4,                     &
                         stress12_1, stress12_2, stress12_3, stress12_4,                 &
                         capping, e_factor, epp2i
    close(lun)
    uvel_init=uvel
    vvel_init=vvel
    binfile = 'input_integer_2d.bin'
    open(lun,file=binfile, form='unformatted', access='stream', action='read', iostat=ios)
    read(lun,iostat=ios) indxti, indxtj,indxui, indxuj, icellt, icellu, navel                                           
    close(lun)
  end subroutine readin_2d
  subroutine writeout_2d()
    use create_nml, only : lconvert21d
    implicit none
    integer(kind=int_kind) :: ios
    character (100)        :: binfile
    binfile = 'output_stress_v0_2d.bin'
    open(lun,file=binfile, form='unformatted', access='stream', action='write', &
         status='replace', iostat=ios)
    write(lun,iostat=ios)                           &
       stressp_1, stressp_2, stressp_3, stressp_4 , &
       stressm_1, stressm_2, stressm_3, stressm_4 , &
       stress12_1,stress12_2,stress12_3,stress12_4
    close(lun)
    binfile = 'output_stepu_v0_2d.bin'
    open(lun,file=binfile, form='unformatted', access='stream', action='write', &
         status='replace', iostat=ios)
    write(lun,iostat=ios) strintx, strinty, uvel, vvel
    close(lun)
    if ( lconvert21d ) then
      call convert_2d_1d()
      call stat_1d()
      binfile = 'output_stress_v0_1d.bin'
      open(lun,file=binfile, form='unformatted', access='stream', action='write', &
           status='replace', iostat=ios)
      write(lun,iostat=ios)                                                      &
            Vstressp_1, Vstressp_2, Vstressp_3, Vstressp_4,                      &
            Vstressm_1, Vstressm_2, Vstressm_3, Vstressm_4,                      &
            Vstress12_1,Vstress12_2,Vstress12_3,Vstress12_4
      close(lun)  
      binfile = 'output_stepu_v0_1d.bin'
      open(lun,file=binfile, form='unformatted', access='stream', action='write', &
           status='replace', iostat=ios)
      write(lun,iostat=ios)                                                      &
           Vstrintx,Vstrinty,Vuvel,Vvvel
      close(lun)
    endif
  end subroutine writeout_2d
  subroutine stat_1d()
    implicit none
    write (*,*) '-----------------------------------------------'
    write (*,*) 'Statistics for output relevant 1D variables    '
    write (*,*) '-----------------------------------------------'
    write(*,*) 'Size of vectors are :',  size(Vstressp_1(:))
    write(*,'(a11,3f25.8)') ' stressp_1:',   minval(Vstressp_1(:)), maxval(Vstressp_1(:)),sum(Vstressp_1(:))
    write(*,'(a11,3f25.8)') ' stressp_2:',   minval(Vstressp_2(:)), maxval(Vstressp_2(:)),sum(Vstressp_2(:))
    write(*,'(a11,3f25.8)') ' stressp_3:',   minval(Vstressp_3(:)), maxval(Vstressp_3(:)),sum(Vstressp_3(:))
    write(*,'(a11,3f25.8)') ' stressp_4:',   minval(Vstressp_4(:)), maxval(Vstressp_4(:)),sum(Vstressp_4(:))
    write(*,'(a11,3f25.8)') ' stressm_1:',   minval(Vstressm_1(:)), maxval(Vstressm_1(:)),sum(Vstressm_1(:))
    write(*,'(a11,3f25.8)') ' stressm_2:',   minval(Vstressm_2(:)), maxval(Vstressm_2(:)),sum(Vstressm_2(:))
    write(*,'(a11,3f25.8)') ' stressm_3:',   minval(Vstressm_3(:)), maxval(Vstressm_3(:)),sum(Vstressm_3(:))
    write(*,'(a11,3f25.8)') ' stressp_4:',   minval(Vstressm_4(:)), maxval(Vstressm_4(:)),sum(Vstressm_4(:))
    write(*,'(a11,3f25.8)') ' stress12_1:',   minval(Vstress12_1(:)), maxval(Vstress12_1(:)),sum(Vstress12_1(:))
    write(*,'(a11,3f25.8)') ' stress12_2:',   minval(Vstress12_2(:)), maxval(Vstress12_2(:)),sum(Vstress12_2(:))
    write(*,'(a11,3f25.8)') ' stress12_3:',   minval(Vstress12_3(:)), maxval(Vstress12_3(:)),sum(Vstress12_3(:))
    write(*,'(a11,3f25.8)') ' stress12_4:',   minval(Vstress12_4(:)), maxval(Vstress12_4(:)),sum(Vstress12_4(:))
    write(*,'(a11,3f25.8)') ' strintx:',   minval(Vstrintx(:)), maxval(Vstrintx(:)),sum(Vstrintx(:))
    write(*,'(a11,3f25.8)') ' strinty:',   minval(Vstrinty(:)), maxval(Vstrinty(:)),sum(Vstrinty(:))
    write(*,'(a11,3f25.8)') ' uvel:',   minval(Vuvel(:)), maxval(Vuvel(:)),sum(Vuvel(:))
    write(*,'(a11,3f25.8)') ' vvel:',   minval(Vvvel(:)), maxval(Vvvel(:)),sum(Vvvel(:))
  end subroutine stat_1d
  subroutine stat_2d()
    implicit none
    write (*,*) '-----------------------------------------------'
    write (*,*) 'Statistics for stress relevant 2D variables'
    write (*,*) '-----------------------------------------------'
    write(*,*) 'Size and shape of real arrays are :',  size(stressp_1), shape(stressp_1)
    write(*,'(a11,3i25)')   '    indxti:',   minval(indxti(:,1)), maxval(indxti(:,1)), sum(indxti(:,1))
    write(*,'(a11,3i25)')   '    indxtj:',   minval(indxtj(:,1)), maxval(indxtj(:,1)), sum(indxtj(:,1))
    write(*,'(a11,3f25.8)') '      uvel:',   minval(uvel(1:nx_block,1:ny_block,1)), &
                                             maxval(uvel(1:nx_block,1:ny_block,1)), &
                                             sum   (uvel(1:nx_block,1:ny_block,1))
    write(*,'(a11,3f25.8)') '      vvel:',   minval(vvel(:,:,1)), maxval(vvel(:,:,1)), sum(vvel(:,:,1))
    write(*,'(a11,3f25.8)') '       dxt:',   minval(dxt(:,:,1)), maxval(dxt(:,:,1)),sum(dxt(:,:,1))
    write(*,'(a11,3f25.8)') '       dyt:',   minval(dyt(:,:,1)), maxval(dyt(:,:,1)),sum(dyt(:,:,1))
    write(*,'(a11,3f25.8)') '      dxhy:',   minval(dxhy(:,:,1)), maxval(dxhy(:,:,1)),sum(dxhy(:,:,1))
    write(*,'(a11,3f25.8)') '      dyhx:',   minval(dyhx(:,:,1)), maxval(dyhx(:,:,1)),sum(dyhx(:,:,1))
    write(*,'(a11,3f25.8)') '       cxp:',   minval(cxp(:,:,1)), maxval(cxp(:,:,1)),sum(cxp(:,:,1))
    write(*,'(a11,3f25.8)') '       cyp:',   minval(cyp(:,:,1)), maxval(cyp(:,:,1)),sum(cyp(:,:,1))
    write(*,'(a11,3f25.8)') '       cxm:',   minval(cxm(:,:,1)), maxval(cxm(:,:,1)),sum(cxm(:,:,1))
    write(*,'(a11,3f25.8)') '       cym:',   minval(cym(:,:,1)), maxval(cym(:,:,1)),sum(cym(:,:,1))
    write(*,'(a11,3f25.8)') 'DminTarea:',   minval(DminTarea(:,:,1)), maxval(DminTarea(:,:,1)),sum(DminTarea(:,:,1))
    write(*,'(a11,3f25.8)') 'strength:',   minval(strength(:,:,1)), maxval(strength(:,:,1)),sum(strength(:,:,1))
    write(*,'(a11,3f25.8)') ' stressp_1:',   minval(stressp_1(:,:,1)), maxval(stressp_1(:,:,1)),sum(stressp_1(:,:,1))
    write(*,'(a11,3f25.8)') ' stressp_2:',   minval(stressp_2(:,:,1)), maxval(stressp_2(:,:,1)),sum(stressp_2(:,:,1))
    write(*,'(a11,3f25.8)') ' stressp_3:',   minval(stressp_3(:,:,1)), maxval(stressp_3(:,:,1)),sum(stressp_3(:,:,1))
    write(*,'(a11,3f25.8)') ' stressp_4:',   minval(stressp_4(:,:,1)), maxval(stressp_4(:,:,1)),sum(stressp_4(:,:,1))
    write(*,'(a11,3f25.8)') ' stressm_1:',   minval(stressm_1(:,:,1)), maxval(stressm_1(:,:,1)),sum(stressm_1(:,:,1))
    write(*,'(a11,3f25.8)') ' stressm_2:',   minval(stressm_2(:,:,1)), maxval(stressm_2(:,:,1)),sum(stressm_2(:,:,1))
    write(*,'(a11,3f25.8)') ' stressm_3:',   minval(stressm_3(:,:,1)), maxval(stressm_3(:,:,1)),sum(stressm_3(:,:,1))
    write(*,'(a11,3f25.8)') ' stressm_4:',   minval(stressm_4(:,:,1)), maxval(stressm_4(:,:,1)),sum(stressm_4(:,:,1))
    write(*,'(a11,3f25.8)') 'stress12_1:',   minval(stress12_1(:,:,1)), maxval(stress12_1(:,:,1)),sum(stress12_1(:,:,1))
    write(*,'(a11,3f25.8)') 'stress12_2:',   minval(stress12_2(:,:,1)), maxval(stress12_2(:,:,1)),sum(stress12_2(:,:,1))
    write(*,'(a11,3f25.8)') 'stress12_3:',   minval(stress12_3(:,:,1)), maxval(stress12_3(:,:,1)),sum(stress12_3(:,:,1))
    write(*,'(a11,3f25.8)') 'stress12_4:',   minval(stress12_4(:,:,1)), maxval(stress12_4(:,:,1)),sum(stress12_4(:,:,1))
    write (*,*) '-----------------------------------------------'
    write (*,*) 'Statistics for stepu relevant 2D variables'
    write (*,*) '-----------------------------------------------'
    write(*,'(a11,3i25)') '      indxui:',   minval(indxui(:,1)), maxval(indxui(:,1)),sum(indxui(:,1))
    write(*,'(a11,3i25)') '      indxuj:',   minval(indxuj(:,1)), maxval(indxuj(:,1)),sum(indxuj(:,1))
    write(*,'(a11,3f25.8)') ' uvel_init:',minval( uvel_init(:,:,1)),maxval( uvel_init(:,:,1)) ,sum( uvel_init(:,:,1))
    write(*,'(a11,3f25.8)') ' vvel_init:',minval( vvel_init(:,:,1)),maxval( vvel_init(:,:,1)) ,sum( vvel_init(:,:,1))
    write(*,'(a11,3f25.8)') '   strintx:',   minval(strintx(:,:,1)), maxval(strintx(:,:,1)),sum(strintx(:,:,1))
    write(*,'(a11,3f25.8)') '   strinty:',   minval(strinty(:,:,1)), maxval(strinty(:,:,1)),sum(strinty(:,:,1))
    write(*,'(a11,3f25.8)') '   cdn_ocn:',   minval( cdn_ocn(:,:,1)),maxval( cdn_ocn(:,:,1)),sum( cdn_ocn(:,:,1))  
    write(*,'(a11,3f25.8)') '       aiX:',      minval( aiX(:,:,1)),      maxval( aiX(:,:,1)),sum( aiX(:,:,1))
    write(*,'(a11,3f25.8)') '      uocn:',     minval( uocn(:,:,1)),     maxval( uocn(:,:,1)),sum( uocn(:,:,1))
    write(*,'(a11,3f25.8)') '      vocn:',     minval( vocn(:,:,1)),     maxval( vocn(:,:,1)),sum( vocn(:,:,1))
    write(*,'(a11,3f25.8)') '    waterx:',   minval( waterx(:,:,1)),   maxval( waterx(:,:,1)),sum( waterx(:,:,1))
    write(*,'(a11,3f25.8)') '    watery:',   minval( watery(:,:,1)),   maxval( watery(:,:,1)),sum( watery(:,:,1))
    write(*,'(a11,3f25.8)') '    forcex:',   minval( forcex(:,:,1)),   maxval( forcex(:,:,1)),sum( forcex(:,:,1))
    write(*,'(a11,3f25.8)') '    forcey:',   minval( forcey(:,:,1)),   maxval( forcey(:,:,1)),sum( forcey(:,:,1))
    write(*,'(a11,3f25.8)') '  umassdti:', minval( umassdti(:,:,1)),maxval( umassdti(:,:,1)),sum( umassdti(:,:,1))
    write(*,'(a11,3f25.8)') '        fm:',       minval( fm(:,:,1)),       maxval( fm(:,:,1)),sum( fm(:,:,1))
    write(*,'(a11,3f25.8)') '    uarear:',   minval( uarear(:,:,1)),maxval( uarear(:,:,1)),sum( uarear(:,:,1))
  end subroutine stat_2d
  subroutine stat_2d_resized()
    ! FIXME verify that the next block is identical to the first - INCOMPLETE
    write(*,'(a11,3i25)')   '    1indxti:', minval(indxti(1:nx_block*ny_block,1)), &
                                            maxval(indxti(1:nx_block*ny_block,1)), &
                                            sum(indxti(1:nx_block*ny_block,1))
    write(*,'(a11,3i25)')   '    2indxti:', minval(indxti(nx_block*ny_block+1:2*nx_block*ny_block,1)), &
                                            maxval(indxti(nx_block*ny_block+1:2*nx_block*ny_block,1)), &
                                            sum(indxti(nx_block*ny_block+1:2*nx_block*ny_block,1))
!    write(*,'(a11,3i25)')   '    indxtj:',   minval(indxtj(:,1)), maxval(indxtj(:,1)), sum(indxtj(:,1))
    write(*,'(a11,3f25.8)') '   1uvel:',   minval(uvel(1:nx_block,1:ny_block,1)), &
                                           maxval(uvel(1:nx_block,1:ny_block,1)), &
                                           sum   (uvel(1:nx_block,1:ny_block,1))
    write(*,'(a11,3f25.8)') '   2uvel:',   minval(uvel(nx_block+1:2*nx_block, 1:ny_block,1)), &
                                           maxval(uvel(nx_block+1:2*nx_block, 1:ny_block,1)), &
                                           sum   (uvel(nx_block+1:2*nx_block, 1:ny_block,1))
    write(*,'(a11,3f25.8)') '   3uvel:',   minval(uvel(1:nx_block, 1+ny_block:2*ny_block,1)), &
                                           maxval(uvel(1:nx_block, 1+ny_block:2*ny_block,1)), &
                                           sum   (uvel(1:nx_block, 1+ny_block:2*ny_block,1))
    write(*,'(a11,3f25.8)') '   4uvel:',   minval(uvel(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1)), &
                                           maxval(uvel(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1)), &
                                           sum   (uvel(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1))
    write(*,'(a11,3f25.8)') '   1vvel:',   minval(vvel(1:nx_block,1:ny_block,1)), &
                                           maxval(vvel(1:nx_block,1:ny_block,1)), &
                                           sum   (vvel(1:nx_block,1:ny_block,1))
    write(*,'(a11,3f25.8)') '   2vvel:',   minval(vvel(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1)), &
                                           maxval(vvel(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1)), &
                                           sum   (vvel(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1))
    write(*,'(a11,3f25.8)') '   1dxt:',    minval(dxt(1:nx_block,1:ny_block,1)), &
                                           maxval(dxt(1:nx_block,1:ny_block,1)), &
                                           sum(dxt(1:nx_block,1:ny_block,1))
    write(*,'(a11,3f25.8)') '   2dxt:',    minval(dxt(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1)), &
                                           maxval(dxt(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1)), &
                                           sum   (dxt(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1))
    write(*,'(a11,3f25.8)') '   1dyt:',    minval(dyt(1:nx_block,1:ny_block,1)), &
                                           maxval(dyt(1:nx_block,1:ny_block,1)), &
                                           sum(dyt(1:nx_block,1:ny_block,1))
    write(*,'(a11,3f25.8)') '   2dyt:',   minval(dyt(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1)), &
                                           maxval(dyt(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1)), &
                                           sum   (dyt(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1))
    write(*,'(a11,3f25.8)') '   1dxhy:',   minval(dxhy(1:nx_block,1:ny_block,1)), &
                                           maxval(dxhy(1:nx_block,1:ny_block,1)), &
                                           sum(dxhy(1:nx_block,1:ny_block,1))
    write(*,'(a11,3f25.8)') '   2dxhy:',   minval(dxhy(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1)), &
                                           maxval(dxhy(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1)), &
                                           sum   (dxhy(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1))
    write(*,'(a11,3f25.8)') '   1dyhx:',   minval(dyhx(1:nx_block,1:ny_block,1)), &
                                           maxval(dyhx(1:nx_block,1:ny_block,1)), &
                                           sum(dyhx(1:nx_block,1:ny_block,1))
    write(*,'(a11,3f25.8)') '   2dyhx:',   minval(dyhx(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1)), &
                                           maxval(dyhx(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1)), &
                                           sum   (dyhx(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1))
    write(*,'(a11,3f25.8)') '   1cxm:',   minval(cxm(1:nx_block,1:ny_block,1)), &
                                           maxval(cxm(1:nx_block,1:ny_block,1)), &
                                           sum(cxm(1:nx_block,1:ny_block,1))
    write(*,'(a11,3f25.8)') '   2cxm:',   minval(cxm(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1)), &
                                           maxval(cxm(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1)), &
                                           sum   (cxm(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1))
    write(*,'(a11,3f25.8)') '   1cym:',   minval(cym(1:nx_block,1:ny_block,1)), &
                                           maxval(cym(1:nx_block,1:ny_block,1)), &
                                           sum(cym(1:nx_block,1:ny_block,1))
    write(*,'(a11,3f25.8)') '   2cym:',   minval(cym(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1)), &
                                           maxval(cym(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1)), &
                                           sum   (cym(nx_block+1:2*nx_block, 1+ny_block:2*ny_block,1))
  end subroutine stat_2d_resized
  subroutine stat_out_2d()
    implicit none
    write (*,*) '-----------------------------------------------'
    write (*,*) 'Statistics for output relevant 2D variables    '
    write (*,*) '-----------------------------------------------'
    write(*,'(a11,3f25.8)') ' stressp_1:',   minval(stressp_1(:,:,1)), maxval(stressp_1(:,:,1)),sum(stressp_1(:,:,1))
    write(*,'(a11,3f25.8)') ' stressp_2:',   minval(stressp_2(:,:,1)), maxval(stressp_2(:,:,1)),sum(stressp_2(:,:,1))
    write(*,'(a11,3f25.8)') ' stressp_3:',   minval(stressp_3(:,:,1)), maxval(stressp_3(:,:,1)),sum(stressp_3(:,:,1))
    write(*,'(a11,3f25.8)') ' stressp_4:',   minval(stressp_4(:,:,1)), maxval(stressp_4(:,:,1)),sum(stressp_4(:,:,1))
    write(*,'(a11,3f25.8)') ' stressm_1:',   minval(stressm_1(:,:,1)), maxval(stressm_1(:,:,1)),sum(stressm_1(:,:,1))
    write(*,'(a11,3f25.8)') ' stressm_2:',   minval(stressm_2(:,:,1)), maxval(stressm_2(:,:,1)),sum(stressm_2(:,:,1))
    write(*,'(a11,3f25.8)') ' stressm_3:',   minval(stressm_3(:,:,1)), maxval(stressm_3(:,:,1)),sum(stressm_3(:,:,1))
    write(*,'(a11,3f25.8)') ' stressm_4:',   minval(stressm_4(:,:,1)), maxval(stressm_4(:,:,1)),sum(stressm_4(:,:,1))
    write(*,'(a11,3f25.8)') 'stress12_1:',    minval(stress12_1(:,:,1)), maxval(stress12_1(:,:,1)),sum(stress12_1(:,:,1))
    write(*,'(a11,3f25.8)') 'stress12_2:',    minval(stress12_2(:,:,1)), maxval(stress12_2(:,:,1)),sum(stress12_2(:,:,1))
    write(*,'(a11,3f25.8)') 'stress12_3:',    minval(stress12_3(:,:,1)), maxval(stress12_3(:,:,1)),sum(stress12_3(:,:,1))
    write(*,'(a11,3f25.8)') 'stress12_4:',    minval(stress12_4(:,:,1)), maxval(stress12_4(:,:,1)),sum(stress12_4(:,:,1))
    write(*,'(a11,3f25.8)') '   strintx:',   minval(strintx(:,:,1)), maxval(strintx(:,:,1)),sum(strintx(:,:,1))
    write(*,'(a11,3f25.8)') '   strinty:',   minval(strinty(:,:,1)), maxval(strinty(:,:,1)),sum(strinty(:,:,1))
    write(*,'(a11,3f25.8)') '      uvel:',   minval(uvel(:,:,1)), maxval(uvel(:,:,1)),sum(uvel(:,:,1))
    write(*,'(a11,3f25.8)') '      vvel:',   minval(vvel(:,:,1)), maxval(vvel(:,:,1)),sum(vvel(:,:,1))
  end subroutine stat_out_2d
  subroutine dump_all_2d()
    implicit none
    integer(kind=int_kind) :: ios
    character (100)        :: binfile
    binfile = 'output_stress_2d.bin'
    open(lun,file=binfile, form='unformatted', access='stream', action='write', status='replace', iostat=ios)
    write(lun,iostat=ios)                           &
       stressp_1, stressp_2, stressp_3, stressp_4 , &
       stressm_1, stressm_2, stressm_3, stressm_4 , &
       stress12_1,stress12_2,stress12_3,stress12_4
    close(lun)
    binfile = 'output_stepu_2d.bin'
    open(lun,file=binfile, form='unformatted', access='stream', action='write', status='replace', iostat=ios)
    write(lun,iostat=ios) strintx, strinty, uvel, vvel 
    close(lun)
  end subroutine dump_all_2d
  subroutine convert_2d_1d()
    use create_nml, only : na, nb, navel
    implicit none
    real (kind=dbl_kind), parameter :: c0   = 0.0_dbl_kind
    integer(kind=int_kind) :: iblk, i, j, ij
    call alloc_1d_v0()
    str1=c0
    str2=c0
    str3=c0
    str4=c0
    str5=c0
    str6=c0
    str7=c0
    str8=c0
    ! matrix2vector
    iblk=1
    do ij=1,na
      i=indxti(ij,iblk)
      j=indxtj(ij,iblk)
      Vstressp_1 (ij)=stressp_1 (i,j,iblk)
      Vstressp_2 (ij)=stressp_2 (i,j,iblk)
      Vstressp_3 (ij)=stressp_3 (i,j,iblk)
      Vstressp_4 (ij)=stressp_4 (i,j,iblk)
      Vstressm_1 (ij)=stressm_1 (i,j,iblk)
      Vstressm_2 (ij)=stressm_2 (i,j,iblk)
      Vstressm_3 (ij)=stressm_3 (i,j,iblk)
      Vstressm_4 (ij)=stressm_4 (i,j,iblk)
      Vstress12_1(ij)=stress12_1(i,j,iblk)
      Vstress12_2(ij)=stress12_2(i,j,iblk)
      Vstress12_3(ij)=stress12_3(i,j,iblk)
      Vstress12_4(ij)=stress12_4(i,j,iblk)
      str1       (ij)=strtmp    (i,j,1)
      str2       (ij)=strtmp    (i,j,2)
      str3       (ij)=strtmp    (i,j,3)
      str4       (ij)=strtmp    (i,j,4)
      str5       (ij)=strtmp    (i,j,5)
      str6       (ij)=strtmp    (i,j,6)
      str7       (ij)=strtmp    (i,j,7)
      str8       (ij)=strtmp    (i,j,8)
    enddo
    iblk=1
    do ij=1,nb
      i=indxui(ij,iblk)
      j=indxuj(ij,iblk)
      Vstrintx   (ij)=strintx   (i,j,iblk)
      Vstrinty   (ij)=strinty   (i,j,iblk)
      Vuvel      (ij)=uvel      (i,j,iblk)
      Vvvel      (ij)=vvel      (i,j,iblk)
    enddo
    write(*,*) 'Kernel 1D parameters na, navel, nb:', na, navel, nb
  end subroutine convert_2d_1d
  
  subroutine dump_all_1d()
    implicit none
    integer(kind=int_kind) :: ios
    character (100)        :: binfile
    binfile = 'output_stress_1d.bin'
    open(lun,file=binfile, form='unformatted', access='stream', action='write', status='replace', iostat=ios)
    write(lun,iostat=ios)                                                          &
      Vstressp_1, Vstressp_2, Vstressp_3, Vstressp_4,                             &
      Vstressm_1, Vstressm_2, Vstressm_3, Vstressm_4,                             &
      Vstress12_1,Vstress12_2,Vstress12_3,Vstress12_4
    close(lun)  
    binfile = 'output_stepu_1d.bin'
    open(lun,file=binfile, form='unformatted', access='stream', action='write', status='replace', iostat=ios)
    write(lun,iostat=ios)                                                          &
      Vstrintx,Vstrinty,Vuvel,Vvvel
    close(lun)
  end subroutine dump_all_1d
  subroutine resize_2d()
    use create_nml, only : nghost
    real(kind=dbl_kind),    dimension(:,:,:), allocatable :: tmpa
    integer(kind=int_kind), dimension(:,:),   allocatable :: tmpi, tmpj
    integer(kind=int_kind)                                :: ierr, len1, len2, len3, dim1, dim2, dim3
    integer(8)                                :: leni1, lenj1, ij, iblk, i, j, ii, jj
    write(*,*) 'Resize active points before:',  sum(icellt(:)), sum(icellu(:))
    len1 = size(strtmp,1)
    len2 = size(strtmp,2)
    len3 = size(strtmp,3)
    dim1 = 2*len1-2*nghost
    dim2 = 2*len2-2*nghost
    dim3 = len3

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = strtmp(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=strtmp)

!    stop 'firt double handled'
    len1 = size(stressp_1,1)
    len2 = size(stressp_1,2)
    len3 = size(stressp_1,3)
    dim1 = 2*len1-2*nghost
    dim2 = 2*len2-2*nghost
    dim3 = len3
    print *, len1, len2, len3
    print *, dim1, dim2, dim3
    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = stressp_1(ii,jj,:)
    enddo
    enddo
!    write(*,*) size(stressp_1)
!    write(*,*) shape(stressp_1)
!    write(*,*) size(tmpa)
    call move_alloc(from=tmpa,to=stressp_1)
!  stop 'second double handled'
!    write(*,*) size(stressp_1)
!    write(*,*) shape(stressp_1)
    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = stressp_2(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=stressp_2)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = stressp_3(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=stressp_3)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = stressp_4(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=stressp_4)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = stressm_1(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=stressm_1)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = stressm_2(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=stressm_2)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = stressm_3(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=stressm_3)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = stressm_4(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=stressm_4)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = stress12_1(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=stress12_1)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = stress12_2(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=stress12_2)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = stress12_3(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=stress12_3)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = stress12_4(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=stress12_4)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = uvel(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=uvel)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = vvel(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=vvel)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = dxT(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=dxT)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = dyT(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=dyT)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = dxhy(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=dxhy)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = dyhx(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=dyhx)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = cxp(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=cxp)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = cyp(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=cyp)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = cxm(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=cxm)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = cym(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=cym)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = DminTarea(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=DminTarea)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = strength(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=strength)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = umassdti(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=umassdti)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = uvel_init(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=uvel_init)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = vvel_init(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=vvel_init)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = strintx(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=strintx)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = strinty(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=strinty)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = forcex(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=forcex)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = forcey(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=forcey)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = waterx(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=waterx)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = watery(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=watery)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = aiX(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=aiX)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = uocn(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=uocn)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = vocn(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=vocn)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = uarear(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=uarear)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = taubx(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=taubx)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = tauby(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=tauby)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = fm(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=fm)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = cdn_ocn(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=cdn_ocn)

    allocate(tmpa(dim1,dim2,dim3), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpa'
    do j=1,dim2
    do i=1,dim1
      ii=min(i/2+1, len1)
      jj=min(j/2+1, len2)
      tmpa(i, j, :)  = Tbu(ii,jj,:)
    enddo
    enddo
    call move_alloc(from=tmpa,to=Tbu)

!    stop 'All doubles handled'

    write(*,*) 'Resizing before indxti :', minval(indxti(:,:)), maxval(indxti(:,:)), size(indxti(:,:))
    leni1 = size(indxti,1) ! nx_block*ny_block
    lenj1 = size(indxtj,1) ! nx_block*ny_block
  !  allocate(tmpi(4*leni1,1), stat=ierr)
  !  allocate(tmpj(4*lenj1,1), stat=ierr) ! fixme hardcoded maxblock=1
    allocate(tmpi(4*leni1,max_blocks), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpi'
    allocate(tmpj(4*lenj1,max_blocks), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpj'
    print *, size(tmpi)
    print *, shape(tmpi)
    print *, size(tmpj)
    print *, shape(tmpj)
    print *, 'leni1, lenj1 and 4x and na, 4xna', leni1, lenj1, 4*leni1, 4*lenj1, na, 4*na
    if (4*leni1 < 4*na) then
        stop 'dimensions wrong'
    endif

    ! quadrant1
    tmpi(1:leni1, 1:max_blocks) = indxti 
    tmpj(1:lenj1, 1:max_blocks) = indxtj
    do iblk = 1, max_blocks
      do ij = 1, icellt(iblk)
        ! quadrant2
        tmpi(ij+icellt(iblk),iblk) = indxti(ij,iblk)+nx_block
        tmpj(ij+icellt(iblk),iblk) = indxtj(ij,iblk)
      enddo
      do ij = 1, icellt(iblk)
        ! quadrant3
        tmpi(ij+2*icellt(iblk),iblk) = indxti(ij,iblk)
        tmpj(ij+2*icellt(iblk),iblk) = indxtj(ij,iblk)+ny_block
      enddo
      do ij = 1, icellt(iblk)
        ! quadrant4
        tmpi(ij+3*icellt(iblk),iblk) = indxti(ij,iblk)+nx_block
        tmpj(ij+3*icellt(iblk),iblk) = indxtj(ij,iblk)+ny_block
      enddo
      icellt(iblk) = icellt(iblk)*4
    enddo

    call move_alloc(from=tmpi,to=indxti)
    call move_alloc(from=tmpj,to=indxtj)
    write(*,*) 'Resizing after indxti debug :', leni1, 1+leni1,  1+na,  icellt(1)+leni1, na+leni1
    write(*,*) 'Resizing after indxti :', minval(indxti(:,:)), maxval(indxti(:,:)), size(indxti(:,:))
    write(*,*) 'Resizing after indxtj :', minval(indxtj(:,:)), maxval(indxtj(:,:)), size(indxtj(:,:))

    write(*,*) 'Resizing before indxui :', minval(indxui(:,:)), maxval(indxui(:,:)), size(indxui(:,:))
    write(*,*) 'si before indxui :', minval(indxui(:,:)), maxval(indxui(:,:)), size(indxui(:,:))
    leni1 = size(indxui,1)
    lenj1 = size(indxuj,1)
    allocate(tmpi(4*leni1,max_blocks), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpi'
    allocate(tmpj(4*lenj1,max_blocks), stat=ierr)
    if (ierr/=0) stop 'Error allocating 2D tmpj'
    ! quadrant1
    tmpi(1:leni1, 1:max_blocks)  = indxui
    tmpj(1:lenj1, 1:max_blocks)  = indxuj
    do iblk = 1, max_blocks
      do ij = 1, icellu(iblk)
        ! quadrant2
        tmpi(ij+icellu(iblk),iblk) = indxui(ij,iblk)+nx_block
        tmpj(ij+icellu(iblk),iblk) = indxuj(ij,iblk)
      enddo
      do ij = 1, icellu(iblk)
        ! quadrant3
        tmpi(ij+2*icellu(iblk),iblk) = indxui(ij,iblk)
        tmpj(ij+2*icellu(iblk),iblk) = indxuj(ij,iblk)+ny_block
      enddo
      do ij = 1, icellu(iblk)
        ! quadrant4
        tmpi(ij+3*icellu(iblk),iblk) = indxui(ij,iblk)+nx_block
        tmpj(ij+3*icellu(iblk),iblk) = indxuj(ij,iblk)+ny_block
      enddo
      icellu(iblk) =  icellu(iblk)*4
    enddo
    call move_alloc(from=tmpi,to=indxui)
    call move_alloc(from=tmpj,to=indxuj)
    write(*,*) 'Resizing after indxui :', minval(indxui(:,:)), maxval(indxui(:,:)), size(indxui(:,:))
    write(*,*) 'Resize active points after:',  sum(icellt(:)), sum(icellu(:))
  end subroutine resize_2d
#endif
#ifdef v1
  subroutine alloc_1d_v1()
    implicit none
    integer(kind=int_kind) :: ierr
    allocate(ee(1:na),ne(1:na),se(1:na),nw(1:na),sw(1:na),sse(1:na),stat=ierr)
    if (ierr/=0) stop 'Error allocating indexes'
    allocate(str1(1:navel), str2(1:navel), str3(1:navel), str4(1:navel), &
             str5(1:navel), str6(1:navel), str7(1:navel), str8(1:navel), stat=ierr)
    if (ierr/=0) stop 'Error allocating str' 
    allocate( uvel      (1:navel), &
              vvel      (1:navel), &
              dxt       (1:na   ), &
              dyt       (1:na   ), &
              dxhy      (1:na   ), &
              dyhx      (1:na   ), &
              cxp       (1:na   ), &
              cyp       (1:na   ), &
              cxm       (1:na   ), &
              cym       (1:na   ), &
              DminTarea (1:na   ), &
              strength  (1:na   ), &
              stressp_1 (1:na   ), &
              stressp_2 (1:na   ), &
              stressp_3 (1:na   ), &
              stressp_4 (1:na   ), &
              stressm_1 (1:na   ), &
              stressm_2 (1:na   ), &
              stressm_3 (1:na   ), &
              stressm_4 (1:na   ), &
              stress12_1(1:na   ), &
              stress12_2(1:na   ), &
              stress12_3(1:na   ), &
              stress12_4(1:na   ),stat=ierr)
    if (ierr/=0) stop 'Error allocating stress'
    allocate(cdn_ocn(1:na), aiX   (1:na),    &
           uocn     (1:na), vocn  (1:na),    &
           waterx   (1:na), watery(1:na),    &
           forcex   (1:na), forcey (1:na),   &
           umassdti (1:na), fm(1:na),        &
           uarear   (1:na),                  &
           strintx  (1:na), strinty(1:na),    &   ! fixme nb
           taubx    (1:na), tauby  (1:na),       &  
           Tbu      (1:na),                    &
           uvel_init(1:na),vvel_init(1:na), &
           skiptcell1d(1:na),skipucell1d(1:na), stat=ierr)
    if (ierr/=0) stop 'Error allocating stepu'
  end subroutine alloc_1d_v1

  subroutine readin_1d()
    implicit none
    integer(kind=int_kind) :: ios
    character (100)        :: binfile
    binfile = 'input_double_1d_v1.bin'
    open(lun,file=binfile, form='unformatted', access='stream', action='read', iostat=ios)
    read(lun,iostat=ios) strength, uvel, vvel, dxt, dyt,                &
      dxhy, dyhx, cxp, cyp, cxm, cym, DminTarea,                      &
      uarear,cdn_ocn,aiX, uocn, vocn, waterx, watery, forcex, forcey, &
      umassdti, fm, strintx, strinty, Tbu,                                &         
      stressp_1, stressp_2, stressp_3, stressp_4,                     &
      stressm_1, stressm_2, stressm_3, stressm_4,                     &
      stress12_1, stress12_2, stress12_3, stress12_4,                 &
      capping, e_factor, epp2i
    close(lun)
    uvel_init=uvel
    vvel_init=vvel
    binfile = 'input_integer_1d.bin'
    open(lun,file=binfile, form='unformatted', access='stream', action='read', iostat=ios)
    read(lun,iostat=ios)  ee, ne, se,  nw, sw, sse
    close(lun)
    write(*,'(a11,4I30)') 'Min/max ee', minval(ee), maxval(ee), sum(ee), size(ee)
    write(*,'(a11,4I30)') 'Min/max ne', minval(ne), maxval(ne), sum(ne), size(ne)
    write(*,'(a11,4I30)') 'Min/max se', minval(se), maxval(se), sum(se), size(se)
    write(*,'(a11,4I30)') 'Min/max nw', minval(nw), maxval(nw), sum(nw), size(nw)
    write(*,'(a11,4I30)') 'Min/max sw', minval(sw), maxval(sw), sum(sw), size(sw)
    write(*,'(a11,4I30)') 'Min/max sse', minval(sse), maxval(sse), sum(sse), size(sse)
    binfile = 'input_logical_1d.bin'
    open(lun,file=binfile, form='unformatted', access='stream', action='read', iostat=ios)
    read(lun,iostat=ios) skipUcell1d, skipTcell1d
    close(lun)
  end subroutine readin_1d
  subroutine stat_1d()
    implicit none
    real   (kind=dbl_kind) :: tmp
    write(*,'(a11,4I30)') 'Min/max ee', minval(ee), maxval(ee), sum(ee), size(ee)
    write(*,'(a11,4I30)') 'Min/max ne', minval(ne), maxval(ne), sum(ne), size(ne)
    write(*,'(a11,4I30)') 'Min/max se', minval(se), maxval(se), sum(se), size(se)
    write(*,'(a11,4I30)') 'Min/max nw', minval(nw), maxval(nw), sum(nw), size(nw)
    write(*,'(a11,4I30)') 'Min/max sw', minval(sw), maxval(sw), sum(sw), size(sw)
    write(*,'(a11,4I30)') 'Min/max sse', minval(sse), maxval(sse), sum(sse), size(sse)
    write(*,'(a11,2I30)') 'count skipu ',  count(skipUcell1d(:)), size(skipUcell1d(:))
    write(*,'(a11,2I30)') 'count skipt ',  count(skipTcell1d(:)), size(skipTcell1d(:))
    tmp=count(skipUcell1d(:))/size(uvel)
    write(*,'(a17,1f25.8)') ' Fraction skipu :', tmp
    tmp=count(skipTcell1d(:))/size(uvel)
    write(*,'(a17,1f25.8)') ' Fraction skipt :', tmp
    write (*,*) '-----------------------------------------------'
    write (*,*) 'Statistics for all relevant 1D variables'
    write (*,*) '-----------------------------------------------'
    write(*,'(a11,4f25.8)') '      uvel:',   minval(uvel(:)), maxval(uvel(:)),sum(uvel(:)), size(uvel)
    write(*,'(a11,4f25.8)') '      vvel:',   minval(vvel(:)), maxval(vvel(:)),sum(vvel(:)), size(vvel)
    write(*,'(a11,4f25.8)') '       dxt:',   minval(dxt(:)), maxval(dxt(:)),sum(dxt(:)), size(dxt)
    write(*,'(a11,4f25.8)') '       dyt:',   minval(dyt(:)), maxval(dyt(:)),sum(dyt(:)), size(dyt)
    ! FIXME incomplete FIXME
    write(*,'(a11,3f25.8)') '      dxhy:',   minval(dxhy(:)), maxval(dxhy(:)),sum(dxhy(:))
    write(*,'(a11,3f25.8)') '      dyhx:',   minval(dyhx(:)), maxval(dyhx(:)),sum(dyhx(:))
    write(*,'(a11,3f25.8)') '       cxp:',   minval(cxp(:)), maxval(cxp(:)),sum(cxp(:))
    write(*,'(a11,3f25.8)') '       cyp:',   minval(cyp(:)), maxval(cyp(:)),sum(cyp(:))
    write(*,'(a11,3f25.8)') '       cxm:',   minval(cxm(:)), maxval(cxm(:)),sum(cxm(:))
    write(*,'(a11,3f25.8)') '       cym:',   minval(cym(:)), maxval(cym(:)),sum(cym(:))
    write(*,'(a11,3f25.8)') 'DminTarea:',   minval(DminTarea(:)), maxval(DminTarea(:)),sum(DminTarea(:))
    write(*,'(a11,3f25.8)') 'strength:',   minval(strength(:)), maxval(strength(:)),sum(strength(:))
    write(*,'(a11,3f25.8)') ' stressp_1:',   minval(stressp_1(:)), maxval(stressp_1(:)),sum(stressp_1(:))
    write(*,'(a11,3f25.8)') ' stressp_2:',   minval(stressp_2(:)), maxval(stressp_2(:)),sum(stressp_2(:))
    write(*,'(a11,3f25.8)') ' stressp_3:',   minval(stressp_3(:)), maxval(stressp_3(:)),sum(stressp_3(:))
    write(*,'(a11,3f25.8)') ' stressp_4:',   minval(stressp_4(:)), maxval(stressp_4(:)),sum(stressp_4(:))
    write(*,'(a11,3f25.8)') ' stressm_1:',   minval(stressm_1(:)), maxval(stressm_1(:)),sum(stressm_1(:))
    write(*,'(a11,3f25.8)') ' stressm_2:',   minval(stressm_2(:)), maxval(stressm_2(:)),sum(stressm_2(:))
    write(*,'(a11,3f25.8)') ' stressm_3:',   minval(stressm_3(:)), maxval(stressm_3(:)),sum(stressm_3(:))
    write(*,'(a11,3f25.8)') ' stressm_4:',   minval(stressm_4(:)), maxval(stressm_4(:)),sum(stressm_4(:))
    write(*,'(a11,3f25.8)') 'stress12_1:',   minval(stress12_1(:)), maxval(stress12_1(:)),sum(stress12_1(:))
    write(*,'(a11,3f25.8)') 'stress12_2:',   minval(stress12_2(:)), maxval(stress12_2(:)),sum(stress12_2(:))
    write(*,'(a11,3f25.8)') 'stress12_3:',   minval(stress12_3(:)), maxval(stress12_3(:)),sum(stress12_3(:))
    write(*,'(a11,3f25.8)') 'stress12_4:',   minval(stress12_4(:)), maxval(stress12_4(:)),sum(stress12_4(:))
    write (*,*) '-----------------------------------------------'
    write (*,*) 'Statistics for stepu relevant 1D variables'
    write (*,*) '-----------------------------------------------'
    write(*,'(a11,3f25.8)') ' uvel_init:',minval( uvel_init(:)),maxval( uvel_init(:)) ,sum(uvel_init(:))
    write(*,'(a11,3f25.8)') ' vvel_init:',minval( vvel_init(:)),maxval( vvel_init(:)) ,sum(vvel_init(:))
    write(*,'(a11,3f25.8)') '   strintxU:',   minval(strintx(:)), maxval(strintx(:)),sum(strintx(:))
    write(*,'(a11,3f25.8)') '   strintyU:',   minval(strinty(:)), maxval(strinty(:)),sum(strinty(:))
    write(*,'(a11,3f25.8)') '   cdn_ocn:',   minval( cdn_ocn(:)),maxval( cdn_ocn(:)),sum( cdn_ocn(:))  
    write(*,'(a11,3f25.8)') '       aiu:',      minval( aiX(:)),      maxval( aiX(:)),sum( aiX(:))
    write(*,'(a11,3f25.8)') '      uocn:',     minval( uocn(:)),     maxval( uocn(:)),sum( uocn(:))
    write(*,'(a11,3f25.8)') '      vocn:',     minval( vocn(:)),     maxval( vocn(:)),sum( vocn(:))
    write(*,'(a11,3f25.8)') '    waterxU:',   minval( waterx(:)),   maxval( waterx(:)),sum( waterx(:))
    write(*,'(a11,3f25.8)') '    wateryU:',   minval( watery(:)),   maxval( watery(:)),sum( watery(:))
    write(*,'(a11,3f25.8)') '    forcexU:',   minval( forcex(:)),   maxval( forcex(:)),sum( forcex(:))
    write(*,'(a11,3f25.8)') '    forceyU:',   minval( forcey(:)),   maxval( forcey(:)),sum( forcey(:))
    write(*,'(a11,3f25.8)') '  umassdti:', minval( umassdti(:)),maxval( umassdti(:)),sum( umassdti(:))
    write(*,'(a11,3f25.8)') '        fmU:',       minval( fm(:)),       maxval( fm(:)),sum( fm(:))
    write(*,'(a11,3f25.8)') '    uarear:',   minval( uarear(:)),maxval( uarear(:)),sum( uarear(:))
  end subroutine stat_1d
  subroutine writeout_1d()
    use create_nml, only : na, nb
    implicit none
    integer(kind=int_kind) :: ios, ierr, i, iw, nbb
    character (100)        :: binfile
    integer (kind=int_kind) :: lun1
    integer (kind=int_kind) :: lun2
    write(*,*) 'Saving output in output_stress_v1_1d.bin and output_stepu_v1_1d.bin'
    binfile = 'output_stress_v1_1d.bin'
    lun1 = io_new_unit()
    open(lun1,file=binfile, form='unformatted', access='stream', action='write', &
         status='replace', iostat=ios)
    if (ios .ne. 0) then
      stop ('Failed open output_stress_v1_1d.bin') 
    endif 
    write(lun1,iostat=ios)                                                      &
          stressp_1, stressp_2, stressp_3, stressp_4,                      &
          stressm_1, stressm_2, stressm_3, stressm_4,                      &
          stress12_1,stress12_2,stress12_3,stress12_4
    if (ios .ne. 0) then
      stop ('Failed write output_stress_v1_1d.bin') 
    endif 
    close(lun1)  
    flush(lun1)
!    binfile = 'output_stepu_v1_1d_all.bin'
!    open(lun,file=binfile, form='unformatted', access='stream', action='write', &
!         status='replace', iostat=ios)
!    write(lun,iostat=ios)                                                      &
!         strintx,strinty,uvel,vvel
!    close(lun)
    nbb=na-count(skipucell1d(:))
    write (*,*) nb, nbb
    allocate(Vstrintx(1:nb),Vstrinty(1:nb),Vuvel(1:nb), Vvvel(1:nb), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D Vstrintx,...'
    i=0
    do iw=1,na
      if (skipucell1d(iw)) cycle
      i=i+1
      Vstrintx   (i)=strintx(iw)
      Vstrinty   (i)=strinty(iw)
      Vuvel      (i)=uvel(iw)
      Vvvel      (i)=vvel(iw)
    enddo
    if (i .ne. nbb) stop 'PROBLEM in generating output_stepu_v1_1d.bin'
    binfile = 'output_stepu_v1_1d.bin'
    lun2 = io_new_unit()
    open(lun2, file=binfile, form='unformatted', access='stream',     &
            action='write', status='replace', iostat=ios)
    if (ios .ne. 0) then
      write(*,*) 'ios is ', ios
      stop ('Failed open output_stepu_v1_1d.bin') 
    endif 
    write(lun2,iostat=ios) Vstrintx,Vstrinty,Vuvel,Vvvel
    if (ios .ne. 0) then
      stop ('Failed write output_stepu_v1_1d.bin') 
    endif 
    close(lun2)
    flush(lun2)
    write(*,*) 'Output files closed and flushed, the rest is up to the FS'
  end subroutine writeout_1d
  subroutine stat_out_1d()
    implicit none
    write (*,*) '-----------------------------------------------'
    write (*,*) 'Statistics for output relevant 1D variables    '
    write (*,*) '-----------------------------------------------'
    write(*,'(a11,3f25.8)') ' stressp_1:',   minval(stressp_1(:)), maxval(stressp_1(:)),sum(stressp_1(:))
    write(*,'(a11,3f25.8)') ' stressp_2:',   minval(stressp_2(:)), maxval(stressp_2(:)),sum(stressp_2(:))
    write(*,'(a11,3f25.8)') ' stressp_3:',   minval(stressp_3(:)), maxval(stressp_3(:)),sum(stressp_3(:))
    write(*,'(a11,3f25.8)') ' stressp_4:',   minval(stressp_4(:)), maxval(stressp_4(:)),sum(stressp_4(:))
    write(*,'(a11,3f25.8)') ' stressm_1:',   minval(stressm_1(:)), maxval(stressm_1(:)),sum(stressm_1(:))
    write(*,'(a11,3f25.8)') ' stressm_2:',   minval(stressm_2(:)), maxval(stressm_2(:)),sum(stressm_2(:))
    write(*,'(a11,3f25.8)') ' stressm_3:',   minval(stressm_3(:)), maxval(stressm_3(:)),sum(stressm_3(:))
    write(*,'(a11,3f25.8)') ' stressm_4:',   minval(stressm_4(:)), maxval(stressm_4(:)),sum(stressm_4(:))
    write(*,'(a11,3f25.8)') 'stress12_1:',    minval(stress12_1(:)), maxval(stress12_1(:)),sum(stress12_1(:))
    write(*,'(a11,3f25.8)') 'stress12_2:',    minval(stress12_2(:)), maxval(stress12_2(:)),sum(stress12_2(:))
    write(*,'(a11,3f25.8)') 'stress12_3:',    minval(stress12_3(:)), maxval(stress12_3(:)),sum(stress12_3(:))
    write(*,'(a11,3f25.8)') 'stress12_4:',    minval(stress12_4(:)), maxval(stress12_4(:)),sum(stress12_4(:))
    write(*,'(a11,3f25.8)') '      str1:',   minval(str1(:)), maxval(str1(:)),sum(str1(:))
    write(*,'(a11,3f25.8)') '      str2:',   minval(str2(:)), maxval(str2(:)),sum(str2(:))
    write(*,'(a11,3f25.8)') '      str3:',   minval(str3(:)), maxval(str3(:)),sum(str3(:))
    write(*,'(a11,3f25.8)') '      str4:',   minval(str4(:)), maxval(str4(:)),sum(str4(:))
    write(*,'(a11,3f25.8)') '      str5:',   minval(str5(:)), maxval(str5(:)),sum(str5(:))
    write(*,'(a11,3f25.8)') '      str6:',   minval(str6(:)), maxval(str6(:)),sum(str6(:))
    write(*,'(a11,3f25.8)') '      str7:',   minval(str7(:)), maxval(str7(:)),sum(str7(:))
    write(*,'(a11,3f25.8)') '      str8:',   minval(str8(:)), maxval(str8(:)),sum(str8(:))
    write(*,'(a11,3f25.8)') '   strintx:',   minval(strintx(:)), maxval(strintx(:)),sum(strintx(:))
    write(*,'(a11,3f25.8)') '   strinty:',   minval(strinty(:)), maxval(strinty(:)),sum(strinty(:))
    write(*,'(a11,3f25.8)') '      uvel:',   minval(uvel(:)), maxval(uvel(:)),sum(uvel(:))
    write(*,'(a11,3f25.8)') '      vvel:',   minval(vvel(:)), maxval(vvel(:)),sum(vvel(:))
  end subroutine stat_out_1d
  subroutine resize_1d()
    ! FIXME incomplete FIXME
    real   (kind=dbl_kind), dimension(:), allocatable :: tmpVa
    integer(kind=int_kind), dimension(:), allocatable :: tmpVi
    logical(kind=log_kind), dimension(:), allocatable :: tmpVl
    integer(kind=int_kind) :: ierr, vlen, ilen, llen
    integer(kind=int_kind), parameter :: scalefactor=2 ! 4 to mimic v0
    vlen = size(str1)
    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = str1
    tmpVa(vlen+1:) = str1
    call move_alloc(from=tmpVa,to=str1)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = str2
    tmpVa(vlen+1:) = str2
    call move_alloc(from=tmpVa,to=str2)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = str3
    tmpVa(vlen+1:) = str3
    call move_alloc(from=tmpVa,to=str3)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = str4
    tmpVa(vlen+1:) = str4
    call move_alloc(from=tmpVa,to=str4)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = str5
    tmpVa(vlen+1:) = str5
    call move_alloc(from=tmpVa,to=str5)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = str6
    tmpVa(vlen+1:) = str6
    call move_alloc(from=tmpVa,to=str6)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = str7
    tmpVa(vlen+1:) = str7
    call move_alloc(from=tmpVa,to=str7)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = str8
    tmpVa(vlen+1:) = str8
    call move_alloc(from=tmpVa,to=str8)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = uvel
    tmpVa(vlen+1:) = uvel
    call move_alloc(from=tmpVa,to=uvel)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = vvel
    tmpVa(vlen+1:) = vvel
    call move_alloc(from=tmpVa,to=vvel)


    vlen = size(stressp_1)
    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = stressp_1
    tmpVa(vlen+1:) = stressp_1
    call move_alloc(from=tmpVa,to=stressp_1)
    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = stressp_2
    tmpVa(vlen+1:) = stressp_2
    call move_alloc(from=tmpVa,to=stressp_2)
    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = stressp_3
    tmpVa(vlen+1:) = stressp_3
    call move_alloc(from=tmpVa,to=stressp_3)
    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = stressp_4
    tmpVa(vlen+1:) = stressp_4
    call move_alloc(from=tmpVa,to=stressp_4)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = stressm_1
    tmpVa(vlen+1:) = stressm_1
    call move_alloc(from=tmpVa,to=stressm_1)
    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = stressm_2
    tmpVa(vlen+1:) = stressm_2
    call move_alloc(from=tmpVa,to=stressm_2)
    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = stressm_3
    tmpVa(vlen+1:) = stressm_3
    call move_alloc(from=tmpVa,to=stressm_3)
    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = stressm_4
    tmpVa(vlen+1:) = stressm_4
    call move_alloc(from=tmpVa,to=stressm_4)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = stress12_1
    tmpVa(vlen+1:) = stress12_1
    call move_alloc(from=tmpVa,to=stress12_1)
    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = stress12_2
    tmpVa(vlen+1:) = stress12_2
    call move_alloc(from=tmpVa,to=stress12_2)
    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = stress12_3
    tmpVa(vlen+1:) = stress12_3
    call move_alloc(from=tmpVa,to=stress12_3)
    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = stress12_4
    tmpVa(vlen+1:) = stress12_4
    call move_alloc(from=tmpVa,to=stress12_4)

    !  uvel, vvel, dxT, dyT, dxhy, dyhx, cxp, cyp, cxm,   &
    !  cym, DminTarea, strength, cdn_ocn,   &
    ! aiX, uocn, vocn, waterx, watery,   & 
    !                        forcex, forcey, umassdti, fm, uarear, strintx,     &
    !                        strinty, taubx, tauby, uvel_init, vvel_init, Tbu,  &
    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = dxT
    tmpVa(vlen+1:) = dxT
    call move_alloc(from=tmpVa,to=dxT)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = dyT
    tmpVa(vlen+1:) = dyT
    call move_alloc(from=tmpVa,to=dyT)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = dxhy
    tmpVa(vlen+1:) = dxhy
    call move_alloc(from=tmpVa,to=dxhy)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = dyhx
    tmpVa(vlen+1:) = dyhx
    call move_alloc(from=tmpVa,to=dyhx)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = cxp
    tmpVa(vlen+1:) = cxp
    call move_alloc(from=tmpVa,to=cxp)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = cyp
    tmpVa(vlen+1:) = cyp
    call move_alloc(from=tmpVa,to=cyp)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = cxm
    tmpVa(vlen+1:) = cxm
    call move_alloc(from=tmpVa,to=cxm)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = cym
    tmpVa(vlen+1:) = cym
    call move_alloc(from=tmpVa,to=cym)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = DminTarea
    tmpVa(vlen+1:) = DminTarea
    call move_alloc(from=tmpVa,to=DminTarea)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = strength
    tmpVa(vlen+1:) = strength
    call move_alloc(from=tmpVa,to=strength)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = cdn_ocn
    tmpVa(vlen+1:) = cdn_ocn
    call move_alloc(from=tmpVa,to=cdn_ocn)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = aiX
    tmpVa(vlen+1:) = aiX
    call move_alloc(from=tmpVa,to=aiX)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = uocn
    tmpVa(vlen+1:) = uocn
    call move_alloc(from=tmpVa,to=uocn)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = vocn
    tmpVa(vlen+1:) = vocn
    call move_alloc(from=tmpVa,to=vocn)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = waterx
    tmpVa(vlen+1:) = waterx
    call move_alloc(from=tmpVa,to=waterx)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = watery
    tmpVa(vlen+1:) = watery
    call move_alloc(from=tmpVa,to=watery)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = forcex
    tmpVa(vlen+1:) = forcex
    call move_alloc(from=tmpVa,to=forcex)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = forcey
    tmpVa(vlen+1:) = forcey
    call move_alloc(from=tmpVa,to=forcey)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = umassdti
    tmpVa(vlen+1:) = umassdti
    call move_alloc(from=tmpVa,to=umassdti)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = fm
    tmpVa(vlen+1:) = fm
    call move_alloc(from=tmpVa,to=fm)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = uarear
    tmpVa(vlen+1:) = uarear
    call move_alloc(from=tmpVa,to=uarear)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = strintx
    tmpVa(vlen+1:) = strintx
    call move_alloc(from=tmpVa,to=strintx)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = strinty
    tmpVa(vlen+1:) = strinty
    call move_alloc(from=tmpVa,to=strinty)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = taubx
    tmpVa(vlen+1:) = taubx
    call move_alloc(from=tmpVa,to=taubx)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = tauby
    tmpVa(vlen+1:) = tauby
    call move_alloc(from=tmpVa,to=tauby)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = uvel_init
    tmpVa(vlen+1:) = uvel_init
    call move_alloc(from=tmpVa,to=uvel_init)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = vvel_init
    tmpVa(vlen+1:) = vvel_init
    call move_alloc(from=tmpVa,to=vvel_init)

    allocate( tmpVa(scalefactor*vlen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVa'
    tmpVa(1:vlen)  = Tbu
    tmpVa(vlen+1:) = Tbu
    call move_alloc(from=tmpVa,to=Tbu)

    ilen = size(ee)
    allocate( tmpVi(scalefactor*ilen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVi'
    tmpVi(1:ilen)  = ee
    tmpVi(ilen+1:) = ee
    call move_alloc(from=tmpVi,to=ee)

    ilen = size(ne)
    allocate( tmpVi(scalefactor*ilen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVi'
    tmpVi(1:ilen)  = ne
    tmpVi(ilen+1:) = ne
    call move_alloc(from=tmpVi,to=ne)

    ilen = size(se)
    allocate( tmpVi(scalefactor*ilen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVi'
    tmpVi(1:ilen)  = se
    tmpVi(ilen+1:) = se
    call move_alloc(from=tmpVi,to=se)

    ilen = size(nw)
    allocate( tmpVi(scalefactor*ilen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVi'
    tmpVi(1:ilen)  = nw
    tmpVi(ilen+1:) = nw
    call move_alloc(from=tmpVi,to=nw)

    ilen = size(sw)
    allocate( tmpVi(scalefactor*ilen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVi'
    tmpVi(1:ilen)  = sw
    tmpVi(ilen+1:) = sw
    call move_alloc(from=tmpVi,to=sw)

    ilen = size(sse)
    allocate( tmpVi(scalefactor*ilen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVi'
    tmpVi(1:ilen)  = sse
    tmpVi(ilen+1:) = sse
    call move_alloc(from=tmpVi,to=sse)

    llen = size(skipTcell1d)
    allocate( tmpVl(scalefactor*llen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVi'
    tmpVl(1:llen)  = skipTcell1d
    tmpVl(llen+1:) = skipTcell1d
    call move_alloc(from=tmpVl,to=skipTcell1d)

    llen = size(skipUcell1d)
    allocate( tmpVl(scalefactor*llen), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D tmpVi'
    tmpVl(1:llen)  = skipUcell1d
    tmpVl(llen+1:) = skipUcell1d
    call move_alloc(from=tmpVl,to=skipUcell1d)
    ! FIXME incomplete FIXME

    ! finally bump size of wetpoints
    write(*,*) 'Resizing debug before na nb navel:', na, nb, navel
!    icellu(1) =  icellu(1)*scalefactor
!    icellt(1) =  icellt(1)*scalefactor
    na=scalefactor*na
    navel=scalefactor*navel
    nb=scalefactor*nb
    write(*,*) 'Resizing debug after na nb navel:', na, nb, navel
  end subroutine resize_1d
#endif

#ifdef v2
  subroutine replicate()
    use create_nml, only : na, nb, navel
    real   (kind=dbl_kind), dimension(:), allocatable :: tmpVa
    integer(kind=int_kind), dimension(:), allocatable :: tmpVi
    logical(kind=log_kind), dimension(:), allocatable :: tmpVl
    integer(kind=int_kind) :: ierr, vlen, ilen, llen, nvlen, nilen, nllen
    integer(kind=int_kind), parameter :: scalefactor=2 ! 4 to mimic v0
    vlen  = navel
    nvlen = scalefactor*vlen
    write(*,*) na, navel, vlen, nvlen
    allocate( tmpVa(nvlen), stat=ierr)
    if (ierr/=0) then
      write (*,*) nvlen, ierr
      stop 'Error allocating 1D tmpVa in replicate'
    endif
    tmpVa(1:vlen)  = str1(1:vlen)
    tmpVa(vlen+1:) = str1(1:vlen)
    str1(1:nvlen)  = tmpVa
    tmpVa(1:vlen)  = str2(1:vlen)
    tmpVa(vlen+1:) = str2(1:vlen)
    str2(1:nvlen)  = tmpVa
    tmpVa(1:vlen)  = str3(1:vlen)
    tmpVa(vlen+1:) = str3(1:vlen)
    str3(1:nvlen)  = tmpVa
    tmpVa(1:vlen)  = str4(1:vlen)
    tmpVa(vlen+1:) = str4(1:vlen)
    str4(1:nvlen)  = tmpVa
    tmpVa(1:vlen)  = str5(1:vlen)
    tmpVa(vlen+1:) = str5(1:vlen)
    str5(1:nvlen)  = tmpVa
    tmpVa(1:vlen)  = str6(1:vlen)
    tmpVa(vlen+1:) = str6(1:vlen)
    str6(1:nvlen)  = tmpVa
    tmpVa(1:vlen)  = str7(1:vlen)
    tmpVa(vlen+1:) = str7(1:vlen)
    str7(1:nvlen)  = tmpVa
    tmpVa(1:vlen)  = str8(1:vlen)
    tmpVa(vlen+1:) = str8(1:vlen)
    str8(1:nvlen)  = tmpVa
    tmpVa(1:vlen)  = uvel(1:vlen)
    tmpVa(vlen+1:) = uvel(1:vlen)
    uvel(1:nvlen)  = tmpVa
    tmpVa(1:vlen)  = vvel(1:vlen)
    tmpVa(vlen+1:) = vvel(1:vlen)
    vvel(1:nvlen)  = tmpVa

    vlen = na
    nvlen = scalefactor*vlen
    deallocate(tmpVa)
    allocate( tmpVa(nvlen), stat=ierr)
    if (ierr/=0) then
      write (*,*) nvlen, ierr
      stop 'Error allocating 1D tmpVa in replicate'
    endif

    tmpVa(1:vlen)  = stressp_1(1:vlen)
    tmpVa(vlen+1:) = stressp_1(1:vlen)
    stressp_1(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = stressp_2(1:vlen)
    tmpVa(vlen+1:) = stressp_2(1:vlen)
    stressp_2(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = stressp_3(1:vlen)
    tmpVa(vlen+1:) = stressp_3(1:vlen)
    stressp_3(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = stressp_4(1:vlen)
    tmpVa(vlen+1:) = stressp_4(1:vlen)
    stressp_4(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = stressm_1(1:vlen)
    tmpVa(vlen+1:) = stressm_1(1:vlen)
    stressm_1(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = stressm_2(1:vlen)
    tmpVa(vlen+1:) = stressm_2(1:vlen)
    stressm_2(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = stressm_3(1:vlen)
    tmpVa(vlen+1:) = stressm_3(1:vlen)
    stressm_3(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = stressm_4(1:vlen)
    tmpVa(vlen+1:) = stressm_4(1:vlen)
    stressm_4(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = stress12_1(1:vlen)
    tmpVa(vlen+1:) = stress12_1(1:vlen)
    stress12_1(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = stress12_2(1:vlen)
    tmpVa(vlen+1:) = stress12_2(1:vlen)
    stress12_2(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = stress12_3(1:vlen)
    tmpVa(vlen+1:) = stress12_3(1:vlen)
    stress12_3(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = stress12_4(1:vlen)
    tmpVa(vlen+1:) = stress12_4(1:vlen)
    stress12_4(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = Cb(1:vlen)
    tmpVa(vlen+1:) = Cb(1:vlen)
    Cb(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = dxT(1:vlen)
    tmpVa(vlen+1:) = dxT(1:vlen)
    dxT(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = dyT(1:vlen)
    tmpVa(vlen+1:) = dyT(1:vlen)
    dyT(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = HTE1d(1:vlen)
    tmpVa(vlen+1:) = HTE1d(1:vlen)
    HTE1d(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = HTN1d(1:vlen)
    tmpVa(vlen+1:) = HTN1d(1:vlen)
    HTN1d(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = HTE1dm1(1:vlen)
    tmpVa(vlen+1:) = HTE1dm1(1:vlen)
    HTE1dm1(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = HTN1dm1(1:vlen)
    tmpVa(vlen+1:) = HTN1dm1(1:vlen)
    HTN1dm1(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = strength(1:vlen)
    tmpVa(vlen+1:) = strength(1:vlen)
    strength(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = cdn_ocn(1:vlen)
    tmpVa(vlen+1:) = cdn_ocn(1:vlen)
    cdn_ocn(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = aiu(1:vlen)
    tmpVa(vlen+1:) = aiu(1:vlen)
    aiu(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = uocn(1:vlen)
    tmpVa(vlen+1:) = uocn(1:vlen)
    uocn(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = vocn(1:vlen)
    tmpVa(vlen+1:) = vocn(1:vlen)
    vocn(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = waterxU(1:vlen)
    tmpVa(vlen+1:) = waterxU(1:vlen)
    waterxU(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = wateryU(1:vlen)
    tmpVa(vlen+1:) = wateryU(1:vlen)
    wateryU(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = forcexU(1:vlen)
    tmpVa(vlen+1:) = forcexU(1:vlen)
    forcexU(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = forceyU(1:vlen)
    tmpVa(vlen+1:) = forceyU(1:vlen)
    forceyU(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = umassdti(1:vlen)
    tmpVa(vlen+1:) = umassdti(1:vlen)
    umassdti(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = fmU(1:vlen)
    tmpVa(vlen+1:) = fmU(1:vlen)
    fmU(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = uarear(1:vlen)
    tmpVa(vlen+1:) = uarear(1:vlen)
    uarear(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = strintxU(1:vlen)
    tmpVa(vlen+1:) = strintxU(1:vlen)
    strintxU(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = strintyU(1:vlen)
    tmpVa(vlen+1:) = strintyU(1:vlen)
    strintyU(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = uvel_init(1:vlen)
    tmpVa(vlen+1:) = uvel_init(1:vlen)
    uvel_init(1:nvlen)  = tmpVa


    tmpVa(1:vlen)  = vvel_init(1:vlen)
    tmpVa(vlen+1:) = vvel_init(1:vlen)
    vvel_init(1:nvlen)  = tmpVa

    tmpVa(1:vlen)  = Tbu(1:vlen)
    tmpVa(vlen+1:) = Tbu(1:vlen)
    Tbu(1:nvlen)  = tmpVa
    deallocate(tmpVa)

    ilen = na
    nilen = scalefactor*ilen
    allocate( tmpVi(nilen), stat=ierr)
    if (ierr/=0) then
      write (*,*) nilen, ierr
      stop 'Error allocating 1D tmpVi in replicate'
    endif

    tmpVi(1:ilen)  = ee(1:ilen)
    tmpVi(ilen+1:) = ee(1:ilen)+ilen
    ee(1:nilen)    = tmpVi

    tmpVi(1:ilen)  = ne(1:ilen)
    tmpVi(ilen+1:) = ne(1:ilen)+ilen
    ne(1:nilen)    = tmpVi

    tmpVi(1:ilen)  = se(1:ilen)
    tmpVi(ilen+1:) = se(1:ilen)+ilen
    se(1:nilen)    = tmpVi

    tmpVi(1:ilen)  = nw(1:ilen)
    tmpVi(ilen+1:) = nw(1:ilen)+ilen
    nw(1:nilen)    = tmpVi

    tmpVi(1:ilen)  = sw(1:ilen)
    tmpVi(ilen+1:) = sw(1:ilen)+ilen
    sw(1:nilen)    = tmpVi

    tmpVi(1:ilen)  = sse(1:ilen)
    tmpVi(ilen+1:) = sse(1:ilen)+ilen
    sse(1:nilen)   = tmpVi

    deallocate(tmpVi)

    llen = na
    nllen = scalefactor*llen
    allocate( tmpVl(nllen), stat=ierr)
    if (ierr/=0) then
      write (*,*) nllen, ierr
      stop 'Error allocating 1D tmpVl in replicate'
    endif

    tmpVl(1:llen)  = skipTcell1d(1:llen)
    tmpVl(llen+1:) = skipTcell1d(1:llen)
    skipTcell1d(1:nllen)  = tmpVl

    tmpVl(1:llen)  = skipUcell1d(1:llen)
    tmpVl(llen+1:) = skipUcell1d(1:llen)
    skipUcell1d(1:nllen)  = tmpVl

    deallocate(tmpVl)

    write(*,*) 'Resizing debug before na nb navel:', na, nb, navel
    write(*,*) 'Resizing debug after na nb navel:',  scalefactor*na, scalefactor*nb, scalefactor*navel
    na=scalefactor*na
    navel=scalefactor*navel
    nb=scalefactor*nb
    write(*,*) 'Resizing debug after na nb navel:', na, nb, navel
  end subroutine replicate
  subroutine alloc_1d_v2(scalefactor)
    implicit none
    integer(kind=int_kind), intent(out) :: scalefactor
    integer(kind=int_kind) :: ierr, sna, snavel
    scalefactor = 2**(testscale-1)
    sna         = na*scalefactor
    snavel      = navel*scalefactor
    allocate(ee(1:sna),ne(1:sna),se(1:sna),nw(1:sna),sw(1:sna),sse(1:sna),     &
             stat=ierr)
    if (ierr/=0) stop 'Error allocating indexes'
    allocate(str1(1:snavel), str2(1:snavel), str3(1:snavel), str4(1:snavel),   &
             str5(1:snavel), str6(1:snavel), str7(1:snavel), str8(1:snavel),   &
             uvel(1:snavel), vvel(1:snavel),                                   &
             stat=ierr)
    if (ierr/=0) stop 'Error allocating navel' 
    allocate( HTE1d     (1:sna   ),                                            &
              HTN1d     (1:sna   ),                                            &
              HTE1dm1   (1:sna   ),                                            &
              HTN1dm1   (1:sna   ),                                            &
              dxt       (1:sna   ),                                            &
              dyt       (1:sna   ),                                            &
              strength  (1:sna   ),                                            &
              stressp_1 (1:sna   ),                                            &
              stressp_2 (1:sna   ),                                            &
              stressp_3 (1:sna   ),                                            &
              stressp_4 (1:sna   ),                                            &
              stressm_1 (1:sna   ),                                            &
              stressm_2 (1:sna   ),                                            &
              stressm_3 (1:sna   ),                                            &
              stressm_4 (1:sna   ),                                            &
              stress12_1(1:sna   ),                                            &
              stress12_2(1:sna   ),                                            &
              stress12_3(1:sna   ),                                            &
              stress12_4(1:sna   ),stat=ierr)
    if (ierr/=0) stop 'Error allocating stress'
    allocate(cdn_ocn    (1:sna), aiu        (1:sna),                           &
             uocn       (1:sna), vocn       (1:sna),                           &
             waterxU    (1:sna), wateryU    (1:sna),                           &
             forcexU    (1:sna), forceyU    (1:sna),                           &
             umassdti   (1:sna), fmU        (1:sna),                           &
             uarear     (1:sna),                                               &
             strintxU   (1:sna), strintyU   (1:sna),                           &   ! fixme nb
             Tbu        (1:sna), Cb         (1:sna),                           &
             uvel_init  (1:sna), vvel_init  (1:sna),                           &
             skiptcell1d(1:sna), skipucell1d(1:sna), stat=ierr)
    if (ierr/=0) stop 'Error allocating stepu'
  end subroutine alloc_1d_v2
  subroutine readin_1d_v2()
    implicit none
    integer(kind=int_kind) :: ios
    character (100)        :: binfile
    binfile = 'input_double_1d_v2.bin'
    open(lun,file=binfile, form='unformatted', access='stream', action='read', &
         iostat=ios)
    read(lun,iostat=ios) strength(1:na), uvel(1:navel), vvel(1:navel),         &
         dxt(1:na), dyt(1:na), uarear(1:na), cdn_ocn(1:na), aiu(1:na),         &
         uocn(1:na), vocn(1:na), waterxU(1:na), wateryU(1:na), forcexU(1:na),  &
         forceyU(1:na), umassdti(1:na), fmU(1:na), strintxU(1:na),             &
         strintyU(1:na), Tbu(1:na), stressp_1(1:na), stressp_2(1:na),          &
         stressp_3(1:na), stressp_4(1:na), stressm_1(1:na), stressm_2(1:na),   &
         stressm_3(1:na), stressm_4(1:na), stress12_1(1:na), stress12_2(1:na), &
         stress12_3(1:na), stress12_4(1:na),                                   &
         capping, e_factor, epp2i, deltaminEVP,                                &
         HTE1d(1:na), HTN1d(1:na), HTE1dm1(1:na), HTN1dm1(1:na)
    close(lun)
    uvel_init=uvel
    vvel_init=vvel
    binfile = 'input_integer_1d.bin'
    open(lun,file=binfile, form='unformatted', access='stream', action='read', iostat=ios)
    read(lun,iostat=ios) ee(1:na), ne(1:na), se(1:na),  nw(1:na), sw(1:na), sse(1:na)
    close(lun)
    write(*,'(a11,4I30)') 'Min/max ee', minval(ee(1:na)), maxval(ee(1:na)), sum(ee(1:na)), size(ee)
    write(*,'(a11,4I30)') 'Min/max ne', minval(ne(1:na)), maxval(ne(1:na)), sum(ne(1:na)), size(ne)
    write(*,'(a11,4I30)') 'Min/max se', minval(se(1:na)), maxval(se(1:na)), sum(se(1:na)), size(se)
    write(*,'(a11,4I30)') 'Min/max nw', minval(nw(1:na)), maxval(nw(1:na)), sum(nw(1:na)), size(nw)
    write(*,'(a11,4I30)') 'Min/max sw', minval(sw(1:na)), maxval(sw(1:na)), sum(sw(1:na)), size(sw)
    write(*,'(a11,4I30)') 'Min/max sse', minval(sse(1:na)), maxval(sse(1:na)), sum(sse(1:na)), size(sse)
    binfile = 'input_logical_1d.bin'
    open(lun,file=binfile, form='unformatted', access='stream', action='read', iostat=ios)
    read(lun,iostat=ios) skipUcell1d(1:na), skipTcell1d(1:na)
    close(lun)
  end subroutine readin_1d_v2
  subroutine stat_1d()
    implicit none
    write(*,'(a11,4I30)') 'Min/max ee', minval(ee(1:na)), maxval(ee(1:na)), sum(ee(1:na)), size(ee)
    write(*,'(a11,4I30)') 'Min/max ne', minval(ne(1:na)), maxval(ne(1:na)), sum(ne(1:na)), size(ne)
    write(*,'(a11,4I30)') 'Min/max se', minval(se(1:na)), maxval(se(1:na)), sum(se(1:na)), size(se)
    write(*,'(a11,4I30)') 'Min/max nw', minval(nw(1:na)), maxval(nw(1:na)), sum(nw(1:na)), size(nw)
    write(*,'(a11,4I30)') 'Min/max sw', minval(sw(1:na)), maxval(sw(1:na)), sum(sw(1:na)), size(sw)
    write(*,'(a11,4I30)') 'Min/max sse', minval(sse(1:na)), maxval(sse(1:na)), sum(sse(1:na)), size(sse)
    write(*,'(a11,1I30)') 'count skipu ',  count(skipUcell1d(1:na))
    write(*,'(a11,1I30)') 'count skipt ',  count(skipTcell1d(1:na))
    write (*,*) '-----------------------------------------------'
    write (*,*) 'Statistics for all relevant 1D variables'
    write (*,*) '-----------------------------------------------'
    write(*,'(a11,3f25.8,1I30)') '     HTE1d:',   minval(HTE1d(1:na)),      maxval(HTE1d(1:na)),     sum(HTE1d(1:na)),   size(HTE1d)
    write(*,'(a11,3f25.8,1I30)') '     HTN1d:',   minval( HTN1d(1:na)),     maxval( HTN1d(1:na)),    sum( HTN1d(1:na)),  size( HTN1d)
    write(*,'(a11,3f25.8,1I30)') '   HTE1dm1:',   minval(HTE1dm1(1:na)),    maxval(HTE1dm1(1:na)),   sum(HTE1dm1(1:na)), size(HTE1dm1)
    write(*,'(a11,3f25.8,1I30)') '   HTN1dm1:',   minval(HTN1dm1(1:na)),    maxval(HTN1dm1(1:na)),   sum(HTN1dm1(1:na)), size(HTN1dm1)
    write(*,'(a11,3f25.8,1I30)') '      uvel:',   minval(uvel(1:navel)),    maxval(uvel(1:navel)),   sum(uvel(1:navel)), size(uvel)
    write(*,'(a11,3f25.8,1I30)') '      vvel:',   minval(vvel(1:navel)),    maxval(vvel(1:navel)),   sum(vvel(1:navel)), size(vvel)
    write(*,'(a11,3f25.8,1I30)') '       dxt:',   minval(dxt(1:na)),        maxval(dxt(1:na)),       sum(dxt(1:na)),     size(dxt)
    write(*,'(a11,3f25.8,1I30)') '       dyt:',   minval(dyt(1:na)),        maxval(dyt(1:na)),       sum(dyt(1:na)),     size(dyt)
    write(*,'(a11,3f25.8)') '  strength:',        minval(strength(1:na)),   maxval(strength(1:na)),  sum(strength(1:na))
    write(*,'(a11,3f25.8)') ' stressp_1:',        minval(stressp_1(1:na)),  maxval(stressp_1(1:na)), sum(stressp_1(1:na))
    write(*,'(a11,3f25.8)') ' stressp_2:',        minval(stressp_2(1:na)),  maxval(stressp_2(1:na)), sum(stressp_2(1:na))
    write(*,'(a11,3f25.8)') ' stressp_3:',        minval(stressp_3(1:na)),  maxval(stressp_3(1:na)), sum(stressp_3(1:na))
    write(*,'(a11,3f25.8)') ' stressp_4:',        minval(stressp_4(1:na)),  maxval(stressp_4(1:na)), sum(stressp_4(1:na))
    write(*,'(a11,3f25.8)') ' stressm_1:',        minval(stressm_1(1:na)),  maxval(stressm_1(1:na)), sum(stressm_1(1:na))
    write(*,'(a11,3f25.8)') ' stressm_2:',        minval(stressm_2(1:na)),  maxval(stressm_2(1:na)), sum(stressm_2(1:na))
    write(*,'(a11,3f25.8)') ' stressm_3:',        minval(stressm_3(1:na)),  maxval(stressm_3(1:na)), sum(stressm_3(1:na))
    write(*,'(a11,3f25.8)') ' stressm_4:',        minval(stressm_4(1:na)),  maxval(stressm_4(1:na)), sum(stressm_4(1:na))
    write(*,'(a11,3f25.8)') 'stress12_1:',        minval(stress12_1(1:na)), maxval(stress12_1(1:na)),sum(stress12_1(1:na))
    write(*,'(a11,3f25.8)') 'stress12_2:',        minval(stress12_2(1:na)), maxval(stress12_2(1:na)),sum(stress12_2(1:na))
    write(*,'(a11,3f25.8)') 'stress12_3:',        minval(stress12_3(1:na)), maxval(stress12_3(1:na)),sum(stress12_3(1:na))
    write(*,'(a11,3f25.8)') 'stress12_4:',        minval(stress12_4(1:na)), maxval(stress12_4(1:na)),sum(stress12_4(1:na))
    write (*,*) '-----------------------------------------------'
    write (*,*) 'Statistics for stepu relevant 1D variables'
    write (*,*) '-----------------------------------------------'
    write(*,'(a11,3f25.8)') ' uvel_init:',   minval( uvel_init(1:na)),maxval( uvel_init(1:na)),sum(uvel_init(1:na))
    write(*,'(a11,3f25.8)') ' vvel_init:',   minval( vvel_init(1:na)),maxval( vvel_init(1:na)),sum(vvel_init(1:na))
    write(*,'(a11,3f25.8)') '  strintxU:',   minval(strintxU(1:na)),  maxval(strintxU(1:na)),  sum(strintxU(1:na))
    write(*,'(a11,3f25.8)') '  strintyU:',   minval(strintyU(1:na)),  maxval(strintyU(1:na)),  sum(strintyU(1:na))
    write(*,'(a11,3f25.8)') '   cdn_ocn:',   minval( cdn_ocn(1:na)),  maxval( cdn_ocn(1:na)),  sum( cdn_ocn(1:na))
    write(*,'(a11,3f25.8)') '       aiu:',   minval( aiu(1:na)),      maxval( aiu(1:na)),      sum( aiu(1:na))
    write(*,'(a11,3f25.8)') '      uocn:',   minval( uocn(1:na)),     maxval( uocn(1:na)),     sum( uocn(1:na))
    write(*,'(a11,3f25.8)') '      vocn:',   minval( vocn(1:na)),     maxval( vocn(1:na)),     sum( vocn(1:na))
    write(*,'(a11,3f25.8)') '   waterxU:',   minval( waterxU(1:na)),  maxval( waterxU(1:na)),  sum( waterxU(1:na))
    write(*,'(a11,3f25.8)') '   wateryU:',   minval( wateryU(1:na)),  maxval( wateryU(1:na)),  sum( wateryU(1:na))
    write(*,'(a11,3f25.8)') '   forcexU:',   minval( forcexU(1:na)),  maxval( forcexU(1:na)),  sum( forcexU(1:na))
    write(*,'(a11,3f25.8)') '   forceyU:',   minval( forceyU(1:na)),  maxval( forceyU(1:na)),  sum( forceyU(1:na))
    write(*,'(a11,3f25.8)') '  umassdti:',   minval( umassdti(1:na)), maxval( umassdti(1:na)), sum( umassdti(1:na))
    write(*,'(a11,3f25.8)') '       fmU:',   minval( fmU(1:na)),      maxval( fmU(1:na)),      sum( fmU(1:na))
    write(*,'(a11,3f25.8)') '    uarear:',   minval( uarear(1:na)),   maxval( uarear(1:na)),   sum( uarear(1:na))
    write(*,'(a11,3f25.8)') '      str1:',   minval(str1(1:na)),      maxval(str1(1:na)),      sum(str1(1:na))
    write(*,'(a11,3f25.8)') '      str2:',   minval(str2(1:na)),      maxval(str2(1:na)),      sum(str2(1:na))
    write(*,'(a11,3f25.8)') '      str3:',   minval(str3(1:na)),      maxval(str3(1:na)),      sum(str3(1:na))
    write(*,'(a11,3f25.8)') '      str4:',   minval(str4(1:na)),      maxval(str4(1:na)),      sum(str4(1:na))
    write(*,'(a11,3f25.8)') '      str5:',   minval(str5(1:na)),      maxval(str5(1:na)),      sum(str5(1:na))
    write(*,'(a11,3f25.8)') '      str6:',   minval(str6(1:na)),      maxval(str6(1:na)),      sum(str6(1:na))
    write(*,'(a11,3f25.8)') '      str7:',   minval(str7(1:na)),      maxval(str7(1:na)),      sum(str7(1:na))
    write(*,'(a11,3f25.8)') '      str8:',   minval(str8(1:na)),      maxval(str8(1:na)),      sum(str8(1:na))
  end subroutine stat_1d
  subroutine stat_out_1d()
    implicit none
    write (*,*) '-----------------------------------------------'
    write (*,*) 'Statistics for output relevant 1D variables    '
    write (*,*) '-----------------------------------------------'
    write(*,'(a11,3f25.8)') ' stressp_1:',   minval(stressp_1(1:na)), maxval(stressp_1(1:na)),sum(stressp_1(1:na))
    write(*,'(a11,3f25.8)') ' stressp_2:',   minval(stressp_2(1:na)), maxval(stressp_2(1:na)),sum(stressp_2(1:na))
    write(*,'(a11,3f25.8)') ' stressp_3:',   minval(stressp_3(1:na)), maxval(stressp_3(1:na)),sum(stressp_3(1:na))
    write(*,'(a11,3f25.8)') ' stressp_4:',   minval(stressp_4(1:na)), maxval(stressp_4(1:na)),sum(stressp_4(1:na))
    write(*,'(a11,3f25.8)') ' stressm_1:',   minval(stressm_1(1:na)), maxval(stressm_1(1:na)),sum(stressm_1(1:na))
    write(*,'(a11,3f25.8)') ' stressm_2:',   minval(stressm_2(1:na)), maxval(stressm_2(1:na)),sum(stressm_2(1:na))
    write(*,'(a11,3f25.8)') ' stressm_3:',   minval(stressm_3(1:na)), maxval(stressm_3(1:na)),sum(stressm_3(1:na))
    write(*,'(a11,3f25.8)') ' stressm_4:',   minval(stressm_4(1:na)), maxval(stressm_4(1:na)),sum(stressm_4(1:na))
    write(*,'(a11,3f25.8)') 'stress12_1:',    minval(stress12_1(1:na)), maxval(stress12_1(1:na)),sum(stress12_1(1:na))
    write(*,'(a11,3f25.8)') 'stress12_2:',    minval(stress12_2(1:na)), maxval(stress12_2(1:na)),sum(stress12_2(1:na))
    write(*,'(a11,3f25.8)') 'stress12_3:',    minval(stress12_3(1:na)), maxval(stress12_3(1:na)),sum(stress12_3(1:na))
    write(*,'(a11,3f25.8)') 'stress12_4:',    minval(stress12_4(1:na)), maxval(stress12_4(1:na)),sum(stress12_4(1:na))
    write(*,'(a11,3f25.8)') '      str1:',   minval(str1(1:navel)), maxval(str1(1:navel)),sum(str1(1:navel))
    write(*,'(a11,3f25.8)') '      str2:',   minval(str2(1:navel)), maxval(str2(1:navel)),sum(str2(1:navel))
    write(*,'(a11,3f25.8)') '      str3:',   minval(str3(1:navel)), maxval(str3(1:navel)),sum(str3(1:navel))
    write(*,'(a11,3f25.8)') '      str4:',   minval(str4(1:navel)), maxval(str4(1:navel)),sum(str4(1:navel))
    write(*,'(a11,3f25.8)') '      str5:',   minval(str5(1:navel)), maxval(str5(1:navel)),sum(str5(1:navel))
    write(*,'(a11,3f25.8)') '      str6:',   minval(str6(1:navel)), maxval(str6(1:navel)),sum(str6(1:navel))
    write(*,'(a11,3f25.8)') '      str7:',   minval(str7(1:navel)), maxval(str7(1:navel)),sum(str7(1:navel))
    write(*,'(a11,3f25.8)') '      str8:',   minval(str8(1:navel)), maxval(str8(1:navel)),sum(str8(1:navel))
    write(*,'(a11,3f25.8)') '   strintx:',   minval(strintxU(:)), maxval(strintxU(:)),sum(strintxU(:))
    write(*,'(a11,3f25.8)') '   strinty:',   minval(strintyU(:)), maxval(strintyU(:)),sum(strintyU(:))
    write(*,'(a11,3f25.8)') '      uvel:',   minval(uvel(:)), maxval(uvel(:)),sum(uvel(:))
    write(*,'(a11,3f25.8)') '      vvel:',   minval(vvel(:)), maxval(vvel(:)),sum(vvel(:))
  end subroutine stat_out_1d
  subroutine writeout_1d(id)
    use create_nml, only : na, nb
    implicit none
    character(len=1), intent(in), optional :: id
    integer(kind=int_kind) :: ios, ierr, i, iw, nbb
    character (100)        :: binfile
    integer (kind=int_kind) :: lun1
    integer (kind=int_kind) :: lun2
    write(*,*) 'Saving output in output_stress_v2_1d.bin and output_stepu_v2_1d.bin'
    if (present(id)) then
      binfile = trim('output_stress_v2'//trim(id)//'_1d.bin')
    else
      binfile = 'output_stress_v2_1d.bin'
    endif
    lun1 = io_new_unit()
    open(lun1,file=binfile, form='unformatted', access='stream', action='write', &
         status='replace', iostat=ios)
    if (ios .ne. 0) then
      write(*,*) ios
      stop ('Failed open output_stress_v2_1d.bin') 
    endif 
    write(lun1,iostat=ios)                                                      &
          stressp_1, stressp_2, stressp_3, stressp_4,                      &
          stressm_1, stressm_2, stressm_3, stressm_4,                      &
          stress12_1,stress12_2,stress12_3,stress12_4
    if (ios .ne. 0) then
      write(*,*) ios
      stop ('Failed write output_stress_v2_1d.bin') 
    endif 
    close(lun1)  
    flush(lun1)
!    binfile = 'output_stepu_v2_1d_all.bin'
!    open(lun,file=binfile, form='unformatted', access='stream', action='write', &
!         status='replace', iostat=ios)
!    write(lun,iostat=ios)                                                      &
!         strintx,strinty,uvel,vvel
!    close(lun)
    nbb=na-count(skipucell1d(:))
    write (*,*) nb, nbb
    allocate(Vstrintx(1:nb),Vstrinty(1:nb),Vuvel(1:nb), Vvvel(1:nb), stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D Vstrintx,...'
    i=0
    do iw=1,na
      if (skipucell1d(iw)) cycle
      i=i+1
      Vstrintx   (i)=strintxU(iw)
      Vstrinty   (i)=strintyU(iw)
      Vuvel      (i)=uvel(iw)
      Vvvel      (i)=vvel(iw)
    enddo
    if (i .ne. nbb) stop 'PROBLEM in generating output_stepu_v1_1d.bin'
    if (present(id)) then
      binfile = trim('output_stepu_v2'//trim(id)//'_1d.bin')
    else
      binfile = 'output_stepu_v2_1d.bin'
    endif
    lun2 = io_new_unit()
    open(lun2, file=binfile, form='unformatted', access='stream',     &
            action='write', status='replace', iostat=ios)
    if (ios .ne. 0) then
      write(*,*) ios
      stop ('Failed open output_stepu_v2_1d.bin') 
    endif 
    write(lun2,iostat=ios) Vstrintx,Vstrinty,Vuvel,Vvvel
    if (ios .ne. 0) then
      write(*,*) ios
      stop ('Failed write output_stepu_v2_1d.bin') 
    endif 
    close(lun2)
    flush(lun2)
    write(*,*) 'Output files closed and flushed, the rest is up to the FS'
  end subroutine writeout_1d
#endif

end module vars

