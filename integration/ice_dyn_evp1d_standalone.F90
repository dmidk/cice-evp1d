! test module
module cicestuff

  use ice_kinds_mod
  use ice_dyn_shared, only: set_evp_parameters, ndte
  implicit none
  public
  integer(kind=int_kind) :: ierr
  integer(kind=int_kind), parameter :: nx_block=100, ny_block=120, max_block=1, nghost=1
  real(kind=dbl_kind), parameter :: dt =300
  integer(kind=int_kind) :: nx_global,ny_global
  logical(kind=log_kind), allocatable, dimension(:,:,:) :: L_iceUmask, L_tmask
  logical(kind=log_kind), allocatable, dimension(:,:,:) :: L_iceTmask
  real(kind=dbl_kind), allocatable, dimension(:,:,:) :: L_dyT, L_dxT, L_uarear
  real(kind=dbl_kind), allocatable, dimension(:,:) :: G_HTE, G_HTN
  real(kind=dbl_kind), allocatable, dimension(:,:,:) ::  &
          L_stressp_1(:,:,:), L_stressp_2(:,:,:), &
          L_stressp_3(:,:,:), L_stressp_4(:,:,:), &
          L_stressm_1(:,:,:), L_stressm_2(:,:,:), &
          L_stressm_3(:,:,:), L_stressm_4(:,:,:), &
         L_stress12_1(:,:,:),L_stress12_2(:,:,:), &
         L_stress12_3(:,:,:),L_stress12_4(:,:,:), &
         L_strength(:,:,:), &
            L_cdn_ocn(:,:,:),       L_aiu(:,:,:), &
               L_uocn(:,:,:),      L_vocn(:,:,:), &
            L_waterxU(:,:,:),   L_wateryU(:,:,:), &
            L_forcexU(:,:,:),   L_forceyU(:,:,:), &
           L_umassdti(:,:,:),       L_fmU(:,:,:), &
           L_strintxU(:,:,:),  L_strintyU(:,:,:), &
                L_Tbu(:,:,:),        L_Cb(:,:,:), &
               L_uvel(:,:,:),      L_vvel(:,:,:)

  contains
  subroutine  ciceinit()
  ! allocate static
   allocate(L_iceTmask(nx_block,ny_block,max_block), L_iceUmask(nx_block,ny_block,max_block), &
            L_tmask(nx_block,ny_block,max_block),                                        &
              L_dxT(nx_block,ny_block,max_block),    L_dyT(nx_block,ny_block,max_block), &
           L_uarear(nx_block,ny_block,max_block),    G_HTE(nx_block,ny_block          ), &           
              G_HTN(nx_block,ny_block          ),stat=ierr)
 ! allocate dynamic
 allocate(L_stressp_1(nx_block,ny_block,max_block), L_stressp_2(nx_block,ny_block,max_block), &
          L_stressp_3(nx_block,ny_block,max_block), L_stressp_4(nx_block,ny_block,max_block), &
          L_stressm_1(nx_block,ny_block,max_block), L_stressm_2(nx_block,ny_block,max_block), &
          L_stressm_3(nx_block,ny_block,max_block), L_stressm_4(nx_block,ny_block,max_block), &
         L_stress12_1(nx_block,ny_block,max_block),L_stress12_2(nx_block,ny_block,max_block), &
         L_stress12_3(nx_block,ny_block,max_block),L_stress12_4(nx_block,ny_block,max_block), &
         L_strength(nx_block,ny_block,max_block), &
            L_cdn_ocn(nx_block,ny_block,max_block),       L_aiu(nx_block,ny_block,max_block), &
               L_uocn(nx_block,ny_block,max_block),      L_vocn(nx_block,ny_block,max_block), &
            L_waterxU(nx_block,ny_block,max_block),   L_wateryU(nx_block,ny_block,max_block), &
            L_forcexU(nx_block,ny_block,max_block),   L_forceyU(nx_block,ny_block,max_block), &
           L_umassdti(nx_block,ny_block,max_block),       L_fmU(nx_block,ny_block,max_block), &
           L_strintxU(nx_block,ny_block,max_block),  L_strintyU(nx_block,ny_block,max_block), &
                L_Tbu(nx_block,ny_block,max_block),        L_Cb(nx_block,ny_block,max_block), &
               L_uvel(nx_block,ny_block,max_block),      L_vvel(nx_block,ny_block,max_block), &
           stat=ierr )
            if (ierr/=0) stop 'Error allocating dynamic'

  L_iceTmask=.false.
!  L_iceTmask=0
  nx_global=nx_block-2*nghost
  ny_global=ny_block-2*nghost
  L_iceUmask=.false.
  L_tmask=.false.
  L_tmask(20:40,20:40,1)=.true.
  L_iceUmask(20:40,20:40,1)=.true.
  L_icetmask(20:40,20:40,1)=.true.
  L_umassdti=1.
  L_aiu=0.9 
  L_dxT=1.
  L_dyT=1.
  L_Uarear=2.
  G_HTE=1.
  G_HTN=1.
  L_uvel=0.3
  L_vvel=-0.3
  ndte=1
  call set_evp_parameters(dt)
  end subroutine  ciceinit
end module cicestuff

! simple test
program cice_evp1dstandalone
  use cicestuff, only : nx_global, ny_global, nx_block, ny_block, max_block,      &
                        ciceinit, nghost, L_dyT, L_dxT, L_Uarear, L_tmask, G_HTE, G_HTN, &
                        L_stressp_1 , L_stressp_2 , L_stressp_3, L_stressp_4,     &
                        L_stressm_1 , L_stressm_2 , L_stressm_3, L_stressm_4,     &
                        L_stress12_1, L_stress12_2, L_stress12_3,L_stress12_4,    &
                        L_strength,                                               &
                        L_cdn_ocn   , L_aiu       , L_uocn     , L_vocn     ,     &
                        L_waterxU   , L_wateryU   , L_forcexU  , L_forceyU  ,     &
                        L_umassdti  , L_fmU       , L_strintxU , L_strintyU ,     &
                        L_Tbu       , L_Cb        , L_uvel     , L_vvel     ,     &
                        L_icetmask , L_iceUmask

  use ice_dyn_evp1d
  character(10),parameter :: debugf='after1d'

  call ciceinit

  ! the evp_1d_init will be called before first timestep and after ALL 2D arrays are allocated and initialized with pre time-step 0 data, will open evp1d.log, will read evp1d.nml, will allocate all 1D arrays and initialize them acc. to content in 2D arrays.
!  call evp_1d_init(10,20,30,cnx,cny)
  call dyn_evp1d_init(nx_global, ny_global, nx_block, ny_block, max_block, nghost, &
                      L_dyT, L_dxT, L_uarear, L_tmask,                             &
                      G_HTE, G_HTN)
  L_stressp_1 = 1.
  L_stressp_2 = 2.
  L_stressp_3 = 3.
  L_stressp_4 = 4.
  L_stressm_1 = 5.
  L_stressm_2 = 6.
  L_stressm_3 = 7.
  L_stressm_4 = 8.
  
  call dyn_evp1d_run(L_stressp_1 , L_stressp_2 , L_stressp_3, L_stressp_4,     &
                     L_stressm_1 , L_stressm_2 , L_stressm_3, L_stressm_4,     &
                     L_stress12_1, L_stress12_2, L_stress12_3,L_stress12_4,    &
                     L_strength,                                               &
                     L_cdn_ocn   , L_aiu       , L_uocn     , L_vocn     ,     &
                     L_waterxU   , L_wateryU   , L_forcexU  , L_forceyU  ,     &
                     L_umassdti  , L_fmU       , L_strintxU , L_strintyU ,     &
                     L_Tbu       , L_Cb        , L_uvel     , L_vvel     ,     &
                     L_icetmask , L_iceUmask)
  write(*,*) 'a'
  call dyn_evp1d_run(L_stressp_1 , L_stressp_2 , L_stressp_3, L_stressp_4,     &
                     L_stressm_1 , L_stressm_2 , L_stressm_3, L_stressm_4,     &
                     L_stress12_1, L_stress12_2, L_stress12_3,L_stress12_4,    &
                     L_strength,                                               &
                     L_cdn_ocn   , L_aiu       , L_uocn     , L_vocn     ,     &
                     L_waterxU   , L_wateryU   , L_forcexU  , L_forceyU  ,     &
                     L_umassdti  , L_fmU       , L_strintxU , L_strintyU ,     &
                     L_Tbu       , L_Cb        , L_uvel     , L_vvel     ,     &
                     L_icetmask , L_iceUmask)
  write(*,*) 'b'
  call dyn_evp2d_dump(L_stressp_1 , L_stressp_2 , L_stressp_3, L_stressp_4,     &
                      L_stressm_1 , L_stressm_2 , L_stressm_3, L_stressm_4,     &
                      L_stress12_1, L_stress12_2, L_stress12_3,L_stress12_4,    &
                      L_strength,                                               &
                      L_cdn_ocn   , L_aiu       , L_uocn     , L_vocn     ,     &
                      L_waterxU   , L_wateryU   , L_forcexU  , L_forceyU  ,     &
                      L_umassdti  , L_fmU       , L_strintxU , L_strintyU ,     &
                      L_Tbu       , L_Cb        , L_uvel     , L_vvel     ,     &
                      L_icetmask , L_iceUmask, debugf)
  ! deallocate all 1D arrays, close files, say goodbye in the logfile.... Will be called after last timestep is completed
  write(*,*) 'c'
  call dyn_evp1d_finalize()
end
