!=======================================================================
!
! Elastic-viscous-plastic sea ice dynamics model
! Computes ice velocity and deformation
!
! See:
!
! Hunke, E. C., and J. K. Dukowicz (1997). An elastic-viscous-plastic model
! for sea ice dynamics. J. Phys. Oceanogr., 27, 1849-1867.
!
! Hunke, E. C. (2001).  Viscous-Plastic Sea Ice Dynamics with the EVP Model:
! Linearization Issues. J. Comput. Phys., 170, 18-38.
!
! Hunke, E. C., and J. K. Dukowicz (2002).  The Elastic-Viscous-Plastic
! Sea Ice Dynamics Model in General Orthogonal Curvilinear Coordinates
! on a Sphere - Incorporation of Metric Terms. Mon. Weather Rev.,
! 130, 1848-1865.
!
! Hunke, E. C., and J. K. Dukowicz (2003).  The sea ice momentum
! equation in the free drift regime.  Los Alamos Tech. Rep. LA-UR-03-2219.
!
! Hibler, W. D. (1979). A dynamic thermodynamic sea ice model. J. Phys.
! Oceanogr., 9, 817-846.
!
! Bouillon, S., T. Fichefet, V. Legat and G. Madec (2013).  The
! elastic-viscous-plastic method revisited.  Ocean Model., 71, 2-12.
!
! author: Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb (LANL)
! 2004: Block structure added by William Lipscomb
! 2005: Removed boundary calls for stress arrays (WHL)
! 2006: Streamlined for efficiency by Elizabeth Hunke
!       Converted to free source form (F90)
!=======================================================================
module bench

   use ice_dyn_shared, only: Ktens, epp2i,capping, e_factor

   use ice_constants, only: c0, c1

contains

!=======================================================================
! Computes the rates of strain and internal stress components for
! each of the four corners on each T-grid cell.
! Computes stress terms for the momentum equation
!
! author: Elizabeth C. Hunke, LANL
subroutine stress (nx_block,   ny_block,   &
                               icellt,     &
                   indxti,     indxtj,     &
                   uvel,       vvel,       & 
                   dxT,        dyT,        & 
                   dxhy,       dyhx,       & 
                   cxp,        cyp,        & 
                   cxm,        cym,        &
                               DminTarea,  &
                   strength,               & 
                   stressp_1,  stressp_2,  & 
                   stressp_3,  stressp_4,  & 
                   stressm_1,  stressm_2,  & 
                   stressm_3,  stressm_4,  & 
                   stress12_1, stress12_2, &
                   stress12_3, stress12_4, &
                   str )

  use ice_kinds_mod
  use ice_constants, only: p027, p055, p111, p166, & 
      p2, p222, p25, p333, p5 
!puny is not used unless eliminate underflow is used (search Elimiate underflow in this function)
  use ice_dyn_shared, only: arlx1i, denom1, revp
  implicit none
  integer (kind=int_kind), intent(in) :: & 
     nx_block, ny_block, & ! block dimensions
     icellt                ! no. of cells where icetmask = 1

  integer (kind=int_kind), dimension (nx_block*ny_block) :: &
     indxti   , & ! compressed index in i-direction
     indxtj       ! compressed index in j-direction

  real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
     strength , & ! ice strength (N/m)
     uvel     , & ! x-component of velocity (m/s)
     vvel     , & ! y-component of velocity (m/s)
     dxT      , & ! width of T-cell through the middle (m)
     dyT      , & ! height of T-cell through the middle (m)
     dxhy     , & ! 0.5*(HTE - HTE)
     dyhx     , & ! 0.5*(HTN - HTN)
     cyp      , & ! 1.5*HTE - 0.5*HTE
     cxp      , & ! 1.5*HTN - 0.5*HTN
     cym      , & ! 0.5*HTE - 1.5*HTE
     cxm      , & ! 0.5*HTN - 1.5*HTN
     DminTarea    ! deltaminEVP*tarea

  real (kind=dbl_kind), dimension (nx_block,ny_block), & 
     intent(inout) :: &
     stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
     stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
     stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12

  real (kind=dbl_kind), dimension(nx_block,ny_block,8), & 
     intent(out) :: &
     str          ! stress combinations

  integer (kind=int_kind) :: &
     i, j, ij

  real (kind=dbl_kind) :: &
    divune, divunw, divuse, divusw            , & ! divergence
    tensionne, tensionnw, tensionse, tensionsw, & ! tension
    shearne, shearnw, shearse, shearsw        , & ! shearing
    Deltane, Deltanw, Deltase, Deltasw        , & ! Delt
    zetax2ne, zetax2nw, zetax2se, zetax2sw    , & ! 2 x zeta (bulk visc)
    etax2ne, etax2nw, etax2se, etax2sw        , & ! 2 x eta (shear visc)
    rep_prsne, rep_prsnw, rep_prsse, rep_prssw, & ! replacement pressure

    ssigpn, ssigps, ssigpe, ssigpw            , &
    ssigmn, ssigms, ssigme, ssigmw            , &
    ssig12n, ssig12s, ssig12e, ssig12w        , &
    ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122, &
    csigpne, csigpnw, csigpse, csigpsw        , &
    csigmne, csigmnw, csigmse, csigmsw        , &
    csig12ne, csig12nw, csig12se, csig12sw    , &
    str12ew, str12we, str12ns, str12sn        , &
    strp_tmp, strm_tmp

  !-----------------------------------------------------------------
  ! Initialize
  !-----------------------------------------------------------------

  str(:,:,:) = c0

  do ij = 1, icellt
     i = indxti(ij)
     j = indxtj(ij)

  !-----------------------------------------------------------------
  ! strain rates
  ! NOTE these are actually strain rates * area  (m^2/s)
  !-----------------------------------------------------------------

         call strain_rates (nx_block,   ny_block,   &
                            i,          j,          &
                            uvel,       vvel,       &
                            dxT,        dyT,        &
                            cxp,        cyp,        &
                            cxm,        cym,        &
                            divune,     divunw,     &
                            divuse,     divusw,     &
                            tensionne,  tensionnw,  &
                            tensionse,  tensionsw,  &
                            shearne,    shearnw,    &
                            shearse,    shearsw,    &
                            Deltane,    Deltanw,    &
                            Deltase,    Deltasw     )

         !-----------------------------------------------------------------
         ! viscosities and replacement pressure
         !-----------------------------------------------------------------

         call visc_replpress (strength(i,j), DminTarea(i,j), Deltane, &
                              zetax2ne, etax2ne, rep_prsne)

         call visc_replpress (strength(i,j), DminTarea(i,j), Deltanw, &
                              zetax2nw, etax2nw, rep_prsnw)

         call visc_replpress (strength(i,j), DminTarea(i,j), Deltasw, &
                              zetax2sw, etax2sw, rep_prssw)

         call visc_replpress (strength(i,j), DminTarea(i,j), Deltase, &
                              zetax2se, etax2se, rep_prsse)

  !-----------------------------------------------------------------
  ! the stresses                            ! kg/s^2
  ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
  !-----------------------------------------------------------------

! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code
         stressp_1 (i,j) = (stressp_1 (i,j)*(c1-arlx1i*revp) &
                           + arlx1i*(zetax2ne*divune - rep_prsne)) * denom1
         stressp_2 (i,j) = (stressp_2 (i,j)*(c1-arlx1i*revp) &
                           + arlx1i*(zetax2nw*divunw - rep_prsnw)) * denom1
         stressp_3 (i,j) = (stressp_3 (i,j)*(c1-arlx1i*revp) &
                           + arlx1i*(zetax2sw*divusw - rep_prssw)) * denom1
         stressp_4 (i,j) = (stressp_4 (i,j)*(c1-arlx1i*revp) &
                           + arlx1i*(zetax2se*divuse - rep_prsse)) * denom1

         stressm_1 (i,j) = (stressm_1 (i,j)*(c1-arlx1i*revp) &
                           + arlx1i*etax2ne*tensionne) * denom1
         stressm_2 (i,j) = (stressm_2 (i,j)*(c1-arlx1i*revp) &
                           + arlx1i*etax2nw*tensionnw) * denom1
         stressm_3 (i,j) = (stressm_3 (i,j)*(c1-arlx1i*revp) &
                           + arlx1i*etax2sw*tensionsw) * denom1
         stressm_4 (i,j) = (stressm_4 (i,j)*(c1-arlx1i*revp) &
                           + arlx1i*etax2se*tensionse) * denom1

         stress12_1(i,j) = (stress12_1(i,j)*(c1-arlx1i*revp) &
                           + arlx1i*p5*etax2ne*shearne) * denom1
         stress12_2(i,j) = (stress12_2(i,j)*(c1-arlx1i*revp) &
                           + arlx1i*p5*etax2nw*shearnw) * denom1
         stress12_3(i,j) = (stress12_3(i,j)*(c1-arlx1i*revp) &
                           + arlx1i*p5*etax2sw*shearsw) * denom1
         stress12_4(i,j) = (stress12_4(i,j)*(c1-arlx1i*revp) &
                           + arlx1i*p5*etax2se*shearse) * denom1

  !-----------------------------------------------------------------
  ! Eliminate underflows.
  ! The following code is commented out because it is relatively 
  ! expensive and most compilers include a flag that accomplishes
  ! the same thing more efficiently.  This code is cheaper than
  ! handling underflows if the compiler lacks a flag; uncomment
  ! it in that case.  The compiler flag is often described with the 
  ! phrase "flush to zero".
  !-----------------------------------------------------------------

!      stressp_1(i,j) = sign(max(abs(stressp_1(i,j)),puny),stressp_1(i,j))
!      stressp_2(i,j) = sign(max(abs(stressp_2(i,j)),puny),stressp_2(i,j))
!      stressp_3(i,j) = sign(max(abs(stressp_3(i,j)),puny),stressp_3(i,j))
!      stressp_4(i,j) = sign(max(abs(stressp_4(i,j)),puny),stressp_4(i,j))

!      stressm_1(i,j) = sign(max(abs(stressm_1(i,j)),puny),stressm_1(i,j))
!      stressm_2(i,j) = sign(max(abs(stressm_2(i,j)),puny),stressm_2(i,j))
!      stressm_3(i,j) = sign(max(abs(stressm_3(i,j)),puny),stressm_3(i,j))
!      stressm_4(i,j) = sign(max(abs(stressm_4(i,j)),puny),stressm_4(i,j))

!      stress12_1(i,j) = sign(max(abs(stress12_1(i,j)),puny),stress12_1(i,j))
!      stress12_2(i,j) = sign(max(abs(stress12_2(i,j)),puny),stress12_2(i,j))
!      stress12_3(i,j) = sign(max(abs(stress12_3(i,j)),puny),stress12_3(i,j))
!      stress12_4(i,j) = sign(max(abs(stress12_4(i,j)),puny),stress12_4(i,j))

  !-----------------------------------------------------------------
  ! combinations of the stresses for the momentum equation ! kg/s^2
  !-----------------------------------------------------------------

     ssigpn  = stressp_1(i,j) + stressp_2(i,j)
     ssigps  = stressp_3(i,j) + stressp_4(i,j)
     ssigpe  = stressp_1(i,j) + stressp_4(i,j)
     ssigpw  = stressp_2(i,j) + stressp_3(i,j)
     ssigp1  =(stressp_1(i,j) + stressp_3(i,j))*p055
     ssigp2  =(stressp_2(i,j) + stressp_4(i,j))*p055

     ssigmn  = stressm_1(i,j) + stressm_2(i,j)
     ssigms  = stressm_3(i,j) + stressm_4(i,j)
     ssigme  = stressm_1(i,j) + stressm_4(i,j)
     ssigmw  = stressm_2(i,j) + stressm_3(i,j)
     ssigm1  =(stressm_1(i,j) + stressm_3(i,j))*p055
     ssigm2  =(stressm_2(i,j) + stressm_4(i,j))*p055

     ssig12n = stress12_1(i,j) + stress12_2(i,j)
     ssig12s = stress12_3(i,j) + stress12_4(i,j)
     ssig12e = stress12_1(i,j) + stress12_4(i,j)
     ssig12w = stress12_2(i,j) + stress12_3(i,j)
     ssig121 =(stress12_1(i,j) + stress12_3(i,j))*p111
     ssig122 =(stress12_2(i,j) + stress12_4(i,j))*p111

     csigpne = p111*stressp_1(i,j) + ssigp2 + p027*stressp_3(i,j)
     csigpnw = p111*stressp_2(i,j) + ssigp1 + p027*stressp_4(i,j)
     csigpsw = p111*stressp_3(i,j) + ssigp2 + p027*stressp_1(i,j)
     csigpse = p111*stressp_4(i,j) + ssigp1 + p027*stressp_2(i,j)
     
     csigmne = p111*stressm_1(i,j) + ssigm2 + p027*stressm_3(i,j)
     csigmnw = p111*stressm_2(i,j) + ssigm1 + p027*stressm_4(i,j)
     csigmsw = p111*stressm_3(i,j) + ssigm2 + p027*stressm_1(i,j)
     csigmse = p111*stressm_4(i,j) + ssigm1 + p027*stressm_2(i,j)
     
     csig12ne = p222*stress12_1(i,j) + ssig122 &
              + p055*stress12_3(i,j)
     csig12nw = p222*stress12_2(i,j) + ssig121 &
              + p055*stress12_4(i,j)
     csig12sw = p222*stress12_3(i,j) + ssig122 &
              + p055*stress12_1(i,j)
     csig12se = p222*stress12_4(i,j) + ssig121 &
              + p055*stress12_2(i,j)

     str12ew = p5*dxt(i,j)*(p333*ssig12e + p166*ssig12w)
     str12we = p5*dxt(i,j)*(p333*ssig12w + p166*ssig12e)
     str12ns = p5*dyt(i,j)*(p333*ssig12n + p166*ssig12s)
     str12sn = p5*dyt(i,j)*(p333*ssig12s + p166*ssig12n)

  !-----------------------------------------------------------------
  ! for dF/dx (u momentum)
  !-----------------------------------------------------------------
     strp_tmp  = p25*dyT(i,j)*(p333*ssigpn  + p166*ssigps)
     strm_tmp  = p25*dyT(i,j)*(p333*ssigmn  + p166*ssigms)

     ! northeast (i,j)
     str(i,j,1) = -strp_tmp - strm_tmp - str12ew &
                  +dxhy(i,j)*(-csigpne + csigmne) + dyhx(i,j)*csig12ne

     ! northwest (i+1,j)
     str(i,j,2) =  strp_tmp + strm_tmp - str12we &
                  +dxhy(i,j)*(-csigpnw + csigmnw) + dyhx(i,j)*csig12nw

     strp_tmp  = p25*dyT(i,j)*(p333*ssigps  + p166*ssigpn)
     strm_tmp  = p25*dyT(i,j)*(p333*ssigms  + p166*ssigmn)

     ! southeast (i,j+1)
     str(i,j,3) = -strp_tmp - strm_tmp + str12ew &
                  +dxhy(i,j)*(-csigpse + csigmse) + dyhx(i,j)*csig12se

     ! southwest (i+1,j+1)
     str(i,j,4) =  strp_tmp + strm_tmp + str12we &
                  +dxhy(i,j)*(-csigpsw + csigmsw) + dyhx(i,j)*csig12sw

  !-----------------------------------------------------------------
  ! for dF/dy (v momentum)
  !-----------------------------------------------------------------
     strp_tmp  = p25*dxT(i,j)*(p333*ssigpe  + p166*ssigpw)
     strm_tmp  = p25*dxT(i,j)*(p333*ssigme  + p166*ssigmw)

     ! northeast (i,j)
     str(i,j,5) = -strp_tmp + strm_tmp - str12ns &
                  -dyhx(i,j)*(csigpne + csigmne) + dxhy(i,j)*csig12ne

     ! southeast (i,j+1)
     str(i,j,6) =  strp_tmp - strm_tmp - str12sn &
                  -dyhx(i,j)*(csigpse + csigmse) + dxhy(i,j)*csig12se

     strp_tmp  = p25*dxT(i,j)*(p333*ssigpw  + p166*ssigpe)
     strm_tmp  = p25*dxT(i,j)*(p333*ssigmw  + p166*ssigme)

     ! northwest (i+1,j)
     str(i,j,7) = -strp_tmp + strm_tmp + str12ns &
                  -dyhx(i,j)*(csigpnw + csigmnw) + dxhy(i,j)*csig12nw

     ! southwest (i+1,j+1)
     str(i,j,8) =  strp_tmp - strm_tmp + str12sn &
                  -dyhx(i,j)*(csigpsw + csigmsw) + dxhy(i,j)*csig12sw

  enddo                     ! ij

end subroutine stress

!=======================================================================
! Compute strain rates
!
! author: Elizabeth C. Hunke, LANL
!
! 2019: subroutine created by Philippe Blain, ECCC

subroutine strain_rates (      nx_block,   ny_block,   &
                               i,          j,          &
                               uvel,       vvel,       &
                               dxT,        dyT,        &
                               cxp,        cyp,        &
                               cxm,        cym,        &
                               divune,     divunw,     &
                               divuse,     divusw,     &
                               tensionne,  tensionnw,  &
                               tensionse,  tensionsw,  &
                               shearne,    shearnw,    &
                               shearse,    shearsw,    &
                               Deltane,    Deltanw,    &
                               Deltase,    Deltasw     )

  use ice_kinds_mod

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block    ! block dimensions

      integer (kind=int_kind), intent(in) :: &
         i, j                  ! indices

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         dxT      , & ! width of T-cell through the middle (m)
         dyT      , & ! height of T-cell through the middle (m)
         cyp      , & ! 1.5*HTE - 0.5*HTW
         cxp      , & ! 1.5*HTN - 0.5*HTS
         cym      , & ! 0.5*HTE - 1.5*HTW
         cxm          ! 0.5*HTN - 1.5*HTS

      real (kind=dbl_kind), intent(out):: &           ! at each corner :
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        Deltane, Deltanw, Deltase, Deltasw            ! Delta

      character(len=*), parameter :: subname = '(strain_rates)'

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------

      ! divergence  =  e_11 + e_22
      divune    = cyp(i,j)*uvel(i  ,j  ) - dyT(i,j)*uvel(i-1,j  ) &
                + cxp(i,j)*vvel(i  ,j  ) - dxT(i,j)*vvel(i  ,j-1)
      divunw    = cym(i,j)*uvel(i-1,j  ) + dyT(i,j)*uvel(i  ,j  ) &
                + cxp(i,j)*vvel(i-1,j  ) - dxT(i,j)*vvel(i-1,j-1)
      divusw    = cym(i,j)*uvel(i-1,j-1) + dyT(i,j)*uvel(i  ,j-1) &
                + cxm(i,j)*vvel(i-1,j-1) + dxT(i,j)*vvel(i-1,j  )
      divuse    = cyp(i,j)*uvel(i  ,j-1) - dyT(i,j)*uvel(i-1,j-1) &
                + cxm(i,j)*vvel(i  ,j-1) + dxT(i,j)*vvel(i  ,j  )

      ! tension strain rate  =  e_11 - e_22
      tensionne = -cym(i,j)*uvel(i  ,j  ) - dyT(i,j)*uvel(i-1,j  ) &
                +  cxm(i,j)*vvel(i  ,j  ) + dxT(i,j)*vvel(i  ,j-1)
      tensionnw = -cyp(i,j)*uvel(i-1,j  ) + dyT(i,j)*uvel(i  ,j  ) &
                +  cxm(i,j)*vvel(i-1,j  ) + dxT(i,j)*vvel(i-1,j-1)
      tensionsw = -cyp(i,j)*uvel(i-1,j-1) + dyT(i,j)*uvel(i  ,j-1) &
                +  cxp(i,j)*vvel(i-1,j-1) - dxT(i,j)*vvel(i-1,j  )
      tensionse = -cym(i,j)*uvel(i  ,j-1) - dyT(i,j)*uvel(i-1,j-1) &
                +  cxp(i,j)*vvel(i  ,j-1) - dxT(i,j)*vvel(i  ,j  )

      ! shearing strain rate  =  2*e_12
      shearne = -cym(i,j)*vvel(i  ,j  ) - dyT(i,j)*vvel(i-1,j  ) &
              -  cxm(i,j)*uvel(i  ,j  ) - dxT(i,j)*uvel(i  ,j-1)
      shearnw = -cyp(i,j)*vvel(i-1,j  ) + dyT(i,j)*vvel(i  ,j  ) &
              -  cxm(i,j)*uvel(i-1,j  ) - dxT(i,j)*uvel(i-1,j-1)
      shearsw = -cyp(i,j)*vvel(i-1,j-1) + dyT(i,j)*vvel(i  ,j-1) &
              -  cxp(i,j)*uvel(i-1,j-1) + dxT(i,j)*uvel(i-1,j  )
      shearse = -cym(i,j)*vvel(i  ,j-1) - dyT(i,j)*vvel(i-1,j-1) &
              -  cxp(i,j)*uvel(i  ,j-1) + dxT(i,j)*uvel(i  ,j  )
      ! Delta (in the denominator of zeta, eta)
      Deltane = sqrt(divune**2 + e_factor*(tensionne**2 + shearne**2))
      Deltanw = sqrt(divunw**2 + e_factor*(tensionnw**2 + shearnw**2))
      Deltasw = sqrt(divusw**2 + e_factor*(tensionsw**2 + shearsw**2))
      Deltase = sqrt(divuse**2 + e_factor*(tensionse**2 + shearse**2))

      end subroutine strain_rates

!===============================================================================
! Computes viscosities and replacement pressure for stress
! calculations. Note that tensile strength is included here.
!
! Hibler, W. D. (1979). A dynamic thermodynamic sea ice model. J. Phys.
! Oceanogr., 9, 817-846.
!
! Konig Beatty, C. and Holland, D. M.  (2010). Modeling landfast ice by
! adding tensile strength. J. Phys. Oceanogr. 40, 185-198.
!
! Lemieux, J. F. et al. (2016). Improving the simulation of landfast ice
! by combining tensile strength and a parameterization for grounded ridges.
! J. Geophys. Res. Oceans, 121, 7354-7368.
!===============================================================================

      subroutine visc_replpress(strength, DminArea, Delta, &
                                zetax2, etax2, rep_prs)
      use ice_kinds_mod

      real (kind=dbl_kind), intent(in)::  &
         strength, & !
         DminArea    !

      real (kind=dbl_kind), intent(in)::  &
         Delta 

      real (kind=dbl_kind), intent(out):: &
         zetax2  , & ! bulk viscosity
         etax2   , & ! shear viscosity
         rep_prs     ! replacement pressure

      ! local variables
      real (kind=dbl_kind) :: &
         tmpcalc     ! temporary

      character(len=*), parameter :: subname = '(visc_replpress)'

      ! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code

      tmpcalc =     capping *(strength/max(Delta,DminArea))+ &
                (c1-capping)*(strength/(Delta + DminArea))
      zetax2  = (c1+Ktens)*tmpcalc
      rep_prs = (c1-Ktens)*tmpcalc*Delta
      etax2   = epp2i*zetax2

      end subroutine visc_replpress

!=======================================================================
! Calculation of the surface stresses
! Integration of the momentum equation to find velocity (u,v)
!
! author: Elizabeth C. Hunke, LANL

      subroutine stepu (nx_block,   ny_block, &
                        icellU,     Cw,       &
                        indxUi,     indxUj,   &
                        aiX,        str,      &
                        uocn,       vocn,     &
                        waterx,     watery,   &
                        forcex,     forcey,   &
                        umassdti,   fm,       &
                        uarear,               &
                        strintx,    strinty,  &
                        taubx,      tauby,    &
                        uvel_init,  vvel_init,&
                        uvel,       vvel,     &
                        Tbu)
  use ice_kinds_mod

  use ice_dyn_shared, only: brlx, revp, u0, cosw, sinw

  implicit none
     integer (kind=int_kind), intent(in) :: &
        nx_block, ny_block, & ! block dimensions
        icellu                ! total count when iceumask is true

     integer (kind=int_kind), dimension (nx_block*ny_block), &
        intent(in) :: &
        indxUi  , & ! compressed index in i-direction
        indxUj      ! compressed index in j-direction

  real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
     Tbu,      & ! coefficient for basal stress (N/m^2)
     uvel_init,& ! x-component of velocity (m/s), beginning of timestep
     vvel_init,& ! y-component of velocity (m/s), beginning of timestep
     aiX,      & ! ice fraction on u-grid
     waterx,   & ! for ocean stress calculation, x (m/s)
     watery,   & ! for ocean stress calculation, y (m/s)
     forcex,   & ! work array: combined atm stress and ocn tilt, x
     forcey,   & ! work array: combined atm stress and ocn tilt, y
     Umassdti, & ! mass of U-cell/dt (kg/m^2 s)
     uocn,     & ! ocean current, x-direction (m/s)
     vocn,     & ! ocean current, y-direction (m/s)
     fm,       & ! Coriolis param. * mass in U-cell (kg/s)
     uarear      ! 1/uarea

  real (kind=dbl_kind), dimension(nx_block,ny_block,8), &
     intent(in) :: &
     str

  real (kind=dbl_kind), dimension (nx_block,ny_block), &
     intent(inout) :: &
     uvel    , & ! x-component of velocity (m/s)
     vvel        ! y-component of velocity (m/s)

  real (kind=dbl_kind), dimension (nx_block,ny_block), &
     intent(inout) ::  & 
     strintx ,  & ! divergence of internal ice stress, x (N/m^2)
     strinty,   & ! divergence of internal ice stress, y (N/m^2)
     taubx,     & ! basal stress, x-direction (N/m^2)
     tauby       ! basal stress, y-direction (N/m^2)

  real (kind=dbl_kind), dimension (nx_block,ny_block), &
     intent(inout) :: &
     Cw                   ! ocean-ice neutral drag coefficient

  real (kind=dbl_kind), parameter ::                                           &
         rhow =  1026._dbl_kind  ! This originally originates from icepack

      ! local variables

      integer (kind=int_kind) :: i, j, ij

      real (kind=dbl_kind) :: &
         uold, vold        , & ! old-time uvel, vvel
         vrel              , & ! relative ice-ocean velocity
         cca,ccb,ab2,cc1,cc2,& ! intermediate variables
         taux, tauy,         & ! part of ocean stress term          
         Cb                    ! complete basal stress coeff


      !-----------------------------------------------------------------
      ! integrate the momentum equation
      !-----------------------------------------------------------------

      do ij =1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         uold = uvel(i,j)
         vold = vvel(i,j)

         ! (magnitude of relative ocean current)*rhow*drag*aice
         vrel = aiX(i,j)*rhow*Cw(i,j)*sqrt((uocn(i,j) - uold)**2 + &
                                           (vocn(i,j) - vold)**2)  ! m/s
         ! ice/ocean stress
         taux = vrel*waterx(i,j) ! NOTE this is not the entire
         tauy = vrel*watery(i,j) ! ocn stress term

         Cb  = Tbu(i,j) / (sqrt(uold**2 + vold**2) + u0) ! for basal stress

         ! revp = 0 for classic evp, 1 for revised evp
         cca = (brlx + revp)*umassdti(i,j) + vrel * cosw +Cb ! kg/m^2 s
         ccb = fm(i,j) + sign(c1,fm(i,j)) * vrel * sinw ! kg/m^2 s

         ab2 = cca**2 + ccb**2

         ! divergence of the internal stress tensor
         strintx(i,j) = uarear(i,j)* &
             (str(i,j,1) + str(i+1,j,2) + str(i,j+1,3) + str(i+1,j+1,4))
         strinty(i,j) = uarear(i,j)* &
             (str(i,j,5) + str(i,j+1,6) + str(i+1,j,7) + str(i+1,j+1,8))

         ! finally, the velocity components
         cc1 = strintx(i,j) + forcex(i,j) + taux &
             + umassdti(i,j)*(brlx*uold + revp*uvel_init(i,j))
         cc2 = strinty(i,j) + forcey(i,j) + tauy &
             + umassdti(i,j)*(brlx*vold + revp*vvel_init(i,j))
         uvel(i,j) = (cca*cc1 + ccb*cc2) / ab2 ! m/s
         vvel(i,j) = (cca*cc2 - ccb*cc1) / ab2

         ! calculate seabed stress component for outputs
         ! only needed on last iteration.
         taubx(i,j) = -uvel(i,j)*Cb
         tauby(i,j) = -vvel(i,j)*Cb

      enddo                     ! ij

      end subroutine stepu
end module bench
