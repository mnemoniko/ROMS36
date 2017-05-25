      MODULE bulk_flux_mod
!
!svn $Id: bulk_flux.F 1484 2012-06-13 19:25:03Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the bulk parameterization of surface wind     !
!  stress and surface net heat fluxes.                                 !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Fairall, C.W., E.F. Bradley, D.P. Rogers, J.B. Edson and G.S.     !
!      Young, 1996:  Bulk parameterization of air-sea fluxes for       !
!      tropical ocean-global atmosphere Coupled-Ocean Atmosphere       !
!      Response Experiment, JGR, 101, 3747-3764.                       !
!                                                                      !
!    Fairall, C.W., E.F. Bradley, J.S. Godfrey, G.A. Wick, J.B.        !
!      Edson, and G.S. Young, 1996:  Cool-skin and warm-layer          !
!      effects on sea surface temperature, JGR, 101, 1295-1308.        !
!                                                                      !
!    Liu, W.T., K.B. Katsaros, and J.A. Businger, 1979:  Bulk          !
!        parameterization of the air-sea exchange of heat and          !
!        water vapor including the molecular constraints at            !
!        the interface, J. Atmos. Sci, 36, 1722-1735.                  !
!                                                                      !
!  Adapted from COARE code written originally by David Rutgers and     !
!  Frank Bradley.                                                      !
!                                                                      !
!   option for equivalent salt fluxes added by Paul Goodman     !
!  (10/2004).                                                          !
!                                                                      !
!  Modified by Kate Hedstrom for COARE version 3.0 (03/2005).          !
!  Modified by Jim Edson to correct specific hunidities.               !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!     Fairall et al., 2003: J. Climate, 16, 571-591.                   !
!                                                                      !
!     Taylor, P. K., and M. A. Yelland, 2001: The dependence of sea    !
!     surface roughness on the height and steepness of the waves.      !
!     J. Phys. Oceanogr., 31, 572-590.                                 !
!                                                                      !
!     Oost, W. A., G. J. Komen, C. M. J. Jacobs, and C. van Oort, 2002:!
!     New evidence for a relation between wind stress and wave age     !
!     from measurements during ASGAMAGE. Bound.-Layer Meteor., 103,    !
!     409-438.                                                         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: bulk_flux, bulk_psiu, bulk_psit
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE bulk_flux (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_ice
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL wclock_on (ng, iNLM, 17)
      CALL bulk_flux_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nrhs(ng),                                    &
     &                     liold(ng),                                   &
     &                     linew(ng),                                   &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
     &                     MIXING(ng) % alpha,                          &
     &                     MIXING(ng) % beta,                           &
     &                     OCEAN(ng) % rho,                             &
     &                     OCEAN(ng) % t,                               &
     &                     FORCES(ng) % Hair,                           &
     &                     FORCES(ng) % Pair,                           &
     &                     FORCES(ng) % Tair,                           &
     &                     FORCES(ng) % Uwind,                          &
     &                     FORCES(ng) % Vwind,                          &
     &                     FORCES(ng) % cloud,                          &
     &                     FORCES(ng) % rain,                           &
     &                     FORCES(ng) % lhflx,                          &
     &                     FORCES(ng) % lrflx,                          &
     &                     FORCES(ng) % shflx,                          &
     &                     FORCES(ng) % srflx,                          &
     &                     FORCES(ng) % stflx,                          &
     &                     ICE(ng) % ai,                                &
     &                     ICE(ng) % hi,                                &
     &                     FORCES(ng) % qai_n,                          &
     &                     ICE(ng) % hsn,                               &
     &                     ICE(ng) % tis,                               &
     &                     ICE(ng) % sfwat,                             &
     &                     ICE(ng) % coef_ice_heat,                     &
     &                     ICE(ng) % rhs_ice_heat,                      &
     &                     FORCES(ng) % snow_n,                         &
     &                     FORCES(ng) % sustr_aw,                       &
     &                     FORCES(ng) % svstr_aw,                       &
     &                     FORCES(ng) % qao_n,                          &
     &                     FORCES(ng) % tau_aix_n,                      &
     &                     FORCES(ng) % tau_aiy_n,                      &
     &                     FORCES(ng) % EminusP,                        &
     &                     FORCES(ng) % evap,                           &
     &                     FORCES(ng) % sustr,                          &
     &                     FORCES(ng) % svstr)
      CALL wclock_off (ng, iNLM, 17)
      RETURN
      END SUBROUTINE bulk_flux
!
!***********************************************************************
      SUBROUTINE bulk_flux_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs,                                  &
     &                           liold, linew,                          &
     &                           rmask, umask, vmask,                   &
     &                           alpha, beta, rho, t,                   &
     &                           Hair, Pair, Tair, Uwind, Vwind,        &
     &                           cloud,                                 &
     &                           rain, lhflx, lrflx, shflx,             &
     &                           srflx, stflx,                          &
     &                           ai, hi,                                &
     &                           qai_n,                                 &
     &                           hsn, tis, sfwat,                       &
     &                           coef_ice_heat, rhs_ice_heat, snow_n,   &
     &                           sustr_aw, svstr_aw, qao_n,             &
     &                           tau_aix_n, tau_aiy_n,                  &
     &                           EminusP, evap,                         &
     &                           sustr, svstr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_ice
!
      USE exchange_2d_mod
      USE mp_exchange_mod, ONLY : mp_exchange2d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
      integer, intent(in) :: liold
      integer, intent(in) :: linew
!
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: alpha(LBi:,LBj:)
      real(r8), intent(in) :: beta(LBi:,LBj:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: Hair(LBi:,LBj:)
      real(r8), intent(in) :: Pair(LBi:,LBj:)
      real(r8), intent(in) :: Tair(LBi:,LBj:)
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
      real(r8), intent(in) :: cloud(LBi:,LBj:)
      real(r8), intent(in) :: rain(LBi:,LBj:)
      real(r8), intent(inout) :: lhflx(LBi:,LBj:)
      real(r8), intent(inout) :: lrflx(LBi:,LBj:)
      real(r8), intent(inout) :: shflx(LBi:,LBj:)
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
      real(r8), intent(out) :: sustr_aw(LBi:,LBj:)
      real(r8), intent(out) :: svstr_aw(LBi:,LBj:)
      real(r8), intent(out) :: qao_n(LBi:,LBj:)
      real(r8), intent(in) :: ai(LBi:,LBj:,:)
      real(r8), intent(in) :: hi(LBi:,LBj:,:)
      real(r8), intent(out) :: qai_n(LBi:,LBj:)
      real(r8), intent(in) :: hsn(LBi:,LBj:,:)
      real(r8), intent(in) :: tis(LBi:,LBj:)
      real(r8), intent(in) :: sfwat(LBi:,LBj:,:)
      real(r8), intent(out) :: coef_ice_heat(LBi:,LBj:)
      real(r8), intent(out) :: rhs_ice_heat(LBi:,LBj:)
      real(r8), intent(out) :: snow_n(LBi:,LBj:)
      real(r8), intent(out) :: tau_aix_n(LBi:,LBj:)
      real(r8), intent(out) :: tau_aiy_n(LBi:,LBj:)
      real(r8), intent(out) :: EminusP(LBi:,LBj:)
      real(r8), intent(out) :: evap(LBi:,LBj:)
      real(r8), intent(out) :: sustr(LBi:,LBj:)
      real(r8), intent(out) :: svstr(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: Iter, i, j, k
      integer, parameter :: IterMax = 3
      real(r8), parameter :: eps = 1.0E-20_r8
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8
      real(r8) :: Bf, Cd, Ce, Ch, Hl, Hlw, Hscale, Hscale2, Hs, Hsr, IER
      real(r8) :: PairM,  RH, Taur
      real(r8) :: Wspeed, ZQoL, ZToL
      real(r8), parameter :: alb_i_dry=0.65_r8
      real(r8), parameter :: alb_i_wet=0.60_r8
      real(r8), parameter :: alb_s_dry=0.85_r8
      real(r8), parameter :: alb_s_wet=0.72_r8
      real(r8), parameter :: alb_w=0.1_r8
      real(r8) :: albs, qswi, qlwi, qlh_i, qsh_i
      real(r8) :: le_i, dq_i,fqlat1, sfc_temp, slp, Qsati
      real(r8) :: vap_p_i
      real(r8) :: albi, albsn, thki_n, thksn_n
      real(r8) :: cff, cff1, cff2, diffh, diffw, oL, upvel
      real(r8) :: twopi_inv, wet_bulb
      real(r8) :: sfwat_n
      real(r8) :: e_sat, vap_p
      real(r8), dimension(IminS:ImaxS) :: CC
      real(r8), dimension(IminS:ImaxS) :: Cd10
      real(r8), dimension(IminS:ImaxS) :: Ch10
      real(r8), dimension(IminS:ImaxS) :: Ct10
      real(r8), dimension(IminS:ImaxS) :: charn
      real(r8), dimension(IminS:ImaxS) :: Ct
      real(r8), dimension(IminS:ImaxS) :: Cwave
      real(r8), dimension(IminS:ImaxS) :: Cwet
      real(r8), dimension(IminS:ImaxS) :: delQ
      real(r8), dimension(IminS:ImaxS) :: delQc
      real(r8), dimension(IminS:ImaxS) :: delT
      real(r8), dimension(IminS:ImaxS) :: delTc
      real(r8), dimension(IminS:ImaxS) :: delW
      real(r8), dimension(IminS:ImaxS) :: L
      real(r8), dimension(IminS:ImaxS) :: L10
      real(r8), dimension(IminS:ImaxS) :: Q
      real(r8), dimension(IminS:ImaxS) :: Qair
      real(r8), dimension(IminS:ImaxS) :: Qpsi
      real(r8), dimension(IminS:ImaxS) :: Qsea
      real(r8), dimension(IminS:ImaxS) :: Qstar
      real(r8), dimension(IminS:ImaxS) :: rhoAir
      real(r8), dimension(IminS:ImaxS) :: rhoSea
      real(r8), dimension(IminS:ImaxS) :: Ri
      real(r8), dimension(IminS:ImaxS) :: Ribcu
      real(r8), dimension(IminS:ImaxS) :: Rr
      real(r8), dimension(IminS:ImaxS) :: Scff
      real(r8), dimension(IminS:ImaxS) :: TairC
      real(r8), dimension(IminS:ImaxS) :: TairK
      real(r8), dimension(IminS:ImaxS) :: Tcff
      real(r8), dimension(IminS:ImaxS) :: Tpsi
      real(r8), dimension(IminS:ImaxS) :: TseaC
      real(r8), dimension(IminS:ImaxS) :: TseaK
      real(r8), dimension(IminS:ImaxS) :: Tstar
      real(r8), dimension(IminS:ImaxS) :: u10
      real(r8), dimension(IminS:ImaxS) :: VisAir
      real(r8), dimension(IminS:ImaxS) :: WaveLength
      real(r8), dimension(IminS:ImaxS) :: Wgus
      real(r8), dimension(IminS:ImaxS) :: Wmag
      real(r8), dimension(IminS:ImaxS) :: Wpsi
      real(r8), dimension(IminS:ImaxS) :: Wstar
      real(r8), dimension(IminS:ImaxS) :: Zetu
      real(r8), dimension(IminS:ImaxS) :: Zo10
      real(r8), dimension(IminS:ImaxS) :: ZoT10
      real(r8), dimension(IminS:ImaxS) :: ZoL
      real(r8), dimension(IminS:ImaxS) :: ZoQ
      real(r8), dimension(IminS:ImaxS) :: ZoT
      real(r8), dimension(IminS:ImaxS) :: ZoW
      real(r8), dimension(IminS:ImaxS) :: ZWoL
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hlv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: LHeat
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: LRad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: SHeat
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: SRad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Taux
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Tauy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Coef_Ice
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: RHS_Ice
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Qai
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Taux_Ice
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Tauy_Ice
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ice_thick
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: snow_thick
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!=======================================================================
!  Atmosphere-Ocean bulk fluxes parameterization.
!=======================================================================
!
      Hscale=rho0*Cp
      Hscale2=1.0_r8/(rho0*Cp)
      twopi_inv=0.5_r8/pi
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
!
!  Input bulk parameterization fields.
!
          Wmag(i)=SQRT(Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j))
          PairM=Pair(i,j)
          TairC(i)=Tair(i,j)
          TairK(i)=TairC(i)+273.16_r8
          TseaC(i)=t(i,j,N(ng),nrhs,itemp)
          TseaK(i)=TseaC(i)+273.16_r8
          rhoSea(i)=rho(i,j,N(ng))+1000.0_r8
          RH=Hair(i,j)
          SRad(i,j)=srflx(i,j)*Hscale
          Tcff(i)=alpha(i,j)
          Scff(i)=beta(i,j)
!
!  Initialize.
!
          delTc(i)=0.0_r8
          delQc(i)=0.0_r8
          LHeat(i,j)=lhflx(i,j)*Hscale
          SHeat(i,j)=shflx(i,j)*Hscale
          Taur=0.0_r8
          Taux(i,j)=0.0_r8
          Tauy(i,j)=0.0_r8
!
!-----------------------------------------------------------------------
!  Compute net longwave radiation (W/m2), LRad.
!-----------------------------------------------------------------------
!
!  Use Berliand (1952) formula to calculate net longwave radiation.
!  The equation for saturation vapor pressure is from Gill (Atmosphere-
!  Ocean Dynamics, pp 606). Here the coefficient in the cloud term
!  is assumed constant, but it is a function of latitude varying from
!  1.0 at poles to 0.5 at the Equator).
!
      IF(RH.lt.2.0_r8) THEN                               !WPB
          cff=(0.7859_r8+0.03477_r8*TairC(i))/                          &
     &        (1.0_r8+0.00412_r8*TairC(i))
          e_sat=10.0_r8**cff   ! saturation vapor pressure (hPa or mbar)
          vap_p=e_sat*RH       ! water vapor pressure (hPa or mbar)
      ELSE                                                !WPB
          vap_p=0.001_r8*PairM*RH/(1.0_r8+0.000378_r8*RH) !WPB
      ENDIF                                               !WPB
          cff2=TairK(i)*TairK(i)*TairK(i)
          cff1=cff2*TairK(i)
          LRad(i,j)=-emmiss*StefBo*                                     &
     &              (cff1*(0.39_r8-0.05_r8*SQRT(vap_p))*                &
!    &              (cff1*(0.39_r8-0.05_r8*SQRT(0.01_r8*vap_p))*        &
     &                    (1.0_r8-0.6823_r8*cloud(i,j)*cloud(i,j))+     &
     &               cff2*4.0_r8*(TseaK(i)-TairK(i)))
          LRad(i,j)=LRad(i,j)*rmask(i,j)
!
!-----------------------------------------------------------------------
!  Compute specific humidities (kg/kg).
!
!    note that Qair(i) is the saturation specific humidity at Tair
!                 Q(i) is the actual specific humidity
!              Qsea(i) is the saturation specific humidity at Tsea
!
!          Saturation vapor pressure in mb is first computed and then
!          converted to specific humidity in kg/kg
!
!          The saturation vapor pressure is computed from Teten formula
!          using the approach of Buck (1981):
!
!          Esat(mb) = (1.0007_r8+3.46E-6_r8*PairM(mb))*6.1121_r8*
!                  EXP(17.502_r8*TairC(C)/(240.97_r8+TairC(C)))
!
!          The ambient vapor is found from the definition of the
!          Relative humidity:
!
!          RH = W/Ws*100 ~ E/Esat*100   E = RH/100*Esat if RH is in %
!                                       E = RH*Esat     if RH fractional
!
!          The specific humidity is then found using the relationship:
!
!          Q = 0.622 E/(P + (0.622-1)e)
!
!          Q(kg/kg) = 0.62197_r8*(E(mb)/(PairM(mb)-0.378_r8*E(mb)))
!
!-----------------------------------------------------------------------
!
!  Compute air saturation vapor pressure (mb), using Teten formula.
!
          cff=(1.0007_r8+3.46E-6_r8*PairM)*6.1121_r8*                   &
     &        EXP(17.502_r8*TairC(i)/(240.97_r8+TairC(i)))
!
!  Compute specific humidity at Saturation, Qair (kg/kg).
!
          Qair(i)=0.62197_r8*(cff/(PairM-0.378_r8*cff))
!
!  Compute specific humidity, Q (kg/kg).
!
          IF (RH.lt.2.0_r8) THEN                       !RH fraction
            cff=cff*RH                                 !Vapor pres (mb)
            Q(i)=0.62197_r8*(cff/(PairM-0.378_r8*cff)) !Spec hum (kg/kg)
          ELSE          !RH input was actually specific humidity in g/kg
            Q(i)=RH/1000.0_r8                          !Spec Hum (kg/kg)
          END IF
!
!  Compute water saturation vapor pressure (mb), using Teten formula.
!
          cff=(1.0007_r8+3.46E-6_r8*PairM)*6.1121_r8*                   &
     &        EXP(17.502_r8*TseaC(i)/(240.97_r8+TseaC(i)))
!
!  Vapor Pressure reduced for salinity (Kraus & Businger, 1994, pp 42).
!
          cff=cff*0.98_r8
!
!  Compute Qsea (kg/kg) from vapor pressure.
!
          Qsea(i)=0.62197_r8*(cff/(PairM-0.378_r8*cff))
!
!-----------------------------------------------------------------------
!  Compute Monin-Obukhov similarity parameters for wind (Wstar),
!  heat (Tstar), and moisture (Qstar), Liu et al. (1979).
!-----------------------------------------------------------------------
!
!  Moist air density (kg/m3).
!
          rhoAir(i)=PairM*100.0_r8/(blk_Rgas*TairK(i)*                  &
     &                              (1.0_r8+0.61_r8*Q(i)))
!
!  Kinematic viscosity of dry air (m2/s), Andreas (1989).
!
          VisAir(i)=1.326E-5_r8*                                        &
     &              (1.0_r8+TairC(i)*(6.542E-3_r8+TairC(i)*             &
     &               (8.301E-6_r8-4.84E-9_r8*TairC(i))))
!
!  Compute latent heat of vaporization (J/kg) at sea surface, Hlv.
!
          Hlv(i,j)=(2.501_r8-0.00237_r8*TseaC(i))*1.0E+6_r8
!
!  Assume that wind is measured relative to sea surface and include
!  gustiness.
!
          Wgus(i)=0.5_r8
          delW(i)=SQRT(Wmag(i)*Wmag(i)+Wgus(i)*Wgus(i))
          delQ(i)=Qsea(i)-Q(i)
          delT(i)=TseaC(i)-TairC(i)
!
!  Neutral coefficients.
!
          ZoW(i)=0.0001_r8
          u10(i)=delW(i)*LOG(10.0_r8/ZoW(i))/LOG(blk_ZW(ng)/ZoW(i))
          Wstar(i)=0.035_r8*u10(i)
          Zo10(i)=0.011_r8*Wstar(i)*Wstar(i)/g+                         &
     &            0.11_r8*VisAir(i)/Wstar(i)
          Cd10(i)=(vonKar/LOG(10.0_r8/Zo10(i)))**2
          Ch10(i)=0.00115_r8
          Ct10(i)=Ch10(i)/sqrt(Cd10(i))
          ZoT10(i)=10.0_r8/EXP(vonKar/Ct10(i))
          Cd=(vonKar/LOG(blk_ZW(ng)/Zo10(i)))**2
!
!  Compute Richardson number.
!
          Ct(i)=vonKar/LOG(blk_ZT(ng)/ZoT10(i))  ! T transfer coefficient
          CC(i)=vonKar*Ct(i)/Cd
          delTc(i)=0.0_r8
          Ribcu(i)=-blk_ZW(ng)/(blk_Zabl*0.004_r8*blk_beta**3)
          Ri(i)=-g*blk_ZW(ng)*((delT(i)-delTc(i))+                      &
     &                          0.61_r8*TairK(i)*delQ(i))/              &
     &          (TairK(i)*delW(i)*delW(i))
          IF (Ri(i).lt.0.0_r8) THEN
            Zetu(i)=CC(i)*Ri(i)/(1.0_r8+Ri(i)/Ribcu(i))       ! Unstable
          ELSE
            Zetu(i)=CC(i)*Ri(i)/(1.0_r8+3.0_r8*Ri(i)/CC(i))   ! Stable
          END IF
          L10(i)=blk_ZW(ng)/Zetu(i)
!
!  First guesses for Monon-Obukhov similarity scales.
!
          Wstar(i)=delW(i)*vonKar/(LOG(blk_ZW(ng)/Zo10(i))-             &
     &                             bulk_psiu(blk_ZW(ng)/L10(i),pi))
          Tstar(i)=-(delT(i)-delTc(i))*vonKar/                          &
     &             (LOG(blk_ZT(ng)/ZoT10(i))-                           &
     &              bulk_psit(blk_ZT(ng)/L10(i),pi))
          Qstar(i)=-(delQ(i)-delQc(i))*vonKar/                          &
     &             (LOG(blk_ZQ(ng)/ZoT10(i))-                           &
     &              bulk_psit(blk_ZQ(ng)/L10(i),pi))
!
!  Modify Charnock for high wind speeds. The 0.125 factor below is for
!  1.0/(18.0-10.0).
!
          IF (delW(i).gt.18.0_r8) THEN
            charn(i)=0.018_r8
          ELSE IF ((10.0_r8.lt.delW(i)).and.(delW(i).le.18.0_r8)) THEN
            charn(i)=0.011_r8+                                          &
     &               0.125_r8*(0.018_r8-0.011_r8)*(delW(i)-10._r8)
          ELSE
            charn(i)=0.011_r8
          END IF
        END DO
!
!  Iterate until convergence. It usually converges within 3 iterations.
!
        DO Iter=1,IterMax
          DO i=Istr-1,IendR
            ZoW(i)=charn(i)*Wstar(i)*Wstar(i)/g+                        &
     &             0.11_r8*VisAir(i)/(Wstar(i)+eps)
            Rr(i)=ZoW(i)*Wstar(i)/VisAir(i)
!
!  Compute Monin-Obukhov stability parameter, Z/L.
!
            ZoQ(i)=MIN(1.15e-4_r8,5.5e-5_r8/Rr(i)**0.6_r8)
            ZoT(i)=ZoQ(i)
            ZoL(i)=vonKar*g*blk_ZW(ng)*                                 &
     &             (Tstar(i)*(1.0_r8+0.61_r8*Q(i))+                     &
     &                        0.61_r8*TairK(i)*Qstar(i))/               &
     &             (TairK(i)*Wstar(i)*Wstar(i)*                         &
     &             (1.0_r8+0.61_r8*Q(i))+eps)
            L(i)=blk_ZW(ng)/(ZoL(i)+eps)
!
!  Evaluate stability functions at Z/L.
!
            Wpsi(i)=bulk_psiu(ZoL(i),pi)
            Tpsi(i)=bulk_psit(blk_ZT(ng)/L(i),pi)
            Qpsi(i)=bulk_psit(blk_ZQ(ng)/L(i),pi)
!
!  Compute wind scaling parameters, Wstar.
!
            Wstar(i)=MAX(eps,delW(i)*vonKar/                            &
     &               (LOG(blk_ZW(ng)/ZoW(i))-Wpsi(i)))
            Tstar(i)=-(delT(i)-delTc(i))*vonKar/                        &
     &               (LOG(blk_ZT(ng)/ZoT(i))-Tpsi(i))
            Qstar(i)=-(delQ(i)-delQc(i))*vonKar/                        &
     &               (LOG(blk_ZQ(ng)/ZoQ(i))-Qpsi(i))
!
!  Compute gustiness in wind speed.
!
            Bf=-g/TairK(i)*                                             &
     &         Wstar(i)*(Tstar(i)+0.61_r8*TairK(i)*Qstar(i))
            IF (Bf.gt.0.0_r8) THEN
              Wgus(i)=blk_beta*(Bf*blk_Zabl)**r3
            ELSE
              Wgus(i)=0.2_r8
            END IF
            delW(i)=SQRT(Wmag(i)*Wmag(i)+Wgus(i)*Wgus(i))
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Compute Atmosphere/Ocean fluxes.
!-----------------------------------------------------------------------
!
        DO i=Istr-1,IendR
!
!  Compute transfer coefficients for momentum (Cd).
!
          Wspeed=SQRT(Wmag(i)*Wmag(i)+Wgus(i)*Wgus(i))
          Cd=Wstar(i)*Wstar(i)/(Wspeed*Wspeed+eps)
          Ch=Wstar(i)*Tstar(i)/(-Wspeed*delT(i)+0.0098_r8*blk_ZT(ng))
          Ce=Wstar(i)*Qstar(i)/(-Wspeed*delQ(i)+eps)
!
!  Compute turbulent sensible heat flux (W/m2), Hs.
!
          Hs=-blk_Cpa*rhoAir(i)*Wstar(i)*Tstar(i)
!
!  Compute sensible heat flux (W/m2) due to rainfall (kg/m2/s), Hsr.
!
          diffw=2.11E-5_r8*(TairK(i)/273.16_r8)**1.94_r8
          diffh=0.02411_r8*(1.0_r8+TairC(i)*                            &
     &                      (3.309E-3_r8-1.44E-6_r8*TairC(i)))/         &
     &          (rhoAir(i)*blk_Cpa)
          cff=Qair(i)*Hlv(i,j)/(blk_Rgas*TairK(i)*TairK(i))
          wet_bulb=1.0_r8/(1.0_r8+0.622_r8*(cff*Hlv(i,j)*diffw)/        &
     &                                     (blk_Cpa*diffh))
          Hsr=rain(i,j)*wet_bulb*blk_Cpw*                               &
     &        ((TseaC(i)-TairC(i))+(Qsea(i)-Q(i))*Hlv(i,j)/blk_Cpa)
          SHeat(i,j)=(Hs+Hsr)
          SHeat(i,j)=SHeat(i,j)*rmask(i,j)
!
!  Compute turbulent latent heat flux (W/m2), Hl.
!
          Hl=-Hlv(i,j)*rhoAir(i)*Wstar(i)*Qstar(i)
!
!  Compute Webb correction (Webb effect) to latent heat flux, Hlw.
!
          upvel=-1.61_r8*Wstar(i)*Qstar(i)-                             &
     &          (1.0_r8+1.61_r8*Q(i))*Wstar(i)*Tstar(i)/TairK(i)
          Hlw=rhoAir(i)*Hlv(i,j)*upvel*Q(i)
          LHeat(i,j)=(Hl+Hlw)
          LHeat(i,j)=LHeat(i,j)*rmask(i,j)
!
!  Compute momentum flux (N/m2) due to rainfall (kg/m2/s).
!
          Taur=0.85_r8*rain(i,j)*Wmag(i)
!
!  Compute wind stress components (N/m2), Tau.
!
          cff=rhoAir(i)*Cd*Wspeed
          Taux(i,j)=(cff*Uwind(i,j)+Taur*SIGN(1.0_r8,Uwind(i,j)))
          Taux(i,j)=Taux(i,j)*rmask(i,j)
          Tauy(i,j)=(cff*Vwind(i,j)+Taur*SIGN(1.0_r8,Vwind(i,j)))
          Tauy(i,j)=Tauy(i,j)*rmask(i,j)
          Taux_Ice(i,j) = rhoAir(i)*cdai(ng)*Uwind(i,j)*Wspeed
          Tauy_Ice(i,j) = rhoAir(i)*cdai(ng)*Vwind(i,j)*Wspeed
!  Correct ocean net short-wave radiation for ice cover and convert to
!  kinematic units
          Coef_Ice(i,j) = 0._r8
          RHS_Ice(i,j) = 0._r8
! Sensible heat flux over ice
          IF (ai(i,j,liold).gt.min_a(ng)) THEN
            sfc_temp = tis(i,j) + 273.16_r8
          ELSE
            sfc_temp = t(i,j,N(ng),nrhs,itemp) + 273.16_r8
          ENDIF
          qsh_i = rhoAir(i)*blk_Cpa*Ch*delW(i)*                         &
     &            (TairK(i)+0.0098*blk_ZT(ng)-sfc_temp)
          Coef_Ice(i,j) = Coef_Ice(i,j) +                               &
     &          rhoAir(i)*blk_Cpa*Ch*delW(i)
          RHS_Ice(i,j) = RHS_Ice(i,j) +                                 &
     &        rhoAir(i)*blk_Cpa*Ch*delW(i)*(TairK(i)+0.0098*blk_ZT(ng))
! Latent heat flux over ice
          le_i=2.834E+6_r8
! Qsati is saturation specific humidity over the ice from
! Parkinson and Washington (1979)
          slp = Pair(i,j)*100._r8 !Convert back to Pascal from millibars
          cff = 611._r8*                                                &
     &      10._r8**(9.5_r8*(sfc_temp-273.16_r8)/(sfc_temp-7.66_r8))
!         Qsati=eps*cff/(slp-(1._r8-.62197_r8)*cff)
          Qsati=0.62197_r8*cff/(slp-(1._r8-.62197_r8)*cff)
          dq_i=Q(i)-Qsati
          fqlat1=rhoAir(i)*Ce*delW(i)
          qlh_i = fqlat1*le_i*dq_i
          RHS_Ice(i,j) = RHS_Ice(i,j) + qlh_i
!
! Calculate the ice/snow albedo
! Ice and snow albedo is calculated from Ebert and Curry,1992b
          albi=0._r8
          albsn=0._r8
          sfwat_n = sfwat(i,j,liold)
          IF (ai(i,j,liold) .gt. min_a(ng)) THEN
            thki_n = hi(i,j,liold)/(ai(i,j,liold)+0.01_r8)
            thksn_n = hsn(i,j,liold)/(ai(i,j,liold)+0.01_r8)
            albi=0.082409_r8*LOG(thki_n)+0.485472_r8
            IF (thki_n.GE.1._r8) albi=0.07616_r8*thki_n+0.414492_r8
            IF (thki_n.GE.2._r8) albi=0.561632_r8
!approximated values for albsn(depends on COSZ, but small variation)
            albsn=0.83_r8
            IF (sfwat_n.GT.0._r8) THEN
              IF (thksn_n.GE.0.1_r8) THEN
                albsn=0.701009_r8
              ELSE
                albsn=albi+(0.701009_r8-albi)*thksn_n/0.1
              ENDIF
            ENDIF
            albs=albi
            IF (hsn(i,j,liold).GT.0._r8) albs=albsn
            IF (sfwat_n.GT.0._r8) albs=0.100737_r8                      &
     &         +0.518_r8*EXP(-8.1_r8 *sfwat_n-0.47_r8)                  &
     &         +0.341_r8*EXP(-31.8_r8*sfwat_n-0.94_r8)                  &
     &         +0.131_r8*EXP(-2.6_r8 *sfwat_n-3.82_r8)
!            ENDIF  !Removed 11/13/2012  !Unremoved 2/22/2013
          ELSE
            albs=alb_w
          ENDIF
!
! Short-wave radiation to ice
          SRad(i,j) = SRad(i,j) * 1.075269_r8
          qswi = (1._r8-albs)*SRad(i,j)
! Already in W/m**2
          RHS_Ice(i,j) = RHS_Ice(i,j) + qswi
! Net long-wave radiation from ice
          cff = Q(i)/(1.0_r8+Q(i))
          vap_p_i = slp*cff/(0.62197 + cff)
          qlwi = emmiss*StefBo*(TairK(i)**4)*                           &
     &         (0.39_r8-0.005_r8*SQRT(vap_p_i))*                        &
     &                (1.0_r8-0.65_r8*cloud(i,j)) +                     &
     &          4._r8*StefBo*emmiss*TairK(i)**3*(sfc_temp-TairK(i))
          Coef_Ice(i,j) = Coef_Ice(i,j) +                               &
     &           4._r8*StefBo*emmiss*TairK(i)**3
          RHS_Ice(i,j) = RHS_Ice(i,j) -                                 &
     &         emmiss*StefBo*(TairK(i)**4)*(                            &
     &    (0.39_r8-0.005_r8*SQRT(vap_p_i))*(1.0_r8-0.65_r8*cloud(i,j))  &
     &    - 4.0_r8)
!      if(i.eq.134.and.j.eq.54) then
!        write(*,*) 'vap_p_i,cloud(i,j),fac1, fac2 :'
!        write(*,*) vap_p_i,cloud(i,j), emmiss*StefBo*(TairK(i)**4)*    &
!     & (0.39_r8-0.005_r8*SQRT(vap_p_i))*(1._r8-0.65_r8*cloud(i,j)),    &
!     & 4._r8*StefBo*emmiss*TairK(i)**3*(sfc_temp-TairK(i))
!      endif
! Correct for Kelvin to Celsius in RHS
           RHS_Ice(i,j) = RHS_Ice(i,j) -                                &
     &            Coef_Ice(i,j)*273.16_r8
! Net heat flux from ice to atmosphere
           Qai(i,j) =-qsh_i - qlh_i - qswi + qlwi
!      if(i.eq.134.and.j.eq.54) then
!         write(*,*) 'TairK, sfc_temp = ',TairK(i),sfc_temp
!         write(*,*) 'i,j,qsh_i,qlh_i,qswi,qlwi,qai_n(i,j) :'
!         write(*,*) i,j,qsh_i,qlh_i,qswi,qlwi,qai_n(i,j)
!         write(*,*) 'rhoAir(i),blk_Cpa,Ch,delW(i),TairK(i),blk_ZT(ng),sfc_temp :'
!         write(*,*) rhoAir(i),blk_Cpa,Ch,delW(i),TairK(i),blk_ZT(ng),sfc_temp
!      endif
!      if(i.eq.49.and.j.eq.149) then
!         write(*,*) 'i,j,qsh_i,qlh_i,qswi,qlwi,qai_n(i,j) :'
!         write(*,*) i,j,qsh_i,qlh_i,qswi,qlwi,qai_n(i,j)
!      endif
        END DO
      END DO
!
!=======================================================================
!  Compute surface net heat flux and surface wind stress.
!=======================================================================
!
!  Compute kinematic, surface, net heat flux (degC m/s).  Notice that
!  the signs of latent and sensible fluxes are reversed because fluxes
!  calculated from the bulk formulations above are positive out of the
!  ocean.
!
!  For  option,  EVAP = LHeat (W/m2) / Hlv (J/kg) = kg/m2/s
!                       PREC = rain = kg/m2/s
!
!  To convert these rates to m/s divide by freshwater density, rhow.
!
!  Note that when the air is undersaturated in water vapor (Q < Qsea)
!  the model will evaporate and LHeat > 0:
!
!                   LHeat positive out of the ocean
!                    evap positive out of the ocean
!
!  Note that if evaporating, the salt flux is positive
!        and if     raining, the salt flux is negative
!
!  Note that fresh water flux is positive out of the ocean and the
!  salt flux (stflx(isalt)) is positive into the ocean. It is converted
!  to (psu m/s) for stflx(isalt) in "set_vbc.F" or "ice_mk.h". The E-P
!  value is saved in variable EminusP for I/O purposes.
!
      cff=1.0_r8/rhow
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          lrflx(i,j)=LRad(i,j)*Hscale2
          lhflx(i,j)=-LHeat(i,j)*Hscale2
          shflx(i,j)=-SHeat(i,j)*Hscale2
          stflx(i,j,itemp)=(srflx(i,j)+lrflx(i,j)+                      &
     &                      lhflx(i,j)+shflx(i,j))
          evap(i,j)=LHeat(i,j)/Hlv(i,j)
          stflx(i,j,isalt)=cff*(evap(i,j)-rain(i,j))
          evap(i,j) = (1.0_r8-ai(i,j,linew))*evap(i,j)
          stflx(i,j,isalt) = (1.0_r8-ai(i,j,linew))*stflx(i,j,isalt)
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask(i,j)
          evap(i,j)=evap(i,j)*rmask(i,j)
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask(i,j)
! Limit shortwave radiation to be used in computing penetrative
! radiation by presence of ice and biology
! Limit shortwave radiation to be used in computing penetrative
! radiation by presence of ice alone
          ice_thick(i,j) = hi(i,j,1)/(ai(i,j,1)+eps)
          snow_thick(i,j) = hsn(i,j,linew)/(ai(i,j,linew)+eps)
          IF (ice_thick(i,j).le.0.1_r8) THEN
            cff1=ice_thick(i,j)*5._r8
          ELSE
            cff1=0.1_r8*5._r8+(ice_thick(i,j)-0.1_r8)*1._r8
          ENDIF
          cff1= cff1+snow_thick(i,j)*20._r8
          srflx(i,j)=(1.0_r8-ai(i,j,linew))*srflx(i,j)                 !&
!     &            +ai(i,j,linew)*srflx(i,j)*exp(-cff1)
! Removed (changed to 3.3 style) 2/22/2013
          qao_n(i,j) = -stflx(i,j,itemp)*Hscale
! Net precipitation - evaporation in m/s
          IF (Tair(i,j).ge.0.0) THEN
            snow_n(i,j) = 0.0_r8
          ELSE
            snow_n(i,j) = rain(i,j)/rhosnow_dry(ng)
          ENDIF
          coef_ice_heat(i,j) = Coef_Ice(i,j)
          rhs_ice_heat(i,j) = RHS_Ice(i,j)
          qai_n(i,j) = Qai(i,j)
          EminusP(i,j)=stflx(i,j,isalt)
        END DO
      END DO
!
!  Compute kinematic, surface wind stress (m2/s2).
!
      cff=0.5_r8/rho0
      DO j=JstrR,JendR
        DO i=Istr,IendR
          sustr(i,j)=cff*(Taux(i-1,j)+Taux(i,j))
          sustr(i,j)=sustr(i,j)*umask(i,j)
          sustr_aw(i,j) = sustr(i,j)
          tau_aix_n(i,j) = Taux_Ice(i,j)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          svstr(i,j)=cff*(Tauy(i,j-1)+Tauy(i,j))
          svstr(i,j)=svstr(i,j)*vmask(i,j)
          svstr_aw(i,j) = svstr(i,j)
          tau_aiy_n(i,j) = Tauy_Ice(i,j)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          lrflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          lhflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          shflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          stflx(:,:,itemp))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        srflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        qao_n)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    qai_n)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    coef_ice_heat)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    rhs_ice_heat)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    snow_n)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    tau_aix_n)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    tau_aiy_n)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          evap)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          EminusP)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          stflx(:,:,isalt))
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        sustr_aw)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        svstr_aw)
      END IF
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    lrflx, lhflx, shflx,                          &
     &                    stflx(:,:,itemp))
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    srflx, qao_n, tau_aix_n, tau_aiy_n )
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    qai_n, coef_ice_heat, rhs_ice_heat, snow_n )
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    evap, EminusP,                                &
     &                    stflx(:,:,isalt))
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr, svstr)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr_aw, svstr_aw)
      RETURN
      END SUBROUTINE bulk_flux_tile
      FUNCTION bulk_psiu (ZoL, pi)
!
!=======================================================================
!                                                                      !
!  This function evaluates the stability function for  wind speed      !
!  by matching Kansas  and free convection forms.  The convective      !
!  form follows Fairall et al. (1996) with profile constants from      !
!  Grachev et al. (2000) BLM.  The  stable  form is from Beljaars      !
!  and Holtslag (1991).                                                !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
!  Function result
!
      real(r8) :: bulk_psiu
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: ZoL, pi
!
!  Local variable declarations.
!
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8
      real(r8) :: Fw, cff, psic, psik, x, y
!
!-----------------------------------------------------------------------
!  Compute stability function, PSI.
!-----------------------------------------------------------------------
!
!  Unstable conditions.
!
      IF (ZoL.lt.0.0_r8) THEN
        x=(1.0_r8-15.0_r8*ZoL)**0.25_r8
        psik=2.0_r8*LOG(0.5_r8*(1.0_r8+x))+                             &
     &       LOG(0.5_r8*(1.0_r8+x*x))-                                  &
     &       2.0_r8*ATAN(x)+0.5_r8*pi
!
!  For very unstable conditions, use free-convection (Fairall).
!
        cff=SQRT(3.0_r8)
        y=(1.0_r8-10.15_r8*ZoL)**r3
        psic=1.5_r8*LOG(r3*(1.0_r8+y+y*y))-                             &
     &       cff*ATAN((1.0_r8+2.0_r8*y)/cff)+pi/cff
!
!  Match Kansas and free-convection forms with weighting Fw.
!
        cff=ZoL*ZoL
        Fw=cff/(1.0_r8+cff)
        bulk_psiu=(1.0_r8-Fw)*psik+Fw*psic
!
!  Stable conditions.
!
      ELSE
        cff=MIN(50.0_r8,0.35_r8*ZoL)
        bulk_psiu=-((1.0_r8+ZoL)+0.6667_r8*(ZoL-14.28_r8)/              &
     &            EXP(cff)+8.525_r8)
      END IF
      RETURN
      END FUNCTION bulk_psiu
      FUNCTION bulk_psit (ZoL, pi)
!
!=======================================================================
!                                                                      !
!  This function evaluates the  stability function  for moisture and   !
!  heat by matching Kansas and free convection forms. The convective   !
!  form follows Fairall et al. (1996) with  profile  constants  from   !
!  Grachev et al. (2000) BLM.  The stable form is from  Beljaars and   !
!  and Holtslag (1991).                                                !
!
!=======================================================================
!                                                                      !
!
      USE mod_kinds
!
!  Function result
!
      real(r8) :: bulk_psit
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: ZoL, pi
!
!  Local variable declarations.
!
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8
      real(r8) :: Fw, cff, psic, psik, x, y
!
!-----------------------------------------------------------------------
!  Compute stability function, PSI.
!-----------------------------------------------------------------------
!
!  Unstable conditions.
!
      IF (ZoL.lt.0.0_r8) THEN
        x=(1.0_r8-15.0_r8*ZoL)**0.5_r8
        psik=2.0_r8*LOG(0.5_r8*(1.0_r8+x))
!
!  For very unstable conditions, use free-convection (Fairall).
!
        cff=SQRT(3.0_r8)
        y=(1.0_r8-34.15_r8*ZoL)**r3
        psic=1.5_r8*LOG(r3*(1.0_r8+y+y*y))-                             &
     &       cff*ATAN((1.0_r8+2.0_r8*y)/cff)+pi/cff
!
!  Match Kansas and free-convection forms with weighting Fw.
!
        cff=ZoL*ZoL
        Fw=cff/(1.0_r8+cff)
        bulk_psit=(1.0_r8-Fw)*psik+Fw*psic
!
!  Stable conditions.
!
      ELSE
        cff=MIN(50.0_r8,0.35_r8*ZoL)
        bulk_psit=-((1.0_r8+2.0_r8*ZoL)**1.5_r8+                        &
     &            0.6667_r8*(ZoL-14.28_r8)/EXP(cff)+8.525_r8)
      END IF
      RETURN
      END FUNCTION bulk_psit
      END MODULE bulk_flux_mod
