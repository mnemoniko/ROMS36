      MODULE ice_vbc_mod
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This module sets the ice-water and ice-air stresses for the
!  ice momentum equation.
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC ice_vbc
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_vbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean
      USE mod_ice
      USE mod_coupling
      USE mod_mixing
      USE mod_stepping
!
      implicit none
!
      integer, intent(in) :: ng, tile
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
      CALL wclock_on (ng, iNLM, 6)
      CALL ice_vbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng),                                      &
     &                   liold(ng), liuol(ng),                          &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   GRID(ng) % zice,                               &
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v,                                 &
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
     &                   OCEAN(ng) % rho,                               &
     &                   COUPLING(ng) % Zt_avg1,                        &
     &                   ICE(ng) % ai,                                  &
     &                   ICE(ng) % hi,                                  &
     &                   ICE(ng) % ui,                                  &
     &                   ICE(ng) % vi,                                  &
     &                   ICE(ng) % tauaiu,                              &
     &                   ICE(ng) % tauaiv,                              &
     &                   ICE(ng) % uwater,                              &
     &                   ICE(ng) % vwater,                              &
     &                   ICE(ng) % sealev,                              &
     &                   FORCES(ng) % sustr_aw,                         &
     &                   FORCES(ng) % svstr_aw,                         &
     &                   FORCES(ng) % tau_aix_n,                        &
     &                   FORCES(ng) % tau_aiy_n,                        &
     &                   ICE(ng) % utau_iw,                             &
     &                   ICE(ng) % chu_iw,                              &
     &                   ICE(ng) % spd_iw                               &
     &                   )
      CALL wclock_off (ng, iNLM, 6)
      RETURN
      END SUBROUTINE ice_vbc
!
!***********************************************************************
      SUBROUTINE ice_vbc_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,                                    &
     &                         liold, liuol,                            &
     &                         z_r, z_w,                                &
     &                         zice,                                    &
     &                         u, v,                                    &
     &                         sustr, svstr, rho,                       &
     &                         Zt_avg1,                                 &
     &                         ai, hi, ui, vi, tauaiu, tauaiv,          &
     &                         uwater, vwater, sealev,                  &
     &                         sustr_aw, svstr_aw,                      &
     &                         tau_aix_n, tau_aiy_n,                    &
     &                         utau_iw, chu_iw, spd_iw)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod
      USE mp_exchange_mod, ONLY : mp_exchange2d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
      integer, intent(in) :: liold, liuol
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(out) :: sustr(LBi:,LBj:)
      real(r8), intent(out) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(in) :: ai(LBi:,LBj:,:)
      real(r8), intent(in) :: hi(LBi:,LBj:,:)
      real(r8), intent(in) :: ui(LBi:,LBj:,:)
      real(r8), intent(in) :: vi(LBi:,LBj:,:)
      real(r8), intent(out) :: tauaiu(LBi:,LBj:)
      real(r8), intent(out) :: tauaiv(LBi:,LBj:)
      real(r8), intent(in) :: uwater(LBi:,LBj:)
      real(r8), intent(in) :: vwater(LBi:,LBj:)
      real(r8), intent(out) :: sealev(LBi:,LBj:)
      real(r8), intent(in) :: sustr_aw(LBi:,LBj:)
      real(r8), intent(in) :: svstr_aw(LBi:,LBj:)
      real(r8), intent(in) :: tau_aix_n(LBi:,LBj:)
      real(r8), intent(in) :: tau_aiy_n(LBi:,LBj:)
      real(r8), intent(out) :: utau_iw(LBi:,LBj:)
      real(r8), intent(inout) :: chu_iw(LBi:,LBj:)
      real(r8), intent(in) :: spd_iw(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: zice(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j
      integer :: k
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: spdiw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: chuiw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: utauiw
      real(r8) :: tauiwu
      real(r8) :: tauiwv
      real(r8) :: tauawu
      real(r8) :: tauawv
      real(r8) :: aix
      real(r8) :: aiy
      real(r8) :: spd
      real(r8) :: hix
      real(r8) :: hiy
      real(r8) :: chux
      real(r8) :: chuy
      real(r8) :: spdu
      real(r8) :: spdv
      real(r8) :: rhoO
      real(r8) :: dztop
      real(r8) :: thic
      real(r8) :: zdz0
      real(r8) :: z0
      real(r8), parameter :: kappa = 0.4_r8
      real(r8), parameter :: z0ii = 0.02_r8
      real(r8), parameter :: eps = 1.e-20
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
! *** Input from ocean model/data
      DO j=Jstr-1,Jend
        DO i=Istr-1,Iend
          rhoO = 1000._r8+rho(i,j,N(ng))
          spd = spd_iw(i,j)
          spd = max(spd,0.15_r8)
          utauiw(i,j) = spd
          chuiw(i,j) = cdiw(ng)*spd
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          rhoO = 1000._r8 + 0.5_r8*(rho(i,j,N(ng))+rho(i-1,j,N(ng)))
          aix = 0.5_r8*(ai(i,j,liold)+ai(i-1,j,liold))
          hix = 0.5_r8*(hi(i,j,liold)+hi(i-1,j,liold))
          chux = 0.5_r8*(chuiw(i,j)+chuiw(i-1,j))
          tauaiu(i,j) = 0.5_r8*aix*(tau_aix_n(i,j)+tau_aix_n(i-1,j))    &
     &                        /rhoice(ng)
          IF (zice(i,j).eq.0.0_r8) THEN
          sustr(i,j) = aix*chux*(ui(i,j,liuol)-uwater(i,j))             &
     &                 + (1.0_r8-aix)*sustr_aw(i,j)
          ENDIF
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          rhoO = 1000._r8 + 0.5_r8*(rho(i,j,N(ng))+rho(i,j-1,N(ng)))
          aiy = 0.5_r8*(ai(i,j,liold)+ai(i,j-1,liold))
          hiy = 0.5_r8*(hi(i,j,liold)+hi(i,j-1,liold))
          chuy = 0.5_r8*(chuiw(i,j)+chuiw(i,j-1))
          tauaiv(i,j) = 0.5_r8*aiy*(tau_aiy_n(i,j)+tau_aiy_n(i,j-1))    &
     &                        /rhoice(ng)
          IF (zice(i,j).eq.0.0_r8) THEN
          svstr(i,j) = aiy*chuy*(vi(i,j,liuol)-vwater(i,j))             &
     &                 + (1.0_r8-aiy)*svstr_aw(i,j)
          ENDIF
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=Istr,Iend
           sealev(i,j) = Zt_avg1(i,j)
           chu_iw(i,j) = chuiw(i,j)
           utau_iw(i,j) = utauiw(i,j)
        END DO
      END DO
!
!  Apply boundary conditions.
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sealev)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          utau_iw)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          chu_iw)
      CALL bc_u2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr)
      CALL bc_u2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tauaiu)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tauaiv)
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    sealev, utau_iw, chu_iw, sustr)
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    svstr, tauaiu, tauaiv)
      RETURN
      END SUBROUTINE ice_vbc_tile
      END MODULE ice_vbc_mod
