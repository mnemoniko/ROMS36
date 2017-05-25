      MODULE ice_elastic_mod
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine timesteps the ice momentum equations.
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC ice_elastic
      CONTAINS
      SUBROUTINE ice_elastic (ng, tile)
      USE mod_param
      USE mod_grid
      USE mod_ice
      USE mod_stepping
      integer, intent(in) :: ng, tile
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
      CALL wclock_on (ng, iNLM, 55)
!
      CALL ice_elastic_tile (ng, tile,                                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold(ng), liuol(ng), liunw(ng),            &
     &                      lieol(ng), lienw(ng),                       &
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
     &                      GRID(ng) % zice,                            &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      GRID(ng) % om_u,                            &
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % on_v,                            &
     &                      GRID(ng) % f,                               &
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      ICE(ng) % uwater,                           &
     &                      ICE(ng) % vwater,                           &
     &                      ICE(ng) % sealev,                           &
     &                      ICE(ng) % uie,                              &
     &                      ICE(ng) % vie,                              &
     &                      ICE(ng) % ai,                               &
     &                      ICE(ng) % hi,                               &
     &                      ICE(ng) % pice,                             &
     &                      ICE(ng) % zetai,                            &
     &                      ICE(ng) % eta,                              &
     &                      ICE(ng) % sig11,                            &
     &                      ICE(ng) % sig22,                            &
     &                      ICE(ng) % sig12,                            &
     &                      ICE(ng) % tauaiu,                           &
     &                      ICE(ng) % tauaiv,                           &
     &                      ICE(ng) % chu_iw                            &
     &                      )
      CALL wclock_off (ng, iNLM, 55)
      RETURN
      END SUBROUTINE ice_elastic
!
!***********************************************************************
      SUBROUTINE ice_elastic_tile (ng, tile,                            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        liold, liuol, liunw, lieol, lienw,        &
     &                        rmask, umask, vmask,                      &
     &                        zice,                                     &
     &                        pm, pn, om_u, on_u, om_v, on_v,           &
     &                        f,                                        &
     &                        ui, vi, uwater, vwater, sealev,           &
     &                        uie, vie,                                 &
     &                        ai, hi, pice, zetai, eta,                 &
     &                        sig11, sig22, sig12,                      &
     &                        tauaiu, tauaiv,                           &
     &                        chu_iw) 
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE uibc_mod, ONLY : uibc_tile
      USE vibc_mod, ONLY : vibc_tile
!
      USE exchange_2d_mod
      USE mp_exchange_mod, ONLY : mp_exchange2d
!
      implicit none
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile                    
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liold, liuol, liunw, lieol, lienw
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: zice(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: f(LBi:,LBj:)
      real(r8), intent(out) :: ui(LBi:,LBj:,:)
      real(r8), intent(out) :: vi(LBi:,LBj:,:)
      real(r8), intent(in) :: uwater(LBi:,LBj:)
      real(r8), intent(in) :: vwater(LBi:,LBj:)
      real(r8), intent(in) :: sealev(LBi:,LBj:)
      real(r8), intent(inout) :: uie(LBi:,LBj:,:)
      real(r8), intent(inout) :: vie(LBi:,LBj:,:)
      real(r8), intent(in) :: ai(LBi:,LBj:,:)
      real(r8), intent(in) :: hi(LBi:,LBj:,:)
      real(r8), intent(in) :: pice(LBi:,LBj:)
      real(r8), intent(in) :: zetai(LBi:,LBj:)
      real(r8), intent(in) :: eta(LBi:,LBj:)
      real(r8), intent(in) :: sig11(LBi:,LBj:,:)
      real(r8), intent(in) :: sig22(LBi:,LBj:,:)
      real(r8), intent(in) :: sig12(LBi:,LBj:,:)
      real(r8), intent(in) :: tauaiu(LBi:,LBj:)
      real(r8), intent(in) :: tauaiv(LBi:,LBj:)
      real(r8), intent(in) :: chu_iw(LBi:,LBj:)
! Local variable definitions
!
      integer :: i, j
      integer :: i1, j1, iterin
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps11
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps22
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps12
      real(r8) :: cff
      real(r8) :: f1
      real(r8) :: f2
      real(r8) :: s1
      real(r8) :: s2
      real(r8) :: ua
      real(r8) :: va
      real(r8) :: fakt
      real(r8) :: uforce
      real(r8) :: vforce
      real(r8) :: alfa
      real(r8) :: beta
      real(r8) :: flx
      real(r8) :: dimax 
      real(r8) :: masu
      real(r8) :: masv 
      real(r8) :: chux
      real(r8) :: chuy
      real(r8) :: auf
      real(r8) :: avf
      real(r8) :: dsum
      real(r8) :: pmu
      real(r8) :: pnu
      real(r8) :: pmv
      real(r8) :: pnv
      real(r8) :: cosstang
      real(r8) :: sinstang
      real(r8) :: mfu11
      real(r8) :: mfu21
      real(r8) :: mfu12
      real(r8) :: mfu22
      real(r8) :: mfv11
      real(r8) :: mfv21
      real(r8) :: mfv12
      real(r8) :: mfv22
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
!----------------
!
      cosstang = COS(deg2rad*stressang(ng))
      sinstang = SIN(deg2rad*stressang(ng))
! 
! Initialize fast-step velocities if beginning of cycle
!
! Update ice velocity
!
      IF (ievp(ng).eq.1) THEN
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            uie(i,j,1) = ui(i,j,liuol)
            uie(i,j,2) = ui(i,j,liuol)
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
             vie(i,j,1) = vi(i,j,liuol)
             vie(i,j,2) = vi(i,j,liuol)
          END DO
        END DO
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
             uie(1,j,1) = ui(1,j,liuol)
             uie(1,j,2) = ui(1,j,liuol)
          END DO
          DO j=JstrV,Jend
             vie(0,j,1) = vi(0,j,liuol)
             vie(0,j,2) = vi(0,j,liuol)
          END DO
        ENDIF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
             uie(Lm(ng)+1,j,1) = ui(Lm(ng)+1,j,liuol)
             uie(Lm(ng)+1,j,2) = ui(Lm(ng)+1,j,liuol)
          END DO
          DO j=JstrV,Jend
             vie(Lm(ng)+1,j,1) = vi(Lm(ng)+1,j,liuol)
             vie(Lm(ng)+1,j,2) = vi(Lm(ng)+1,j,liuol)
          END DO
        ENDIF
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=IstrU,Iend
             uie(i,0,1) = ui(i,0,liuol)
             uie(i,0,2) = ui(i,0,liuol)
          END DO
          DO i=Istr,Iend
             vie(i,1,1) = vi(i,1,liuol)
             vie(i,1,2) = vi(i,1,liuol)
          END DO
        ENDIF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=IstrU,Iend
             uie(i,Mm(ng)+1,1) = ui(i,Mm(ng)+1,liuol)
             uie(i,Mm(ng)+1,2) = ui(i,Mm(ng)+1,liuol)
          END DO
          DO i=Istr,Iend
             vie(i,Mm(ng)+1,1) = vi(i,Mm(ng)+1,liuol)
             vie(i,Mm(ng)+1,2) = vi(i,Mm(ng)+1,liuol)
          END DO
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          uie(1,0,1) = ui(1,0,liuol)
          uie(1,0,2) = ui(1,0,liuol)
          vie(0,1,1) = vi(0,1,liuol)
          vie(0,1,2) = vi(0,1,liuol)
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          uie(Lm(ng)+1,0,1) = ui(Lm(ng)+1,0,liuol)
          uie(Lm(ng)+1,0,2) = ui(Lm(ng)+1,0,liuol)
          vi(Lm(ng)+1,1,1) = vie(Lm(ng)+1,1,lienw)
          vi(Lm(ng)+1,1,2) = vie(Lm(ng)+1,1,lienw)
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          uie(1,Mm(ng)+1,1) = ui(1,Mm(ng)+1,liuol)
          uie(1,Mm(ng)+1,2) = ui(1,Mm(ng)+1,liuol)
          vie(0,Mm(ng)+1,1) = vi(0,Mm(ng)+1,liuol)
          vie(0,Mm(ng)+1,2) = vi(0,Mm(ng)+1,liuol)
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          uie(Lm(ng)+1,Mm(ng)+1,1) = ui(Lm(ng)+1,Mm(ng)+1,liuol)
          uie(Lm(ng)+1,Mm(ng)+1,2) = ui(Lm(ng)+1,Mm(ng)+1,liuol)
          vie(Lm(ng)+1,Mm(ng)+1,1) = vi(Lm(ng)+1,Mm(ng)+1,liuol)
          vie(Lm(ng)+1,Mm(ng)+1,2) = vi(Lm(ng)+1,Mm(ng)+1,liuol)
        END IF
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    uie(:,:,1), vie(:,:,1))
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    uie(:,:,2), vie(:,:,2))
      ENDIF
! ...................................
!
! *** momentum equation 
!
! u-eqn
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          uforce = 0._r8
          chux = 0.5_r8*(chu_iw(i-1,j)+chu_iw(i,j))
          masu = 0.5_r8*(hi(i,j,liold)+hi(i-1,j,liold))
          masu = max(masu,0.1_r8)
          masu = masu*rhoice(ng)
          auf = max(0.5_r8*(ai(i-1,j,liold)+ai(i,j,liold)),0.01_r8)
          pmu = 0.5_r8*(pm(i,j) + pm(i-1,j)) 
          pnu = 0.5_r8*(pn(i,j) + pn(i-1,j)) 
!
! *** forces from ice rheology x-direction 
          s1 = (sig11(i,j,lienw)-sig11(i-1,j,lienw))*pmu 
          IF (umask(i,j).gt.0.0_r8.and.umask(i,j+1).lt.1.0_r8) THEN 
            f1 = 0.5_r8*(sig12(i,j,lienw)+sig12(i-1,j,lienw))
          ELSE IF (zice(i-1,j).eq.0.0_r8.and.zice(i,j).eq.0.0_r8.and.   &
     &             zice(i-1,j+1)+zice(i,j+1).ne.0.0_r8) THEN
            f1 = 0.5_r8*(sig12(i,j,lienw)+sig12(i-1,j,lienw))
          ELSE
            f1 = 0.25_r8*(sig12(i,j,lienw)+sig12(i,j+1,lienw)           &
     &           + sig12(i-1,j+1,lienw)+sig12(i-1,j,lienw))
          END IF
          IF (umask(i,j).gt.0.0_r8.and.umask(i,j-1).lt.1.0_r8) THEN
             f2 = 0.5_r8*(sig12(i,j,lienw)+sig12(i-1,j,lienw))
          ELSE IF (zice(i-1,j).eq.0.0_r8.and.zice(i,j).eq.0.0_r8.and.   &
     &             zice(i-1,j-1)+zice(i,j-1).ne.0.0_r8) then
             f2 = 0.5_r8*(sig12(i,j,lienw)+sig12(i-1,j,lienw))
          ELSE
             f2 = 0.25_r8*(sig12(i,j,lienw)+sig12(i-1,j,lienw)          &
     &           + sig12(i-1,j-1,lienw)+sig12(i,j-1,lienw))
          END IF
!
           s2 = (f1-f2)*pnu
           uforce = (s1 + s2)*om_u(i,j)*on_u(i,j)
! *** wind force 
           uforce = uforce + rhoice(ng)*tauaiu(i,j)*om_u(i,j)*on_u(i,j)
! *** pressure from tilting ocean surface 
           uforce = uforce - g*(masu)*on_u(i,j)                         &
     &            *(sealev(i,j)-sealev(i-1,j))
           fakt = 0.0_r8
           dsum = vmask(i-1,j)+vmask(i,j)+vmask(i-1,j+1)+vmask(i,j+1)
           if (dsum.gt.0.0_r8) fakt = 1.0_r8/dsum
           mfv11 = 0.5_r8*(hi(i-1,j-1,liold)*f(i-1,j-1)+                &
     &                     hi(i-1,j,liold)*f(i-1,j))*                   &
     &                      vie(i-1,j,lieol)*vmask(i-1,j)
           mfv21 = 0.5_r8*(hi(i,j-1,liold)*f(i,j-1)+                    &
     &                     hi(i,j,liold)*f(i,j))*                       &
     &                      vie(i,j,lieol)*vmask(i,j)
           mfv12 = 0.5_r8*(hi(i-1,j,liold)*f(i-1,j)+                    &
     &                     hi(i-1,j+1,liold)*f(i-1,j+1))*               &
     &                      vie(i-1,j+1,lieol)*vmask(i-1,j+1)
           mfv22 = 0.5_r8*(hi(i,j,liold)*f(i,j)+                        &
     &                     hi(i,j+1,liold)*f(i,j+1))*                   &
     &                      vie(i,j+1,lieol)*vmask(i,j+1)
! *** coriolis force u
          uforce = uforce + fakt*rhoice(ng)*om_u(i,j)*on_u(i,j)*        &
     &                            (mfv11 + mfv21 + mfv12 + mfv22)
!
! *** stress from ocean current 
          uforce = uforce/(om_u(i,j)*on_u(i,j)) +                       &
     &                         auf*chux*rho0*uwater(i,j)
!
          alfa = masu + dte(ng)*auf*rho0*chux
!
!  solving the momentum equation for u
          uie(i,j,lienw) = (masu*uie(i,j,lieol) +                       &
     &                       dte(ng)*uforce)/alfa 
!
          uie(i,j,lienw) = umask(i,j)*uie(i,j,lienw)
          IF(zice(i-1,j)+zice(i,j).ne.0.0_r8) THEN
              uie(i,j,lienw) = 0.0_r8
          END IF
        END DO
      END DO
!
!
! *** momentum equation 
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
          vforce = 0._r8
          masv = 0.5_r8*(hi(i,j,liold)+hi(i,j-1,liold))
          masv = max(masv,0.1_r8)
          masv = masv*rhoice(ng) 
          avf = max(0.5_r8*(ai(i,j-1,liold)+ai(i,j,liold)),0.01_r8)
          chuy = 0.5_r8*(chu_iw(i,j-1)+chu_iw(i,j))
          pmv = 0.5_r8*(pm(i,j) + pm(i,j-1)) 
          pnv = 0.5_r8*(pn(i,j) + pn(i,j-1)) 
!
! *** forces from ice rheology y-direction 
          s1 = (sig22(i,j,lienw) - sig22(i,j-1,lienw))*pnv
          IF (vmask(i,j).gt.0.0_r8.and.vmask(i+1,j).lt.1.0_r8) THEN
             f1 = 0.5_r8*(sig12(i,j,lienw)+sig12(i,j-1,lienw))
          ELSE IF (zice(i,j-1).eq.0.0_r8.and.zice(i,j).eq.0.0_r8.and.   &
     &             zice(i+1,j-1)+zice(i+1,j).ne.0.0_r8) then
            f1 = 0.5_r8*(sig12(i,j,lienw)+sig12(i,j-1,lienw))
          ELSE
            f1 = 0.25_r8*(sig12(i,j,lienw)+sig12(i+1,j,lienw)           &
     &           + sig12(i+1,j-1,lienw)+sig12(i,j-1,lienw))
          END IF
          IF (vmask(i,j).gt.0.0_r8.and.vmask(i-1,j).lt.1.0_r8 ) THEN
              f2 = 0.5_r8*(sig12(i,j,lienw)+sig12(i,j-1,lienw))
          ELSE IF (zice(i,j-1).eq.0.0_r8.and.zice(i,j).eq.0.0_r8.and.   &
     &          zice(i-1,j-1)+zice(i-1,j).ne.0.0_r8) then
            f2 = 0.5_r8*(sig12(i,j,lienw)+sig12(i,j-1,lienw))
          ELSE
            f2 = 0.25_r8*(sig12(i,j,lienw)+sig12(i-1,j,lienw)           &
     &           + sig12(i-1,j-1,lienw)+sig12(i,j-1,lienw))
          END IF
!
          s2 = (f1-f2)*pmv 
          vforce = (s1+s2)*om_v(i,j)*on_v(i,j)
! *** wind force
          vforce = vforce + rhoice(ng)*tauaiv(i,j)*om_v(i,j)*on_v(i,j)
! *** pressure from tilting ocean surface 
          vforce = vforce - g*(masv)*om_v(i,j)                          &
     &            *(sealev(i,j)-sealev(i,j-1))
!
          fakt = 0.0_r8
          dsum = umask(i,j-1)+umask(i+1,j-1)+umask(i,j)+umask(i+1,j)
          IF (dsum.gt.0.0_r8) fakt = 1.0_r8/dsum
          mfu11 = 0.5_r8*(hi(i-1,j-1,liold)*f(i-1,j-1)+                 &
     &                     hi(i,j-1,liold)*f(i,j-1))*                   &
     &                      uie(i,j-1,lieol)*umask(i,j-1)
          mfu21 = 0.5_r8*(hi(i,j-1,liold)*f(i,j-1)+                     &
     &                     hi(i+1,j-1,liold)*f(i+1,j-1))*               &
     &                      uie(i+1,j-1,lieol)*umask(i+1,j-1)
          mfu12 = 0.5_r8*(hi(i-1,j,liold)*f(i-1,j)+                     &
     &                     hi(i,j,liold)*f(i,j))*                       &
     &                      uie(i,j,lieol)*umask(i,j)
          mfu22 = 0.5_r8*(hi(i,j,liold)*f(i,j)+                         &
     &                     hi(i+1,j,liold)*f(i+1,j))*                   &
     &                      uie(i+1,j,lieol)*umask(i+1,j)
! *** coriolis force v
          vforce = vforce - fakt*rhoice(ng)*om_v(i,j)*on_v(i,j)*        &
     &                            (mfu11 + mfu21 + mfu12 + mfu22)
!
! *** stress from ocean current 
          vforce = vforce/(om_v(i,j)*on_v(i,j)) +                       &
     &                         avf*chuy*rho0*vwater(i,j)
!
          alfa = masv + dte(ng)*avf*rho0*chuy
!  solving the momentum equation for v
          vie(i,j,lienw) = (masv*vie(i,j,lieol) +                       &
     &                       dte(ng)*vforce)/alfa 
!
          vie(i,j,lienw) = vmask(i,j)*vie(i,j,lienw)
          IF (zice(i,j-1)+zice(i,j).ne.0.0_r8) THEN
            vie(i,j,lienw) = 0.0_r8
          ENDIF
        END DO 
      END DO 
      CALL uibc_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                IminS, ImaxS, JminS, JmaxS,                       &
     &                lieol, lienw, uie)
      CALL vibc_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                IminS, ImaxS, JminS, JmaxS,                       &
     &                lieol, lienw, vie)
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          uie(:,:,lienw))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vie(:,:,lienw))
      END IF
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    uie(:,:,lienw), vie(:,:,lienw))
!
! Update ice velocity
!
      IF(ievp(ng).eq.nevp(ng)) THEN
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            ui(i,j,liunw) = uie(i,j,lienw)
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            vi(i,j,liunw) = vie(i,j,lienw)
          END DO
        END DO
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            ui(1,j,liunw) = uie(1,j,lienw)
          END DO
          DO j=JstrV,Jend
            vi(0,j,liunw) = vie(0,j,lienw)
          END DO
        ENDIF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            ui(Lm(ng)+1,j,liunw) = uie(Lm(ng)+1,j,lienw)
          END DO
          DO j=JstrV,Jend
            vi(Lm(ng)+1,j,liunw) = vie(Lm(ng)+1,j,lienw)
          END DO
        ENDIF
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=IstrU,Iend
            ui(i,0,liunw) = uie(i,0,lienw)
          END DO
          DO i=Istr,Iend
            vi(i,1,liunw) = vie(i,1,lienw)
          END DO
        ENDIF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=IstrU,Iend
            ui(i,Mm(ng)+1,liunw) = uie(i,Mm(ng)+1,lienw)
          END DO
          DO i=Istr,Iend
            vi(i,Mm(ng)+1,liunw) = vie(i,Mm(ng)+1,lienw)
          END DO
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          ui(1,0,liunw) = uie(1,0,lienw)
          vi(0,1,liunw) = vie(0,1,lienw)
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          ui(Lm(ng)+1,0,liunw) = uie(Lm(ng)+1,0,lienw)
          vi(Lm(ng)+1,1,liunw) = vie(Lm(ng)+1,1,lienw)
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          ui(1,Mm(ng)+1,liunw) = uie(1,Mm(ng)+1,lienw)
          vi(0,Mm(ng)+1,liunw) = vie(0,Mm(ng)+1,lienw)
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          ui(Lm(ng)+1,Mm(ng)+1,liunw) = uie(Lm(ng)+1,Mm(ng)+1,lienw)
          vi(Lm(ng)+1,Mm(ng)+1,liunw) = vie(Lm(ng)+1,Mm(ng)+1,lienw)
        END IF
!      CALL uibc_tile (ng, tile,                                        &
!     &                LBi, UBi, LBj, UBj,                              &
!     &                IminS, ImaxS, JminS, JmaxS,                      &
!     &                liuol, liunw, ui)
!      CALL vibc_tile (ng, tile,                                        &
!     &                LBi, UBi, LBj, UBj,                              &
!     &                IminS, ImaxS, JminS, JmaxS,                      &
!     &                liuol, liunw, vi)
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ui(:,:,liunw))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vi(:,:,liunw))
      END IF
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ui(:,:,liunw), vi(:,:,liunw))
      ENDIF
!
      RETURN
      END SUBROUTINE ice_elastic_tile
      END MODULE ice_elastic_mod
