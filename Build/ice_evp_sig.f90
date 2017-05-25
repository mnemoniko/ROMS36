      MODULE ice_evp_sig_mod
!
!================================================ W. Paul Budgell ======
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine time steps the EVP stresses                            !
!  term.                                                               !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC ice_evp_sig
      CONTAINS
      SUBROUTINE ice_evp_sig (ng, tile)
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
      CALL ice_evp_sig_tile (ng, tile,                                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold(ng), lieol(ng), lienw(ng),            &
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      ICE(ng) % uie,                              &
     &                      ICE(ng) % vie,                              &
     &                      ICE(ng) % hi,                               &
     &                      ICE(ng) % pice,                             &
     &                      ICE(ng) % zetai,                            &
     &                      ICE(ng) % eta,                              &
     &                      ICE(ng) % sig11,                            &
     &                      ICE(ng) % sig22,                            &
     &                      ICE(ng) % sig12                             &
     &                      )
      CALL wclock_off (ng, iNLM, 55)
      RETURN
      END SUBROUTINE ice_evp_sig
!
!***********************************************************************
      SUBROUTINE ice_evp_sig_tile (ng, tile,                            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        liold, lieol, lienw,                      &
     &                        rmask,                                    &
     &                        pm, pn,                                   &
     &                        uie, vie,                                 &
     &                        hi, pice, zetai, eta,                     &
     &                        sig11, sig22, sig12                       &
     &                        ) 
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_boundary
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE mp_exchange_mod, ONLY : mp_exchange2d
!
      USE i2d_bc_mod, ONLY : i2d_bc_tile
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liold, lieol, lienw
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: uie(LBi:,LBj:,:)
      real(r8), intent(in) :: vie(LBi:,LBj:,:)
      real(r8), intent(in) :: hi(LBi:,LBj:,:)
      real(r8), intent(in) :: pice(LBi:,LBj:)
      real(r8), intent(in) :: zetai(LBi:,LBj:)
      real(r8), intent(in) :: eta(LBi:,LBj:)
      real(r8), intent(inout) :: sig11(LBi:,LBj:,:)
      real(r8), intent(inout) :: sig22(LBi:,LBj:,:)
      real(r8), intent(inout) :: sig12(LBi:,LBj:,:)
! Local variable definitions
!
      integer :: i, j
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps11
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps22
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps12
      real(r8) :: e
      real(r8) :: e0
      real(r8) :: ep
      real(r8) :: ee
      real(r8) :: ees 
      real(r8) :: gamma
      real(r8) :: f1
      real(r8) :: f2
      real(r8) :: f3 
      real(r8) :: s1
      real(r8) :: s2
      real(r8) :: alfa
      real(r8) :: beta
      real(r8) :: pmu
      real(r8) :: pnu
      real(r8) :: pmv
      real(r8) :: pnv
      real(r8) :: epx
      real(r8) :: epy
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
!.....initial value for Youngs modulus (between 0 and 1) 
!
      e0=0.25_r8
!
!...........stress tensor 
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          eps11(i,j) = (uie(i+1,j,lieol)-uie(i,j,lieol))*pm(i,j)
          eps22(i,j) = (vie(i,j+1,lieol)-vie(i,j,lieol))*pn(i,j)
          epx = 0.25_r8*( vie(i+1,j+1,lieol)+vie(i+1,j,lieol)           &
     &             - vie(i-1,j+1,lieol)-vie(i-1,j,lieol) )*pm(i,j)
          epy = 0.25_r8*( uie(i+1,j+1,lieol)+uie(i,j+1,lieol)           &
     &             - uie(i+1,j-1,lieol)-uie(i,j-1,lieol) )*pn(i,j)
          eps12(i,j) = 0.5_r8*(epx + epy)
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
!
          IF (hi(i,j,liold).gt.0.01_r8) THEN 
            e = 2.0_r8*e0*rhoice(ng)*hi(i,j,liold)/(pm(i,j)*dte(ng))**2 
            ep = e*pice(i,j)/(4.0_r8*zetai(i,j)+1.0E-8_r8) 
            ee = e/(2.0_r8*eta(i,j)+1.0E-8_r8) 
            ees= e*(eta(i,j)-zetai(i,j))/                               &
     &             (4.0_r8*eta(i,j)*zetai(i,j)+1.0E-8_r8) 
!
            alfa = 1.0_r8/dte(ng) + ee + ees 
            beta = ees 
            gamma = 1/dte(ng) + ee 
            f1 = e*eps11(i,j) - ep + sig11(i,j,lieol)/dte(ng) 
            f2 = e*eps22(i,j) - ep + sig22(i,j,lieol)/dte(ng) 
            f3 = e*eps12(i,j) + sig12(i,j,lieol)/dte(ng) 
            sig11(i,j,lienw) = (alfa*f1 - beta*f2)/(alfa**2 - beta**2) 
            sig12(i,j,lienw) = f3/gamma 
            sig22(i,j,lienw) = (alfa*f2 - beta*f1)/(alfa**2 - beta**2) 
!
          ELSE
! 
            sig11(i,j,lienw) = 2.0_r8*eta(i,j)*eps11(i,j)+              &
     &        (zetai(i,j)-eta(i,j))*                                    &
     &        (eps11(i,j)+eps22(i,j))  - pice(i,j)/2.0_r8
!     
            sig22(i,j,lienw) = 2.0_r8*eta(i,j)*eps22(i,j)+              &
     &        (zetai(i,j)-eta(i,j))*                                    &
     &        (eps11(i,j)+eps22(i,j))  - pice(i,j)/2.0_r8 
!
            sig12(i,j,lienw) = 2.0_r8*eta(i,j)*eps12(i,j)
!
          END IF 
        END DO 
      END DO
      CALL i2d_bc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  lieol, lienw,                                   &
     &                  BOUNDARY(ng)%sig11_west(LBj:UBj),               &
     &                  BOUNDARY(ng)%sig11_east(LBj:UBj),               &
     &                  BOUNDARY(ng)%sig11_north(LBi:UBi),              &
     &                  BOUNDARY(ng)%sig11_south(LBi:UBi),              &
     &                  uie, vie, sig11, LBC(:,isSig11,ng),0)
      CALL i2d_bc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  lieol, lienw,                                   &
     &                  BOUNDARY(ng)%sig22_west(LBj:UBj),               &
     &                  BOUNDARY(ng)%sig22_east(LBj:UBj),               &
     &                  BOUNDARY(ng)%sig22_north(LBi:UBi),              &
     &                  BOUNDARY(ng)%sig22_south(LBi:UBi),              &
     &                  uie, vie, sig22, LBC(:,isSig22,ng),0)
      CALL i2d_bc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  lieol, lienw,                                   &
     &                  BOUNDARY(ng)%sig12_west(LBj:UBj),               &
     &                  BOUNDARY(ng)%sig12_east(LBj:UBj),               &
     &                  BOUNDARY(ng)%sig12_north(LBi:UBi),              &
     &                  BOUNDARY(ng)%sig12_south(LBi:UBi),              &
     &                  uie, vie, sig12, LBC(:,isSig12,ng),0)
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sig11(:,:,lienw))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sig22(:,:,lienw))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sig12(:,:,lienw))
      END IF
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    sig11(:,:,lienw), sig22(:,:,lienw),           &
     &                    sig12(:,:,lienw))
!
      RETURN
      END SUBROUTINE ice_evp_sig_tile
      END MODULE ice_evp_sig_mod
