      MODULE ice_evp_mod
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the parameters for elastic-viscous-plastic    !
!  rheology.                                                           !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC ice_evp
      CONTAINS
      SUBROUTINE ice_evp (ng, tile)
      USE mod_param
      USE mod_grid
      USE mod_ice
      USE mod_stepping
!
      implicit none
!
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
      CALL wclock_on (ng, iNLM, 52)
!
      CALL ice_evp_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   liold(ng), lieol(ng),                          &
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % pm,                                 &
     &                   GRID(ng) % pn,                                 &
     &                   ICE(ng) % uie,                                 &
     &                   ICE(ng) % vie,                                 &
     &                   ICE(ng) % ai,                                  &
     &                   ICE(ng) % hi,                                  &
     &                   ICE(ng) % pice,                                &
     &                   ICE(ng) % zetai,                               &
     &                   ICE(ng) % eta                                  &
     &                   )
      CALL wclock_off (ng, iNLM, 52)
      RETURN
      END SUBROUTINE ice_evp
!
!***********************************************************************
      SUBROUTINE ice_evp_tile (ng, tile,                                &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        liold, lieol,                             &
     &                        rmask,                                    &
     &                        pm, pn, uie, vie,                         &
     &                        ai, hi, pice, zetai, eta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod, ONLY : bc_r2d_tile
      USE mp_exchange_mod, ONLY : mp_exchange2d
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liold, lieol
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: uie(LBi:,LBj:,:)
      real(r8), intent(in) :: vie(LBi:,LBj:,:)
      real(r8), intent(in) :: ai(LBi:,LBj:,:)
      real(r8), intent(in) :: hi(LBi:,LBj:,:)
      real(r8), intent(out) :: pice(LBi:,LBj:)
      real(r8), intent(out) :: zetai(LBi:,LBj:)
      real(r8), intent(out) :: eta(LBi:,LBj:)
! Local variable definitions
!
      integer :: i, j
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps11
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps22
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps12
      real(r8) :: eone
      real(r8) :: etwos
      real(r8) :: epx
      real(r8) :: epy
      real(r8) :: e2r
      real(r8) :: delta
      real(r8) :: zmax
      real(r8), parameter :: epso = 1.E-12_r8
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
! ------------------------------------------------------------
!
      e2r = 1.0_r8/(ellip_sq(ng))
!
! *** Compute strain rates
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
!
          eone=eps11(i,j)+eps22(i,j) 
          etwos=(eps11(i,j)-eps22(i,j))*(eps11(i,j)-eps22(i,j))+        &
     &         4.0_r8*eps12(i,j)*eps12(i,j)
! 
          delta=abs(eone**2+e2r*etwos)
          delta=max(sqrt(delta),epso)
          pice(i,j)=pstar(ng)*hi(i,j,liold)                             &
     &                *exp(-astren(ng)*(1.0_r8-ai(i,j,liold)))
          zetai(i,j)=pice(i,j)/(2.0_r8*delta)
          zmax = 2.5E+8_r8*pice(i,j)
          zetai(i,j)= min(zetai(i,j),zetamax(ng))
          zetai(i,j)= max(zetai(i,j),zetamin(ng))
          eta(i,j)=e2r*zetai(i,j)
        ENDDO
      ENDDO
      CALL bc_r2d_tile (ng, tile,                                       &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    pice)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    zetai)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    eta)
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    pice, zetai, eta)
      DO j=JstrR,JendR
        DO i=IstrR,IendR
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE ice_evp_tile
      END MODULE ice_evp_mod
