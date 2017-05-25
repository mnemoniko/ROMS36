      MODULE ice_frazil_mod
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the frazil ice growth in the water when the
!  water temperature gets below freezing. It adjusts the water
!  temperature and salinity accordingly.
!
!  Reference: Steele et al. (1989). JPO, 19, 139-147.
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC ice_frazil
      CONTAINS
      SUBROUTINE ice_frazil (ng, tile)
      USE mod_param
      USE mod_grid
      USE mod_ocean
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
      CALL wclock_on (ng, iNLM, 51)
!
      CALL ice_frazil_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nnew(ng),                                   &
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % z_r,                             &
     &                      OCEAN(ng) % rho,                            &
     &                      OCEAN(ng) % t,                              &
     &                      ICE(ng) % wfr)
      CALL wclock_off (ng, iNLM, 51)
      RETURN
      END SUBROUTINE ice_frazil
!
!***********************************************************************
      subroutine ice_frazil_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nnew,                                 &
     &                            rmask,                                &
     &                            Hz, z_r, rho, t, wfr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod, ONLY : bc_r2d_tile
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
      USE mp_exchange_mod
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nnew
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(out) :: wfr(LBi:,LBj:)
!
! Local variable definitions
!
      integer :: i, j, k, itrc
      real(r8), parameter :: Lhat = 79.2_r8
      real(r8), parameter :: r = 0.5_r8
      real(r8) :: t_freeze
      real(r8) :: s1
      real(r8) :: z1
      real(r8) :: gamma_k
      real(r8) :: t_fr
!  Inline functions
!  Freezing temperature (Gill, 1982)
!     t_freeze(s1,z1) = -0.0575*s1 + 1.710523d-3*sqrt(s1)**3
!    &       - 2.154996d-4*s1*s1 + 0.000753*z1
!  Freezing temperature (Steele et al. 1989)
      t_freeze(s1,z1) = -0.0543*s1 + 0.000759*z1
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
      DO j=Jstr,Jend
        DO i=Istr,Iend
          wfr(i,j) = 0.0_r8
        ENDDO
      ENDDO
      DO j=Jstr,Jend
        DO i=Istr,Iend
          DO k=1,N(ng)
            IF (rmask(i,j) .ne. 0.0_r8) THEN
              t_fr = t_freeze(t(i,j,k,nnew,isalt),z_r(i,j,k))
              IF (t(i,j,k,nnew,itemp) .lt. t_fr) THEN
                gamma_k = (t_fr - t(i,j,k,nnew,itemp)) /                &
     &                     (Lhat + t(i,j,k,nnew,itemp)*(1.0_r8 - r)     &
     &                         + 0.0543_r8 * t(i,j,k,nnew,isalt))
                IF (gamma_k .lt. 0.0_r8) THEN
                  print *, 'trouble in ice_frazil1', i, j, k,           &
     &             t(i,j,k,nnew,itemp), t(i,j,k,nnew,isalt),            &
     &             t_fr, wfr(i,j), gamma_k, Hz(i,j,k)
                END IF
                wfr(i,j) = wfr(i,j) + gamma_k * Hz(i,j,k) *             &
     &                    (1000.0_r8 + rho(i,j,k) ) / rhoice(ng)
                t(i,j,k,nnew,itemp) = t(i,j,k,nnew,itemp) + gamma_k *   &
     &                 (Lhat + t(i,j,k,nnew,itemp)*(1.0_r8 - r))
                t(i,j,k,nnew,isalt) = t(i,j,k,nnew,isalt) *             &
     &                                  (1.0_r8 + gamma_k)
              END IF
            END IF
          END DO
          wfr(i,j) = wfr(i,j)/dt(ng)
          IF (wfr(i,j) .lt. 0.0_r8) THEN
            print *, 'trouble in ice_frazil', i, j,                     &
     &         t(i,j,N(ng),nnew,itemp), t(i,j,N(ng),nnew,isalt),        &
     &         wfr(i,j), gamma_k, Hz(i,j,N(ng))
          END IF
        END DO
      END DO
      CALL bc_r2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          wfr)
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    wfr)
      RETURN
      END SUBROUTINE ice_frazil_tile
      END MODULE ice_frazil_mod
