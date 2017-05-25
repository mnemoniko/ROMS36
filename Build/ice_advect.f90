      MODULE ice_advect_mod
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the advection of the ice tracer fields.       !
!  Currently, the only option is to use the MPDATA advection scheme.   !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC ice_advect
      CONTAINS
       SUBROUTINE ice_advect (ng, tile)
!
!*************************************************** W. Paul Budgell ***
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!************************************************** Hernan G. Arango ***
!                                                                      !
!  This subroutine performs advection of ice scalars using the         !
!  Smolarkiewicz second-order upwind scheme.                           !
!  Reference:                                                          !
!  Smolarkiewicz and Grabowski (1990)                                  !
!***********************************************************************
!
      USE mod_param
      implicit none
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
      CALL wclock_on (ng, iNLM, 56)
! ---------------------------------------------------------------------
!  Advect the ice concentration.
! ---------------------------------------------------------------------
      CALL ice_advect_all_tile (ng, tile,                               &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS                  &
     &                      )
      CALL wclock_off (ng, iNLM, 56)
      RETURN
      END SUBROUTINE ice_advect
      SUBROUTINE ice_advect_all_tile (ng, tile,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS                  &
     &                      )
      USE mod_param
      USE mod_ncparam
      USE mod_grid
      USE mod_ocean
      USE mod_ice
      USE mod_forces
      USE mod_scalars
      USE mod_stepping
      USE mod_boundary
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE mp_exchange_mod, ONLY : mp_exchange2d
      USE i2d_bc_mod
      USE tibc_mod, ONLY : tibc_tile
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variables
!
      integer :: i, j
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
      CALL ice_advect_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), linew(ng), liold(ng), liunw(ng),  &
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % zice,                            &
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      ICE(ng) % ai                                &
     &                      )
!
      CALL i2d_bc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  liold(ng), linew(ng),                           &
     &                  BOUNDARY(ng)%ai_west(LBj:UBj),                  &
     &                  BOUNDARY(ng)%ai_east(LBj:UBj),                  &
     &                  BOUNDARY(ng)%ai_north(LBi:UBi),                 &
     &                  BOUNDARY(ng)%ai_south(LBi:UBi),                 &
     &                  ICE(ng)%ui,                                     &
     &                  ICE(ng)%vi,                                     &
     &                  ICE(ng)%ai, LBC(:,isAice,ng),0)
!
! ---------------------------------------------------------------------
!  Advect the ice thickness.
! ---------------------------------------------------------------------
      CALL ice_advect_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), linew(ng), liold(ng), liunw(ng),  &
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % zice,                            &
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      ICE(ng) % hi                                &
     &                      )
!
      CALL i2d_bc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  liold(ng), linew(ng),                           &
     &                  BOUNDARY(ng)%hi_west(LBj:UBj),                  &
     &                  BOUNDARY(ng)%hi_east(LBj:UBj),                  &
     &                  BOUNDARY(ng)%hi_north(LBi:UBi),                 &
     &                  BOUNDARY(ng)%hi_south(LBi:UBi),                 &
     &                  ICE(ng)%ui,                                     &
     &                  ICE(ng)%vi,                                     &
     &                  ICE(ng)%hi, LBC(:,isHice,ng),1)
!
! ---------------------------------------------------------------------
!  Advect the snow thickness.
! ---------------------------------------------------------------------
      CALL ice_advect_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), linew(ng), liold(ng), liunw(ng),  &
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % zice,                            &
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      ICE(ng) % hsn                               &
     &                      )
!
      CALL i2d_bc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  liold(ng), linew(ng),                           &
     &                  BOUNDARY(ng)%hsn_west(LBj:UBj),                 &
     &                  BOUNDARY(ng)%hsn_east(LBj:UBj),                 &
     &                  BOUNDARY(ng)%hsn_north(LBi:UBi),                &
     &                  BOUNDARY(ng)%hsn_south(LBi:UBi),                &
     &                  ICE(ng)%ui,                                     &
     &                  ICE(ng)%vi,                                     &
     &                  ICE(ng)%hsn, LBC(:,isHsno,ng),1)
!
! ---------------------------------------------------------------------
!  Advect the surface melt water.
! ---------------------------------------------------------------------
      CALL ice_advect_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), linew(ng), liold(ng), liunw(ng),  &
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % zice,                            &
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      ICE(ng) % sfwat                             &
     &                      )
!
      CALL i2d_bc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  liold(ng), linew(ng),                           &
     &                  BOUNDARY(ng)%sfwat_west(LBj:UBj),               &
     &                  BOUNDARY(ng)%sfwat_east(LBj:UBj),               &
     &                  BOUNDARY(ng)%sfwat_north(LBi:UBi),              &
     &                  BOUNDARY(ng)%sfwat_south(LBi:UBi),              &
     &                  ICE(ng)%ui,                                     &
     &                  ICE(ng)%vi,                                     &
     &                  ICE(ng)%sfwat, LBC(:,isSfwat,ng),0)
!
! ---------------------------------------------------------------------
!  Advect the interior ice temperature.
! ---------------------------------------------------------------------
!
      CALL ice_advect_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), linew(ng), liold(ng), liunw(ng),  &
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % zice,                            &
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      ICE(ng) % enthalpi                          &
     &                      )
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ICE(ng)%ti(i,j,linew(ng)) = ICE(ng)%enthalpi(i,j,linew(ng))/  &
       &                  MAX(ICE(ng)%hi(i,j,linew(ng)),1.0E-6_r8)
          IF (ICE(ng)%hi(i,j,linew(ng)).LE.min_h(ng)) THEN
            ICE(ng)%enthalpi(i,j,linew(ng)) = 0.0_r8
            ICE(ng)%ti(i,j,linew(ng)) = 0.0_r8
          END IF
        ENDDO
      ENDDO
!
      CALL tibc_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                liold(ng), linew(ng),                             &
     &                ICE(ng)%ui,                                       &
     &                ICE(ng)%vi,                                       &
     &                ICE(ng)%hi,                                       &
     &                ICE(ng)%ti,                                       &
     &                ICE(ng)%enthalpi)
!
! ---------------------------------------------------------------------
!  Advect the ice age.
! ---------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ICE(ng)%hage(i,j,liold(ng)) = ICE(ng)%hi(i,j,liold(ng))*      &
     &                       ICE(ng)%ageice(i,j,liold(ng))
          ICE(ng)%hage(i,j,linew(ng)) = ICE(ng)%hi(i,j,linew(ng))*      &
     &                       ICE(ng)%ageice(i,j,linew(ng))
          IF (ICE(ng)%hi(i,j,liold(ng)).LE.min_h(ng)) THEN
            ICE(ng)%hage(i,j,liold(ng)) = 0.0_r8
            ICE(ng)%ageice(i,j,liold(ng)) = 0.0_r8
          END IF
        ENDDO
      ENDDO
      CALL ice_advect_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), linew(ng), liold(ng), liunw(ng),  &
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % zice,                            &
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      ICE(ng) % hage                              &
     &                      )
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ICE(ng)%ageice(i,j,linew(ng)) = ICE(ng)%hage(i,j,linew(ng))/  &
     &                  MAX(ICE(ng)%hi(i,j,linew(ng)),1.0E-6_r8)
          IF (ICE(ng)%hi(i,j,linew(ng)).LE.min_h(ng)) THEN
            ICE(ng)%hage(i,j,linew(ng)) = 0.0_r8
            ICE(ng)%ageice(i,j,linew(ng)) = 0.0_r8
          END IF
        ENDDO
      ENDDO
!
!        CALL i2d_bc_tile (ng, tile,                                    &
!     &                    LBi, UBi, LBj, UBj,                          &
!     &                    IminS, ImaxS, JminS, JmaxS,                  &
!     &                    liold(ng), linew(ng),                        &
!     &                    BOUNDARY(ng)%ageice_west(LBj:UBj),           &
!     &                    BOUNDARY(ng)%ageice_east(LBj:UBj),           &
!     &                    BOUNDARY(ng)%ageice_north(LBi:UBi),          &
!     &                    BOUNDARY(ng)%ageice_south(LBi:UBi),          &
!     &                    ICE(ng)%ui,                                  &
!     &                    ICE(ng)%vi,                                  &
!     &                    ICE(ng)%ageice, LBC(:,isHice,ng),0)
!      CALL ageicebc_tile (ng, tile,                                    &
!     &                          LBi, UBi, LBj, UBj,                    &
!     &                          liold(ng), linew(ng),                  &
!     &                          ICE(ng)%ui,                            &
!     &                          ICE(ng)%vi,                            &
!     &                          ICE(ng)%hi,                            &
!     &                          ICE(ng)%ageice,                        &
!     &                          ICE(ng)%hage)
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%ai(:,:,linew(ng)))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%hi(:,:,linew(ng)))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%hsn(:,:,linew(ng)))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%sfwat(:,:,linew(ng)))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%ti(:,:,linew(ng)))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%enthalpi(:,:,linew(ng)))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%ageice(:,:,linew(ng)))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%hage(:,:,linew(ng)))
      END IF
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ICE(ng)%ai(:,:,linew(ng)),                    &
     &                    ICE(ng)%hi(:,:,linew(ng)),                    &
     &                    ICE(ng)%hsn(:,:,linew(ng)),                   &
     &                    ICE(ng)%sfwat(:,:,linew(ng)))
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ICE(ng)%ti(:,:,linew(ng)),                    &
     &                    ICE(ng)%enthalpi(:,:,linew(ng)),              &
     &                    ICE(ng)%ageice(:,:,linew(ng)),                &
     &                    ICE(ng)%hage(:,:,linew(ng)))
      RETURN
      END SUBROUTINE ice_advect_all_tile
!
!==========================================================================!
      SUBROUTINE ice_advect_tile (ng, tile,                             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nrhs, linew, liold, liunw,                &
     &                        rmask,                                    &
     &                        zice,                                     &
     &                        on_u, om_v, omn,                          &
     &                        ui, vi, scr)
!==========================================================================!
      USE mod_param
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, linew, liold, liunw
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: zice(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: ui(LBi:,LBj:,:)
      real(r8), intent(in) :: vi(LBi:,LBj:,:)
      real(r8), intent(inout) :: scr(LBi:,LBj:,:)
!
! Local variable definitions
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ar
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: aflxu
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: aflxv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: aif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), parameter :: epsil = 1.0E-15_r8
      real(r8), parameter :: add = 3.0E+3_r8
      real(r8) :: fakt1
      real(r8) :: fakt2
      real(r8) :: aim1
      real(r8) :: ajm1
      real(r8) :: aim1jm1u
      real(r8) :: aim1jm1v
      real(r8) :: ajp1
      real(r8) :: aim1jp1
      real(r8) :: aip1
      real(r8) :: aip1jm1
      real(r8) :: Cu_crss, Cu
      real(r8) :: rateu
      real(r8) :: ratev
      real(r8) :: rateyiu
      real(r8) :: ratexiv
      real(r8) :: ratiou
      real(r8) :: ratiov
      real(r8) :: ratioyiu
      real(r8) :: ratioxiv
      real(r8) :: uiv
      real(r8) :: viu
      real(r8) :: uspeed
      real(r8) :: vspeed
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
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
!
! upstream:
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax+1
          aflxu(i,j)=on_u(i,j)*                                         &
     &          (max(0.0_r8,ui(i,j,liunw))*scr(i-1,j,liold)             &
     &          +min(0.0_r8,ui(i,j,liunw))*scr(i,j,liold))
        END DO
      END DO
      DO j=Jmin,Jmax+1
        DO i=Imin,Imax
          aflxv(i,j)=om_v(i,j)*                                         &
     &          (max(0.0_r8,vi(i,j,liunw))*scr(i,j-1,liold)             &
     &          +min(0.0_r8,vi(i,j,liunw))*scr(i,j,liold))
!
        END DO
      END DO
!
! step number 1 in mpdata:
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
!
          ar(i,j)=1.0_r8/omn(i,j)
          aif(i,j)=(scr(i,j,liold)-dtice(ng)*(aflxu(i+1,j)-aflxu(i,j)   &
     &        +aflxv(i,j+1)-aflxv(i,j))*ar(i,j))
          aif(i,j) = aif(i,j)*rmask(i,j)
          IF (zice(i,j).ne.0.0_r8) aif(i,j) = 0.0_r8
        END DO
      END DO
!
! set values at the open boundaries
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jmin,Jmax
            aif(Istr-1,j)=aif(Istr,j)   !? scr(Istr-1,j,liold)
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jmin,Jmax
            aif(Iend+1,j)=aif(Iend,j)  !? scr(Iend+1,j,liold)
          END DO
        END IF
      END IF
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Imin,Imax
             aif(i,Jstr-1)=aif(i,Jstr)  !??? scr(i,Jstr-1,liold)
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Imin,Imax
             aif(i,Jend+1)=aif(i,Jend)  !??? scr(i,Jend+1,liold)
          END DO
        END IF
        IF (.not.EWperiodic(ng)) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            aif(Istr-1,Jstr-1)=aif(Istr,Jstr)
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            aif(Istr-1,Jend+1)=aif(Istr,Jend)
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            aif(Iend+1,Jstr-1)=aif(Iend,Jstr)
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            aif(Iend+1,Jend+1)=aif(Iend,Jend)
          END IF
        END IF
      END IF
!
! mask ???
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
            scr(i,j,linew) = aif(i,j)
        END DO
      END DO
!
      RETURN
      END SUBROUTINE ice_advect_tile
      END MODULE ice_advect_mod
