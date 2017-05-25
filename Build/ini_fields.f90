      MODULE ini_fields_mod
!
!svn $Id: ini_fields.F 1451 2012-02-02 20:56:14Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes other time levels for 2D fields. It also   !
!  couples 3D and 2D momentum equations:  it initializes 2D momentum   !
!  (ubar,vbar) to the vertical integral of initial 3D momentum (u,v).  !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC :: ini_fields
      PUBLIC :: ini_zeta
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ini_fields (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_coupling
      USE mod_ocean
      USE mod_ice
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
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
      CALL wclock_on (ng, iNLM, 2)
      CALL ini_fields_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      kstp(ng), krhs(ng), knew(ng),               &
     &                      nstp(ng), nnew(ng),                         &
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
     &                      GRID(ng) % Hz,                              &
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % u,                              &
     &                      OCEAN(ng) % v,                              &
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      ICE(ng) % ai,                               &
     &                      ICE(ng) % hi,                               &
     &                      ICE(ng) % hsn,                              &
     &                      ICE(ng) % ti,                               &
     &                      ICE(ng) % sfwat,                            &
     &                      ICE(ng) % ageice,                           &
     &                      OCEAN(ng) % ubar,                           &
     &                      OCEAN(ng) % vbar,                           &
     &                      OCEAN(ng) % zeta)
      CALL wclock_off (ng, iNLM, 2)
      RETURN
      END SUBROUTINE ini_fields
!
!***********************************************************************
      SUBROUTINE ini_fields_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            kstp, krhs, knew,                     &
     &                            nstp, nnew,                           &
     &                            rmask, umask, vmask,                  &
     &                            Hz,                                   &
     &                            t, u, v,                              &
     &                            ui, vi, ai, hi, hsn,                  &
     &                            ti, sfwat, ageice,                    &
     &                            ubar, vbar, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_2d_mod
      USE exchange_3d_mod
      USE mp_exchange_mod, ONLY : mp_exchange2d
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
      USE t3dbc_mod, ONLY : t3dbc_tile
      USE pt3dbc_mod, ONLY : pt3dbc_tile
      USE u3dbc_mod, ONLY : u3dbc_tile
      USE v3dbc_mod, ONLY : v3dbc_tile
      USE u2dbc_mod, ONLY : u2dbc_tile
      USE v2dbc_mod, ONLY : v2dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kstp, krhs, knew
      integer, intent(in) :: nstp, nnew
!
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ui(LBi:,LBj:,:)
      real(r8), intent(inout) :: vi(LBi:,LBj:,:)
      real(r8), intent(inout) :: ai(LBi:,LBj:,:)
      real(r8), intent(inout) :: hi(LBi:,LBj:,:)
      real(r8), intent(inout) :: hsn(LBi:,LBj:,:)
      real(r8), intent(inout) :: ti(LBi:,LBj:,:)
      real(r8), intent(inout) :: sfwat(LBi:,LBj:,:)
      real(r8), intent(inout) :: ageice(LBi:,LBj:,:)
      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k
      real(r8) :: cff1, cff2, cff3, cff4, cff5, cff6
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
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
!-----------------------------------------------------------------------
!  If not perfect restart, initialize other time levels for 3D momentum.
!-----------------------------------------------------------------------
!
      IF (.not.PerfectRST(ng)) THEN
        DO j=JstrB,JendB
          DO k=1,N(ng)
            DO i=IstrM,IendB
              cff1=u(i,j,k,nstp)
              cff1=cff1*umask(i,j)
              u(i,j,k,nstp)=cff1
              u(i,j,k,nnew)=cff1
            END DO
          END DO
!
          IF (j.ge.JstrM) THEN
            DO k=1,N(ng)
              DO i=IstrB,IendB
                cff2=v(i,j,k,nstp)
                cff2=cff2*vmask(i,j)
                v(i,j,k,nstp)=cff2
                v(i,j,k,nnew)=cff2
              END DO
            END DO
          END IF
        END DO
!
!  Apply boundary conditions.
!
        CALL u3dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng),                     &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, nstp,                                    &
     &                   u)
        CALL v3dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng),                     &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, nstp,                                    &
     &                   v)
        CALL u3dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng),                     &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, nnew,                                    &
     &                   u)
        CALL v3dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng),                     &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, nnew,                                    &
     &                   v)
      END IF
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          u(:,:,:,nstp))
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          v(:,:,:,nstp))
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          u(:,:,:,nnew))
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          v(:,:,:,nnew))
      END IF
!
      CALL mp_exchange3d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    u(:,:,:,nstp), v(:,:,:,nstp),                 &
     &                    u(:,:,:,nnew), v(:,:,:,nnew))
!
!-----------------------------------------------------------------------
!  If not perfect restart, compute vertically-integrated momentum
!  (ubar, vbar) from initial 3D momentum (u, v).
!-----------------------------------------------------------------------
!
!  Here DC(i,1:N) are the grid cell thicknesses, DC(i,0) is the total
!  depth of the water column, and CF(i,0) is the vertical integral.
!
      IF (.not.PerfectRST(ng)) THEN
        DO j=JstrB,JendB
          DO i=IstrM,IendB
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=IstrM,IendB
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+DC(i,k)*u(i,j,k,nstp)
            END DO
          END DO
          DO i=IstrM,IendB
            cff1=1.0_r8/DC(i,0)
            cff2=CF(i,0)*cff1
            cff2=cff2*umask(i,j)
            ubar(i,j,kstp)=cff2
            ubar(i,j,knew)=cff2
          END DO
!
          IF (j.ge.JstrM) THEN
            DO i=IstrB,IendB
              DC(i,0)=0.0_r8
              CF(i,0)=0.0_r8
            END DO
            DO k=1,N(ng)
              DO i=IstrB,IendB
                DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
                DC(i,0)=DC(i,0)+DC(i,k)
                CF(i,0)=CF(i,0)+DC(i,k)*v(i,j,k,nstp)
              END DO
            END DO
            DO i=IstrB,IendB
              cff1=1.0_r8/DC(i,0)
              cff2=CF(i,0)*cff1
              cff2=cff2*vmask(i,j)
              vbar(i,j,kstp)=cff2
              vbar(i,j,knew)=cff2
            END DO
          END IF
        END DO
!
!  Apply boundary conditions.
!
        IF (.not.(ANY(LBC(:,isUbar,ng)%radiation).or.                   &
     &            ANY(LBC(:,isVbar,ng)%radiation).or.                   &
     &            ANY(LBC(:,isUbar,ng)%Flather).or.                     &
     &            ANY(LBC(:,isVbar,ng)%Flather))) THEN
          CALL u2dbc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs, kstp, kstp,                            &
     &                     ubar, vbar, zeta)
          CALL v2dbc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs, kstp, kstp,                            &
     &                     ubar, vbar, zeta)
          CALL u2dbc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs, kstp, knew,                            &
     &                     ubar, vbar, zeta)
          CALL v2dbc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs, kstp, knew,                            &
     &                     ubar, vbar, zeta)
        END IF
      END IF
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ubar(:,:,kstp))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vbar(:,:,kstp))
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ubar(:,:,knew))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vbar(:,:,knew))
      END IF
!
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    ubar(:,:,kstp), vbar(:,:,kstp),               &
     &                    ubar(:,:,knew), vbar(:,:,knew))
!
!-----------------------------------------------------------------------
!  If not perfect restart, initialize other time levels for tracers.
!-----------------------------------------------------------------------
!
      IF (.not.PerfectRST(ng)) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrB,JendB
              DO i=IstrB,IendB
                cff1=t(i,j,k,nstp,itrc)
                cff1=cff1*rmask(i,j)
                t(i,j,k,nstp,itrc)=cff1
                t(i,j,k,nnew,itrc)=cff1
              END DO
            END DO
          END DO
!
!  Apply boundary conditions.
!
!       IF (itrc.le.NAT) THEN
          CALL t3dbc_tile (ng, tile, itrc,                              &
     &                     LBi, UBi, LBj, UBj, N(ng), NT(ng),           &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nstp,                                  &
     &                     t)
          CALL t3dbc_tile (ng, tile, itrc,                              &
     &                     LBi, UBi, LBj, UBj, N(ng), NT(ng),           &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
     &                     t)
!         ELSE
!          CALL pt3dbc_tile (ng, tile, itrc,                            &
!     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),         &
!     &                      IminS, ImaxS, JminS, JmaxS,                &
!     &                      nstp, nstp,                                &
!     &                      t)
!          CALL pt3dbc_tile (ng, tile, itrc,                            &
!     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),         &
!     &                      IminS, ImaxS, JminS, JmaxS,                &
!     &                      nstp, nnew,                                &
!     &                      t)
!         END IF
        END DO
      END IF
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NT(ng)
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            t(:,:,:,nstp,itrc))
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            t(:,:,:,nnew,itrc))
        END DO
      END IF
!
      CALL mp_exchange4d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    t(:,:,:,nstp,:),                              &
     &                    t(:,:,:,nnew,:))
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff1=ui(i,j,nstp)
          cff1=cff1*umask(i,j)
          ui(i,j,nstp)=cff1
          ui(i,j,nnew)=cff1
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff1=vi(i,j,nstp)
          cff1=cff1*vmask(i,j)
          vi(i,j,nstp)=cff1
          vi(i,j,nnew)=cff1
        END DO
      END DO
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          cff1=ai(i,j,nstp)
          cff2=hi(i,j,nstp)
          cff3=hsn(i,j,nstp)
          cff4=ti(i,j,nstp)
          cff5=sfwat(i,j,nstp)
          cff6=ageice(i,j,nstp)
          cff1=cff1*rmask(i,j)
          cff2=cff2*rmask(i,j)
          cff3=cff3*rmask(i,j)
          cff4=cff4*rmask(i,j)
          cff5=cff5*rmask(i,j)
          cff6=cff6*rmask(i,j)
          ai(i,j,nstp)=cff1
          ai(i,j,nnew)=cff1
          hi(i,j,nstp)=cff2
          hi(i,j,nnew)=cff2
          hsn(i,j,nstp)=cff3
          hsn(i,j,nnew)=cff3
          ti(i,j,nstp)=cff4
          ti(i,j,nnew)=cff4
          sfwat(i,j,nstp)=cff5
          sfwat(i,j,nnew)=cff5
          ageice(i,j,nstp)=cff6
          ageice(i,j,nnew)=cff6
        END DO
      END DO
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO i=1,2
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            ui(:,:,i))
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            vi(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            ai(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            hi(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            hsn(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            ti(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            sfwat(:,:,i))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            ageice(:,:,i))
        END DO
      END IF
        CALL mp_exchange3d (ng, tile, model, 4,                         &
     &                      LBi, UBi, LBj, UBj, 1, 2,                   &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      ui, vi, ai, hi)
        CALL mp_exchange3d (ng, tile, model, 4,                         &
     &                      LBi, UBi, LBj, UBj, 1, 2,                   &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      hsn, ti, sfwat, ageice)
      RETURN
      END SUBROUTINE ini_fields_tile
!
!***********************************************************************
      SUBROUTINE ini_zeta (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_coupling
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
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
      CALL wclock_on (ng, iNLM, 2)
      CALL ini_zeta_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    kstp(ng), krhs(ng), knew(ng),                 &
     &                    GRID(ng) % rmask,                             &
     &                    COUPLING(ng) % Zt_avg1,                       &
     &                    OCEAN(ng) % zeta)
      CALL wclock_off (ng, iNLM, 2)
      RETURN
      END SUBROUTINE ini_zeta
!
!***********************************************************************
      SUBROUTINE ini_zeta_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          kstp, krhs, knew,                       &
     &                          rmask,                                  &
     &                          Zt_avg1,                                &
     &                          zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE mp_exchange_mod, ONLY : mp_exchange2d
      USE zetabc_mod, ONLY : zetabc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kstp, krhs, knew
!
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(inout) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: zeta(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, kbed
      real(r8) :: cff1
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
!-----------------------------------------------------------------------
!  Initialize other time levels for free-surface.
!-----------------------------------------------------------------------
!
      IF (.not.PerfectRST(ng)) THEN
        IF (.not.(ANY(LBC(:,isFsur,ng)%radiation).or.                   &
     &            ANY(LBC(:,isFsur,ng)%Chapman))) THEN
          Imin=IstrB
          Imax=IendB
          Jmin=JstrB
          Jmax=JendB
        ELSE
          Imin=IstrT
          Imax=IendT
          Jmin=JstrT
          Jmax=JendT
        END IF
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            cff1=zeta(i,j,kstp)
            cff1=cff1*rmask(i,j)
            zeta(i,j,kstp)=cff1
            zeta(i,j,knew)=cff1
          END DO
        END DO
!
!  Apply boundary conditions.
!
        IF (.not.(ANY(LBC(:,isFsur,ng)%radiation).or.                   &
     &            ANY(LBC(:,isFsur,ng)%Chapman))) THEN
          CALL zetabc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, kstp,                           &
     &                      zeta)
          CALL zetabc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, knew,                           &
     &                      zeta)
        END IF
      END IF
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          zeta(:,:,kstp))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          zeta(:,:,knew))
        IF (PerfectRST(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            zeta(:,:,krhs))
        END IF
      END IF
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    zeta(:,:,kstp),                               &
     &                    zeta(:,:,knew))
      IF (PerfectRST(ng)) THEN
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      zeta(:,:,krhs))
      END IF
!
!-----------------------------------------------------------------------
!  Initialize fast-time averaged free-surface (Zt_avg1) with the inital
!  free-surface
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Zt_avg1(i,j)=zeta(i,j,kstp)
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Zt_avg1)
      END IF
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Zt_avg1)
      RETURN
      END SUBROUTINE ini_zeta_tile
      END MODULE ini_fields_mod
