      MODULE step3d_t_mod
!
!svn $Id: step3d_t.F 1488 2012-07-17 15:49:19Z kate $
!=======================================================================
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine time-steps tracer equations.  Notice that advective    !
!  and diffusive terms are time-stepped differently. It applies the    !
!  corrector time-step for horizontal/vertical advection,  vertical    !
!  diffusion, nudging if necessary, and lateral boundary conditions.   !
!                                                                      !
!  Notice that at input the tracer arrays have:                        !
!                                                                      !
!    t(:,:,:,nnew,:)   m Tunits  n+1     horizontal/vertical diffusion !
!                                        terms plus source/sink terms  !
!                                        (biology, sediment), if any   !
!                                                                      !
!    t(:,:,:,3   ,:)   Tunits    n+1/2   advective terms and vertical  !
!                                        diffusion predictor step      !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: step3d_t
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE step3d_t (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_clima
      USE mod_diags
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
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
      CALL wclock_on (ng, iNLM, 35)
      CALL step3d_t_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nrhs(ng), nstp(ng), nnew(ng),                 &
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng) % umask,                             &
     &                    GRID(ng) % vmask,                             &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % Huon,                              &
     &                    GRID(ng) % Hvom,                              &
     &                    GRID(ng) % z_r,                               &
     &                    CLIMA(ng) % Tnudgcof,                         &
     &                    CLIMA(ng) % tclm,                             &
     &                    MIXING(ng) % Akt,                             &
     &                    OCEAN(ng) % W,                                &
     &                    DIAGS(ng) % DiaTwrk,                          &
     &                    OCEAN(ng) % t                                 &
     &                               )
      CALL wclock_off (ng, iNLM, 35)
      RETURN
      END SUBROUTINE step3d_t
!
!***********************************************************************
      SUBROUTINE step3d_t_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nrhs, nstp, nnew,                       &
     &                          rmask, umask, vmask,                    &
     &                          pm, pn,                                 &
     &                          Hz, Huon, Hvom,                         &
     &                          z_r,                                    &
     &                          Tnudgcof, tclm,                         &
     &                          Akt,                                    &
     &                          W,                                      &
     &                          DiaTwrk,                                &
     &                          t                                       &
     &                           )
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE ice_frazil_mod, ONLY : ice_frazil
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
      USE mp_exchange_mod, ONLY : mp_exchange4d
      USE t3dbc_mod, ONLY : t3dbc_tile
      USE pt3dbc_mod, ONLY : pt3dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nstp, nnew
!
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Tnudgcof(LBi:,LBj:,:)
      real(r8), intent(in) :: tclm(LBi:,LBj:,:,:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
      real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
!
!  Local variable declarations.
!
      integer :: i, ibt, is, itrc, j, k, ltrc
      integer :: idiag
      real(r8), parameter :: eps = 1.0E-16_r8
      real(r8) :: cff, cff1, cff2, cff3
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: oHz
      real(r8) :: my_maxbio(15)
!
!  Define glacial dFe solubility in dye units = (dFe_soly / end member
!   concentration) * end member dye value = (0.6 nM/ 29 nM) * 10000.0
!   (see set_vbc.F and Johnson et al., 1997)
!
!  Trying new values of solubility (0.4 nM, see Fig. 3b, Taglaibue et
!   al., 2012) and time constant (0.03/yr, residence time of Th...see
!   note from Pete Sedwick) - msd/slm (6/18/14)
!
!     real(r8), parameter :: Glac_Fe_Soly = 206.90_r8
!      real(r8), parameter :: Glac_Fe_Soly = 137.93_r8
!     real(r8), parameter :: kFe = 1.585E-10_r8  ! 0.005/yr in 1/s
!      real(r8), parameter :: kFe = 9.513E-10_r8  ! 0.03/yr in 1/s
!
!
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
!  Time-step horizontal advection term.
!-----------------------------------------------------------------------
!
!  Compute inverse thickness.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            oHz(i,j,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
      END DO
!
!  Compute horizontal tracer advection fluxes.
!
      T_LOOP : DO itrc=1,NT(ng)
        K_LOOP : DO k=1,N(ng)
!
!  Third-order, uptream-biased horizontal advective fluxes.
!
          DO j=Jstr,Jend
            DO i=Istrm1,Iendp2
              FX(i,j)=t(i  ,j,k,3,itrc)-                                &
     &                t(i-1,j,k,3,itrc)
              FX(i,j)=FX(i,j)*umask(i,j)
            END DO
          END DO
          IF (.not.ComposedGrid(ng)) THEN
            IF (.not.EWperiodic(ng)) THEN
              IF (DOMAIN(ng)%Western_Edge(tile)) THEN
                DO j=Jstr,Jend
                  FX(Istr-1,j)=FX(Istr,j)
                END DO
              END IF
              IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
                DO j=Jstr,Jend
                  FX(Iend+2,j)=FX(Iend+1,j)
                END DO
              END IF
            END IF
          END IF
!
          DO j=Jstr,Jend
            DO i=Istr-1,Iend+1
              curv(i,j)=FX(i+1,j)-FX(i,j)
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              FX(i,j)=Huon(i,j,k)*0.5_r8*                               &
     &                (t(i-1,j,k,3,itrc)+                               &
     &                 t(i  ,j,k,3,itrc))-                              &
     &                cff1*(curv(i-1,j)*MAX(Huon(i,j,k),0.0_r8)+        &
     &                      curv(i  ,j)*MIN(Huon(i,j,k),0.0_r8))
            END DO
          END DO
!
          DO j=Jstrm1,Jendp2
            DO i=Istr,Iend
              FE(i,j)=t(i,j  ,k,3,itrc)-                                &
     &                t(i,j-1,k,3,itrc)
              FE(i,j)=FE(i,j)*vmask(i,j)
            END DO
          END DO
          IF (.not.ComposedGrid(ng)) THEN
            IF (.not.NSperiodic(ng)) THEN
              IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
                DO i=Istr,Iend
                  FE(i,Jstr-1)=FE(i,Jstr)
                END DO
              END IF
              IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
                DO i=Istr,Iend
                  FE(i,Jend+2)=FE(i,Jend+1)
                END DO
              END IF
            END IF
          END IF
!
          DO j=Jstr-1,Jend+1
            DO i=Istr,Iend
              curv(i,j)=FE(i,j+1)-FE(i,j)
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              FE(i,j)=Hvom(i,j,k)*0.5_r8*                               &
     &                (t(i,j-1,k,3,itrc)+                               &
     &                 t(i,j  ,k,3,itrc))-                              &
     &                cff1*(curv(i,j-1)*MAX(Hvom(i,j,k),0.0_r8)+        &
     &                      curv(i,j  )*MIN(Hvom(i,j,k),0.0_r8))
            END DO
          END DO
!
!  Time-step horizontal advection term.  Advective fluxes have units
!  of Tunits m3/s.  The new tracer has units of m Tunits.
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
              cff=dt(ng)*pm(i,j)*pn(i,j)
              cff1=cff*(FX(i+1,j)-FX(i,j))
              cff2=cff*(FE(i,j+1)-FE(i,j))
              cff3=cff1+cff2
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff3
              DiaTwrk(i,j,k,itrc,iTxadv)=-cff1
              DiaTwrk(i,j,k,itrc,iTyadv)=-cff2
              DiaTwrk(i,j,k,itrc,iThadv)=-cff3
            END DO
          END DO
        END DO K_LOOP
      END DO T_LOOP
!
!-----------------------------------------------------------------------
!  Time-step vertical advection term.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO itrc=1,NT(ng)
!
!  Fourth-order, central differences vertical advective flux.
!
          cff1=0.5_r8
          cff2=7.0_r8/12.0_r8
          cff3=1.0_r8/12.0_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                (cff2*(t(i,j,k  ,3,itrc)+                         &
     &                       t(i,j,k+1,3,itrc))-                        &
     &                 cff3*(t(i,j,k-1,3,itrc)+                         &
     &                       t(i,j,k+2,3,itrc)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,1)=W(i,j,1)*                                           &
     &              (cff1*t(i,j,1,3,itrc)+                              &
     &               cff2*t(i,j,2,3,itrc)-                              &
     &               cff3*t(i,j,3,3,itrc))
            FC(i,N(ng)-1)=W(i,j,N(ng)-1)*                               &
     &                    (cff1*t(i,j,N(ng)  ,3,itrc)+                  &
     &                     cff2*t(i,j,N(ng)-1,3,itrc)-                  &
     &                     cff3*t(i,j,N(ng)-2,3,itrc))
            FC(i,N(ng))=0.0_r8
          END DO
!
!  Time-step vertical advection term (Tunits).
!  Convert units of tracer diagnostic terms to Tunits.
!
          DO i=Istr,Iend
            CF(i,0)=dt(ng)*pm(i,j)*pn(i,j)
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=CF(i,0)*(FC(i,k)-FC(i,k-1))
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff1
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)*oHz(i,j,k)
              DiaTwrk(i,j,k,itrc,iTvadv)=-cff1
              DO idiag=1,NDT
                DiaTwrk(i,j,k,itrc,idiag)=DiaTwrk(i,j,k,itrc,idiag)*    &
     &                                    oHz(i,j,k)
              END DO
            END DO
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Time-step vertical diffusion term.
!-----------------------------------------------------------------------
!
        DO itrc=1,NT(ng)
          ltrc=MIN(NAT,itrc)
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff1*Hz(i,j,k  )-                                 &
     &                dt(ng)*Akt(i,j,k-1,ltrc)*oHz(i,j,k  )
              CF(i,k)=cff1*Hz(i,j,k+1)-                                 &
     &                dt(ng)*Akt(i,j,k+1,ltrc)*oHz(i,j,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
            DC(i,0)=0.0_r8
          END DO
!
!  LU decomposition and forward substitution.
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(Hz(i,j,k)+Hz(i,j,k+1))+                     &
     &                dt(ng)*Akt(i,j,k,ltrc)*(oHz(i,j,k)+oHz(i,j,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
              DC(i,k)=cff*(t(i,j,k+1,nnew,itrc)-t(i,j,k,nnew,itrc)-     &
     &                     FC(i,k)*DC(i,k-1))
            END DO
          END DO
!
!  Backward substitution.
!
          DO i=Istr,Iend
            DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            END DO
          END DO
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)*Akt(i,j,k,ltrc)
              cff1=dt(ng)*oHz(i,j,k)*(DC(i,k)-DC(i,k-1))
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+cff1
              DiaTwrk(i,j,k,itrc,iTvdif)=DiaTwrk(i,j,k,itrc,iTvdif)+    &
     &                                   cff1
            END DO
          END DO
        END DO
      END DO
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions and, if appropriate, nudge
!  to tracer data and apply Land/Sea mask.
!-----------------------------------------------------------------------
!
      DO itrc=1,NT(ng)
!
!  Set lateral boundary conditions.
!
!       IF (itrc.le.NAT) THEN
        CALL t3dbc_tile (ng, tile, itrc,                                &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, nnew,                                    &
     &                   t)
!       ELSE
!        CALL pt3dbc_tile (ng, tile, itrc,                              &
!     &                    LBi, UBi, LBj, UBj, N(ng), NT(ng),           &
!     &                    IminS, ImaxS, JminS, JmaxS,                  &
!     &                    nstp, nnew,                                  &
!     &                    t)
!       END IF
!
!  Nudge towards tracer climatology.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              IF (itrc.eq.inert(4) .and. k.eq.1) THEN
                 IF (tclm(i,j,k,itrc).gt.0.01) THEN
                    t(i,j,k,nnew,itrc)=tclm(i,j,k,itrc)
                 END IF
              END IF
            END DO
          END DO
        END DO
!
!  Apply Land/Sea mask.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)*rmask(i,j)
            END DO
          END DO
        END DO
!
!  Apply clamp to salt
!
      IF (itrc .ge. isalt) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              t(i,j,k,nnew,itrc) = max(t(i,j,k,nnew,itrc),0.0_r8)
            END DO
          END DO
        END DO
      END IF
!
!  Apply scavenging term (as defined in Eq. 5 of Johnson et al., 1997) to
!   dye_03 (glacial melt dye)
!
!      IF (itrc == inert(3)) THEN
!        DO k=1,N(ng)
!          DO j=JstrR,JendR
!            DO i=IstrR,IendR
!              IF (t(i,j,k,nnew,itrc) .gt. Glac_Fe_Soly) THEN
!                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+                 &
!     &                             dt(ng)*kFe*                         &
!     &                             (Glac_Fe_Soly-t(i,j,k,nnew,itrc))
!              END IF
!            END DO
!          END DO
!        END DO
!      END IF
!
!  Compute time-rate-of-change diagnostic term.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              DiaTwrk(i,j,k,itrc,iTrate)=t(i,j,k,nnew,itrc)-            &
     &                                   DiaTwrk(i,j,k,itrc,iTrate)
            END DO
          END DO
        END DO
        IF (itrc == isalt) THEN
          CALL ice_frazil(ng, tile)
        END IF
!
!  Apply periodic boundary conditions.
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            t(:,:,:,nnew,itrc))
        END IF
      END DO
!
!  Exchange boundary data.
!
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    t(:,:,:,nnew,:))
      RETURN
      END SUBROUTINE step3d_t_tile
      END MODULE step3d_t_mod
