      MODULE set_vbc_mod
!
!svn $Id: set_vbc.F 1474 2012-05-29 15:15:08Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module sets vertical boundary conditons for momentum and       !
!  tracers.                                                            !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: set_vbc
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_vbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
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
      CALL wclock_on (ng, iNLM, 6)
      CALL set_vbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng),                                      &
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % rdrag2,                             &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % zice,                               &
     &                   GRID(ng) % f,                                  &
     &                   OCEAN(ng) % t,                                 &
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v,                                 &
     &                   OCEAN(ng) % rho,                               &
     &                   FORCES(ng) % srflx,                            &
     &                   FORCES(ng) % Tair,                             &
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
     &                   FORCES(ng) % bustr,                            &
     &                   FORCES(ng) % bvstr,                            &
     &                   FORCES(ng) % stflx,                            &
     &                   FORCES(ng) % btflx)
      CALL wclock_off (ng, iNLM, 6)
      RETURN
      END SUBROUTINE set_vbc
!
!***********************************************************************
      SUBROUTINE set_vbc_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,                                    &
     &                         Hz,                                      &
     &                         rdrag2,                                  &
     &                         z_r, z_w,                                &
     &                         rmask,                                   &
     &                         zice,                                    &
     &                         f,                                       &
     &                         t,                                       &
     &                         u, v,                                    &
     &                         rho,                                     &
     &                         srflx,                                   &
     &                         Tair,                                    &
     &                         sustr, svstr,                            &
     &                         bustr, bvstr,                            &
     &                         stflx, btflx)
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
!
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: rdrag2(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: zice(LBi:,LBj:)
      real(r8), intent(in) :: f(LBi:,LBj:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
      real(r8), intent(in)    :: Tair(LBi:,LBj:)
      real(r8), intent(inout) :: sustr(LBi:,LBj:)
      real(r8), intent(inout) :: svstr(LBi:,LBj:)
      real(r8), intent(inout) :: bustr(LBi:,LBj:)
      real(r8), intent(inout) :: bvstr(LBi:,LBj:)
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
      real(r8), intent(inout) :: btflx(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, j, itrc
      real(r8) :: Hscale, Hlice, uStar, gturb, cff
      real(r8) :: salt_b, temp_b, melt_rate
      real(r8) :: cff1, cff2, cff3
      real(r8), parameter :: eps = 1.e-20_r8
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
!  Multiply fresh water flux with surface salinity. If appropriate,
!  apply correction.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          stflx(i,j,isalt)=stflx(i,j,isalt)
          stflx(i,j,isalt) = rmask(i,j)*stflx(i,j,isalt)
          btflx(i,j,isalt)=btflx(i,j,isalt)*t(i,j,1,nrhs,isalt)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  If ice shelf cavities, replace surface wind stress with ice shelf
!  cavity stress (m2/s2).
!-----------------------------------------------------------------------
!
!  Set quadratic ice shelf cavity stress.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
            cff1=0.25_r8*(v(i  ,j  ,N(ng),nrhs)+                        &
     &                    v(i  ,j+1,N(ng),nrhs)+                        &
     &                    v(i-1,j  ,N(ng),nrhs)+                        &
     &                    v(i-1,j+1,N(ng),nrhs))
            cff2=SQRT(u(i,j,N(ng),nrhs)*u(i,j,N(ng),nrhs)+cff1*cff1)
            sustr(i,j)=-0.5_r8*(rdrag2(i-1,j)+rdrag2(i,j))*             &
     &                 u(i,j,N(ng),nrhs)*cff2
          END IF
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
            cff1=0.25_r8*(u(i  ,j  ,N(ng),nrhs)+                        &
     &                    u(i+1,j  ,N(ng),nrhs)+                        &
     &                    u(i  ,j-1,N(ng),nrhs)+                        &
     &                    u(i+1,j-1,N(ng),nrhs))
            cff2=SQRT(cff1*cff1+v(i,j,N(ng),nrhs)*v(i,j,N(ng),nrhs))
            svstr(i,j)=-0.5_r8*(rdrag2(i,j-1)+rdrag2(i,j))*             &
     &                 v(i,j,N(ng),nrhs)*cff2
          END IF
        END DO
      END DO
!
!  Apply periodic or gradient boundary conditions for output
!  purposes only.
!
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  sustr)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  svstr)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr, svstr)
!
!-----------------------------------------------------------------------
!
!  Two sets of thermodynamics are currently available for the ice shelf
!   cavities now:
!
!  1) ISOMIP: surface tracer fluxes underneath ice are computed in the
!  manner defined by the ISOMIP test specification.  Note that the units
!  for stflx for heat are degC m/s (which means the rho*C_p part is
!  removed here) and for salt are psu m/s (unless sea ice is defined).
!  This version is setup for the ISOMIP 2.01 case where the heat and salt
!  fluxes for the open water/atmosphere interface match the relaxation
!  calculation in sections 3.2.6 and 3.2.7 of John Hunter's document.
!
!  2) : 3 equation formulation as per Holland and
!  Jenkins (JPO, 1999) with the transfer coefficients computed from the
!  friction velocity as in Schmidt et al. (OM, 2004). (Note that the
!  the Gamma_turb equation and check value for gamma_t in Schmidt use
!  log base 10 when it should be natural log (McPhee et al. 1987 JGR).)
!
!-----------------------------------------------------------------------
!
!
!  If there is an iceshelf, recalculate the total heat and salt flux
!   at any point below the iceshelf.
!
      Hscale = 1.0_r8/(rho0*Cp)
      DO j=Jstr,Jend      ! Have to index this way because of need
        DO i=Istr,Iend    !  for surface stresses below
          IF (zice(i,j) .ne. 0.0_r8) THEN
!
!  First, calculate the salinity in the boundary layer from solving the
!   3 equations simultaneously
!
            IF (Tair(i,j) .lt. -1.95_r8) THEN
              Hlice = Hlfreeze - blk_Cpi*(Tair(i,j)+1.95_r8)
            ELSE
              Hlice = Hlfreeze
            END IF
!  Compute friction velocity from surface stresses
            uStar = SQRT(SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+    &
     &                        (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2))
            uStar = uStar*rmask(i,j)
!  Compute transfer coefficients
            IF (uStar .gt. 1.0E-4_r8 .and. ABS(f(i,j)) .gt. 1.0E-8) THEN
              gturb = 2.5*LOG(5300.0_r8*uStar*uStar/ABS(f(i,j))) +      &
     &                7.12  !Changed to LOG from LOG10   ~2013-11-25 SM
            ELSE
              gturb = 0.0_r8
            ENDIF
            gamma_t = uStar/(gturb+65.9_r8)
            gamma_s = uStar/(gturb+2255.0_r8)
!           gamma_t = 9.00E-3_r8*uStar
!           gamma_s = 0.025_r8*gamma_t
!  Solve for the boundary layer salinity
            cff=blk_Cpw*gamma_t/Hlice
            cff1=cff*((0.0939_r8-t(i,j,N(ng),nrhs,itemp))+              &
     &                 7.6410E-4_r8*zice(i,j))-gamma_s
            cff2=cff1*cff1+4.0_r8*                                      &
     &           cff*0.057_r8*gamma_s*t(i,j,N(ng),nrhs,isalt)
            IF (cff2 .gt. 0.0_r8 .and. cff .ne. 0.0_r8) THEN
              cff3=-cff1+sqrt(cff2)
              IF (cff3 .lt. 0.0_r8) THEN
                salt_b=cff3/(-2.0_r8*0.057_r8*cff)
              ELSE
                salt_b=(-cff1-sqrt(cff2))/(-2.0_r8*0.057_r8*cff)
              END IF
            ELSE
              salt_b=5.0_r8
            END IF
!
!  Then, calculate the actual heat (degC m/s) and salt (psu m/s) fluxes
!   into the top model layer.  Note that this formulation includes the
!   calculation and addition of the melt rate term for the proper flux
!   condition through a material surface (see Jenkins et al. 2001 JPO).
!
!  Note: The salt flux has to be m/s if the SEA ICE model is being used.
!
            temp_b=0.0939_r8-0.057_r8*salt_b+7.6410E-4_r8*zice(i,j)
            IF (salt_b .gt. 5.0_r8) THEN
              melt_rate=-gamma_s*                                       &
     &                   (1.0_r8-(t(i,j,N(ng),nrhs,isalt)/salt_b))
            ELSE
              melt_rate=0.0_r8
            END IF
            stflx(i,j,itemp)=(rho(i,j,N(ng))+1000.0_r8)*blk_Cpw*        &
     &                       (gamma_t+melt_rate)*                       &
     &                       (temp_b-t(i,j,N(ng),nrhs,itemp))*Hscale
            IF (t(i,j,N(ng),nrhs,isalt) .gt. 0.001_r8) THEN
              stflx(i,j,isalt)=(gamma_s+melt_rate)*                     &
     &                         (salt_b-t(i,j,N(ng),nrhs,isalt))/        &
     &                         t(i,j,N(ng),nrhs,isalt)
            ELSE
              stflx(i,j,isalt)=0.0_r8
            END IF
            IF (melt_rate .gt. 0.0_r8) THEN
               stflx(i,j,inert(3))=10000.0_r8*melt_rate
!           ELSE  !SM 2015-09-16 Fixed bug, but not turned on yet
!              stflx(i,j,inert(3))=0.0_r8
            END IF
!
!   Diagnostic print checks
!
!            IF (i.eq.180 .and. (j.eq.6 .or. j.eq.15)) THEN
!              print*,'Ice-shelf checks: '
!              print*,sustr(i,j),sustr(i+1,j),svstr(i,j),svstr(i,j+1)
!              print*,uStar,f(i,j),gturb,gamma_t,gamma_s
!              print*,' '
!              print*,Hlice,Hlfreeze,blk_Cpi,Tair(i,j)
!              print*,rho(i,j,N(ng)),blk_Cpw,gamma_t,cff
!              print*,t(i,j,N(ng),nrhs,itemp),zice(i,j),gamma_s,cff1
!              print*,t(i,j,N(ng),nrhs,isalt),cff2,cff3
!              print*,salt_b
!              print*,' '
!              print*,rho_ice,(rho(i,j,N(ng))+1000.0_r8),melt_rate
!              print*,zice(i,j),salt_b,temp_b
!              print*,Hscale,gamma_t,blk_Cpw,stflx(i,j,itemp)
!              print*,stflx(i,j,isalt)
!              print*,j,rmask(i,j)
!              print*,' '
!            END IF
          END IF
        END DO
      END DO
!
!  Apply gradient or periodic boundary conditions for the two fluxes
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  stflx(:,:,itemp))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  stflx(:,:,isalt))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  stflx(:,:,inert(3)))
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    stflx(:,:,itemp))
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    stflx(:,:,isalt))
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    stflx(:,:,inert(3)))
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          IF (zice(i,j).ne.0.0_r8) THEN
            srflx(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Set kinematic bottom momentum flux (m2/s2).
!-----------------------------------------------------------------------
!
!  Set quadratic bottom stress.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff1=0.25_r8*(v(i  ,j  ,1,nrhs)+                              &
     &                  v(i  ,j+1,1,nrhs)+                              &
     &                  v(i-1,j  ,1,nrhs)+                              &
     &                  v(i-1,j+1,1,nrhs))
          cff2=SQRT(u(i,j,1,nrhs)*u(i,j,1,nrhs)+cff1*cff1)
          cff1=0.5_r8*(rdrag2(i,j)+rdrag2(i-1,j))
          bustr(i,j)=cff1*u(i,j,1,nrhs)*cff2
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff1=0.25_r8*(u(i  ,j  ,1,nrhs)+                              &
     &                  u(i+1,j  ,1,nrhs)+                              &
     &                  u(i  ,j-1,1,nrhs)+                              &
     &                  u(i+1,j-1,1,nrhs))
          cff2=SQRT(cff1*cff1+v(i,j,1,nrhs)*v(i,j,1,nrhs))
          cff1=0.5_r8*(rdrag2(i,j)+rdrag2(i,j-1))
          bvstr(i,j)=cff1*v(i,j,1,nrhs)*cff2
        END DO
      END DO
!
!  Apply boundary conditions.
!
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bustr)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bvstr)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bustr, bvstr)
      RETURN
      END SUBROUTINE set_vbc_tile
      END MODULE set_vbc_mod
