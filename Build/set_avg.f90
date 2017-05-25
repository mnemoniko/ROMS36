      MODULE set_avg_mod
!
!svn $Id: set_avg.F 1490 2012-08-23 19:29:48Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine accumulates and computes output time-averaged       !
!  fields.  Due to synchronization, the time-averaged fields are       !
!  computed in delayed mode. All averages are accumulated at the       !
!  beggining of the next time-step.                                    !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC :: set_avg
      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_avg (ng, tile)
!***********************************************************************
!
      USE mod_param
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
      CALL wclock_on (ng, iNLM, 5)
      CALL set_avg_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   linew(ng), liunw(ng), lienw(ng),               &
     &                   nrhs(ng),                                      &
     &                   kstp(ng))
      CALL wclock_off (ng, iNLM, 5)
      RETURN
      END SUBROUTINE set_avg
!
!***********************************************************************
      SUBROUTINE set_avg_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Iout, Iuout, Ieout,                      &
     &                         Nout,                                    &
     &                         Kout)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_average
      USE mod_coupling
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ice
      USE mod_ocean
      USE mod_scalars
!
      USE exchange_2d_mod
      USE exchange_3d_mod
      USE mp_exchange_mod, ONLY : mp_exchange2d
      USE mp_exchange_mod, ONLY : mp_exchange3d
      USE uv_rotate_mod, ONLY : uv_rotate2d
      USE uv_rotate_mod, ONLY : uv_rotate3d
      USE vorticity_mod, ONLY : vorticity_tile
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Kout
      integer, intent(in) :: Iout, Iuout, Ieout
      integer, intent(in) :: Nout
!
!  Local variable declarations.
!
      integer :: i, it, itrc, j, k
      real(r8) :: fac
      real(r8) :: pfac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: rfac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: ufac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: vfac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: potvor(LBi:UBi,LBj:UBj,N(ng))
      real(r8) :: relvor(LBi:UBi,LBj:UBj,N(ng))
      real(r8) :: potvor_bar(LBi:UBi,LBj:UBj)
      real(r8) :: relvor_bar(LBi:UBi,LBj:UBj)
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
!  Return if time-averaging window is zero.
!-----------------------------------------------------------------------
!
      IF (nAVG(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
!  Compute vorticity fields.
!-----------------------------------------------------------------------
!
      IF (Aout(id2dPV,ng).or.Aout(id2dRV,ng).or.                        &
     &    Aout(id3dPV,ng).or.Aout(id3dRV,ng)) THEN
        CALL vorticity_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       Kout, Nout,                                &
     &                       GRID(ng) % pmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
     &                       GRID(ng) % fomn,                           &
     &                       GRID(ng) % h,                              &
     &                       GRID(ng) % om_u,                           &
     &                       GRID(ng) % on_v,                           &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
     &                       GRID(ng) % z_r,                            &
     &                       OCEAN(ng) % pden,                          &
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
     &                       OCEAN(ng) % ubar,                          &
     &                       OCEAN(ng) % vbar,                          &
     &                       OCEAN(ng) % zeta,                          &
     &                       potvor, relvor,                            &
     &                       potvor_bar, relvor_bar)
      END IF
!
!-----------------------------------------------------------------------
!  Initialize time-averaged arrays when appropriate.  Notice that
!  fields are initilized twice during re-start.  However, the time-
!  averaged fields are computed correctly.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsAVG(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nAVG(ng)).eq.1)).or.                          &
     &    ((iic(ng).ge.ntsAVG(ng)).and.(nAVG(ng).eq.1)).or.             &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
!
!  Initialize state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=OCEAN(ng)%zeta(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2d(i,j)=OCEAN(ng)%ubar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2d(i,j)=OCEAN(ng)%vbar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idu2dE,ng).and.Aout(idv2dN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .FALSE., .FALSE.,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
     &                      GRID(ng)%rmask_io,                          &
     &                      OCEAN(ng) % ubar(:,:,Kout),                 &
     &                      OCEAN(ng) % vbar(:,:,Kout),                 &
     &                      AVERAGE(ng)%avgu2dE,                        &
     &                      AVERAGE(ng)%avgv2dN)
        END IF
        IF (Aout(idUvel,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3d(i,j,k)=OCEAN(ng)%u(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVvel,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3d(i,j,k)=OCEAN(ng)%v(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idu3dE,ng).and.Aout(idv3dN,ng)) THEN
          CALL uv_rotate3d (ng, tile, .FALSE., .FALSE.,                 &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
     &                      GRID(ng)%rmask_io,                          &
     &                      OCEAN(ng) % u(:,:,:,Nout),                  &
     &                      OCEAN(ng) % v(:,:,:,Nout),                  &
     &                      AVERAGE(ng)%avgu3dE,                        &
     &                      AVERAGE(ng)%avgv3dN)
        END IF
        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=OCEAN(ng)%W(i,j,k)*           &
     &                                    GRID(ng)%pm(i,j)*             &
     &                                    GRID(ng)%pn(i,j)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idWvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgwvel(i,j,k)=OCEAN(ng)%wvel(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idDano,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgrho(i,j,k)=OCEAN(ng)%rho(i,j,k)
              END DO
            END DO
          END DO
        END IF
        DO it=1,NT(ng)
          IF (Aout(idTvar(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,it)=OCEAN(ng)%t(i,j,k,Nout,it)
                END DO
              END DO
            END DO
          END IF
        END DO
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=MIXING(ng)%Akv(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idTdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKt(i,j,k)=MIXING(ng)%Akt(i,j,k,itemp)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idSdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKs(i,j,k)=MIXING(ng)%Akt(i,j,k,isalt)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idHsbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghsbl(i,j)=MIXING(ng)%hsbl(i,j)
            END DO
          END DO
        END IF
!
!  Initialize 2D/3D coupling terms.
!
        IF (Aout(idUfx1,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg1(i,j)=COUPLING(ng)%DU_avg1(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUfx2,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg2(i,j)=COUPLING(ng)%DU_avg2(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVfx1,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg1(i,j)=COUPLING(ng)%DV_avg1(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVfx2,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg2(i,j)=COUPLING(ng)%DV_avg2(i,j)
            END DO
          END DO
        END IF
!
!  Initialize surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=FORCES(ng)%sustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsvs(i,j)=FORCES(ng)%svstr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgbus(i,j)=FORCES(ng)%bustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgbvs(i,j)=FORCES(ng)%bvstr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idPair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPair(i,j)=FORCES(ng)%Pair(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUwind(i,j)=FORCES(ng)%Uwind(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVwind(i,j)=FORCES(ng)%Vwind(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUairE,ng).and.Aout(idVairN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .FALSE., .FALSE.,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
     &                      GRID(ng)%rmask_io,                          &
     &                      FORCES(ng) % Uwind,                         &
     &                      FORCES(ng) % Vwind,                         &
     &                      AVERAGE(ng)%avgUwindE,                      &
     &                      AVERAGE(ng)%avgVwindN)
        END IF
        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=FORCES(ng)%stflx(i,j,itemp)
            END DO
          END DO
        END IF
        IF (Aout(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgswf(i,j)=FORCES(ng)%stflx(i,j,isalt)
            END DO
          END DO
        END IF
        IF (Aout(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsrf(i,j)=FORCES(ng)%srflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglhf(i,j)=FORCES(ng)%lhflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglrf(i,j)=FORCES(ng)%lrflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgshf(i,j)=FORCES(ng)%shflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgevap(i,j)=FORCES(ng)%evap(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idrain,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgrain(i,j)=FORCES(ng)%rain(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avguice(i,j)=ICE(ng)%ui(i,j,Iuout)
            END DO
          END DO
        END IF
        IF (Aout(idUice,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgvice(i,j)=ICE(ng)%vi(i,j,Iuout)
            END DO
          END DO
        END IF
        IF (Aout(idUiceE,ng).and.Aout(idViceN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .FALSE., .FALSE.,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
     &                      GRID(ng)%rmask_io,                          &
     &                      ICE(ng) % ui(:,:,Kout),                     &
     &                      ICE(ng) % vi(:,:,Kout),                     &
     &                      AVERAGE(ng)%avguiceE,                       &
     &                      AVERAGE(ng)%avgviceN)
        END IF
        IF (Aout(idAice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgaice(i,j)=ICE(ng)%ai(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idHice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghice(i,j)=ICE(ng)%hi(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idHsno,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghsno(i,j)=ICE(ng)%hsn(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idTice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgtice(i,j)=ICE(ng)%tis(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idTimid,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgtimid(i,j)=ICE(ng)%ti(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idSfwat,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsfwat(i,j)=ICE(ng)%sfwat(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idAgeice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgageice(i,j)=ICE(ng)%ageice(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idIomflx,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgiomflx(i,j)=ICE(ng)%io_mflux(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idSig11,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsig11(i,j)=ICE(ng)%sig11(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout(idSig12,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsig12(i,j)=ICE(ng)%sig12(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout(idSig22,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsig22(i,j)=ICE(ng)%sig22(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout(idT0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgT0mk(i,j)=ICE(ng)%t0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idS0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgS0mk(i,j)=ICE(ng)%s0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWfr,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWfr(i,j)=ICE(ng)%wfr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWai,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWai(i,j)=ICE(ng)%wai(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWao,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWao(i,j)=ICE(ng)%wao(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWio,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWio(i,j)=ICE(ng)%wio(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWro,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWro(i,j)=ICE(ng)%wro(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idTauiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgutau_iw(i,j)=ICE(ng)%utau_iw(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idChuiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgchu_iw(i,j)=ICE(ng)%chu_iw(i,j)
            END DO
          END DO
        END IF
!
!  Initialize vorticity fields.
!
        IF (Aout(id2dPV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgpvor2d(i,j)=potvor_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(id2dRV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgrvor2d(i,j)=relvor_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(id3dPV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgpvor3d(i,j,k)=potvor(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(id3dRV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgrvor3d(i,j,k)=relvor(i,j,k)
              END DO
            END DO
          END DO
        END IF
!
!  Initialize quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=OCEAN(ng)%zeta(i,j,Kout)*          &
     &                               OCEAN(ng)%zeta(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idU2av,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgU2(i,j)=OCEAN(ng)%ubar(i,j,Kout)*          &
     &                               OCEAN(ng)%ubar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idV2av,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgV2(i,j)=OCEAN(ng)%vbar(i,j,Kout)*          &
     &                               OCEAN(ng)%vbar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idUUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgUU(i,j,k)=OCEAN(ng)%u(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%u(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgVV(i,j,k)=OCEAN(ng)%v(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%v(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idUVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgUV(i,j,k)=0.25_r8*                       &
     &                                   (OCEAN(ng)%u(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%u(i+1,j  ,k,Nout))* &
     &                                   (OCEAN(ng)%v(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%v(i  ,j+1,k,Nout))
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idHUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgHuon(i,j,k)=GRID(ng)%Huon(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idHVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgHvom(i,j,k)=GRID(ng)%Hvom(i,j,k)
              END DO
            END DO
          END DO
        END IF
        DO it=1,NAT
          IF (Aout(idTTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgTT(i,j,k,it)=OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)*     &
     &                                        OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgUT(i,j,k,it)=0.5_r8*                   &
     &                                        OCEAN(ng)%u(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i-1,j,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i  ,j,k,     &
     &                                                     Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgVT(i,j,k,it)=0.5_r8*                   &
     &                                        OCEAN(ng)%v(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i,j-1,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i,j  ,k,     &
     &                                                     Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(iHUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgHuonT(i,j,k,it)=0.5_r8*                &
     &                                           GRID(ng)%Huon(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i-1,j,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i  ,j,k,  &
     &                                                        Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(iHVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgHvomT(i,j,k,it)=0.5_r8*                &
     &                                           GRID(ng)%Hvom(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i,j-1,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i,j  ,k,  &
     &                                                        Nout,it))
                END DO
              END DO
            END DO
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsAVG(ng)) THEN
!
!  Accumulate state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=AVERAGE(ng)%avgzeta(i,j)+        &
     &                                 OCEAN(ng)%zeta(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2d(i,j)=AVERAGE(ng)%avgu2d(i,j)+          &
     &                                OCEAN(ng)%ubar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2d(i,j)=AVERAGE(ng)%avgv2d(i,j)+          &
     &                                OCEAN(ng)%vbar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idu2dE,ng).and.Aout(idv2dN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .TRUE., .FALSE.,                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
     &                      GRID(ng)%rmask_io,                          &
     &                      OCEAN(ng) % ubar(:,:,Kout),                 &
     &                      OCEAN(ng) % vbar(:,:,Kout),                 &
     &                      AVERAGE(ng)%avgu2dE,                        &
     &                      AVERAGE(ng)%avgv2dN)
        END IF
        IF (Aout(idUvel,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3d(i,j,k)=AVERAGE(ng)%avgu3d(i,j,k)+    &
     &                                    OCEAN(ng)%u(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVvel,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3d(i,j,k)=AVERAGE(ng)%avgv3d(i,j,k)+    &
     &                                    OCEAN(ng)%v(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idu3dE,ng).and.Aout(idv3dN,ng)) THEN
          CALL uv_rotate3d (ng, tile, .TRUE., .FALSE.,                  &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
     &                      GRID(ng)%rmask_io,                          &
     &                      OCEAN(ng) % u(:,:,:,Nout),                  &
     &                      OCEAN(ng) % v(:,:,:,Nout),                  &
     &                      AVERAGE(ng)%avgu3dE,                        &
     &                      AVERAGE(ng)%avgv3dN)
        END IF
        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=AVERAGE(ng)%avgw3d(i,j,k)+    &
     &                                    OCEAN(ng)%W(i,j,k)*           &
     &                                    GRID(ng)%pm(i,j)*             &
     &                                    GRID(ng)%pn(i,j)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idWvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgwvel(i,j,k)=AVERAGE(ng)%avgwvel(i,j,k)+  &
     &                                     OCEAN(ng)%wvel(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idDano,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgrho(i,j,k)=AVERAGE(ng)%avgrho(i,j,k)+    &
     &                                    OCEAN(ng)%rho(i,j,k)
              END DO
            END DO
          END DO
        END IF
        DO it=1,NT(ng)
          IF (Aout(idTvar(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,it)=AVERAGE(ng)%avgt(i,j,k,it)+&
     &                                       OCEAN(ng)%t(i,j,k,Nout,it)
                END DO
              END DO
            END DO
          END IF
        END DO
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=AVERAGE(ng)%avgAKv(i,j,k)+    &
     &                                    MIXING(ng)%Akv(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idTdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKt(i,j,k)=AVERAGE(ng)%avgAKt(i,j,k)+    &
     &                                    MIXING(ng)%Akt(i,j,k,itemp)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idSdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKs(i,j,k)=AVERAGE(ng)%avgAKs(i,j,k)+    &
     &                                    MIXING(ng)%Akt(i,j,k,isalt)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idHsbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghsbl(i,j)=AVERAGE(ng)%avghsbl(i,j)+        &
     &                                 MIXING(ng)%hsbl(i,j)
            END DO
          END DO
        END IF
!
!  Accumulate 2D/3D coupling terms.
!
        IF (Aout(idUfx1,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg1(i,j)=AVERAGE(ng)%avgDU_avg1(i,j)+  &
     &                                    COUPLING(ng)%DU_avg1(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUfx2,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg2(i,j)=AVERAGE(ng)%avgDU_avg2(i,j)+  &
     &                                    COUPLING(ng)%DU_avg2(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVfx1,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg1(i,j)=AVERAGE(ng)%avgDV_avg1(i,j)+  &
     &                                    COUPLING(ng)%DV_avg1(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVfx2,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg2(i,j)=AVERAGE(ng)%avgDV_avg2(i,j)+  &
     &                                    COUPLING(ng)%DV_avg2(i,j)
            END DO
          END DO
        END IF
!
!  Accumulate surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=AVERAGE(ng)%avgsus(i,j)+          &
     &                                FORCES(ng)%sustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsvs(i,j)=AVERAGE(ng)%avgsvs(i,j)+          &
     &                                FORCES(ng)%svstr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgbus(i,j)=AVERAGE(ng)%avgbus(i,j)+          &
     &                                FORCES(ng)%bustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgbvs(i,j)=AVERAGE(ng)%avgbvs(i,j)+          &
     &                                FORCES(ng)%bvstr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idPair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPair(i,j)=AVERAGE(ng)%avgPair(i,j)+        &
     &                                 FORCES(ng)%Pair(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUwind(i,j)=AVERAGE(ng)%avgUwind(i,j)+      &
     &                                  FORCES(ng)%Uwind(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVwind(i,j)=AVERAGE(ng)%avgVwind(i,j)+      &
     &                                  FORCES(ng)%Vwind(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUairE,ng).and.Aout(idVairN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .TRUE., .FALSE.,                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
     &                      GRID(ng)%rmask_io,                          &
     &                      FORCES(ng) % Uwind,                         &
     &                      FORCES(ng) % Vwind,                         &
     &                      AVERAGE(ng)%avgUwindE,                      &
     &                      AVERAGE(ng)%avgVwindN)
        END IF
        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=AVERAGE(ng)%avgstf(i,j)+          &
     &                                FORCES(ng)%stflx(i,j,itemp)
            END DO
          END DO
        END IF
        IF (Aout(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgswf(i,j)=AVERAGE(ng)%avgswf(i,j)+          &
     &                                FORCES(ng)%stflx(i,j,isalt)
            END DO
          END DO
        END IF
        IF (Aout(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsrf(i,j)=AVERAGE(ng)%avgsrf(i,j)+          &
     &                                FORCES(ng)%srflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglhf(i,j)=AVERAGE(ng)%avglhf(i,j)+          &
     &                                FORCES(ng)%lhflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgshf(i,j)=AVERAGE(ng)%avgshf(i,j)+          &
     &                                FORCES(ng)%shflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglrf(i,j)=AVERAGE(ng)%avglrf(i,j)+          &
     &                                FORCES(ng)%lrflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgevap(i,j)=AVERAGE(ng)%avgevap(i,j)+        &
     &                                 FORCES(ng)%evap(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idrain,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgrain(i,j)=AVERAGE(ng)%avgrain(i,j)+        &
     &                                 FORCES(ng)%rain(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avguice(i,j)=AVERAGE(ng)%avguice(i,j)+        &
     &                                 ICE(ng)%ui(i,j,Iuout)
            END DO
          END DO
        END IF
        IF (Aout(idVice,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgvice(i,j)=AVERAGE(ng)%avgvice(i,j)+        &
     &                                 ICE(ng)%vi(i,j,Iuout)
            END DO
          END DO
        END IF
        IF (Aout(idUiceE,ng).and.Aout(idViceN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .TRUE., .FALSE.,                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
     &                      GRID(ng)%rmask_io,                          &
     &                      ICE(ng) % ui(:,:,Kout),                     &
     &                      ICE(ng) % vi(:,:,Kout),                     &
     &                      AVERAGE(ng)%avgUiceE,                       &
     &                      AVERAGE(ng)%avgViceN)
        END IF
        IF (Aout(idAice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgaice(i,j)=AVERAGE(ng)%avgaice(i,j)+        &
     &                                 ICE(ng)%ai(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idHice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghice(i,j)=AVERAGE(ng)%avghice(i,j)+        &
     &                                 ICE(ng)%hi(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idHsno,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghsno(i,j)=AVERAGE(ng)%avghsno(i,j)+        &
     &                                 ICE(ng)%hsn(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idTice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgtice(i,j)=AVERAGE(ng)%avgtice(i,j)+        &
     &                                 ICE(ng)%tis(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idTimid,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgtimid(i,j)=AVERAGE(ng)%avgtimid(i,j)+      &
     &                                 ICE(ng)%ti(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idSfwat,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsfwat(i,j)=AVERAGE(ng)%avgsfwat(i,j)+      &
     &                                 ICE(ng)%sfwat(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idAgeice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgageice(i,j)=AVERAGE(ng)%avgageice(i,j)+    &
     &                                 ICE(ng)%ageice(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idIomflx,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgiomflx(i,j)=AVERAGE(ng)%avgiomflx(i,j)+    &
     &                                 ICE(ng)%io_mflux(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idSig11,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsig11(i,j)=AVERAGE(ng)%avgsig11(i,j)+      &
     &                                 ICE(ng)%sig11(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout(idSig12,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsig12(i,j)=AVERAGE(ng)%avgsig12(i,j)+      &
     &                                 ICE(ng)%sig12(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout(idSig22,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsig22(i,j)=AVERAGE(ng)%avgsig22(i,j)+      &
     &                                 ICE(ng)%sig22(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout(idT0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgT0mk(i,j)=AVERAGE(ng)%avgT0mk(i,j)+        &
     &                                 ICE(ng)%t0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idS0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgS0mk(i,j)=AVERAGE(ng)%avgS0mk(i,j)+        &
     &                                 ICE(ng)%s0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWfr,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWfr(i,j)=AVERAGE(ng)%avgWfr(i,j)+          &
     &                                 ICE(ng)%wfr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWai,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWai(i,j)=AVERAGE(ng)%avgWai(i,j)+          &
     &                                 ICE(ng)%wai(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWao,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWao(i,j)=AVERAGE(ng)%avgWao(i,j)+          &
     &                                 ICE(ng)%wao(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWio,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWio(i,j)=AVERAGE(ng)%avgWio(i,j)+          &
     &                                 ICE(ng)%wio(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWro,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWro(i,j)=AVERAGE(ng)%avgWro(i,j)+          &
     &                                 ICE(ng)%wro(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idTauiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgutau_iw(i,j)=AVERAGE(ng)%avgutau_iw(i,j)+  &
     &                                 ICE(ng)%utau_iw(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idChuiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgchu_iw(i,j)=AVERAGE(ng)%avgchu_iw(i,j)+    &
     &                                 ICE(ng)%chu_iw(i,j)
            END DO
          END DO
        END IF
!
!  Accumulate vorticity fields.
!
        IF (Aout(id2dPV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgpvor2d(i,j)=AVERAGE(ng)%avgpvor2d(i,j)+    &
     &                                   potvor_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(id2dRV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgrvor2d(i,j)=AVERAGE(ng)%avgrvor2d(i,j)+    &
     &                                   relvor_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(id3dPV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgpvor3d(i,j,k)=AVERAGE(ng)%avgpvor3d(i,j, &
     &                                                             k)+  &
     &                                       potvor(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(id3dRV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgrvor3d(i,j,k)=AVERAGE(ng)%avgrvor3d(i,j, &
     &                                                             k)+  &
     &                                       relvor(i,j,k)
              END DO
            END DO
          END DO
        END IF
!
!  Accumulate quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=AVERAGE(ng)%avgZZ(i,j)+            &
     &                               OCEAN(ng)%zeta(i,j,Kout)*          &
     &                               OCEAN(ng)%zeta(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idU2av,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgU2(i,j)=AVERAGE(ng)%avgU2(i,j)+            &
     &                               OCEAN(ng)%ubar(i,j,Kout)*          &
     &                               OCEAN(ng)%ubar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idV2av,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgV2(i,j)=AVERAGE(ng)%avgV2(i,j)+            &
     &                               OCEAN(ng)%vbar(i,j,Kout)*          &
     &                               OCEAN(ng)%vbar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idUUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgUU(i,j,k)=AVERAGE(ng)%avgUU(i,j,k)+      &
     &                                   OCEAN(ng)%u(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%u(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgVV(i,j,k)=AVERAGE(ng)%avgVV(i,j,k)+      &
     &                                   OCEAN(ng)%v(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%v(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idUVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgUV(i,j,k)=AVERAGE(ng)%avgUV(i,j,k)+      &
     &                                   0.25_r8*                       &
     &                                   (OCEAN(ng)%u(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%u(i+1,j  ,k,Nout))* &
     &                                   (OCEAN(ng)%v(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%v(i  ,j+1,k,Nout))
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idHUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgHuon(i,j,k)=AVERAGE(ng)%avgHuon(i,j,k)+  &
     &                                     GRID(ng)%Huon(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idHVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgHvom(i,j,k)=AVERAGE(ng)%avgHvom(i,j,k)+  &
     &                                     GRID(ng)%Hvom(i,j,k)
              END DO
            END DO
          END DO
        END IF
        DO it=1,NAT
          IF (Aout(idTTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgTT(i,j,k,it)=AVERAGE(ng)%avgTT(i,j,k,  &
     &                                                          it)+    &
     &                                        OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)*     &
     &                                        OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgUT(i,j,k,it)=AVERAGE(ng)%avgUT(i,j,k,  &
     &                                                          it)+    &
     &                                        0.5_r8*                   &
     &                                        OCEAN(ng)%u(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i-1,j,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i  ,j,k,     &
     &                                                     Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgVT(i,j,k,it)=AVERAGE(ng)%avgVT(i,j,k,  &
     &                                                          it)+    &
     &                                        0.5_r8*                   &
     &                                        OCEAN(ng)%v(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i,j-1,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i,j  ,k,     &
     &                                                     Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(iHUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgHuonT(i,j,k,it)=AVERAGE(ng)%avgHuonT(i,&
     &                                                       j,k,it)+   &
     &                                           0.5_r8*                &
     &                                           GRID(ng)%Huon(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i-1,j,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i  ,j,k,  &
     &                                                        Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(iHVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgHvomT(i,j,k,it)=AVERAGE(ng)%avgHvomT(i,&
     &                                                       j,k,it)+   &
     &                                           0.5_r8*                &
     &                                           GRID(ng)%Hvom(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i,j-1,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i,j  ,k,  &
     &                                                        Nout,it))
                END DO
              END DO
            END DO
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Convert accumulated sums into time-averages, if appropriate.
!  Notice that we need to apply periodic conditions, if any, since
!  the full I- and J-ranges are different.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsAVG(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nAVG(ng)).eq.0).and.                          &
     &     ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))).or.          &
     &    ((iic(ng).ge.ntsAVG(ng)).and.(nAVG(ng).eq.1))) THEN
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (nAVG(ng).eq.1) THEN
            AVGtime(ng)=time(ng)
          ELSE
            AVGtime(ng)=AVGtime(ng)+REAL(nAVG(ng),r8)*dt(ng)
          END IF
        END IF
!
!  Set time-averaged factors for each C-grid variable type. Notice that
!  the I- and J-ranges are all grid types are the same for convinience.
!
        fac=1.0_r8/REAL(nAVG(ng),r8)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            pfac(i,j)=fac
            rfac(i,j)=fac
            ufac(i,j)=fac
            vfac(i,j)=fac
          END DO
        END DO
!
!  Process state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgzeta(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgzeta)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgzeta)
          END IF
        END IF
        IF (Aout(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2d(i,j)=ufac(i,j)*                        &
     &                                AVERAGE(ng)%avgu2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgu2d)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu2d)
          END IF
        END IF
        IF (Aout(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2d(i,j)=vfac(i,j)*                        &
     &                                AVERAGE(ng)%avgv2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgv2d)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv2d)
          END IF
        END IF
        IF (Aout(idu2dE,ng).and.Aout(idv2dN,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgu2dE(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgu2dE(i,j)
              AVERAGE(ng)%avgv2dN(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgv2dN(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgu2dE)
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgv2dN)
            CALL mp_exchange2d (ng, tile, iNLM, 2,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu2dE,                    &
     &                          AVERAGE(ng)%avgv2dN)
          END IF
        END IF
        IF (Aout(idUvel,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3d(i,j,k)=ufac(i,j)*                    &
     &                                    AVERAGE(ng)%avgu3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgu3d)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu3d)
          END IF
        END IF
        IF (Aout(idVvel,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3d(i,j,k)=vfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgv3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgv3d)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv3d)
          END IF
        END IF
        IF (Aout(idu3dE,ng).and.Aout(idv3dN,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgu3dE(i,j,k)=rfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgu3dE(i,j,k)
                AVERAGE(ng)%avgv3dN(i,j,k)=rfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgv3dN(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgu3dE)
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgv3dN)
            CALL mp_exchange3d (ng, tile, iNLM, 2,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu3dE,                    &
     &                          AVERAGE(ng)%avgv3dN)
          END IF
        END IF
        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgw3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              AVERAGE(ng)%avgw3d)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgw3d)
          END IF
        END IF
        IF (Aout(idWvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgwvel(i,j,k)=rfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgwvel(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              AVERAGE(ng)%avgwvel)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgwvel)
          END IF
        END IF
        IF (Aout(idDano,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgrho(i,j,k)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgrho(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgrho)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgrho)
          END IF
        END IF
        DO it=1,NT(ng)
          IF (Aout(idTvar(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,it)=rfac(i,j)*                 &
     &                                       AVERAGE(ng)%avgt(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgt(:,:,:,it))
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgt(:,:,:,it))
            END IF
          END IF
        END DO
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgAKv(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              AVERAGE(ng)%avgAKv)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgAKv)
          END IF
        END IF
        IF (Aout(idTdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKt(i,j,k)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgAKt(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              AVERAGE(ng)%avgAKt)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgAKt)
          END IF
        END IF
        IF (Aout(idSdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKs(i,j,k)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgAKs(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              AVERAGE(ng)%avgAKs)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgAKs)
          END IF
        END IF
        IF (Aout(idHsbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghsbl(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avghsbl(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avghsbl)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avghsbl)
          END IF
        END IF
!
!  Process 2D/3D coupling terms.
!
        IF (Aout(idUfx1,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg1(i,j)=ufac(i,j)*                    &
     &                                    AVERAGE(ng)%avgDU_avg1(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgDU_avg1)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgDU_avg1)
          END IF
        END IF
        IF (Aout(idUfx2,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg2(i,j)=ufac(i,j)*                    &
     &                                    AVERAGE(ng)%avgDU_avg2(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgDU_avg2)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgDU_avg2)
          END IF
        END IF
        IF (Aout(idVfx1,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg1(i,j)=vfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgDV_avg1(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgDV_avg1)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgDV_avg1)
          END IF
        END IF
        IF (Aout(idVfx2,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg2(i,j)=vfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgDV_avg2(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgDV_avg2)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgDV_avg2)
          END IF
        END IF
!
!  Process surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=ufac(i,j)*                        &
     &                                AVERAGE(ng)%avgsus(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgsus)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgsus)
          END IF
        END IF
        IF (Aout(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsvs(i,j)=vfac(i,j)*                        &
     &                                AVERAGE(ng)%avgsvs(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgsvs)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgsvs)
          END IF
        END IF
        IF (Aout(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgbus(i,j)=ufac(i,j)*                        &
     &                                AVERAGE(ng)%avgbus(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgbus)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgbus)
          END IF
        END IF
        IF (Aout(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgbvs(i,j)=vfac(i,j)*                        &
     &                                AVERAGE(ng)%avgbvs(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgbvs)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgbvs)
          END IF
        END IF
        IF (Aout(idPair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPair(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgPair(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgPair)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgPair)
          END IF
        END IF
        IF (Aout(idUair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUwind(i,j)=rfac(i,j)*                      &
     &                                  AVERAGE(ng)%avgUwind(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgUwind)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUwind)
          END IF
        END IF
        IF (Aout(idVair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVwind(i,j)=rfac(i,j)*                      &
     &                                  AVERAGE(ng)%avgVwind(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgVwind)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgVwind)
          END IF
        END IF
        IF (Aout(idUairE,ng).and.Aout(idVairN,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgUwindE(i,j)=rfac(i,j)*                     &
     &                                 AVERAGE(ng)%avgUwindE(i,j)
              AVERAGE(ng)%avgVwindN(i,j)=rfac(i,j)*                     &
     &                                 AVERAGE(ng)%avgVwindN(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgUwindE)
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgVwindN)
            CALL mp_exchange2d (ng, tile, iNLM, 2,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUwindE,                  &
     &                          AVERAGE(ng)%avgVwindN)
          END IF
        END IF
        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avgstf(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgstf)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgstf)
          END IF
        END IF
        IF (Aout(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgswf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avgswf(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgswf)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgswf)
          END IF
        END IF
        IF (Aout(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsrf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avgsrf(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgsrf)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgsrf)
          END IF
        END IF
        IF (Aout(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglhf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avglhf(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avglhf)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avglhf)
          END IF
        END IF
        IF (Aout(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgshf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avgshf(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgshf)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgshf)
          END IF
        END IF
        IF (Aout(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglrf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avglrf(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avglrf)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avglrf)
          END IF
        END IF
        IF (Aout(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgevap(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgevap(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgevap)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgevap)
          END IF
        END IF
        IF (Aout(idrain,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgrain(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgrain(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgrain)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgrain)
          END IF
        END IF
        IF (Aout(idUice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avguice(i,j)=ufac(i,j)*                       &
     &                AVERAGE(ng)%avguice(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avguice)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avguice)
          END IF
        END IF
        IF (Aout(idVice,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgvice(i,j)=vfac(i,j)*                       &
     &                    AVERAGE(ng)%avgvice(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgvice)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgvice)
          END IF
        END IF
        IF (Aout(idUiceE,ng).and.Aout(idViceN,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avguiceE(i,j)=rfac(i,j)*                      &
     &                                 AVERAGE(ng)%avguiceE(i,j)
              AVERAGE(ng)%avgviceN(i,j)=rfac(i,j)*                      &
     &                                 AVERAGE(ng)%avgviceN(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avguiceE)
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgviceN)
            CALL mp_exchange2d (ng, tile, iNLM, 2,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avguiceE,                   &
     &                          AVERAGE(ng)%avgviceN)
          END IF
        END IF
        IF (Aout(idAice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgaice(i,j)=rfac(i,j)*                       &
     &                       AVERAGE(ng)%avgaice(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgaice)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgaice)
          END IF
        END IF
        IF (Aout(idHice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghice(i,j)=rfac(i,j)*                       &
     &                          AVERAGE(ng)%avghice(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avghice)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avghice)
          END IF
        END IF
        IF (Aout(idHsno,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghsno(i,j)=rfac(i,j)*                       &
     &                        AVERAGE(ng)%avghsno(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avghsno)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avghsno)
          END IF
        END IF
        IF (Aout(idTice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgtice(i,j)=rfac(i,j)*                       &
     &                       AVERAGE(ng)%avgtice(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgtice)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgtice)
          END IF
        END IF
        IF (Aout(idTimid,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgtimid(i,j)=rfac(i,j)*                      &
     &                    AVERAGE(ng)%avgtimid(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgtimid)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgtimid)
          END IF
        END IF
        IF (Aout(idSfwat,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsfwat(i,j)=rfac(i,j)*                      &
     &                   AVERAGE(ng)%avgsfwat(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgsfwat)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgsfwat)
          END IF
        END IF
        IF (Aout(idAgeice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgageice(i,j)=rfac(i,j)*                     &
     &                     AVERAGE(ng)%avgageice(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgageice)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgageice)
          END IF
        END IF
        IF (Aout(idIomflx,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgiomflx(i,j)=rfac(i,j)*                     &
     &                     AVERAGE(ng)%avgiomflx(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgiomflx)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgiomflx)
          END IF
        END IF
        IF (Aout(idSig11,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsig11(i,j)=rfac(i,j)*                      &
     &                    AVERAGE(ng)%avgsig11(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgsig11)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgsig11)
          END IF
        END IF
        IF (Aout(idSig12,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsig12(i,j)=rfac(i,j)*                      &
     &                   AVERAGE(ng)%avgsig12(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgsig12)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgsig12)
          END IF
        END IF
        IF (Aout(idSig22,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsig22(i,j)=rfac(i,j)*                      &
     &                   AVERAGE(ng)%avgsig22(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgsig22)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgsig22)
          END IF
        END IF
        IF (Aout(idT0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgT0mk(i,j)=rfac(i,j)*                       &
     &                   AVERAGE(ng)%avgT0mk(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgT0mk)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgT0mk)
          END IF
        END IF
        IF (Aout(idS0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgS0mk(i,j)=rfac(i,j)*                       &
     &                   AVERAGE(ng)%avgS0mk(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgS0mk)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgS0mk)
          END IF
        END IF
        IF (Aout(idWfr,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWfr(i,j)=rfac(i,j)*                        &
     &                   AVERAGE(ng)%avgWfr(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgWfr)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgWfr)
          END IF
        END IF
        IF (Aout(idWai,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWai(i,j)=rfac(i,j)*                        &
     &                   AVERAGE(ng)%avgWai(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            AVERAGE(ng)%avgWai)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        AVERAGE(ng)%avgWai)
          END IF
        END IF
        IF (Aout(idWao,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWao(i,j)=rfac(i,j)*                        &
     &                   AVERAGE(ng)%avgWao(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWao)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWao)
          END IF
        END IF
        IF (Aout(idWio,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWio(i,j)=rfac(i,j)*                        &
     &                   AVERAGE(ng)%avgWio(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWio)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWio)
          END IF
        END IF
        IF (Aout(idWro,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWro(i,j)=rfac(i,j)*                        &
     &                   AVERAGE(ng)%avgWro(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWro)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWro)
          END IF
        END IF
        IF (Aout(idTauiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgutau_iw(i,j)=                              &
     &                      fac*AVERAGE(ng)%avgutau_iw(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgutau_iw)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgutau_iw)
          END IF
        END IF
        IF (Aout(idChuiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgchu_iw(i,j)=rfac(i,j)*                     &
     &                   AVERAGE(ng)%avgchu_iw(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgchu_iw)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgchu_iw)
          END IF
        END IF
!
!  Process vorticity fields.
!
        IF (Aout(id2dPV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgpvor2d(i,j)=pfac(i,j)*                     &
     &                                   AVERAGE(ng)%avgpvor2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgpvor2d)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgpvor2d)
          END IF
        END IF
        IF (Aout(id2dRV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgrvor2d(i,j)=pfac(i,j)*                     &
     &                                   AVERAGE(ng)%avgrvor2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgrvor2d)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgrvor2d)
          END IF
        END IF
        IF (Aout(id3dPV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgpvor3d(i,j,k)=pfac(i,j)*                 &
     &                                      AVERAGE(ng)%avgpvor3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgpvor3d)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgpvor3d)
          END IF
        END IF
        IF (Aout(id3dRV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgrvor3d(i,j,k)=pfac(i,j)*                 &
     &                                      AVERAGE(ng)%avgrvor3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgrvor3d)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgrvor3d)
          END IF
        END IF
!
!  Process quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=rfac(i,j)*                         &
     &                               AVERAGE(ng)%avgZZ(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgZZ)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgZZ)
          END IF
        END IF
        IF (Aout(idU2av,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgU2(i,j)=ufac(i,j)*                         &
     &                               AVERAGE(ng)%avgU2(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgU2)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgU2)
          END IF
        END IF
        IF (Aout(idV2av,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgV2(i,j)=vfac(i,j)*                         &
     &                               AVERAGE(ng)%avgV2(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgV2)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgV2)
          END IF
        END IF
        IF (Aout(idUUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgUU(i,j,k)=ufac(i,j)*                     &
     &                                   AVERAGE(ng)%avgUU(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgUU)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUU)
          END IF
        END IF
        IF (Aout(idVVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgVV(i,j,k)=vfac(i,j)*                     &
     &                                   AVERAGE(ng)%avgVV(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgVV)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgVV)
          END IF
        END IF
        IF (Aout(idUVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgUV(i,j,k)=rfac(i,j)*                     &
     &                                   AVERAGE(ng)%avgUV(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgUV)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUV)
          END IF
        END IF
        IF (Aout(idHUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgHuon(i,j,k)=ufac(i,j)*                   &
     &                                     AVERAGE(ng)%avgHuon(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgHuon)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgHuon)
          END IF
        END IF
        IF (Aout(idHVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgHvom(i,j,k)=vfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgHvom(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgHvom)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgHvom)
          END IF
        END IF
        DO it=1,NAT
          IF (Aout(idTTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgTT(i,j,k,it)=rfac(i,j)*                &
     &                                       AVERAGE(ng)%avgTT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgTT(:,:,:,it))
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgTT(:,:,:,it))
            END IF
          END IF
          IF (Aout(idUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgUT(i,j,k,it)=ufac(i,j)*                &
     &                                       AVERAGE(ng)%avgUT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgUT(:,:,:,it))
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgUT(:,:,:,it))
            END IF
          END IF
          IF (Aout(idVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgVT(i,j,k,it)=vfac(i,j)*                &
     &                                       AVERAGE(ng)%avgVT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgVT(:,:,:,it))
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgVT(:,:,:,it))
            END IF
          END IF
          IF (Aout(iHUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgHuonT(i,j,k,it)=ufac(i,j)*             &
     &                                    AVERAGE(ng)%avgHuonT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgHuonT(:,:,:,it))
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgHuonT(:,:,:,it))
            END IF
          END IF
          IF (Aout(iHVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgHvomT(i,j,k,it)=vfac(i,j)*             &
     &                                    AVERAGE(ng)%avgHvomT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgHvomT(:,:,:,it))
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgHvomT(:,:,:,it))
            END IF
          END IF
        END DO
      END IF
      RETURN
      END SUBROUTINE set_avg_tile
      END MODULE set_avg_mod
