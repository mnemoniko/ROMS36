      SUBROUTINE wrt_station (ng)
!
!svn $Id: wrt_station.F 1484 2012-06-13 19:25:03Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes out data into stations NetCDF file.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ice
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE extract_sta_mod, ONLY : extract_sta2d, extract_sta3d
      USE uv_rotate_mod, ONLY : uv_rotate2d
      USE uv_rotate_mod, ONLY : uv_rotate3d
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical :: Cgrid
      integer :: NposB, NposR, NposW, LBi, UBi, LBj, UBj
      integer :: Fcount, i, ifield, k, np, status, tile
      real(r8) :: scale
      real(r8), dimension(Nstation(ng)) :: Xpos, Ypos, Zpos, psta
      real(r8), dimension(Nstation(ng)*(N(ng))) :: XposR, YposR, ZposR
      real(r8), dimension(Nstation(ng)*(N(ng)+1)) :: XposW, YposW, ZposW
      real(r8), dimension(Nstation(ng)*(N(ng)+1)) :: rsta
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
      real(r8), allocatable :: Ur3d(:,:,:)
      real(r8), allocatable :: Vr3d(:,:,:)
!
      SourceFile='wrt_station.F'
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out station data at RHO-points.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set time record index.
!
      STA(ng)%Rindex=STA(ng)%Rindex+1
      Fcount=STA(ng)%Fcount
      STA(ng)%Nrec(Fcount)=STA(ng)%Nrec(Fcount)+1
!
!  Set switch to extract station data at native C-grid position (TRUE)
!  or at RHO-points (FALSE).
!
      Cgrid=.FALSE.
!
!  Set positions for generic extraction routine.
!
      NposB=Nstation(ng)*Nbed
      NposR=Nstation(ng)*N(ng)
      NposW=Nstation(ng)*(N(ng)+1)
      DO i=1,Nstation(ng)
        Xpos(i)=SCALARS(ng)%SposX(i)
        Ypos(i)=SCALARS(ng)%SposY(i)
        Zpos(i)=1.0_r8
        DO k=1,N(ng)
          np=k+(i-1)*N(ng)
          XposR(np)=SCALARS(ng)%SposX(i)
          YposR(np)=SCALARS(ng)%SposY(i)
          ZposR(np)=REAL(k,r8)
        END DO
        DO k=0,N(ng)
          np=k+1+(i-1)*(N(ng)+1)
          XposW(np)=SCALARS(ng)%SposX(i)
          YposW(np)=SCALARS(ng)%SposY(i)
          ZposW(np)=REAL(k,r8)
        END DO
      END DO
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                     &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/STA(ng)%Rindex/), (/1/),                  &
     &                      ncid = STA(ng)%ncid,                        &
     &                      varid = STA(ng)%Vid(idtime))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out free-surface (m).
!
      IF (Sout(idFsur,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idFsur, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%zeta(:,:,kstp(ng)),        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idFsur)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idFsur))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUbar,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbar, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%ubar(:,:,kstp(ng)),        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUbar)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbar))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Sout(idVbar,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbar, v2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%vbar(:,:,kstp(ng)),        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVbar)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbar))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out 2D Eastward and Northward momentum components (m/s) at
!  RHO-points
!
      IF (Sout(idu2dE,ng).and.Sout(idv2dN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
            Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
            Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        tile=MyRank
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
     &                    GRID(ng) % rmask_io,                          &
     &                    OCEAN(ng) % ubar(:,:,kstp(ng)),               &
     &                    OCEAN(ng) % vbar(:,:,kstp(ng)),               &
     &                    Ur2d, Vr2d)
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idu2dE, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Ur2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idu2dE)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idu2dE))
        IF (exit_flag.ne.NoError) RETURN
        CALL extract_sta2d (ng, iNLM, Cgrid, idv2dN, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Vr2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idv2dN)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idv2dN))
        IF (exit_flag.ne.NoError) RETURN
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUvel,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idUvel, u3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%u(:,:,:,nrhs(ng)),         &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUvel)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUvel))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Sout(idVvel,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idVvel, v3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%v(:,:,:,nrhs(ng)),         &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVvel)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVvel))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out 3D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Sout(idu3dE,ng).and.Sout(idv3dN,ng)) THEN
        IF (.not.allocated(Ur3d)) THEN
          allocate (Ur3d(LBi:UBi,LBj:UBj,N(ng)))
          Ur3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        IF (.not.allocated(Vr3d)) THEN
          allocate (Vr3d(LBi:UBi,LBj:UBj,N(ng)))
          Vr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        tile=MyRank
        CALL uv_rotate3d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
     &                    GRID(ng) % rmask_io,                          &
     &                    OCEAN(ng) % u(:,:,:,nrhs(ng)),                &
     &                    OCEAN(ng) % v(:,:,:,nrhs(ng)),                &
     &                    Ur3d, Vr3d)
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idu3dE, r3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Ur3d,                                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idu3dE)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idu3dE))
        IF (exit_flag.ne.NoError) RETURN
        CALL extract_sta3d (ng, iNLM, Cgrid, idv3dN, r3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Vr3d,                                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idv3dN)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idv3dN))
        IF (exit_flag.ne.NoError) RETURN
        deallocate (Ur3d)
        deallocate (Vr3d)
      END IF
!
!  Write out vertical velocity (m/s).
!
      IF (Sout(idWvel,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idWvel, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng)%wvel,                      &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWvel)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWvel))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write S-coordinate "omega" vertical velocity (m3/s).
!
      IF (Sout(idOvel,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idOvel, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng)%W,                         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idOvel)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idOvel))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out tracer type variables.
!
      DO i=1,NT(ng)
        ifield=idTvar(i)
        IF (Sout(ifield,ng)) THEN
          scale=1.0_r8
          CALL extract_sta3d (ng, iNLM, Cgrid, ifield, r3dvar,          &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        scale, OCEAN(ng)%t(:,:,:,nrhs(ng),i),     &
     &                        NposR, XposR, YposR, ZposR, rsta)
          CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                 &
     &                          TRIM(Vname(1,idTvar(i))), rsta,         &
     &                          (/1,1,STA(ng)%Rindex/),                 &
     &                          (/N(ng),Nstation(ng),1/),               &
     &                          ncid = STA(ng)%ncid,                    &
     &                          varid = STA(ng)%Tid(i))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Sout(idDano,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idDano, r3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%rho,                       &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idDano)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idDano))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out depth of surface boundary layer.
!
      IF (Sout(idHsbl,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idHsbl, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng)%hsbl,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idHsbl)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idHsbl))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out vertical viscosity coefficient.
!
      IF (Sout(idVvis,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idVvis, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akv,                      &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVvis)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVvis))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Sout(idTdif,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idTdif, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akt(:,:,:,itemp),         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idTdif)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idTdif))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Sout(idSdif,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idSdif, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akt(:,:,:,isalt),         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idSdif)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idSdif))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out surface air pressure.
!
      IF (Sout(idPair,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idPair, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%Pair,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idPair)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idPair))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out surface winds.
!
      IF (Sout(idUair,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idUair, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%Uwind,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUair)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUair))
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (Sout(idVair,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idVair, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%Vwind,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVair)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVair))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out 2D Eastward and Northward surface winds (m/s) at
!  RHO-points
!
      IF (Sout(idUairE,ng).and.Sout(idVairN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
            Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
            Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        tile=MyRank
        CALL uv_rotate2d (ng, TILE, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
     &                    GRID(ng) % rmask_io,                          &
     &                    FORCES(ng) % Uwind,                           &
     &                    FORCES(ng) % Vwind,                           &
     &                    Ur2d, Vr2d)
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idUairE, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Ur2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUairE)), psta,             &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUairE))
        IF (exit_flag.ne.NoError) RETURN
        CALL extract_sta2d (ng, iNLM, Cgrid, idVairN, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Vr2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVairN)), psta,             &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVairN))
        IF (exit_flag.ne.NoError) RETURN
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
!
!  Write out surface net heat flux.
!
      IF (Sout(idTsur(itemp),ng)) THEN
        ifield=idTsur(itemp)
        scale=rho0*Cp
        CALL extract_sta2d (ng, iNLM, Cgrid, ifield, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflx(:,:,itemp),         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,ifield)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(ifield))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out surface salt flux.
!
      IF (Sout(idTsur(isalt),ng)) THEN
        ifield=idTsur(isalt)
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, ifield, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflx(:,:,isalt),         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,ifield)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(ifield))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out latent heat flux.
!
      IF (Sout(idLhea,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, iNLM, Cgrid, idLhea, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%lhflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idLhea)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idLhea))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out sensible heat flux.
!
      IF (Sout(idShea,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, iNLM, Cgrid, idShea, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%shflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idShea)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idShea))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Sout(idLrad,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, iNLM, Cgrid, idLrad, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%lrflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idLrad)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idLrad))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out shortwave radiation flux.
!
      IF (Sout(idSrad,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, iNLM, Cgrid, idSrad, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%srflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idSrad)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idSrad))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out E-P (m/s).
!
      IF (Sout(idEmPf,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idEmPf, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%EminusP,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idEmPf)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idEmPf))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out evaportaion rate (kg/m2/s).
!
      IF (Sout(idevap,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idevap, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%evap,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idevap)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idevap))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Sout(idrain,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idrain, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%rain,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idrain)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idrain))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out surface U-momentum stress.
!
      IF (Sout(idUsms,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idUsms, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%sustr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUsms)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUsms))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Sout(idVsms,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idVsms, v2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%svstr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVsms)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVsms))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Sout(idUbms,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbms, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%bustr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUbms)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbms))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Sout(idVbms,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbms, v2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%bvstr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVbms)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbms))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out ice 2D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUice,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idUice, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%ui(:,:,liunw(ng)),           &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUice)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUice))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out ice 2D momentum component (m/s) in the ETA-direction.
!
      IF (Sout(idVice,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idVice, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%vi(:,:,liunw(ng)),           &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVice)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVice))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out 2D Eastward and Northward ice momentum components (m/s) at
!  RHO-points
!
      IF (Sout(idUiceE,ng).and.Sout(idViceN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
            Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
            Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        tile=MyRank
        CALL uv_rotate2d (ng, TILE, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
     &                    GRID(ng) % rmask_io,                          &
     &                    ICE(ng) % ui(:,:,liunw(ng)),                  &
     &                    ICE(ng) % vi(:,:,liunw(ng)),                  &
     &                    Ur2d, Vr2d)
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idUiceE, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Ur2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUiceE)), psta,             &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUiceE))
        IF (exit_flag.ne.NoError) RETURN
        CALL extract_sta2d (ng, iNLM, Cgrid, idViceN, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Vr2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idViceN)), psta,             &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idViceN))
        IF (exit_flag.ne.NoError) RETURN
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
!
!  Write out ice concentration
!
      IF (Sout(idAice,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idAice, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%ai(:,:,linew(ng)),           &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idAice)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idAice))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out ice average thickness
!
      IF (Sout(idHice,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idHice, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%hi(:,:,linew(ng)),           &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idHice)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idHice))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out snow average thickness
!
      IF (Sout(idHsno,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idHsno, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%hsn(:,:,linew(ng)),          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idHsno)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idHsno))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out surface water thickness (on ice)
!
      IF (Sout(idSfwat,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idSfwat, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%sfwat(:,:,linew(ng)),        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idSfwat)), psta,             &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idSfwat))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out ice age.
!
      IF (Sout(idAgeice,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idAgeice, r2dvar,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%ageice(:,:,linew(ng)),       &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idAgeice)), psta,            &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idAgeice))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out ice-ocean mass flux
!
      IF (Sout(idIomflx,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idIomflx, r2dvar,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%io_mflux(:,:),               &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idIomflx)), psta,            &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idIomflx))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out ice/snow surface temperature
!
      IF (Sout(idTice,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idTice, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%tis,                         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idTice)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idTice))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out ice interior temperature
!
      IF (Sout(idTimid,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idTimid, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%ti(:,:,linew(ng)),           &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idTimid)), psta,             &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idTimid))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out ice-water friction velocity
!
      IF (Sout(idTauiw,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idTauiw, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%utau_iw,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idTauiw)), psta,             &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idTauiw))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out ice-water momentum transfer coefficient
!
      IF (Sout(idChuiw,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idChuiw, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%chu_iw,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idChuiw)), psta,             &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idChuiw))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out under-ice temperature
!
      IF (Sout(idT0mk,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idT0mk, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%t0mk,                        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idT0mk)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idT0mk))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out under-ice salinity
!
      IF (Sout(idS0mk,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idS0mk, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%s0mk,                        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idS0mk)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idS0mk))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out frazil ice growth rate
!
      IF (Sout(idWfr,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idWfr, r2dvar,             &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%wfr,                         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWfr)), psta,               &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWfr))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out ice growth/melt rate
!
      IF (Sout(idWai,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idWfr, r2dvar,             &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%wai,                         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWai)), psta,               &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWai))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out ice growth/melt rate
!
      IF (Sout(idWao,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idWao, r2dvar,             &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%wao,                         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWao)), psta,               &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWao))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out ice growth/melt rate
!
      IF (Sout(idWio,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idWio, r2dvar,             &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%wio,                         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWio)), psta,               &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWio))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out ice melt runoff rate
!
      IF (Sout(idWro,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idWro, r2dvar,             &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%wro,                         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWro)), psta,               &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWro))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out internal ice stress sig11
!
      IF (Sout(idSig11,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idSig11, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%sig11(:,:,lienw(ng)),        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idSig11)), psta,             &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idSig11))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out internal ice stress sig12
!
      IF (Sout(idSig12,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idSig12, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%sig12(:,:,lienw(ng)),        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idSig12)), psta,             &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idSig12))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out internal ice stress sig22
!
      IF (Sout(idSig22,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idSig22, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%sig22(:,:,lienw(ng)),        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idSig22)), psta,             &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idSig22))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out freezing water wfr
!
      IF (Sout(idWfr,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idWfr, r2dvar,             &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, ICE(ng)%wfr(:,:),                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWfr)), psta,               &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWfr))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize stations NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, STA(ng)%name, STA(ng)%ncid)
      RETURN
      END SUBROUTINE wrt_station
