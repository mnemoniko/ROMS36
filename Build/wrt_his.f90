      SUBROUTINE wrt_his (ng)
!
!svn $Id: wrt_his.F 1484 2012-06-13 19:25:03Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes requested model fields at requested levels      !
!  into history NetCDF file.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_coupling
      USE mod_forces
      USE mod_grid
      USE mod_ice
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
      USE omega_mod, ONLY : scale_omega
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
      integer :: LBi, UBi, LBj, UBj
      integer :: Fcount, gfactor, gtype, status, tile
      integer :: i, itrc, j, k
      real(r8) :: scale
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
      real(r8), allocatable :: Ur3d(:,:,:)
      real(r8), allocatable :: Vr3d(:,:,:)
      real(r8), allocatable :: Wr3d(:,:,:)
!
      SourceFile='wrt_his.F'
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out history fields.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
      gfactor=1
!
!  Set time record index.
!
      HIS(ng)%Rindex=HIS(ng)%Rindex+1
      Fcount=HIS(ng)%Fcount
      HIS(ng)%Nrec(Fcount)=HIS(ng)%Nrec(Fcount)+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, HIS(ng)%name,                     &
     &                      TRIM(Vname(idtime,ng)), time(ng:),          &
     &                      (/HIS(ng)%Rindex/), (/1/),                  &
     &                      ncid = HIS(ng)%ncid,                        &
     &                      varid = HIS(ng)%Vid(idtime))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out free-surface (m)
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idFsur), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     OCEAN(ng) % zeta(:,:,kstp(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D U-momentum component (m/s).
!
      IF (Hout(idUbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUbar), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask_io,                         &
     &                     OCEAN(ng) % ubar(:,:,kstp(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUfx1), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask_io,                         &
     &                     COUPLING(ng) % DU_avg1)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUfx1)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUfx2), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask_io,                         &
     &                     COUPLING(ng) % DU_avg2)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUfx2)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D V-momentum component (m/s).
!
      IF (Hout(idVbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVbar), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask_io,                         &
     &                     OCEAN(ng) % vbar(:,:,kstp(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVfx1), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask_io,                         &
     &                     COUPLING(ng) % DV_avg1)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVfx1)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVfx2), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask_io,                         &
     &                     COUPLING(ng) % DV_avg2)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVfx2)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Hout(idu2dE,ng).and.Hout(idv2dN,ng)) THEN
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
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idu2dE), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     Ur2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu2dE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idv2dN), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     Vr2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv2dN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
!
!  Write out 3D U-momentum component (m/s).
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUvel), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % umask_io,                         &
     &                     OCEAN(ng) % u(:,:,:,nrhs(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D V-momentum component (m/s).
!
      IF (Hout(idVvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVvel), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % vmask_io,                         &
     &                     OCEAN(ng) % v(:,:,:,nrhs(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Hout(idu3dE,ng).and.Hout(idv3dN,ng)) THEN
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
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idu3dE), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask_io,                         &
     &                     Ur3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu3dE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idv3dN), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask_io,                         &
     &                     Vr3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv3dN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur3d)
        deallocate (Vr3d)
      END IF
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Hout(idOvel,ng)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate (Wr3d(LBi:UBi,LBj:UBj,0:N(ng)))
          Wr3d(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_r8
        gtype=gfactor*w3dvar
        tile=MyRank
        CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    OCEAN(ng) % W,                                &
     &                    Wr3d)
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idOvel), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask_io,                         &
     &                     Wr3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out vertical velocity (m/s).
!
      IF (Hout(idWvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idWvel), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask_io,                         &
     &                     OCEAN(ng) % wvel)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (Hout(idTvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Tid(itrc), &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % rmask_io,                       &
     &                       OCEAN(ng) % t(:,:,:,nrhs(ng),itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),            &
     &                          HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!----------------------------
!  Write out density anomaly.
!----------------------------
      IF (Hout(idDano,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idDano), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask_io,                         &
     &                     OCEAN(ng) % rho)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out depth surface boundary layer.
!
      IF (Hout(idHsbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idHsbl), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     MIXING(ng) % hsbl)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsbl)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out KPP nonlocal transport.
!
      DO i=1,NAT
        IF (Hout(idGhat(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*w3dvar
          status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid,                    &
     &                       HIS(ng)%Vid(idGhat(i)),                    &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 0, N(ng), scale,       &
     &                       GRID(ng) % rmask_io,                       &
     &                       MIXING(ng) % ghats(:,:,:,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idGhat(i))), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out vertical viscosity coefficient.
!
      IF (Hout(idVvis,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVvis), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask_io,                         &
     &                     MIXING(ng) % Akv,                            &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Hout(idTdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idTdif), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask_io,                         &
     &                     MIXING(ng) % Akt(:,:,:,itemp),               &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Hout(idSdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idSdif), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask_io,                         &
     &                     MIXING(ng) % Akt(:,:,:,isalt),               &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice 2D momentum component (m/s) in the XI-direction.
!
      IF (Hout(idUice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idUice),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask_io,                         &
     &                     ICE(ng) % ui(:,:,liunw(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUice)),                    &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice 2D momentum component (m/s) in the ETA-direction.
!
      IF (Hout(idVice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idVice),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask_io,                         &
     &                     ICE(ng) % vi(:,:,liunw(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVice)),                    &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Eastward and Northward ice momentum components (m/s) at
!  RHO-points.
!
      IF (Hout(idUiceE,ng).and.Hout(idViceN,ng)) THEN
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
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUiceE),&
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     Ur2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUiceE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idViceN),&
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     Vr2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idViceN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
!
!  Write out ice concentration
!
      IF (Hout(idAice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idAice),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % ai(:,:,liunw(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idAice)),                    &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice average thickness
!
      IF (Hout(idHice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idHice),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % hi(:,:,liunw(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHice)),                    &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out snow average thickness
!
      IF (Hout(idHsno,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idHsno),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % hsn(:,:,liunw(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsno)),                    &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface water thickness (on ice)
!
      IF (Hout(idSfwat,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idSfwat),                        &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % sfwat(:,:,liunw(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSfwat)),                   &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice age.
!
      IF (Hout(idAgeice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idAgeice),                       &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % ageice(:,:,liunw(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idAgeice)),                  &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice-ocean mass flux
!
      IF (Hout(idIomflx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idIomflx),                       &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % io_mflux)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idIomflx)),                  &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice/snow surface temperature
!
      IF (Hout(idTice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idTice),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % tis)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTice)),                    &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice interior temperature
!
      IF (Hout(idTimid,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idTimid),                        &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % ti(:,:,liunw(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTimid)),                   &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out internal ice stress component 11
!
      IF (Hout(idSig11,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idSig11),                        &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % sig11(:,:,liunw(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig11)),                   &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out internal ice stress component 12
!
      IF (Hout(idSig12,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idSig12),                        &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % sig12(:,:,liunw(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig12)),                   &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out internal ice stress component 22
!
      IF (Hout(idSig22,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idSig22),                        &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % sig22(:,:,liunw(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig22)),                   &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice-ocean friction velocity
!
      IF (Hout(idTauiw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idTauiw),                        &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % utau_iw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTauiw)),                   &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice-ocean momentum transfer coefficient
!
      IF (Hout(idChuiw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idChuiw),                        &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % chu_iw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idChuiw)),                   &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out temperature of molecular sublayer under ice
!
      IF (Hout(idT0mk,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idT0mk),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % t0mk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idT0mk)),                    &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out salinity of molecular sublayer under ice
!
      IF (Hout(idS0mk,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idS0mk),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % s0mk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idS0mk)),                    &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice freeze Wfr
!
      IF (Hout(idWfr,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idWfr),                          &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % wfr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWfr)),                     &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice melt/freeze wai
!
      IF (Hout(idWai,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idWai),                          &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % wai)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWai)),                     &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice melt/freeze Wao
!
      IF (Hout(idWao,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idWao),                          &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % wao)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWao)),                     &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice melt/freeze wio
!
      IF (Hout(idWio,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idWio),                          &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % wio)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWio)),                     &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice melt/freeze wro
!
      IF (Hout(idWro,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idWro),                          &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     ICE(ng) % wro)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWro)),                     &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface air pressure.
!
      IF (Hout(idPair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idPair), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     FORCES(ng) % Pair)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idPair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface winds.
!
      IF (Hout(idUair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUair), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     FORCES(ng) % Uwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
      IF (Hout(idVair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVair), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     FORCES(ng) % Vwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Eastward and Northward surface winds (m/s) at
!  RHO-points.
!
      IF (Hout(idUairE,ng).and.Hout(idVairN,ng)) THEN
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
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUairE),&
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     Ur2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUairE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVairN),&
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     Vr2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVairN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
!
!  Write out surface active traces fluxes.
!
      DO itrc=1,NAT
        IF (Hout(idTsur(itrc),ng)) THEN
          IF (itrc.eq.itemp) THEN
            scale=rho0*Cp                   ! Celsius m/s to W/m2
          ELSE IF (itrc.eq.isalt) THEN
            scale=1.0_r8
          END IF
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                    &
     &                       HIS(ng)%Vid(idTsur(itrc)),                 &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask_io,                       &
     &                       FORCES(ng) % stflx(:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out latent heat flux.
!
      IF (Hout(idLhea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idLhea), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     FORCES(ng) % lhflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLhea)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Hout(idShea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idShea), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     FORCES(ng) % shflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idShea)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Hout(idLrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idLrad), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     FORCES(ng) % lrflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLrad)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out E-P (m/s).
!
      IF (Hout(idEmPf,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idEmPf), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     FORCES(ng) % EminusP)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idEmPf)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out evaporation rate (kg/m2/s).
!
      IF (Hout(idevap,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idevap), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     FORCES(ng) % evap)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idevap)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Hout(idrain,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idrain), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     FORCES(ng) % rain)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idrain)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out shortwave radiation flux.
!
      IF (Hout(idSrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idSrad), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_io,                         &
     &                     FORCES(ng) % srflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSrad)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface U-momentum stress.
!
      IF (Hout(idUsms,ng)) THEN
        scale=rho0                          ! m2/s2 to Pa
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUsms), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask_io,                         &
     &                     FORCES(ng) % sustr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Hout(idVsms,ng)) THEN
        scale=rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVsms), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask_io,                         &
     &                     FORCES(ng) % svstr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Hout(idUbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idUbms), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask_io,                         &
     &                     FORCES(ng) % bustr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Hout(idVbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, HIS(ng)%Vid(idVbms), &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask_io,                         &
     &                     FORCES(ng) % bvstr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize history NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, HIS(ng)%name, HIS(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
      IF (Master) WRITE (stdout,20) kstp(ng), nrhs(ng), HIS(ng)%Rindex
!
  10  FORMAT (/,' WRT_HIS - error while writing variable: ',a,/,11x,    &
     &        'into history NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_HIS   - wrote history  fields (Index=', i1,       &
     &        ',',i1,') into time record = ',i7.7)
      RETURN
      END SUBROUTINE wrt_his
