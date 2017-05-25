      SUBROUTINE wrt_rst (ng)
!
!svn $Id: wrt_rst.F 1451 2012-02-02 20:56:14Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes fields into restart NetCDF file.                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ice
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
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
      integer :: Fcount, gfactor, gtype, i, itrc, status, varid
      integer :: ntmp(1)
      real(r8) :: scale
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
      SourceFile='wrt_rst.F'
!
!-----------------------------------------------------------------------
!  Write out restart fields.
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
      RST(ng)%Rindex=RST(ng)%Rindex+1
      Fcount=RST(ng)%Fcount
      RST(ng)%Nrec(Fcount)=RST(ng)%Nrec(Fcount)+1
!
!  If requested, set time index to recycle time records in restart
!  file.
!
      IF (LcycleRST(ng)) THEN
        RST(ng)%Rindex=MOD(RST(ng)%Rindex-1,2)+1
      END IF
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, RST(ng)%name,                     &
     &                      TRIM(Vname(idtime,ng)), time(ng:),          &
     &                      (/RST(ng)%Rindex/), (/1/),                  &
     &                      ncid = RST(ng)%ncid,                        &
     &                      varid = RST(ng)%Vid(idtime))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out free-surface (m).
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idFsur),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % rmask_io,                           &
     &                   OCEAN(ng) % zeta(:,:,kstp(ng)))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUbar),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % umask_io,                           &
     &                   OCEAN(ng) % ubar(:,:,kstp(ng)))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVbar),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % vmask_io,                           &
     &                   OCEAN(ng) % vbar(:,:,kstp(ng)))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUvel),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
     &                   GRID(ng) % umask_io,                           &
     &                   OCEAN(ng) % u(:,:,:,nrhs(ng)))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVvel),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
     &                   GRID(ng) % vmask_io,                           &
     &                   OCEAN(ng) % v(:,:,:,nrhs(ng)))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Tid(itrc),   &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask_io,                         &
     &                     OCEAN(ng) % t(:,:,:,nrhs(ng),itrc))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
!
!  Write out density anomaly.
!
      scale=1.0_r8
      gtype=gfactor*r3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idDano),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
     &                   GRID(ng) % rmask_io,                           &
     &                   OCEAN(ng) % rho)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idDano)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out depth of surface boundary layer.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idHsbl),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % rmask_io,                           &
     &                   MIXING(ng) % hsbl)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idHsbl)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out vertical viscosity coefficient.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVvis),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
     &                   GRID(ng) % rmask_io,                           &
     &                   MIXING(ng) % Akv,                              &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvis)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTdif),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
     &                   GRID(ng) % rmask_io,                           &
     &                   MIXING(ng) % Akt(:,:,:,itemp),                 &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTdif)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out vertical diffusion coefficient for salinity.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idSdif),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
     &                   GRID(ng) % rmask_io,                           &
     &                   MIXING(ng) % Akt(:,:,:,isalt),                 &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idSdif)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice 2D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUice),   &
     &                 RST(ng)%Rindex, gtype,                           &
     &                 LBi, UBi, LBj, UBj, scale,                       &
     &                 GRID(ng) % umask,                                &
     &                 ICE (ng) % ui(:,:,liunw(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idUice)), RST(ng)%Rindex
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice 2D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idvice),   &
     &                 RST(ng)%Rindex, gtype,                           &
     &                 LBi, UBi, LBj, UBj, scale,                       &
     &                 GRID(ng) % vmask,                                &
     &                 ICE (ng) % vi(:,:,liunw(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVice)), RST(ng)%Rindex
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice concentration
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idAice),   &
     &                 RST(ng)%Rindex, gtype,                           &
     &                 LBi, UBi, LBj, UBj, scale,                       &
     &                 GRID(ng) % rmask_io,                             &
     &                 ICE (ng) % ai(:,:,linew(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idAice)), RST(ng)%Rindex
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice average thickness
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idHice),   &
     &                 RST(ng)%Rindex, gtype,                           &
     &                 LBi, UBi, LBj, UBj, scale,                       &
     &                 GRID(ng) % rmask_io,                             &
     &                 ICE (ng) % hi(:,:,linew(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idHice)), RST(ng)%Rindex
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out snow average thickness
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idHsno),   &
     &                 RST(ng)%Rindex, gtype,                           &
     &                 LBi, UBi, LBj, UBj, scale,                       &
     &                 GRID(ng) % rmask_io,                             &
     &                 ICE (ng) % hsn(:,:,linew(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idHsno)), RST(ng)%Rindex
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface water thickness (on ice)
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idSfwat),  &
     &                 RST(ng)%Rindex, gtype,                           &
     &                 LBi, UBi, LBj, UBj, scale,                       &
     &                 GRID(ng) % rmask_io,                             &
     &                 ICE (ng) % sfwat(:,:,linew(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idSfwat)), RST(ng)%Rindex
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice age
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idAgeice), &
     &                 RST(ng)%Rindex, gtype,                           &
     &                 LBi, UBi, LBj, UBj, scale,                       &
     &                 GRID(ng) % rmask_io,                             &
     &                 ICE (ng) % ageice(:,:,linew(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idAgeice)), RST(ng)%Rindex
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice/snow surface temperature
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTice),   &
     &                 RST(ng)%Rindex, gtype,                           &
     &                 LBi, UBi, LBj, UBj, scale,                       &
     &                 GRID(ng) % rmask_io,                             &
     &                 ICE (ng) % tis)
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idTice)), RST(ng)%Rindex
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice interior temperature
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTimid),  &
     &                 RST(ng)%Rindex, gtype,                           &
     &                 LBi, UBi, LBj, UBj, scale,                       &
     &                 GRID(ng) % rmask_io,                             &
     &                 ICE (ng) % ti(:,:,linew(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idTimid)), RST(ng)%Rindex
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out internal ice stress component 11
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idSig11),  &
     &                 RST(ng)%Rindex, gtype,                           &
     &                 LBi, UBi, LBj, UBj, scale,                       &
     &                 GRID(ng) % rmask_io,                             &
     &                 ICE (ng) % sig11(:,:,lienw(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idSig11)), RST(ng)%Rindex
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out internal ice stress component 22
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idSig22),  &
     &                 RST(ng)%Rindex, gtype,                           &
     &                 LBi, UBi, LBj, UBj, scale,                       &
     &                 GRID(ng) % rmask_io,                             &
     &                 ICE (ng) % sig22(:,:,lienw(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idSig22)), RST(ng)%Rindex
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out internal ice stress component 12
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idSig12),  &
     &                 RST(ng)%Rindex, gtype,                           &
     &                 LBi, UBi, LBj, UBj, scale,                       &
     &                 GRID(ng) % rmask_io,                             &
     &                 ICE (ng) % sig12(:,:,lienw(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idSig12)), RST(ng)%Rindex
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice-ocean friction velocity
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTauiw),  &
     &                 RST(ng)%Rindex, gtype,                           &
     &                 LBi, UBi, LBj, UBj, scale,                       &
     &                 GRID(ng) % rmask_io,                             &
     &                 ICE (ng) % utau_iw)
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idTauiw)), RST(ng)%Rindex
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice-ocean momentum transfer coefficient
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idChuiw),  &
     &                 RST(ng)%Rindex, gtype,                           &
     &                 LBi, UBi, LBj, UBj, scale,                       &
     &                 GRID(ng) % rmask_io,                             &
     &                 ICE (ng) % chu_iw)
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idChuiw)), RST(ng)%Rindex
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out salinity in molecular sub-layer under ice
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idS0mk),   &
     &                 RST(ng)%Rindex, gtype,                           &
     &                 LBi, UBi, LBj, UBj, scale,                       &
     &                 GRID(ng) % rmask_io,                             &
     &                 ICE (ng) % s0mk)
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idS0mk)), RST(ng)%Rindex
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out temperature in molecular sub-layer under ice
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idT0mk),   &
     &                 RST(ng)%Rindex, gtype,                           &
     &                 LBi, UBi, LBj, UBj, scale,                       &
     &                 GRID(ng) % rmask_io,                             &
     &                 ICE (ng) % t0mk)
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idT0mk)), RST(ng)%Rindex
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize restart NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, RST(ng)%name, RST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
      IF (Master) WRITE (stdout,20) kstp(ng), nrhs(ng), RST(ng)%Rindex
!
  10  FORMAT (/,' WRT_RST - error while writing variable: ',a,/,11x,    &
     &        'into restart NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_RST   - wrote re-start fields (Index=', i1,       &
     &        ',',i1,') into time record = ',i7.7)
      RETURN
      END SUBROUTINE wrt_rst
