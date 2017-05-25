      SUBROUTINE def_station (ng, ldef)
!
!svn $Id: def_station.F 1484 2012-06-13 19:25:03Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates station data NetCDF file, it defines its       !
!  dimensions, attributes, and variables.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE def_var_mod, ONLY : def_var
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
      logical, intent(in) :: ldef
!
!  Local variable declarations.
!
      integer, parameter :: Natt = 25
      logical :: got_var(NV)
      integer :: i, j, recdim, stadim
      integer :: status
      integer :: DimIDs(32), pgrd(2)
      integer :: Vsize(4)
      integer :: def_dim
      integer :: itrc
      integer :: bgrd(3), rgrd(3), wgrd(3)
      real(r8) :: Aval(6)
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      SourceFile='def_station.F'
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=STA(ng)%name
!
      IF (Master) THEN
        IF (ldef) THEN
          WRITE (stdout,10) TRIM(ncname)
        ELSE
          WRITE (stdout,20) TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new station data file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, iNLM, TRIM(ncname), STA(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          IF (Master) WRITE (stdout,30) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define file dimensions.
!-----------------------------------------------------------------------
!
        DimIDs=0
!
        status=def_dim(ng, iNLM, STA(ng)%ncid, ncname, 's_rho',         &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, STA(ng)%ncid, ncname, 's_w',           &
     &                 N(ng)+1, DimIDs(10))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, STA(ng)%ncid, ncname, 'tracer',        &
     &                 NT(ng), DimIDs(11))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, STA(ng)%ncid, ncname, 'station' ,      &
     &                 Nstation(ng), DimIDs(13))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, STA(ng)%ncid, ncname, 'boundary',      &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, STA(ng)%ncid, ncname,                  &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (exit_flag.ne.NoError) RETURN
        recdim=DimIDs(12)
        stadim=DimIDs(13)
!
!  Define dimension vector for point variables.
!
        pgrd(1)=DimIDs(13)
        pgrd(2)=DimIDs(12)
!
!  Define dimension vector for cast variables at vertical RHO-points.
!
        rgrd(1)=DimIDs( 9)
        rgrd(2)=DimIDs(13)
        rgrd(3)=DimIDs(12)
!
!  Define dimension vector for cast variables at vertical W-points.
!
        wgrd(1)=DimIDs(10)
        wgrd(2)=DimIDs(13)
        wgrd(3)=DimIDs(12)
!
!  Define dimension vector for sediment bed layer type variables.
!
        bgrd(1)=DimIDs(16)
        bgrd(2)=DimIDs(13)
        bgrd(3)=DimIDs(12)
!
!  Initialize unlimited time record dimension.
!
        STA(ng)%Rindex=0
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, iNLM, STA(ng)%ncid, ncname, DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        IF (INT(time_ref).eq.-2) THEN
          Vinfo( 3)='seconds since 1968-05-23 00:00:00 GMT'
          Vinfo( 4)='gregorian'
        ELSE IF (INT(time_ref).eq.-1) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='360_day'
        ELSE IF (INT(time_ref).eq.0) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='julian'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          WRITE (Vinfo( 3),'(a,1x,a)') 'seconds since', TRIM(r_text)
          Vinfo( 4)='gregorian'
        END IF
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idtime),     &
     &                 NF_TYPE, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define free-surface.
!
        IF (Sout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          Vinfo( 2)=Vname(2,idFsur)
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idFsur),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D momentum in the XI-direction.
!
        IF (Sout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          Vinfo( 2)=Vname(2,idUbar)
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idUbar),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D momentum in the ETA-direction.
!
        IF (Sout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          Vinfo( 2)=Vname(2,idVbar)
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idVbar),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D Eastward momentum component at RHO-points.
!
        IF (Sout(idu2dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu2dE)
          Vinfo( 2)=Vname(2,idu2dE)
          Vinfo( 3)=Vname(3,idu2dE)
          Vinfo(14)=Vname(4,idu2dE)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idu2dE),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D Northward momentum component at RHO-points.
!
        IF (Sout(idv2dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv2dN)
          Vinfo( 2)=Vname(2,idv2dN)
          Vinfo( 3)=Vname(3,idv2dN)
          Vinfo(14)=Vname(4,idv2dN)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idv2dN),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D momentum component in the XI-direction.
!
        IF (Sout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          Vinfo( 2)=Vname(2,idUvel)
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idUvel),   &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D momentum component in the ETA-direction.
!
        IF (Sout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          Vinfo( 2)=Vname(2,idVvel)
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idVvel),   &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D Eastward momentum component at RHO-points.
!
        IF (Sout(idu3dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu3dE)
          Vinfo( 2)=Vname(2,idu3dE)
          Vinfo( 3)=Vname(3,idu3dE)
          Vinfo(14)=Vname(4,idu3dE)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idu3dE),   &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D Northward momentum component at RHO-points.
!
        IF (Sout(idv3dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv3dN)
          Vinfo( 2)=Vname(2,idv3dN)
          Vinfo( 3)=Vname(3,idv3dN)
          Vinfo(14)=Vname(4,idv3dN)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idv3dN),   &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D momentum component in the S-direction.
!
        IF (Sout(idWvel,ng)) THEN
          Vinfo( 1)=Vname(1,idWvel)
          Vinfo( 2)=Vname(2,idWvel)
          Vinfo( 3)=Vname(3,idWvel)
          Vinfo(14)=Vname(4,idWvel)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idWvel),   &
     &                   NF_FOUT, 3, wgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define S-coordinate vertical "omega" momentum component (m3/s).
!
        IF (Sout(idOvel,ng)) THEN
          Vinfo( 1)=Vname(1,idOvel)
          Vinfo( 2)=Vname(2,idOvel)
          Vinfo( 3)='meter3 second-1'
          Vinfo(14)=Vname(4,idOvel)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idOvel),   &
     &                   NF_FOUT, 3, wgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Sout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            Vinfo( 2)=Vname(2,idTvar(itrc))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
            status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Tid(itrc),   &
     &                     NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,       &
     &                     SetFillVal = .TRUE.,                         &
     &                     SetParAccess = .FALSE.)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
!
!  Define density anomaly.
!
        IF (Sout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          Vinfo( 2)=Vname(2,idDano)
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idDano),   &
     &                   NF_FOUT, 3, rgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define depth of surface boundary layer.
!
        IF (Sout(idHsbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHsbl)
          Vinfo( 2)=Vname(2,idHsbl)
          Vinfo( 3)=Vname(3,idHsbl)
          Vinfo(14)=Vname(4,idHsbl)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idHsbl),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define vertical viscosity coefficient.
!
        IF (Sout(idVvis,ng)) THEN
          Vinfo( 1)=Vname(1,idVvis)
          Vinfo( 2)=Vname(2,idVvis)
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idVvis),   &
     &                   NF_FOUT, 3, wgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Sout(idTdif,ng)) THEN
          Vinfo( 1)=Vname(1,idTdif)
          Vinfo( 2)=Vname(2,idTdif)
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idTdif),   &
     &                   NF_FOUT, 3, wgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Sout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          Vinfo( 2)=Vname(2,idSdif)
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idSdif),   &
     &                   NF_FOUT, 3, wgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface air pressure.
!
        IF (Sout(idPair,ng)) THEN
          Vinfo( 1)=Vname(1,idPair)
          Vinfo( 2)=Vname(2,idPair)
          Vinfo( 3)=Vname(3,idPair)
          Vinfo(14)=Vname(4,idPair)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idPair),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface winds.
!
        IF (Sout(idUair,ng)) THEN
          Vinfo( 1)=Vname(1,idUair)
          Vinfo( 2)=Vname(2,idUair)
          Vinfo( 3)=Vname(3,idUair)
          Vinfo(14)=Vname(4,idUair)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idUair),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (Sout(idVair,ng)) THEN
          Vinfo( 1)=Vname(1,idVair)
          Vinfo( 2)=Vname(2,idVair)
          Vinfo( 3)=Vname(3,idVair)
          Vinfo(14)=Vname(4,idVair)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idVair),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface East winds.
!
        IF (Sout(idUairE,ng)) THEN
          Vinfo( 1)=Vname(1,idUairE)
          Vinfo( 2)=Vname(2,idUairE)
          Vinfo( 3)=Vname(3,idUairE)
          Vinfo(14)=Vname(4,idUairE)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idUairE),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (Sout(idVairN,ng)) THEN
          Vinfo( 1)=Vname(1,idVairN)
          Vinfo( 2)=Vname(2,idVairN)
          Vinfo( 3)=Vname(3,idVairN)
          Vinfo(14)=Vname(4,idVairN)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idVairN),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface net heat flux.
!
        IF (Sout(idTsur(itemp),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(itemp))
          Vinfo( 2)=Vname(2,idTsur(itemp))
          Vinfo( 3)=Vname(3,idTsur(itemp))
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idTsur(itemp))
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid,                        &
     &                   STA(ng)%Vid(idTsur(itemp)), NF_FOUT,           &
     &                   2, pgrd, Aval, Vinfo, ncname,                  &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface net salt flux.
!
        IF (Sout(idTsur(isalt),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(isalt))
          Vinfo( 2)=Vname(2,idTsur(isalt))
          Vinfo( 3)=Vname(3,idTsur(isalt))
          Vinfo(11)='upward flux, freshening (net precipitation)'
          Vinfo(12)='downward flux, salting (net evaporation)'
          Vinfo(14)=Vname(4,idTsur(isalt))
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid,                        &
     &                   STA(ng)%Vid(idTsur(isalt)), NF_FOUT,           &
     &                   2, pgrd, Aval, Vinfo, ncname,                  &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define latent heat flux.
!
        IF (Sout(idLhea,ng)) THEN
          Vinfo( 1)=Vname(1,idLhea)
          Vinfo( 2)=Vname(2,idLhea)
          Vinfo( 3)=Vname(3,idLhea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLhea)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idLhea),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define sensible heat flux.
!
        IF (Sout(idShea,ng)) THEN
          Vinfo( 1)=Vname(1,idShea)
          Vinfo( 2)=Vname(2,idShea)
          Vinfo( 3)=Vname(3,idShea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idShea)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idShea),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define longwave radiation flux.
!
        IF (Sout(idLrad,ng)) THEN
          Vinfo( 1)=Vname(1,idLrad)
          Vinfo( 2)=Vname(2,idLrad)
          Vinfo( 3)=Vname(3,idLrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLrad)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idLrad),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define shortwave radiation flux.
!
        IF (Sout(idSrad,ng)) THEN
          Vinfo( 1)=Vname(1,idSrad)
          Vinfo( 2)=Vname(2,idSrad)
          Vinfo( 3)=Vname(3,idSrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idSrad)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idSrad),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define E-P flux (as computed by bulk_flux.F).
!
        IF (Sout(idEmPf,ng)) THEN
          Vinfo( 1)=Vname(1,idEmPf)
          Vinfo( 2)=Vname(2,idEmPf)
          Vinfo( 3)=Vname(3,idEmPf)
          Vinfo(11)='upward flux, freshening (net precipitation)'
          Vinfo(12)='downward flux, salting (net evaporation)'
          Vinfo(14)=Vname(4,idEmPf)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idEmPf),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define evaporation rate.
!
        IF (Sout(idevap,ng)) THEN
          Vinfo( 1)=Vname(1,idevap)
          Vinfo( 2)=Vname(2,idevap)
          Vinfo( 3)=Vname(3,idevap)
          Vinfo(11)='downward flux, freshening (condensation)'
          Vinfo(12)='upward flux, salting (evaporation)'
          Vinfo(14)=Vname(4,idevap)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idevap),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define precipitation rate.
!
        IF (Sout(idrain,ng)) THEN
          Vinfo( 1)=Vname(1,idrain)
          Vinfo( 2)=Vname(2,idrain)
          Vinfo( 3)=Vname(3,idrain)
          Vinfo(11)='upward flux, salting (NOT POSSIBLE)'
          Vinfo(12)='downward flux, freshening (precipitation)'
          Vinfo(14)=Vname(4,idrain)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idrain),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface U-momentum stress.
!
        IF (Sout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          Vinfo( 2)=Vname(2,idUsms)
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idUsms),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface V-momentum stress.
!
        IF (Sout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idVsms),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom U-momentum stress.
!
        IF (Sout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          Vinfo( 2)=Vname(2,idUbms)
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idUbms),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom V-momentum stress.
!
        IF (Sout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          Vinfo( 2)=Vname(2,idVbms)
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idVbms),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D ice momentum in the XI-direction.
!
        IF (Sout(idUice,ng)) THEN
          Vinfo( 1)=Vname(1,idUice)
          Vinfo( 2)=Vname(2,idUice)
          Vinfo( 3)=Vname(3,idUice)
          Vinfo(14)=Vname(4,idUice)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idUice),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D ice momentum in the ETA-direction.
!
        IF (Sout(idVice,ng)) THEN
          Vinfo( 1)=Vname(1,idVice)
          Vinfo( 2)=Vname(2,idVice)
          Vinfo( 3)=Vname(3,idVice)
          Vinfo(14)=Vname(4,idVice)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idVice),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D ice momentum in the East-direction.
!
        IF (Sout(idUiceE,ng)) THEN
          Vinfo( 1)=Vname(1,idUiceE)
          Vinfo( 2)=Vname(2,idUiceE)
          Vinfo( 3)=Vname(3,idUiceE)
          Vinfo(14)=Vname(4,idUiceE)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idUiceE),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D ice momentum in the North-direction.
!
        IF (Sout(idViceN,ng)) THEN
          Vinfo( 1)=Vname(1,idViceN)
          Vinfo( 2)=Vname(2,idViceN)
          Vinfo( 3)=Vname(3,idViceN)
          Vinfo(14)=Vname(4,idViceN)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idViceN),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice concentration.
!
        IF (Sout(idAice,ng)) THEN
          Vinfo( 1)=Vname(1,idAice)
          Vinfo( 2)=Vname(2,idAice)
          Vinfo( 3)=Vname(3,idAice)
          Vinfo(14)=Vname(4,idAice)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idAice),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice average thickness.
!
        IF (Sout(idHice,ng)) THEN
          Vinfo( 1)=Vname(1,idHice)
          Vinfo( 2)=Vname(2,idHice)
          Vinfo( 3)=Vname(3,idHice)
          Vinfo(14)=Vname(4,idHice)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idHice),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice/snow surface temperature.
!
        IF (Sout(idTice,ng)) THEN
          Vinfo( 1)=Vname(1,idTice)
          Vinfo( 2)=Vname(2,idTice)
          Vinfo( 3)=Vname(3,idTice)
          Vinfo(14)=Vname(4,idTice)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idTice),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define snow thickness.
!
        IF (Sout(idHsno,ng)) THEN
          Vinfo( 1)=Vname(1,idHsno)
          Vinfo( 2)=Vname(2,idHsno)
          Vinfo( 3)=Vname(3,idHsno)
          Vinfo(14)=Vname(4,idHsno)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idHsno),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface water thickness (on ice).
!
        IF (Sout(idSfwat,ng)) THEN
          Vinfo( 1)=Vname(1,idSfwat)
          Vinfo( 2)=Vname(2,idSfwat)
          Vinfo( 3)=Vname(3,idSfwat)
          Vinfo(14)=Vname(4,idSfwat)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idSfwat),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice age
!
        IF (Sout(idAgeice,ng)) THEN
          Vinfo( 1)=Vname(1,idAgeice)
          Vinfo( 2)=Vname(2,idAgeice)
          Vinfo( 3)=Vname(3,idAgeice)
          Vinfo(14)=Vname(4,idAgeice)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idAgeice), &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice-ocean mass flux
!
        IF (Sout(idIomflx,ng)) THEN
          Vinfo( 1)=Vname(1,idIomflx)
          Vinfo( 2)=Vname(2,idIomflx)
          Vinfo( 3)=Vname(3,idIomflx)
          Vinfo(14)=Vname(4,idIomflx)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idIomflx), &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice internal temperature.
!
        IF (Sout(idTimid,ng)) THEN
          Vinfo( 1)=Vname(1,idTimid)
          Vinfo( 2)=Vname(2,idTimid)
          Vinfo( 3)=Vname(3,idTimid)
          Vinfo(14)=Vname(4,idTimid)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idTimid),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice-water friction velocity.
!
        IF (Sout(idTauiw,ng)) THEN
          Vinfo( 1)=Vname(1,idTauiw)
          Vinfo( 2)=Vname(2,idTauiw)
          Vinfo( 3)=Vname(3,idTauiw)
          Vinfo(14)=Vname(4,idTauiw)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idTauiw),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice-water momentum transfer coefficient.
!
        IF (Sout(idChuiw,ng)) THEN
          Vinfo( 1)=Vname(1,idChuiw)
          Vinfo( 2)=Vname(2,idChuiw)
          Vinfo( 3)=Vname(3,idChuiw)
          Vinfo(14)=Vname(4,idChuiw)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idChuiw),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define under-ice temperature
!
        IF (Sout(idT0mk,ng)) THEN
          Vinfo( 1)=Vname(1,idT0mk)
          Vinfo( 2)=Vname(2,idT0mk)
          Vinfo( 3)=Vname(3,idT0mk)
          Vinfo(14)=Vname(4,idT0mk)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idT0mk),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define under-ice salinity
!
        IF (Sout(idS0mk,ng)) THEN
          Vinfo( 1)=Vname(1,idS0mk)
          Vinfo( 2)=Vname(2,idS0mk)
          Vinfo( 3)=Vname(3,idS0mk)
          Vinfo(14)=Vname(4,idS0mk)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idS0mk),   &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define frazil ice growth rate
!
        IF (Sout(idWfr,ng)) THEN
          Vinfo( 1)=Vname(1,idWfr)
          Vinfo( 2)=Vname(2,idWfr)
          Vinfo( 3)=Vname(3,idWfr)
          Vinfo(14)=Vname(4,idWfr)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idWfr),    &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice growth/melt rate
!
        IF (Sout(idWai,ng)) THEN
          Vinfo( 1)=Vname(1,idWai)
          Vinfo( 2)=Vname(2,idWai)
          Vinfo( 3)=Vname(3,idWai)
          Vinfo(14)=Vname(4,idWai)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idWai),    &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!
!  Define ice growth/melt rate
!
        IF (Sout(idWao,ng)) THEN
          Vinfo( 1)=Vname(1,idWao)
          Vinfo( 2)=Vname(2,idWao)
          Vinfo( 3)=Vname(3,idWao)
          Vinfo(14)=Vname(4,idWao)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idWao),    &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice growth/melt rate
!
        IF (Sout(idWio,ng)) THEN
          Vinfo( 1)=Vname(1,idWio)
          Vinfo( 2)=Vname(2,idWio)
          Vinfo( 3)=Vname(3,idWio)
          Vinfo(14)=Vname(4,idWio)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idWio),    &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice melt runoff rate
!
        IF (Sout(idWro,ng)) THEN
          Vinfo( 1)=Vname(1,idWro)
          Vinfo( 2)=Vname(2,idWro)
          Vinfo( 3)=Vname(3,idWro)
          Vinfo(14)=Vname(4,idWro)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idWro),    &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define internal ice stress sig11
!
        IF (Sout(idSig11,ng)) THEN
          Vinfo( 1)=Vname(1,idSig11)
          Vinfo( 2)=Vname(2,idSig11)
          Vinfo( 3)=Vname(3,idSig11)
          Vinfo(14)=Vname(4,idSig11)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idSig11),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define internal ice stress sig12
!
        IF (Sout(idSig12,ng)) THEN
          Vinfo( 1)=Vname(1,idSig12)
          Vinfo( 2)=Vname(2,idSig12)
          Vinfo( 3)=Vname(3,idSig12)
          Vinfo(14)=Vname(4,idSig12)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idSig12),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define internal ice stress sig22
!
        IF (Sout(idSig22,ng)) THEN
          Vinfo( 1)=Vname(1,idSig22)
          Vinfo( 2)=Vname(2,idSig22)
          Vinfo( 3)=Vname(3,idSig22)
          Vinfo(14)=Vname(4,idSig22)
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ng, iNLM, STA(ng)%ncid, STA(ng)%Vid(idSig22),  &
     &                   NF_FOUT, 2, pgrd, Aval, Vinfo, ncname,         &
     &                   SetFillVal = .TRUE.,                           &
     &                   SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iNLM, ncname, STA(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iNLM, STA(ng)%ncid, ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing stations file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=STA(ng)%name
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open stations file for read/write.
!
        CALL netcdf_open (ng, iNLM, ncname, 1, STA(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,50) TRIM(ncname)
          RETURN
        END IF
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  stations variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            STA(ng)%Vid(idtime)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            STA(ng)%Vid(idFsur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            STA(ng)%Vid(idUbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            STA(ng)%Vid(idVbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dE))) THEN
            got_var(idu2dE)=.TRUE.
            STA(ng)%Vid(idu2dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dN))) THEN
            got_var(idv2dN)=.TRUE.
            STA(ng)%Vid(idv2dN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            STA(ng)%Vid(idUvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            STA(ng)%Vid(idVvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dE))) THEN
            got_var(idu3dE)=.TRUE.
            STA(ng)%Vid(idu3dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dN))) THEN
            got_var(idv3dN)=.TRUE.
            STA(ng)%Vid(idv3dN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            STA(ng)%Vid(idWvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            STA(ng)%Vid(idOvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            STA(ng)%Vid(idDano)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            STA(ng)%Vid(idHsbl)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            STA(ng)%Vid(idVvis)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            STA(ng)%Vid(idTdif)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            STA(ng)%Vid(idSdif)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPair))) THEN
            got_var(idPair)=.TRUE.
            STA(ng)%Vid(idPair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUair))) THEN
            got_var(idUair)=.TRUE.
            STA(ng)%Vid(idUair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVair))) THEN
            got_var(idVair)=.TRUE.
            STA(ng)%Vid(idVair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUairE))) THEN
            got_var(idUairE)=.TRUE.
            STA(ng)%Vid(idUairE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVairN))) THEN
            got_var(idVairN)=.TRUE.
            STA(ng)%Vid(idVairN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(itemp)))) THEN
            got_var(idTsur(itemp))=.TRUE.
            STA(ng)%Vid(idTsur(itemp))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(isalt)))) THEN
            got_var(idTsur(isalt))=.TRUE.
            STA(ng)%Vid(idTsur(isalt))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            STA(ng)%Vid(idLhea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            STA(ng)%Vid(idShea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            STA(ng)%Vid(idLrad)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            STA(ng)%Vid(idSrad)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idEmPf))) THEN
            got_var(idEmPf)=.TRUE.
            STA(ng)%Vid(idEmPf)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            STA(ng)%Vid(idevap)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            STA(ng)%Vid(idrain)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            STA(ng)%Vid(idUsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            STA(ng)%Vid(idVsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            STA(ng)%Vid(idUbms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            STA(ng)%Vid(idVbms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUice))) THEN
            got_var(idUice)=.TRUE.
            STA(ng)%Vid(idUice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVice))) THEN
            got_var(idVice)=.TRUE.
            STA(ng)%Vid(idVice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUiceE))) THEN
            got_var(idUiceE)=.TRUE.
            STA(ng)%Vid(idUiceE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idViceN))) THEN
            got_var(idViceN)=.TRUE.
            STA(ng)%Vid(idViceN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAice))) THEN
            got_var(idAice)=.TRUE.
            STA(ng)%Vid(idAice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHice))) THEN
            got_var(idHice)=.TRUE.
            STA(ng)%Vid(idHice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTice))) THEN
            got_var(idTice)=.TRUE.
            STA(ng)%Vid(idTice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsno))) THEN
            got_var(idHsno)=.TRUE.
            STA(ng)%Vid(idHsno)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSfwat))) THEN
            got_var(idSfwat)=.TRUE.
            STA(ng)%Vid(idSfwat)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAgeice))) THEN
            got_var(idAgeice)=.TRUE.
            STA(ng)%Vid(idAgeice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIomflx))) THEN
            got_var(idIomflx)=.TRUE.
            STA(ng)%Vid(idIomflx)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTimid))) THEN
            got_var(idTimid)=.TRUE.
            STA(ng)%Vid(idTimid)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTauiw))) THEN
            got_var(idTauiw)=.TRUE.
            STA(ng)%Vid(idTauiw)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idChuiw))) THEN
            got_var(idChuiw)=.TRUE.
            STA(ng)%Vid(idChuiw)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idT0mk))) THEN
            got_var(idT0mk)=.TRUE.
            STA(ng)%Vid(idT0mk)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idS0mk))) THEN
            got_var(idS0mk)=.TRUE.
            STA(ng)%Vid(idS0mk)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWfr))) THEN
            got_var(idWfr)=.TRUE.
            STA(ng)%Vid(idWfr)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWai))) THEN
            got_var(idWai)=.TRUE.
            STA(ng)%Vid(idWai)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWao))) THEN
            got_var(idWao)=.TRUE.
            STA(ng)%Vid(idWao)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWio))) THEN
            got_var(idWio)=.TRUE.
            STA(ng)%Vid(idWio)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWro))) THEN
            got_var(idWro)=.TRUE.
            STA(ng)%Vid(idWro)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig11))) THEN
            got_var(idSig11)=.TRUE.
            STA(ng)%Vid(idSig11)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig12))) THEN
            got_var(idSig12)=.TRUE.
            STA(ng)%Vid(idSig12)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig22))) THEN
            got_var(idSig22)=.TRUE.
            STA(ng)%Vid(idSig22)=var_id(i)
          END IF
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              STA(ng)%Tid(itrc)=var_id(i)
            END IF
          END DO
        END DO
!
!  Check if station variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur).and.Sout(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Sout(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Sout(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu2dE).and.Sout(idu2dE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu2dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dN).and.Sout(idv2dN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv2dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUvel).and.Sout(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Sout(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu3dE).and.Sout(idu3dE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu3dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dN).and.Sout(idv3dN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv3dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Sout(idWvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvel).and.Sout(idOvel,ng)) THEN
          IF (Master) WRITE(stdout,60) TRIM(Vname(1,idOvel)),           &
     &                                 TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano).and.Sout(idDano,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHsbl).and.Sout(idHsbl,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHsbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvis).and.Sout(idVvis,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif).and.Sout(idTdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSdif).and.Sout(idSdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPair).and.Sout(idPair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idPair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUair).and.Sout(idUair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVair).and.Sout(idVair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUairE).and.Sout(idUairE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUairE)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVairN).and.Sout(idVairN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVairN)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTsur(itemp)).and.Sout(idTsur(itemp),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(itemp))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTsur(isalt)).and.Sout(idTsur(isalt),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(isalt))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLhea).and.Sout(idLhea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLhea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Sout(idShea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idShea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Sout(idLrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSrad).and.Sout(idSrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idEmPf).and.Sout(idEmPf,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idEmPf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idevap).and.Sout(idevap,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idevap)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Sout(idrain,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idrain)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUsms).and.Sout(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Sout(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Sout(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Sout(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Sout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        IF (.not.got_var(idUice).and.Sout(idUice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVice).and.Sout(idVice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUiceE).and.Sout(idUiceE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUiceE)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idViceN).and.Sout(idViceN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idViceN)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idAice).and.Sout(idAice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idAice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHice).and.Sout(idHice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTice).and.Sout(idTice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHsno).and.Sout(idHsno,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHsno)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSfwat).and.Sout(idSfwat,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSfwat)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idAgeice).and.Sout(idAgeice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idAgeice)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idIomflx).and.Sout(idIomflx,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idIomflx)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTimid).and.Sout(idTimid,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTimid)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTauiw).and.Sout(idTauiw,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTauiw)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idChuiw).and.Sout(idChuiw,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idChuiw)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idT0mk).and.Sout(idT0mk,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idT0mk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idS0mk).and.Sout(idS0mk,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idS0mk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWfr).and.Sout(idWfr,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWfr)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWai).and.Sout(idWai,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWai)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWao).and.Sout(idWao,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWao)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWio).and.Sout(idWio,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWio)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWro).and.Sout(idWro,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWro)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig11).and.Sout(idSig11,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSig11)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig12).and.Sout(idSig12,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSig12)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig22).and.Sout(idSig22,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSig22)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
!
!  Set unlimited time record dimension to the appropriate value.
!
        STA(ng)%Rindex=(ntstart(ng)-1)/nSTA(ng)
      END IF QUERY
!
  10  FORMAT (6x,'DEF_STATION - creating stations file: ',a)
  20  FORMAT (6x,'DEF_STATION - inquiring stations file: ',a)
  30  FORMAT (/,' DEF_STATION - unable to create stations NetCDF ',     &
     &        'file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_STATION - unable to open stations NetCDF file: ', &
     &        a)
  60  FORMAT (/,' DEF_STATION - unable to find variable: ',a,2x,        &
     &        ' in stations NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_station
