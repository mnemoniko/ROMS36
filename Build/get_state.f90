      SUBROUTINE get_state (ng, model, msg, ncname, IniRec, Tindex)
!
!svn $Id: get_state.F 1460 2012-03-23 23:12:41Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in requested model state from specified NetCDF   !
!  file. It is usually used to read initial conditions.                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     msg        Message index for StateMsg.                           !
!     ncname     Nonlinear initial conditions NetCDF file name.        !
!     IniRec     Nonlinear initial conditions time record to read.     !
!     Tindex     State variable time index to intialize.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_ice
      USE mod_stepping
      USE mod_strings
      USE mp_exchange_mod, ONLY : mp_exchange2d
      USE mp_exchange_mod, ONLY : mp_exchange3d
      USE nf_fread2d_mod, ONLY : nf_fread2d
      USE nf_fread3d_mod, ONLY : nf_fread3d
      USE nf_fread4d_mod, ONLY : nf_fread4d
      USE strings_mod, ONLY : find_string
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, msg, Tindex
      integer, intent(inout) :: IniRec
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      logical :: Perfect2D, Perfect3D, foundit
      logical, dimension(NV) :: get_var, have_var
      integer :: LBi, UBi, LBj, UBj
      integer :: IDmod, InpRec, gtype, i, ifield, itrc, lstr, lend
      integer :: Nrec, mySize, ncINPid, nvatts, nvdim, status, varid
      integer :: Vsize(4), start(4), total(4)
      real(r8), parameter :: Fscl = 1.0_r8
      real(r8) :: Fmax, Fmin, INPtime, Tmax, scale, time_scale
      real(r8), allocatable :: myTime(:)
      character (len=6 ) :: string
      character (len=14) :: t_code
      character (len=15) :: attnam, tvarnam
      character (len=40) :: tunits
!
      SourceFile='get_state.F'
!
!-----------------------------------------------------------------------
!  Determine variables to read and their availability.
!-----------------------------------------------------------------------
!
!  Set model identification string.
!
      IF (model.eq.iNLM.or.(model.eq.0)) THEN
        string=' NLM: '                    ! nonlinear model, restart
        IDmod=iNLM
      ELSE IF (model.eq.iTLM) THEN
        string=' TLM: '                    ! tangent linear model
        IDmod=iTLM
      ELSE IF (model.eq.iRPM) THEN
        string=' RPM: '                    ! representer model
        IDmod=iRPM
      ELSE IF (model.eq.iADM) THEN
        string=' ADM: '                    ! adjoint model
        IDmod=iADM
      ELSE IF (model.eq.5) THEN
        string=' NRM: '                    ! normalization factor
        IDmod=iNLM                         ! model or initial conditions
      ELSE IF (model.eq.6) THEN
        string=' STD: '                    ! standard deviation
        IDmod=iNLM                         ! model or initial conditions
      ELSE IF (model.eq.7) THEN
        string=' FRC: '                    ! impulse forcing
        IDmod=iNLM
      ELSE IF (model.eq.8) THEN
        string=' STD: '                    ! standard deviation
        IDmod=iNLM                         ! boundary conditions
      ELSE IF (model.eq.9) THEN
        string=' STD: '                    ! standard deviation
        IDmod=iNLM                         ! surface forcing
      ELSE IF (model.eq.10) THEN
        string=' NRM: '                    ! normalization factor
        IDmod=iNLM                         ! boundary conditions
      ELSE IF (model.eq.11) THEN
        string=' NRM: '                    ! normalization factor
        IDmod=iNLM                         ! surface forcing
      ELSE IF (model.eq.12) THEN
        string=' NLM: '                    ! tangent linear forcing and
        IDmod=iNLM                         ! obc increments
      END IF
!
!  Set switch to process variables for nonlinear model perfect restart.
!
      Perfect2D=.FALSE.
      Perfect3D=.FALSE.
      PerfectRST(ng)=Perfect2d.or.Perfect3d
!
!  Determine variables to read.
!
      CALL checkvars (ng, model, ncname, string, Nrec, NV, tvarnam,     &
     &                get_var, have_var)
      IF (exit_flag.ne.NoError) RETURN
!
!  Set Vsize to zero to deactivate interpolation of input data to model
!  grid in "nf_fread2d" and "nf_fread3d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!-----------------------------------------------------------------------
!  Open input NetCDF file and check time variable.
!-----------------------------------------------------------------------
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!  Open input NetCDF file.
!
      CALL netcdf_open (ng, IDmod, ncname, 0, ncINPid)
      IF (exit_flag.ne.NoError) THEN
        IF (Master) WRITE (stdout,20) TRIM(ncname)
        RETURN
      END IF
!
!  If restart, read in lateral boundary conditions global attribute
!  from restart file and check keyword strings with structure vlues
!  for consistency.
!
      IF (((model.eq.0).or.(model.eq.iNLM)).and.(nrrec(ng).ne.0)) THEN
        CALL lbc_getatt (ng, model, ncINPid, ncname, 'NLM_LBC', LBC)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Inquire about the input time variable.
!
      CALL netcdf_inq_var (ng, IDmod, ncname,                           &
     &                     MyVarName = TRIM(tvarnam),                   &
     &                     VarID = varid,                               &
     &                     nVarDim =  nvdim,                            &
     &                     nVarAtt = nvatts)
      IF (exit_flag.ne.NoError) RETURN
!
!  Allocate input time variable and read its value(s).  Recall that
!  input time variable is a one-dimensional array with one or several
!  values.
!
      mySize=var_Dsize(1)
      IF (.not.allocated(myTime)) allocate (myTime(mySize))
      CALL netcdf_get_fvar (ng, IDmod, ncname, TRIM(tvarnam), myTime)
      IF (exit_flag.ne.NoError) RETURN
!
!  If using the latest time record from input NetCDF file as the
!  initialization record, assign input time.
!
      IF (LastRec(ng)) THEN
        Tmax=-1.0_r8
        DO i=1,mySize
          IF (myTime(i).gt.Tmax) THEN
            Tmax=myTime(i)
            IniRec=i
          END IF
        END DO
        INPtime=Tmax
        InpRec=IniRec
      ELSE
        IF ((IniRec.ne.0).and.(IniRec.gt.mySize)) THEN
          IF (Master)  WRITE (stdout,40) string, IniRec, TRIM(ncname),  &
     &                                   mySize
          exit_flag=2
          RETURN
        END IF
        IF (IniRec.ne.0) THEN
          InpRec=IniRec
        ELSE
          InpRec=1
        END IF
        INPtime=myTime(InpRec)
      END IF
      IF (allocated(myTime)) deallocate ( myTime )
!
!  Set input time scale by looking at the "units" attribute.
!
      time_scale=0.0_r8
      DO i=1,nvatts
        IF (TRIM(var_Aname(i)).eq.'units') THEN
          IF (INDEX(TRIM(var_Achar(i)),'day').ne.0) THEN
            time_scale=day2sec
          ELSE IF (INDEX(TRIM(var_Achar(i)),'second').ne.0) THEN
            time_scale=1.0_r8
          END IF
        END IF
      END DO
      IF (time_scale.gt.0.0_r8) THEN
        INPtime=INPtime*time_scale
      END IF
!
!  Set starting time index and time clock in days.  Notice that the
!  global time variables and indices are only over-written when
!  processing initial conditions (msg = 1).
!
      IF ((model.eq.0).or.(model.eq.iNLM).or.                           &
     &    (model.eq.iTLM).or.(model.eq.iRPM)) THEN
        IF (((model.eq.iTLM).or.(model.eq.iRPM)).and.(msg.eq.1).and.    &
     &      (INPtime.ne.(dstart*day2sec))) THEN
          INPtime=dstart*day2sec
        END IF
        IF (msg.eq.1) THEN            ! processing initial conditions
          time(ng)=INPtime
          tdays(ng)=time(ng)*sec2day
          ntstart(ng)=NINT((time(ng)-dstart*day2sec)/dt(ng))+1
          IF (ntstart(ng).lt.1) ntstart(ng)=1
          IF (PerfectRST(ng)) THEN
            ntfirst(ng)=1
          ELSE
            ntfirst(ng)=ntstart(ng)
          END IF
        END IF
      ELSE IF (model.eq.iADM) THEN
        IF ((msg.eq.1).and.(INPtime.eq.0.0_r8)) THEN
          INPtime=time(ng)
        ELSE IF (msg.ne.1) THEN
          time(ng)=INPtime
          tdays(ng)=time(ng)*sec2day
        END IF
        ntstart(ng)=ntimes(ng)+1
        ntend(ng)=1
        ntfirst(ng)=ntend(ng)
      END IF
      CALL time_string (time(ng), time_code(ng))
!
!  Over-write "IniRec" to the actual initial record processed.
!
      IF (model.eq.iNLM) THEN
        IniRec=InpRec
      END IF
!
!  Report information.
!
      lstr=SCAN(ncname,'/',BACK=.TRUE.)+1
      lend=LEN_TRIM(ncname)
      IF (Master) THEN
        CALL time_string (INPtime, t_code)
        IF (ERend.gt.ERstr) THEN
          WRITE (stdout,80) string, TRIM(StateMsg(msg)), t_code, Nrun,  &
     &                      ncname(lstr:lend), InpRec, Tindex
        ELSE
          WRITE (stdout,90) string, TRIM(StateMsg(msg)), t_code,        &
     &                      ncname(lstr:lend), InpRec, Tindex
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Read in nonlinear state variables. If applicable, read in perfect
!  restart variables.
!-----------------------------------------------------------------------
!
      NLM_STATE: IF ((model.eq.iNLM).or.(model.eq.0).or.                &
     &               (model.eq.13)) THEN
!
!  Read in nonlinear free-surface (m).
!
        IF (get_var(idFsur)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idFsur)),   &
     &                        varid)
          IF (Perfect2d) THEN
            gtype=var_flag(varid)*r3dvar
          ELSE
            gtype=var_flag(varid)*r2dvar
          END IF
          IF (Perfect2d) THEN
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idFsur), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 3,                 &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        OCEAN(ng) % zeta)
          ELSE
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idFsur), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        OCEAN(ng) % zeta(:,:,Tindex))
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idFsur)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idFsur)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear RHS of free-surface.
!
        IF (get_var(idRzet).and.Perfect2D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRzet)),   &
     &                        varid)
          gtype=var_flag(varid)*r3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idRzet), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, 2,                   &
     &                      Fscl, Fmin, Fmax,                           &
     &                      GRID(ng) % rmask,                           &
     &                      OCEAN(ng) % rzeta)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idRzet)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idRzet)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear 2D U-momentum component (m/s).
!
        IF (get_var(idUbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUbar)),   &
     &                        varid)
          IF (Perfect2D) THEN
            gtype=var_flag(varid)*u3dvar
          ELSE
            gtype=var_flag(varid)*u2dvar
          END IF
          IF (Perfect2D) THEN
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idUbar), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 3,                 &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % umask,                         &
     &                        OCEAN(ng) % ubar)
          ELSE
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idUbar), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % umask,                         &
     &                        OCEAN(ng) % ubar(:,:,Tindex))
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idUbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idUbar)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear RHS of 2D U-momentum component.
!
        IF (get_var(idRu2d).and.Perfect2D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRu2d)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idRu2d), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, 2,                   &
     &                      Fscl, Fmin, Fmax,                           &
     &                      GRID(ng) % umask,                           &
     &                      OCEAN(ng) % rubar)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idRu2d)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idRu2d)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear 2D U-momentum component (m/s).
!
        IF (get_var(idVbar)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVbar)),   &
     &                        varid)
          IF (Perfect2D) THEN
            gtype=var_flag(varid)*v3dvar
          ELSE
            gtype=var_flag(varid)*v2dvar
          END IF
          IF (Perfect2D) THEN
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVbar), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, 3,                 &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % vmask,                         &
     &                        OCEAN(ng) % vbar)
          ELSE
            status=nf_fread2d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVbar), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % vmask,                         &
     &                        OCEAN(ng) % vbar(:,:,Tindex))
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVbar)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVbar)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear RHS 2D U-momentum component.
!
        IF (get_var(idRv2d).and.Perfect2D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRv2d)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idRv2d), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, 2,                   &
     &                      Fscl, Fmin, Fmax,                           &
     &                      GRID(ng) % vmask,                           &
     &                      OCEAN(ng) % rvbar)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idRv2d)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idRv2d)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear 3D U-momentum component (m/s).
!
        IF (get_var(idUvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idUvel)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          IF (Perfect3D) THEN
            status=nf_fread4d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idUvel), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng), 1, 2,       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % umask,                         &
     &                        OCEAN(ng) % u)
          ELSE
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idUvel), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % umask,                         &
     &                        OCEAN(ng) % u(:,:,:,Tindex))
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idUvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idUvel)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear RHS of 3D U-momentum component.
!
        IF (get_var(idRu3d).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRu3d)),   &
     &                        varid)
          gtype=var_flag(varid)*u3dvar
          status=nf_fread4d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idRu3d), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 0, N(ng), 1, 2,         &
     &                      Fscl, Fmin, Fmax,                           &
     &                      GRID(ng) % umask,                           &
     &                      OCEAN(ng) % ru)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idRu3d)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idRu3d)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear 3D V-momentum component (m/s).
!
        IF (get_var(idVvel)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvel)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          IF (Perfect3D) THEN
            status=nf_fread4d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVvel), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng), 1, 2,       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % vmask,                         &
     &                        OCEAN(ng) % v)
          ELSE
            status=nf_fread3d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idVvel), varid,                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % vmask,                         &
     &                        OCEAN(ng) % v(:,:,:,Tindex))
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVvel)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVvel)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear RHS of 3D V-momentum component.
!
        IF (get_var(idRv3d).and.Perfect3D) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idRv3d)),   &
     &                        varid)
          gtype=var_flag(varid)*v3dvar
          status=nf_fread4d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idRv3d), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 0, N(ng), 1, 2,         &
     &                      Fscl, Fmin, Fmax,                           &
     &                      GRID(ng) % vmask,                           &
     &                      OCEAN(ng) % rv)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idRv3d)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idRv3d)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in nonlinear tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (get_var(idTvar(itrc))) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idTvar(itrc))), varid)
            gtype=var_flag(varid)*r3dvar
            IF (Perfect3D) THEN
              status=nf_fread4d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idTvar(itrc)), varid,           &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 1, N(ng), 1, 2,     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % rmask,                       &
     &                          OCEAN(ng) % t(:,:,:,:,itrc))
            ELSE
              status=nf_fread3d(ng, IDmod, ncname, ncINPid,             &
     &                          Vname(1,idTvar(itrc)), varid,           &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % rmask,                       &
     &                          OCEAN(ng) % t(:,:,:,Tindex,itrc))
            END IF
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,30) string, TRIM(Vname(1,idTvar(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,100) TRIM(Vname(2,idTvar(itrc))),         &
     &                             Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
!
!  Read in vertical viscosity.
!
        IF (have_var(idVvis)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idVvis)),   &
     &                        varid)
          gtype=var_flag(varid)*w3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idVvis), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      Fscl, Fmin,Fmax,                            &
     &                      GRID(ng) % rmask,                           &
     &                      MIXING(ng) % AKv)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idVvis)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idVvis)), Fmin, Fmax
            END IF
          END IF
          CALL mp_exchange3d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        MIXING(ng) % AKv)
        END IF
!
!  Read in temperature vertical diffusion.
!
        IF (have_var(idTdif)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idTdif)),   &
     &                        varid)
          gtype=var_flag(varid)*w3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idTdif), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      Fscl, Fmin,Fmax,                            &
     &                      GRID(ng) % rmask,                           &
     &                      MIXING(ng) % AKt(:,:,:,itemp))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idTdif)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idTdif)), Fmin, Fmax
            END IF
          END IF
          CALL mp_exchange3d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        MIXING(ng) % AKt(:,:,:,itemp))
        END IF
!
!  Read in salinity vertical diffusion.
!
        IF (have_var(idSdif)) THEN
          foundit=find_string(var_name, n_var, TRIM(Vname(1,idSdif)),   &
     &                        varid)
          gtype=var_flag(varid)*w3dvar
          status=nf_fread3d(ng, IDmod, ncname, ncINPid,                 &
     &                      Vname(1,idSdif), varid,                     &
     &                      InpRec, gtype, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      Fscl, Fmin,Fmax,                            &
     &                      GRID(ng) % rmask,                           &
     &                      MIXING(ng) % AKt(:,:,:,isalt))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) string, TRIM(Vname(1,idSdif)), InpRec,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,100) TRIM(Vname(2,idSdif)), Fmin, Fmax
            END IF
          END IF
          CALL mp_exchange3d (ng, MyRank, IDmod, 1,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        MIXING(ng) % AKt(:,:,:,isalt))
        END IF
!
!  Read in 2D ice momentum component (m/s) in the XI-direction.
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idUice)),     &
     &                        varid)
        gtype=var_flag(varid)*u2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idUice), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % umask,                                 &
     &                ICE(ng) % ui(:,:,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) string, TRIM(Vname(1,idUice)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,100) TRIM(Vname(2,idUice)), Fmin, Fmax
          END IF
        END IF
!
!  Read in 2D ice momentum component (m/s) in the ETA-direction.
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idVice)),     &
     &                        varid)
        gtype=var_flag(varid)*v2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idVice), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % vmask,                                 &
     &                ICE(ng) % vi(:,:,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) string, TRIM(Vname(1,idVice)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,100) TRIM(Vname(2,idVice)), Fmin, Fmax
          END IF
        END IF
!
!  Read in ice concentration.
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idAice)),     &
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idAice), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % rmask,                                 &
     &                ICE(ng) % ai(:,:,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) string, TRIM(Vname(1,idAice)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,100) TRIM(Vname(2,idAice)), Fmin, Fmax
          END IF
        END IF
!
!  Read in ice average thickness.
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idHice)),     &
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idHice), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % rmask,                                 &
     &                ICE(ng) % hi(:,:,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) string, TRIM(Vname(1,idHice)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,100) TRIM(Vname(2,idHice)), Fmin, Fmax
          END IF
        END IF
!
!  Read in snow average thickness.
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idHsno)),     &
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idHsno), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % rmask,                                 &
     &                ICE(ng) % hsn(:,:,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) string, TRIM(Vname(1,idHsno)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,100) TRIM(Vname(2,idHsno)), Fmin, Fmax
          END IF
        END IF
!
!  Read in ice/snow surface temperature.
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idTice)),     &
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idTice), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % rmask,                                 &
     &                ICE(ng) % tis)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) string, TRIM(Vname(1,idTice)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,100) TRIM(Vname(2,idTice)), Fmin, Fmax
          END IF
        END IF
!
!  Read in surface water thickness (on ice).
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idSfwat)),    &
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idSfwat), varid, InpRec, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % rmask,                                 &
     &                ICE(ng) % sfwat(:,:,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) string, TRIM(Vname(1,idSfwat)),           &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,100) TRIM(Vname(2,idSfwat)), Fmin, Fmax
          END IF
        END IF
!
!  Read in ice age.
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idAgeice)),   &
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idAgeice), varid, InpRec, gtype,          &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % rmask,                                 &
     &                ICE(ng) % ageice(:,:,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) string, TRIM(Vname(1,idAgeice)),          &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,100) TRIM(Vname(2,idAgeice)), Fmin, Fmax
          END IF
        END IF
!
!  Read in ice internal temperature.
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idTimid)),    &
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idTimid), varid, InpRec, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % rmask,                                 &
     &                ICE(ng) % ti(:,:,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) string, TRIM(Vname(1,idTimid)),           &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,100) TRIM(Vname(2,idTimid)), Fmin, Fmax
          END IF
        END IF
!
!  Read in internal ice stress component 11
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idSig11)),    &
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idSig11), varid, InpRec, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % rmask,                                 &
     &                ICE(ng) % sig11(:,:,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) string, TRIM(Vname(1,idSig11)),           &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,100) TRIM(Vname(2,idSig11)), Fmin, Fmax
          END IF
        END IF
!
!  Read in internal ice stress component 22
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idSig22)),    &
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idSig22), varid, InpRec, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % rmask,                                 &
     &                ICE(ng) % sig22(:,:,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) string, TRIM(Vname(1,idSig22)),           &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,100) TRIM(Vname(2,idSig22)), Fmin, Fmax
          END IF
        END IF
!
!  Read in internal ice stress component 12
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idSig12)),    &
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idSig12), varid, InpRec, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % rmask,                                 &
     &                ICE(ng) % sig12(:,:,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) string, TRIM(Vname(1,idSig12)),           &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,100) TRIM(Vname(2,idSig12)), Fmin, Fmax
          END IF
        END IF
!
!  Read in ice-water friction velocity.
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idTauiw)),    &
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idTauiw), varid, InpRec, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % rmask,                                 &
       &                ICE(ng) % utau_iw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) string, TRIM(Vname(1,idTauiw)),           &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,100) TRIM(Vname(2,idTauiw)), Fmin, Fmax
          END IF
        END IF
!
!  Read in ice-water momentum transfer coefficient.
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idChuiw)),    &
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idChuiw), varid, InpRec, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % rmask,                                 &
     &                ICE(ng) % chu_iw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) string, TRIM(Vname(1,idChuiw)),           &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,100) TRIM(Vname(2,idChuiw)), Fmin, Fmax
          END IF
        END IF
!
!  Read in salinity of molecular sub-layer under ice
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idS0mk)),     &
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idS0mk), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % rmask,                                 &
     &                ICE(ng) % s0mk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) string, TRIM(Vname(1,idS0mk)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,100) TRIM(Vname(2,idS0mk)), Fmin, Fmax
          END IF
        END IF
!
!  Read in temperature of molecular sub-layer under ice
!
        foundit=find_string(var_name, n_var, TRIM(Vname(1,idT0mk)),     &
     &                        varid)
        gtype=var_flag(varid)*r2dvar
        status=nf_fread2d(ng, IDmod, ncname, ncINPid,                   &
     &                Vname(1,idT0mk), varid, InpRec, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % rmask,                                 &
     &                ICE(ng) % t0mk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) string, TRIM(Vname(1,idT0mk)),            &
     &                        InpRec, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,100) TRIM(Vname(2,idT0mk)), Fmin, Fmax
          END IF
        END IF
!
      END IF NLM_STATE
!
!-----------------------------------------------------------------------
!  Close input NetCDF file.
!-----------------------------------------------------------------------
!
      CALL netcdf_close (ng, IDmod, ncINPid, ncname, .FALSE.)
!
  20  FORMAT (/,a,'GET_STATE - unable to open input NetCDF file: ',a)
  30  FORMAT (/,a,'GET_STATE - error while reading variable: ',a,2x,    &
     &        'at time record = ',i3,/,18x,'in input NetCDF file: ',a)
  40  FORMAT (/,a,'GET_STATE - requested input time record = ',i3,/,    &
     &        18x,'not found in input NetCDF: ',a,/,                    &
     &        18x,'number of available records = ',i3)
  50  FORMAT (/,a,'GET_STATE - error while reading attribute: ',a)
  60  FORMAT (/,a,'GET_STATE - error while inquiring attributes',       &
     &          ' for variable: ',a)
  70  FORMAT (/,a,'GET_STATE - cannot inquire about time variable',/,   &
     &        18x,'in input NetCDF file: ',a)
  80  FORMAT (/,a,'GET_STATE - ',a,t62,'t = ',a,                        &
     &        /,19x,'(Iter=',i4.4,', File: ',a,', Rec=',i4.4,           &
     &        ', Index=',i1,')')
  90  FORMAT (/,a,'GET_STATE - ',a,t62,'t = ',a,                        &
     &        /,19x,'(File: ',a,', Rec=',i4.4,', Index=',i1,')')
 100  FORMAT (16x,'- ',a,/,19x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
      RETURN
      END SUBROUTINE get_state
