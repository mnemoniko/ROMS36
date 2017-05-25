      SUBROUTINE checkvars (ng, model, ncname, string, Nrec, Nvar,      &
     &                      tvarnam, get_var, have_var)
!
!svn $Id: checkvars.F 1460 2012-03-23 23:12:41Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine checks if needed state variables are available in      !
!  requested NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     ncname     NetCDF file name.                                     !
!     string     Identification string.                                !
!     Nvar       Size of logical switches arrays.                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Nrec       Number of time records available.                     !
!     tvarnam    Name of time record variable.                         !
!     get_var    Logical switches (T/F), in terms of variable ID,      !
!                  indicating state variables needed by the model.     !
!     have_var   Logical switches (T/F), in terms of variable ID,      !
!                  indicating state variables available in NetCDF      !
!                  file.                                               !
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
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Nvar
      integer, intent(inout) :: Nrec
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: string
      character (len=*), intent(inout) :: tvarnam
      logical, dimension(Nvar), intent(out) :: get_var
      logical, dimension(Nvar), intent(out) :: have_var
!
!  Local variable declarations.
!
      integer :: IDmod, i, itrc
!
      SourceFile='checkvars.F'
!
!-----------------------------------------------------------------------
!  Determine state variables needed and check if they are available in
!  requested NetCDF file.
!-----------------------------------------------------------------------
!
!  Limit model identifier. The profiling is limited to iNLM, iTLM, iRPM,
!  and iADM.
!
      IF ((model.lt.1).or.(model.gt.4)) THEN
        IDmod=iNLM
      ELSE
        IDmod=model
      END IF
!
!  Inquire about the dimensions and check for consistency.
!
      CALL netcdf_check_dim (ng, IDmod, ncname)
      IF (exit_flag.ne.NoError) RETURN
      Nrec=rec_size
!
!  Inquire about the variables.
!
      CALL netcdf_check_var (ng, model, ncname)
      IF (exit_flag.ne.NoError) RETURN
!
!  Initialize logical switches.
!
      DO i=1,Nvar
        get_var(i)=.FALSE.
        have_var(i)=.FALSE.
      END DO
!
!  Determine state variables to read from input NetCDF file.  Notice
!  that these state variable are only assigned if input flag model < 8.
!  Higher values are used to process only boundary condition or
!  surface forcing in the variational data assimilation algorithms.
!
      IF (model.lt.8) THEN
        get_var(idFsur)=.TRUE.
        get_var(idUbar)=.TRUE.
        get_var(idVbar)=.TRUE.
        get_var(idUvel)=.TRUE.
        get_var(idVvel)=.TRUE.
        DO itrc=1,NAT
          get_var(idTvar(itrc))=.TRUE.
        END DO
        DO itrc=1,NPT
          get_var(idTvar(inert(itrc)))=.TRUE.
        END DO
        get_var(idUice)=.true.
        get_var(idVice)=.true.
        get_var(idAice)=.true.
        get_var(idHice)=.true.
        get_var(idTice)=.true.
        get_var(idHsno)=.true.
        get_var(idSfwat)=.true.
        get_var(idAgeice)=.true.
        get_var(idTimid)=.true.
        get_var(idSig11)=.true.
        get_var(idSig22)=.true.
        get_var(idSig12)=.true.
        get_var(idTauiw)=.true.
        get_var(idChuiw)=.true.
        get_var(idS0mk)=.true.
        get_var(idT0mk)=.true.
      END IF
!
!  Scan variable list from input NetCDF and activate switches for
!  model state variables.
!
      DO i=1,n_var
        IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
          tvarnam=TRIM(var_name(i))
          have_var(idtime)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idKver))) THEN
          have_var(idKver)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
          have_var(idUvel)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRu3d))) THEN
          have_var(idRu3d)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
          have_var(idVvel)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRv3d))) THEN
          have_var(idRv3d)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
          have_var(idVvis)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
          have_var(idTdif)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
          have_var(idSdif)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
          have_var(idFsur)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRzet))) THEN
          have_var(idRzet)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
          have_var(idUbar)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRu2d))) THEN
          have_var(idRu2d)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
          have_var(idVbar)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRv2d))) THEN
          have_var(idRv2d)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idKhor))) THEN
          have_var(idKhor)=.TRUE.
        END IF
        DO itrc=1,NT(ng)
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
            have_var(idTvar(itrc))=.TRUE.
          END IF
        END DO
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUice))) THEN
            have_var(idUice)=.true.
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVice))) THEN
            have_var(idVice)=.true.
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAice))) THEN
            have_var(idAice)=.true.
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHice))) THEN
            have_var(idHice)=.true.
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTice))) THEN
            have_var(idTice)=.true.
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsno))) THEN
            have_var(idHsno)=.true.
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSfwat))) THEN
            have_var(idSfwat)=.true.
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAgeice))) THEN
            have_var(idAgeice)=.true.
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTimid))) THEN
            have_var(idTimid)=.true.
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig11))) THEN
            have_var(idSig11)=.true.
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig22))) THEN
            have_var(idSig22)=.true.
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig12))) THEN
            have_var(idSig12)=.true.
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTauiw))) THEN
            have_var(idTauiw)=.true.
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idChuiw))) THEN
            have_var(idChuiw)=.true.
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idS0mk))) THEN
            have_var(idS0mk)=.true.
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idT0mk))) THEN
            have_var(idT0mk)=.true.
          END IF
      END DO
!
!  Check if model state variables are available in input NetCDF file.
!
      IF (.not.have_var(idtime)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idtime)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idFsur).and.get_var(idFsur)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idFsur)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idRzet).and.get_var(idRzet)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idRzet)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idUbar).and.get_var(idUbar)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idUbar)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idRu2d).and.get_var(idRu2d)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idRu2d)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idVbar).and.get_var(idVbar)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idVbar)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idRv2d).and.get_var(idRv2d)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idRv2d)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idUvel).and.get_var(idUvel)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idUvel)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idRu3d).and.get_var(idRu3d)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idRu3d)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idVvel).and.get_var(idVvel)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idVvel)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idRv3d).and.get_var(idRv3d)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idRv3d)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      DO itrc=1,NT(ng)
        IF (.not.have_var(idTvar(itrc)).and.                            &
     &      get_var(idTvar(itrc))) THEN
          IF (Master) WRITE (stdout,10) string,                         &
     &                                  TRIM(Vname(1,idTvar(itrc))),    &
     &                                  TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END DO
      IF (.not.have_var(idUice)) THEN
        IF (Master) WRITE (stdout,10) string,                           &
     &                                TRIM(Vname(1,idUice)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idVice)) THEN
        IF (Master) WRITE (stdout,10) string,                           &
     &                                TRIM(Vname(1,idVice)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idAice)) THEN
        IF (Master) WRITE (stdout,10) string,                           &
     &                                TRIM(Vname(1,idAice)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idHice)) THEN
        IF (Master) WRITE (stdout,10) string,                           &
     &                                TRIM(Vname(1,idHice)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idTice)) THEN
        IF (Master) WRITE (stdout,10) string,                           &
     &                                TRIM(Vname(1,idTice)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idHsno)) THEN
        IF (Master) WRITE (stdout,10) string,                           &
     &                                TRIM(Vname(1,idHsno)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idSfwat)) THEN
        IF (Master) WRITE (stdout,10) string,                           &
     &                                TRIM(Vname(1,idSfwat)),           &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idAgeice)) THEN
        IF (Master) WRITE (stdout,10) string,                           &
     &                                TRIM(Vname(1,idAgeice)),          &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idTimid)) THEN
        IF (Master) WRITE (stdout,10) string,                           &
     &                                TRIM(Vname(1,idTimid)),           &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idSig11)) THEN
        IF (Master) WRITE (stdout,10) string,                           &
     &                                TRIM(Vname(1,idSig11)),           &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idSig22)) THEN
        IF (Master) WRITE (stdout,10) string,                           &
     &                                TRIM(Vname(1,idSig22)),           &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idSig12)) THEN
        IF (Master) WRITE (stdout,10) string,                           &
     &                                TRIM(Vname(1,idSig12)),           &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idTauiw)) THEN
        IF (Master) WRITE (stdout,10) string,                           &
     &                                TRIM(Vname(1,idTauiw)),           &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idChuiw)) THEN
        IF (Master) WRITE (stdout,10) string,                           &
     &                                TRIM(Vname(1,idChuiw)),           &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idS0mk)) THEN
        IF (Master) WRITE (stdout,10) string,                           &
     &                                TRIM(Vname(1,idS0mk)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idT0mk)) THEN
        IF (Master) WRITE (stdout,10) string,                           &
     &                                TRIM(Vname(1,idT0mk)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
  10  FORMAT (/,a,'CHECKVARS - unable to find model variable: ',a,      &
     &        /,18x,'in file: ',a)
      RETURN
      END SUBROUTINE checkvars
