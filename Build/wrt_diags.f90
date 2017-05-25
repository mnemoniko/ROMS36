      SUBROUTINE wrt_diags (ng)
!
!svn $Id: wrt_diags.F 1451 2012-02-02 20:56:14Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine writes model time-averaged diagnostic fields into   !
!  diagnostics NetCDF file.                                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_diags
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
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
      integer :: Fcount, gfactor, gtype, ifield, itrc, ivar, status
      real(r8) :: scale
!
      SourceFile='wrt_diags.F'
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out time-averaged diagnostic fields when appropriate.
!-----------------------------------------------------------------------
!
      if (exit_flag.ne.NoError) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
        gfactor=1
!
!  Set time and time-record index.
!
      DIA(ng)%Rindex=DIA(ng)%Rindex+1
      Fcount=DIA(ng)%Fcount
      DIA(ng)%Nrec(Fcount)=DIA(ng)%Nrec(Fcount)+1
!
!  Write out averaged time.
!
      CALL netcdf_put_fvar (ng, iNLM, DIA(ng)%name,                     &
     &                      TRIM(Vname(idtime,ng)), DIAtime(ng:),       &
     &                      (/DIA(ng)%Rindex/), (/1/),                  &
     &                      ncid = DIA(ng)%ncid,                        &
     &                      varid = DIA(ng)%Vid(idtime))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out tracer diagnostic fields.
!
      DO itrc=1,3,1                   ! Temperature, Salt, and CDW
        DO ivar=1,NDT
          ifield=idDtrc(itrc,ivar)
          IF (Dout(ifield,ng)) THEN
            scale=1.0_r8/dt(ng)
            gtype=gfactor*r3dvar
            status=nf_fwrite3d(ng, iNLM, DIA(ng)%ncid,                  &
     &                         DIA(ng)%Vid(ifield),                     &
     &                         DIA(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
     &                         GRID(ng) % rmask_dia,                    &
     &                         DIAGS(ng) % DiaTrc(:,:,:,itrc,ivar))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vname(1,ifield)), DIA(ng)%Rindex
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END DO
      END DO
!
!  Synchronize time-average NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!
      CALL netcdf_sync (ng, iNLM, DIA(ng)%name, DIA(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
      IF (Master) WRITE (stdout,20) DIA(ng)%Rindex
!
  10  FORMAT (/,' WRT_DIAGS - error while writing variable: ',a,/,13x,  &
     &        'into diagnostics NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_DIAGS - wrote diagnostics fields into time ',     &
     &        'record = ',t72,i7.7)
      RETURN
      END SUBROUTINE wrt_diags
