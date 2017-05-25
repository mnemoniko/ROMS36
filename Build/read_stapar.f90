      SUBROUTINE read_StaPar (model, inp, out, Lwrite)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads and reports stations input parameters.           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      integer :: Mstation, Npts, Nval
      integer :: flag, i, igrid, ista, itrc, ng, status
      integer :: decode_line, load_i, load_l, load_r
      real(r8) :: Xpos, Ypos
      logical, dimension(MT,Ngrids) :: Lsta
      integer, dimension(Ngrids) :: is
      real(r8), dimension(200) :: Rval
      character (len=40 ) :: KeyWord
      character (len=256) :: line
      character (len=256), dimension(200) :: Cval
!
!-----------------------------------------------------------------------
!  Read in stations parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=20,END=30) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('Lstations')
              Npts=load_l(Nval, Cval, Ngrids, Lstations)
            CASE ('Sout(idUvel)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idUvel,:))
            CASE ('Sout(idVvel)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idVvel,:))
            CASE ('Sout(idu3dE)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idu3dE,:))
            CASE ('Sout(idv3dN)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idv3dN,:))
            CASE ('Sout(idWvel)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idWvel,:))
            CASE ('Sout(idOvel)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idOvel,:))
            CASE ('Sout(idUbar)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idUbar,:))
            CASE ('Sout(idVbar)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idVbar,:))
            CASE ('Sout(idu2dE)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idu2dE,:))
            CASE ('Sout(idv2dN)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idv2dN,:))
            CASE ('Sout(idFsur)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idFsur,:))
            CASE ('Sout(idTvar)')
              Npts=load_l(Nval, Cval, MT*Ngrids, Lsta)
              DO ng=1,Ngrids
                DO itrc=1,NT(ng)
                  Sout(idTvar(itrc),ng)=Lsta(itrc,ng)
                END DO
              END DO
            CASE ('Sout(idUsms)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idUsms,:))
            CASE ('Sout(idVsms)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idVsms,:))
            CASE ('Sout(idUbms)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idUbms,:))
            CASE ('Sout(idVbms)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idVbms,:))
            CASE ('Sout(idUbrs)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idUbrs,:))
            CASE ('Sout(idVbrs)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idVbrs,:))
            CASE ('Sout(idUbws)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idUbws,:))
            CASE ('Sout(idVbws)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idVbws,:))
            CASE ('Sout(idUbcs)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idUbcs,:))
            CASE ('Sout(idVbcs)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idVbcs,:))
            CASE ('Sout(idUbot)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idUbot,:))
            CASE ('Sout(idVbot)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idVbot,:))
            CASE ('Sout(idUbur)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idUbur,:))
            CASE ('Sout(idVbvr)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idVbvr,:))
            CASE ('Sout(idW2xx)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idW2xx,:))
            CASE ('Sout(idW2xy)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idW2xy,:))
            CASE ('Sout(idW2yy)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idW2yy,:))
            CASE ('Sout(idU2rs)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idU2rs,:))
            CASE ('Sout(idV2rs)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idV2rs,:))
            CASE ('Sout(idU2Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idU2Sd,:))
            CASE ('Sout(idV2Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idV2Sd,:))
            CASE ('Sout(idW3xx)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idW3xx,:))
            CASE ('Sout(idW3xy)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idW3xy,:))
            CASE ('Sout(idW3yy)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idW3yy,:))
            CASE ('Sout(idW3zx)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idW3zx,:))
            CASE ('Sout(idW3zy)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idW3zy,:))
            CASE ('Sout(idU3rs)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idU3rs,:))
            CASE ('Sout(idV3rs)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idV3rs,:))
            CASE ('Sout(idU3Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idU3Sd,:))
            CASE ('Sout(idV3Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idV3Sd,:))
            CASE ('Sout(idWamp)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idWamp,:))
            CASE ('Sout(idWlen)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idWlen,:))
            CASE ('Sout(idWdir)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idWdir,:))
            CASE ('Sout(idPair)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idPair,:))
            CASE ('Sout(idUair)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idUair,:))
            CASE ('Sout(idVair)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idVair,:))
            CASE ('Sout(idUairE)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idUairE,:))
            CASE ('Sout(idVairN)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idVairN,:))
            CASE ('Sout(idTsur)')
              Npts=load_l(Nval, Cval, NAT*Ngrids, Lsta)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Sout(idTsur(itrc),ng)=Lsta(itrc,ng)
                END DO
              END DO
            CASE ('Sout(idLhea)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idLhea,:))
            CASE ('Sout(idShea)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idShea,:))
            CASE ('Sout(idLrad)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idLrad,:))
            CASE ('Sout(idSrad)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idSrad,:))
            CASE ('Sout(idEmPf)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idEmPf,:))
            CASE ('Sout(idevap)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idevap,:))
            CASE ('Sout(idrain)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idrain,:))
            CASE ('Sout(idDano)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idDano,:))
            CASE ('Sout(idVvis)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idVvis,:))
            CASE ('Sout(idTdif)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idTdif,:))
            CASE ('Sout(idSdif)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idSdif,:))
            CASE ('Sout(idHsbl)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idHsbl,:))
            CASE ('Sout(idHbbl)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idHbbl,:))
            CASE ('Sout(idMtke)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idMtke,:))
            CASE ('Sout(idMtls)')
              Npts=load_l(Nval, Cval, Ngrids, Sout(idMtls,:))
            CASE ('NSTATION')
              Npts=load_i(Nval, Rval, Ngrids, Nstation)
            CASE ('POS')
              DO ng=1,Ngrids
                allocate ( SCALARS(ng) % Sflag(Nstation(ng)) )
                allocate ( SCALARS(ng) % SposX(Nstation(ng)) )
                allocate ( SCALARS(ng) % SposY(Nstation(ng)) )
              END DO
              is(1:Ngrids)=0
              DO WHILE (.TRUE.)
                READ (inp,*,ERR=10,END=10) igrid, flag, Xpos, Ypos
                ng=MAX(1,ABS(igrid))
                is(ng)=is(ng)+1
                SCALARS(ng)%Sflag(is(ng))=flag
                SCALARS(ng)%SposX(is(ng))=Xpos
                SCALARS(ng)%SposY(is(ng))=Ypos
              END DO
 10           DO ng=1,Ngrids
                IF (Nstation(ng).ne.is(ng)) THEN
                  IF (Master) WRITE (out,40) Nstation(ng), is(ng)
                  exit_flag=4
                  RETURN
                END IF
              END DO
          END SELECT
        END IF
      END DO
 20   IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
 30   CONTINUE
!
!-----------------------------------------------------------------------
!  Process input parameters.
!-----------------------------------------------------------------------
!
!  Turn off the processing of stations if not running long enough to
!  create a stations file (LdefSTA=.FALSE. because nSTA < ntimes or
!  nSTA = 0 when nrrec = 0).
!
      DO ng=1,Ngrids
        IF (.not.LdefSTA(ng).and.Lstations(ng)) THEN
          Lstations(ng)=.FALSE.
        END IF
      END DO
!
!  Make sure that both component switches are activated when processing
!  (Eastward,Northward) momentum components at RHO-points.
!
      DO ng=1,Ngrids
        IF (.not.Sout(idu2dE,ng).and.Sout(idv2dN,ng)) THEN
          Sout(idu2dE,ng)=.TRUE.
        END IF
        IF (.not.Sout(idv2dN,ng).and.Sout(idu2dE,ng)) THEN
          Sout(idv2dN,ng)=.TRUE.
        END IF
        IF (.not.Sout(idu3dE,ng).and.Sout(idv3dN,ng)) THEN
          Sout(idu3dE,ng)=.TRUE.
        END IF
        IF (.not.Sout(idv3dN,ng).and.Sout(idu3dE,ng)) THEN
          Sout(idv3dN,ng)=.TRUE.
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lstations(ng)) THEN
            WRITE (out,60) ng
            WRITE (out,70) Nstation(ng), 'Nstation',                    &
     &            'Number of stations to write out into stations file.'
            IF (Sout(idFsur,ng)) WRITE (out,80) Sout(idFsur,ng),        &
     &          'Sout(idFsur)',                                         &
     &          'Write out free-surface.'
            IF (Sout(idUbar,ng)) WRITE (out,80) Sout(idUbar,ng),        &
     &          'Sout(idUbar)',                                         &
     &          'Write out 2D U-momentum component.'
            IF (Sout(idVbar,ng)) WRITE (out,80) Sout(idVbar,ng),        &
     &          'Sout(idVbar)',                                         &
     &          'Write out 2D V-momentum component.'
            IF (Sout(idu2dE,ng)) WRITE (out,80) Sout(idu2dE,ng),        &
     &          'Sout(idu2dE)',                                         &
     &          'Write out 2D U-eastward  at RHO-points.'
            IF (Sout(idv2dN,ng)) WRITE (out,80) Sout(idv2dN,ng),        &
     &          'Sout(idv2dN)',                                         &
     &          'Write out 2D V-northward at RHO-points.'
            IF (Sout(idUvel,ng)) WRITE (out,80) Sout(idUvel,ng),        &
     &          'Sout(idUvel)',                                         &
     &          'Write out 3D U-momentum component.'
            IF (Sout(idVvel,ng)) WRITE (out,80) Sout(idVvel,ng),        &
     &          'Sout(idVvel)',                                         &
     &          'Write out 3D V-momentum component.'
            IF (Sout(idu3dE,ng)) WRITE (out,80) Sout(idu3dE,ng),        &
     &          'Sout(idu3dE)',                                         &
     &          'Write out 3D U-eastward  at RHO-points.'
            IF (Sout(idv3dN,ng)) WRITE (out,80) Sout(idv3dN,ng),        &
     &          'Sout(idv3dN)',                                         &
     &          'Write out 3D V-northward at RHO-points.'
            IF (Sout(idWvel,ng)) WRITE (out,80) Sout(idWvel,ng),        &
     &          'Sout(idWvel)',                                         &
     &          'Write out W-momentum component.'
            IF (Sout(idOvel,ng)) WRITE (out,80) Sout(idOvel,ng),        &
     &          'Sout(idOvel)',                                         &
     &          'Write out omega vertical velocity.'
            DO itrc=1,NT(ng)
              IF (Sout(idTvar(itrc),ng)) WRITE (out,90)                 &
     &            Sout(idTvar(itrc),ng), 'Sout(idTvar)',                &
     &            'Write out tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
            END DO
            IF (Sout(idUsms,ng)) WRITE (out,80) Sout(idUsms,ng),        &
     &          'Sout(idUsms)',                                         &
     &          'Write out surface U-momentum stress.'
            IF (Sout(idVsms,ng)) WRITE (out,80) Sout(idVsms,ng),        &
     &          'Sout(idVsms)',                                         &
     &          'Write out surface V-momentum stress.'
            IF (Sout(idUbms,ng)) WRITE (out,80) Sout(idUbms,ng),        &
     &          'Sout(idUbms)',                                         &
     &          'Write out bottom U-momentum stress.'
            IF (Sout(idVbms,ng)) WRITE (out,80) Sout(idVbms,ng),        &
     &          'Sout(idVbms)',                                         &
     &          'Write out bottom V-momentum stress.'
            IF (Sout(idPair,ng)) WRITE (out,80) Sout(idPair,ng),        &
     &          'Sout(idPair)',                                         &
     &          'Write out surface air pressure.'
            IF (Sout(idUair,ng)) WRITE (out,80) Sout(idUair,ng),        &
     &          'Sout(idUair)',                                         &
     &          'Write out surface U-wind component.'
            IF (Sout(idVair,ng)) WRITE (out,80) Sout(idVair,ng),        &
     &          'Sout(idVair)',                                         &
     &          'Write out surface V-wind component.'
            IF (Sout(idUairE,ng)) WRITE (out,80) Sout(idUairE,ng),      &
     &          'Sout(idUairE)',                                        &
     &          'Write out surface Eastward U-wind component.'
            IF (Sout(idVairN,ng)) WRITE (out,80) Sout(idVairN,ng),      &
     &          'Sout(idVairN)',                                        &
     &          'Write out surface Northward V-wind component.'
            IF (Sout(idTsur(itemp),ng)) WRITE (out,80)                  &
     &          Sout(idTsur(itemp),ng), 'Sout(idTsur)',                 &
     &          'Write out surface net heat flux.'
            IF (Sout(idTsur(isalt),ng)) WRITE (out,80)                  &
     &          Sout(idTsur(isalt),ng), 'Sout(idTsur)',                 &
     &          'Write out surface net salt flux.'
            IF (Sout(idSrad,ng)) WRITE (out,80) Sout(idSrad,ng),        &
     &          'Sout(idSrad)',                                         &
     &          'Write out shortwave radiation flux.'
            IF (Sout(idLrad,ng)) WRITE (out,80) Sout(idLrad,ng),        &
     &          'Sout(idLrad)',                                         &
     &          'Write out longwave radiation flux.'
            IF (Sout(idLhea,ng)) WRITE (out,80) Sout(idLhea,ng),        &
     &          'Sout(idLhea)',                                         &
     &          'Write out latent heat flux.'
            IF (Sout(idShea,ng)) WRITE (out,80) Sout(idShea,ng),        &
     &          'Sout(idShea)',                                         &
     &          'Write out sensible heat flux.'
            IF (Sout(idEmPf,ng)) WRITE (out,80) Sout(idEmPf,ng),        &
     &         'Sout(idEmPf)',                                          &
     &         'Write out E-P flux.'
            IF (Sout(idevap,ng)) WRITE (out,80) Sout(idevap,ng),        &
     &         'Sout(idevap)',                                          &
     &         'Write out evaporation rate.'
            IF (Sout(idrain,ng)) WRITE (out,80) Sout(idrain,ng),        &
     &         'Sout(idrain)',                                          &
     &         'Write out rain rate.'
            IF (Sout(idDano,ng)) WRITE (out,80) Sout(idDano,ng),        &
     &          'Sout(idDano)',                                         &
     &          'Write out density anomaly.'
            IF (Sout(idVvis,ng)) WRITE (out,80) Sout(idVvis,ng),        &
     &          'Sout(idVvis)',                                         &
     &          'Write out vertical viscosity coefficient.'
            IF (Sout(idTdif,ng)) WRITE (out,80) Sout(idTdif,ng),        &
     &          'Sout(idTdif)',                                         &
     &          'Write out vertical T-diffusion coefficient.'
            IF (Sout(idSdif,ng)) WRITE (out,80) Sout(idSdif,ng),        &
     &          'Sout(idSdif)',                                         &
     &          'Write out vertical S-diffusion coefficient.'
            IF (Sout(idHsbl,ng)) WRITE (out,80) Sout(idHsbl,ng),        &
     &          'Sout(idHsbl)',                                         &
     &          'Write out depth of surface boundary layer.'
            WRITE (out,*)
            DO i=1,Nstation(ng)
              WRITE (out,100) i, SCALARS(ng)%Sflag(i),                  &
     &                           SCALARS(ng)%SposX(i),                  &
     &                           SCALARS(ng)%SposY(i)
            END DO
          END IF
        END DO
      END IF
  40  FORMAT (/,' READ_StaPar - Inconsistent number of stations, ',     &
     &        'Nstation = ',2i8,/,15x,'change input script values.')
  50  FORMAT (/,' READ_StaPar - Error while processing line: ',/,a)
  60  FORMAT (/,/,' Stations Parameters, Grid: ',i2.2,                  &
     &        /,  ' =============================',/)
  70  FORMAT (1x,i10,2x,a,t30,a)
  80  FORMAT (10x,l1,2x,a,t30,a)
  90  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
 100  FORMAT (13x,'Flag and positions for station ',i4.4,':',           &
     &        i3,1x,2f10.4)
 110  FORMAT (/,' READ_StaPAR - variable info not yet loaded, ', a)
      RETURN
      END SUBROUTINE read_StaPar
