      SUBROUTINE get_data (ng)
!
!svn $Id: get_data.F 1478 2012-06-05 23:30:15Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in forcing, climatology and other data from      !
!  NetCDF files.  If there is more than one time-record,  data is      !
!  loaded into global  two-time  record arrays. The interpolation      !
!  is carried elsewhere.                                               !
!                                                                      !
!  Currently, this routine is only executed in serial mode by the      !
!  main thread.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
      USE mod_clima
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_stepping
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical, dimension(3) :: update =                                 &
     &         (/ .FALSE., .FALSE., .FALSE. /)
      integer :: ILB, IUB, JLB, JUB
      integer :: LBi, UBi, LBj, UBj
      integer :: i, my_tile
!
!  Lower and upper bounds for nontiled (global values) boundary arrays.
!
      my_tile=-1                           ! for global values
      ILB=BOUNDS(ng)%LBi(my_tile)
      IUB=BOUNDS(ng)%UBi(my_tile)
      JLB=BOUNDS(ng)%LBj(my_tile)
      JUB=BOUNDS(ng)%UBj(my_tile)
!
!  Lower and upper bounds for tiled arrays.
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 3)
!
!=======================================================================
!  Read in forcing data from FORCING NetCDF file.
!=======================================================================
!
!
!-----------------------------------------------------------------------
!  Surface wind components.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idUair, ncFRCid(idUair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                GRID(ng) % rmask,                                 &
     &                FORCES(ng) % UwindG)
      IF (exit_flag.ne.NoError) RETURN
      CALL get_2dfld (ng , iNLM, idVair, ncFRCid(idVair,ng),            &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                GRID(ng) % rmask,                                 &
     &                FORCES(ng) % VwindG)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Surface air pressure.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idPair, ncFRCid(idPair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                GRID(ng) % rmask,                                 &
     &                FORCES(ng) % PairG)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Cloud fraction.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idCfra, ncFRCid(idCfra,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                GRID(ng) % rmask,                                 &
     &                FORCES(ng) % cloudG)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Surface air temperature.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idTair, ncFRCid(idTair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                GRID(ng) % rmask,                                 &
     &                FORCES(ng) % TairG)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Surface air humidity.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idQair, ncFRCid(idQair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                GRID(ng) % rmask,                                 &
     &                FORCES(ng) % HairG)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Rain fall rate.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idrain, ncFRCid(idrain,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                GRID(ng) % rmask,                                 &
     &                FORCES(ng) % rainG)
      IF (exit_flag.ne.NoError) RETURN
!
!=======================================================================
!  Read in open boundary conditions from BOUNDARY NetCDF file.
!=======================================================================
!
      IF (LBC(iwest,isFsur,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idZbry(iwest),                        &
     &                  ncBRYid(idZbry(iwest),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % zetaG_west)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(ieast,isFsur,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idZbry(ieast),                        &
     &                  ncBRYid(idZbry(ieast),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % zetaG_east)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(isouth,isFsur,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idZbry(isouth),                       &
     &                  ncBRYid(idZbry(isouth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % zetaG_south)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(inorth,isFsur,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idZbry(inorth),                       &
     &                  ncBRYid(idZbry(inorth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % zetaG_north)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
      IF (LBC(iwest,isUbar,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idU2bc(iwest),                        &
     &                  ncBRYid(idU2bc(iwest),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % ubarG_west)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(iwest,isVbar,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idV2bc(iwest),                        &
     &                  ncBRYid(idV2bc(iwest),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % vbarG_west)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(ieast,isUbar,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idU2bc(ieast),                        &
     &                  ncBRYid(idU2bc(ieast),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % ubarG_east)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(ieast,isVbar,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idV2bc(ieast),                        &
     &                  ncBRYid(idV2bc(ieast),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % vbarG_east)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(isouth,isUbar,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idU2bc(isouth),                       &
     &                  ncBRYid(idU2bc(isouth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % ubarG_south)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(isouth,isVbar,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idV2bc(isouth),                       &
     &                  ncBRYid(idV2bc(isouth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % vbarG_south)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(inorth,isUbar,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idU2bc(inorth),                       &
     &                  ncBRYid(idU2bc(inorth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % ubarG_north)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(inorth,isVbar,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idV2bc(inorth),                       &
     &                  ncBRYid(idV2bc(inorth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % vbarG_north)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
      IF (LBC(iwest,isUvel,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idU3bc(iwest),                        &
     &                  ncBRYid(idU3bc(iwest),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % uG_west)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(iwest,isVvel,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idV3bc(iwest),                        &
     &                  ncBRYid(idV3bc(iwest),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % vG_west)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(ieast,isUvel,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idU3bc(ieast),                        &
     &                  ncBRYid(idU3bc(ieast),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % uG_east)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(ieast,isVvel,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idV3bc(ieast),                        &
     &                  ncBRYid(idV3bc(ieast),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % vG_east)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(isouth,isUvel,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idU3bc(isouth),                       &
     &                  ncBRYid(idU3bc(isouth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % uG_south)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(isouth,isVvel,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idV3bc(isouth),                       &
     &                  ncBRYid(idV3bc(isouth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % vG_south)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(inorth,isUvel,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idU3bc(inorth),                       &
     &                  ncBRYid(idU3bc(inorth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % uG_north)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(inorth,isVvel,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idV3bc(inorth),                       &
     &                  ncBRYid(idV3bc(inorth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % vG_north)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
      DO i=1,NT(ng)
        IF (LBC(iwest,isTvar(i),ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idTbry(iwest,i),                    &
     &                    ncBRYid(idTbry(iwest,i),ng),                  &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % tG_west(:,:,:,i))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
      DO i=1,NT(ng)
        IF (LBC(ieast,isTvar(i),ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idTbry(ieast,i),                    &
     &                    ncBRYid(idTbry(ieast,i),ng),                  &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % tG_east(:,:,:,i))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
      DO i=1,NT(ng)
        IF (LBC(isouth,isTvar(i),ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idTbry(isouth,i),                   &
     &                    ncBRYid(idTbry(isouth,i),ng),                 &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % tG_south(:,:,:,i))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
      DO i=1,NT(ng)
        IF (LBC(inorth,isTvar(i),ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idTbry(inorth,i),                   &
     &                    ncBRYid(idTbry(inorth,i),ng),                 &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % tG_north(:,:,:,i))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
      IF (LBC(iwest,isUice,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idUibc(iwest),                        &
     &                  ncBRYid(idUibc(iwest),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % uiG_west)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(iwest,isVice,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idVibc(iwest),                        &
     &                  ncBRYid(idVibc(iwest),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % viG_west)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(ieast,isUice,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idUibc(ieast),                        &
     &                  ncBRYid(idUibc(ieast),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % uiG_east)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(ieast,isVice,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idVibc(ieast),                        &
     &                  ncBRYid(idVibc(ieast),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % viG_east)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(isouth,isUice,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idUibc(isouth),                       &
     &                  ncBRYid(idUibc(isouth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % uiG_south)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(isouth,isVice,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idVibc(isouth),                       &
     &                  ncBRYid(idVibc(isouth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % viG_south)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(inorth,isUice,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idUibc(inorth),                       &
     &                  ncBRYid(idUibc(inorth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % uiG_north)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(inorth,isVice,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idVibc(inorth),                       &
     &                  ncBRYid(idVibc(inorth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % viG_north)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(iwest,isAice,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idAibc(iwest),                        &
     &                  ncBRYid(idAibc(iwest),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % aiG_west)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(ieast,isAice,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idAibc(ieast),                        &
     &                  ncBRYid(idAibc(ieast),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % aiG_east)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(isouth,isAice,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idAibc(isouth),                       &
     &                  ncBRYid(idAibc(isouth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % aiG_south)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(inorth,isAice,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idAibc(inorth),                       &
     &                  ncBRYid(idAibc(inorth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % aiG_north)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(iwest,isTice,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idTibc(iwest),                        &
     &                  ncBRYid(idTibc(iwest),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % tiG_west)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(ieast,isTice,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idTibc(ieast),                        &
     &                  ncBRYid(idTibc(ieast),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % tiG_east)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(isouth,isTice,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idTibc(isouth),                       &
     &                  ncBRYid(idTibc(isouth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % tiG_south)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(inorth,isTice,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idTibc(inorth),                       &
     &                  ncBRYid(idTibc(inorth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % tiG_north)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(iwest,isSfwat,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idSfwbc(iwest),                       &
     &                  ncBRYid(idSfwbc(iwest),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % sfwatG_west)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(ieast,isSfwat,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idSfwbc(ieast),                       &
     &                  ncBRYid(idSfwbc(ieast),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % sfwatG_east)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(isouth,isSfwat,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idSfwbc(isouth),                      &
     &                  ncBRYid(idSfwbc(isouth),ng),                    &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % sfwatG_south)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(inorth,isSfwat,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idSfwbc(inorth),                      &
     &                  ncBRYid(idSfwbc(inorth),ng),                    &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % sfwatG_north)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(iwest,isSig11,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idS11bc(iwest),                       &
     &                  ncBRYid(idS11bc(iwest),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % sig11G_west)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(ieast,isSig11,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idS11bc(ieast),                       &
     &                  ncBRYid(idS11bc(ieast),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % sig11G_east)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(isouth,isSig11,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idS11bc(isouth),                      &
     &                  ncBRYid(idS11bc(isouth),ng),                    &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % sig11G_south)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(inorth,isSig11,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idS11bc(inorth),                      &
     &                  ncBRYid(idS11bc(inorth),ng),                    &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % sig11G_north)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(iwest,isSig22,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idS22bc(iwest),                       &
     &                  ncBRYid(idS22bc(iwest),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % sig22G_west)
        IF (exit_flag.ne.NoError) RETURN
      IF (LBC(ieast,isSig22,ng)%acquire) THEN
      END IF
        CALL get_ngfld (ng, iNLM, idS22bc(ieast),                       &
     &                  ncBRYid(idS22bc(ieast),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % sig22G_east)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(isouth,isSig22,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idS22bc(isouth),                      &
     &                  ncBRYid(idS22bc(isouth),ng),                    &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % sig22G_south)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(inorth,isSig22,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idS22bc(inorth),                      &
     &                  ncBRYid(idS22bc(inorth),ng),                    &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % sig22G_north)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(iwest,isSig12,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idS12bc(iwest),                       &
     &                  ncBRYid(idS12bc(iwest),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % sig12G_west)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(ieast,isSig12,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idS12bc(ieast),                       &
     &                  ncBRYid(idS12bc(ieast),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % sig12G_east)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(isouth,isSig12,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idS12bc(isouth),                      &
     &                  ncBRYid(idS12bc(isouth),ng),                    &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % sig12G_south)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (LBC(inorth,isSig12,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idS12bc(inorth),                      &
     &                  ncBRYid(idS12bc(inorth),ng),                    &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % sig12G_north)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!=======================================================================
!  Read in climatology data from  NetCDF file.
!=======================================================================
!
!  KLUDGE - how many tracers to read? It depends...
      DO i=1,NT(ng)
        CALL get_3dfld (ng, iNLM, idTclm(i), ncCLMid(idTclm(i),ng),     &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
     &                  GRID(ng) % rmask,                               &
     &                  CLIMA(ng) % tclmG(:,:,:,:,i))
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 3)
      RETURN
      END SUBROUTINE get_data
