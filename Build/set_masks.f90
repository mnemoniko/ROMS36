      MODULE set_masks_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines set internal Land/Sea masking arrays that are used   !
!  to process fields into output NetCDF files.  The Land grid points   !
!  are replaced by the _FillValue in the output files to  facilitate   !
!  post-processing with generic tools.                                 !
!                                                                      !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC :: set_masks
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_masks (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
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
      CALL wclock_on (ng, model, 2)
      CALL set_masks_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     GRID(ng) % pmask,                            &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
     &                     GRID(ng) % pmask_avg,                        &
     &                     GRID(ng) % rmask_avg,                        &
     &                     GRID(ng) % umask_avg,                        &
     &                     GRID(ng) % vmask_avg,                        &
     &                     GRID(ng) % pmask_dia,                        &
     &                     GRID(ng) % rmask_dia,                        &
     &                     GRID(ng) % umask_dia,                        &
     &                     GRID(ng) % vmask_dia,                        &
     &                     GRID(ng) % pmask_io,                         &
     &                     GRID(ng) % rmask_io,                         &
     &                     GRID(ng) % umask_io,                         &
     &                     GRID(ng) % vmask_io)
      CALL wclock_off (ng, model, 2)
      RETURN
      END SUBROUTINE set_masks
!
!***********************************************************************
      SUBROUTINE set_masks_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           pmask, rmask,                          &
     &                           umask, vmask,                          &
     &                           pmask_avg, rmask_avg,                  &
     &                           umask_avg, vmask_avg,                  &
     &                           pmask_dia, rmask_dia,                  &
     &                           umask_dia, vmask_dia,                  &
     &                           pmask_io, rmask_io,                    &
     &                           umask_io, vmask_io)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
      USE mp_exchange_mod, ONLY : mp_exchange2d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(inout) :: pmask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: umask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: pmask_dia(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_dia(LBi:,LBj:)
      real(r8), intent(inout) :: umask_dia(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_dia(LBi:,LBj:)
      real(r8), intent(inout) :: pmask_io(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_io(LBi:,LBj:)
      real(r8), intent(inout) :: umask_io(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_io(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j
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
!  Initialize internal history files Land/Sea masks with its respective
!  application grid mask.
!-----------------------------------------------------------------------
!
      DO j=Jstr,JendR
        DO i=Istr,IendR
          pmask_io(i,j)=pmask(i,j)
        END DO
      END DO
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          rmask_io(i,j)=rmask(i,j)
        END DO
      END DO
      DO j=JstrR,JendR
        DO i=Istr,IendR
          umask_io(i,j)=umask(i,j)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vmask_io(i,j)=vmask(i,j)
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask_io)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_io)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_io)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_io)
      END IF
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmask_io, rmask_io, umask_io, vmask_io)
!
!-----------------------------------------------------------------------
!  Initialize average file Land/Sea masks for time-averaged fields.
!-----------------------------------------------------------------------
!
      DO j=Jstr,JendR
        DO i=Istr,IendR
          pmask_avg(i,j)=pmask_io(i,j)
        END DO
      END DO
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          rmask_avg(i,j)=rmask_io(i,j)
        END DO
      END DO
      DO j=JstrR,JendR
        DO i=Istr,IendR
          umask_avg(i,j)=umask_io(i,j)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vmask_avg(i,j)=vmask_io(i,j)
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask_avg)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_avg)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_avg)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_avg)
      END IF
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmask_avg, rmask_avg, umask_avg, vmask_avg)
!
!-----------------------------------------------------------------------
!  Initialize diagnostic file Land/Sea masks for time-averaged fields.
!-----------------------------------------------------------------------
!
      DO j=Jstr,JendR
        DO i=Istr,IendR
          pmask_dia(i,j)=pmask_io(i,j)
        END DO
      END DO
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          rmask_dia(i,j)=rmask_io(i,j)
        END DO
      END DO
      DO j=JstrR,JendR
        DO i=Istr,IendR
          umask_dia(i,j)=umask_io(i,j)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vmask_dia(i,j)=vmask_io(i,j)
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask_dia)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_dia)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_dia)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_dia)
      END IF
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmask_dia, rmask_dia, umask_dia, vmask_dia)
      RETURN
      END SUBROUTINE set_masks_tile
      END MODULE set_masks_mod
