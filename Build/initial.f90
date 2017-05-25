      SUBROUTINE initial (ng)
!
!svn $Id: initial.F 1474 2012-05-29 15:15:08Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes all model variables.                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE analytical_mod
      USE distribute_mod, ONLY : mp_bcasti
      USE ini_hmixcoef_mod, ONLY : ini_hmixcoef
      USE metrics_mod, ONLY : metrics
      USE set_depth_mod, ONLY : set_depth
      USE omega_mod, ONLY : omega
      USE rho_eos_mod, ONLY : rho_eos
      USE set_massflux_mod, ONLY : set_massflux
      USE set_masks_mod, ONLY : set_masks
      USE stiffness_mod, ONLY : stiffness
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical, save :: First = .TRUE.
      logical :: update = .FALSE.
      integer :: LBi, UBi, LBj, UBj
      integer :: Fcount, IniRec, Tindex
      integer :: thread, tile
!
!=======================================================================
!   Initialize model variables.
!=======================================================================
!
      IF (Master) THEN
        WRITE (stdout,20) 'INITIAL: Configuring and initializing ',     &
     &                    'forward nonlinear model ...'
 20     FORMAT (/,1x,a,a,/)
      END IF
!
!-----------------------------------------------------------------------
!  Initialize time stepping indices and counters.
!-----------------------------------------------------------------------
!
      iif(ng)=1
      indx1(ng)=1
      kstp(ng)=1
      krhs(ng)=1
      knew(ng)=1
      PREDICTOR_2D_STEP(ng)=.FALSE.
!
      iic(ng)=0
      nstp(ng)=1
      nrhs(ng)=1
      nnew(ng)=1
!
      synchro_flag(ng)=.TRUE.
      first_time(ng)=0
      tdays(ng)=dstart
      time(ng)=tdays(ng)*day2sec
      ntstart(ng)=INT((time(ng)-dstart*day2sec)/dt(ng))+1
      ntend(ng)=ntimes(ng)
      ntfirst(ng)=ntstart(ng)
      CALL time_string (time(ng), time_code(ng))
      IniRec=nrrec(ng)
      Tindex=1
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Start time wall clocks.
!-----------------------------------------------------------------------
!
      DO thread=MyRank,MyRank
        CALL wclock_on (ng, iNLM, 2)
      END DO
!
!=======================================================================
!  On first pass of ensemble/perturbation/iteration loop, initialize
!  model configuration.
!=======================================================================
!
      IF (Nrun.eq.ERstr) THEN
!
!-----------------------------------------------------------------------
!  Set horizontal grid, bathymetry, and Land/Sea masking (if any).
!  Use analytical functions or read in from a grid NetCDF.
!-----------------------------------------------------------------------
!
        CALL get_grid (ng, iNLM)
        CALL mp_bcasti (ng, iNLM, exit_flag)
        if (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Set vertical S-coordinate transformation function.
!-----------------------------------------------------------------------
!
        CALL set_scoord (ng)
!
!-----------------------------------------------------------------------
!  Set barotropic time-steps average weighting function.
!-----------------------------------------------------------------------
!
        CALL set_weights (ng)
!
!-----------------------------------------------------------------------
!  Compute various metric term combinations.
!-----------------------------------------------------------------------
!
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL metrics (ng, tile, iNLM)
        END DO
!
!-----------------------------------------------------------------------
!  If appropriate, set nudging coefficiests time scales.
!-----------------------------------------------------------------------
!
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_nudgcoef (ng, tile, iNLM)
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Initialize horizontal mixing coefficients.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL ini_hmixcoef (ng, tile, iNLM)
      END DO
!
!=======================================================================
!  Initialize model state variables and forcing.  This part is
!  executed for each ensemble/perturbation/iteration run.
!=======================================================================
!
!-----------------------------------------------------------------------
!  Set primitive variables initial conditions.
!-----------------------------------------------------------------------
!
!  Read in initial conditions from initial NetCDF file.
!
      CALL get_state (ng, iNLM, 1, INI(ng)%name, IniRec, Tindex)
      CALL mp_bcasti (ng, iNLM, exit_flag)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Compute initial time-evolving depths.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_depth (ng, tile)
      END DO
!
!-----------------------------------------------------------------------
!  Compute initial horizontal mass fluxes, Hz*u/n and Hz*v/m.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_massflux (ng, tile)
      END DO
!
!-----------------------------------------------------------------------
!  Compute initial S-coordinates vertical velocity. Compute initial
!  density anomaly from potential temperature and salinity via equation
!  of state for seawater.  Also compute other equation of state related
!  quatities.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL omega (ng, tile)
        CALL rho_eos (ng, tile)
      END DO
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------
!
!
!  If applicable, read in input data.
!
      CALL close_inp (ng, iNLM)
      CALL get_data (ng)
      CALL mp_bcasti (ng, iNLM, exit_flag)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Set internal I/O mask arrays.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_masks (ng, tile, iNLM)
      END DO
!
!-----------------------------------------------------------------------
!  Compute grid stiffness.
!-----------------------------------------------------------------------
!
      IF (Lstiffness) THEN
        Lstiffness=.FALSE.
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL stiffness (ng, tile, iNLM)
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  If applicable, convert initial locations to fractional grid
!  coordinates.
!-----------------------------------------------------------------------
!
      CALL grid_coords (ng, iNLM)
!
!-----------------------------------------------------------------------
!  Initialize time-stepping counter and clock.
!-----------------------------------------------------------------------
!
!  Subsract one time unit to avoid special case due to initialization
!  in the main time-stepping routine.
!
      iic(ng)=ntstart(ng)-1
      time(ng)=time(ng)-dt(ng)
!
!-----------------------------------------------------------------------
!  Turn off initialization time wall clock.
!-----------------------------------------------------------------------
!
      DO thread=MyRank,MyRank
        CALL wclock_off (ng, iNLM, 2)
      END DO
      RETURN
      END SUBROUTINE initial
