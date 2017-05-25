      SUBROUTINE main3d (RunInterval)
!
!svn $Id: main3d.F 1474 2012-05-29 15:15:08Z kate $
!=======================================================================
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine is the main driver for nonlinear ROMS/TOMS when     !
!  configurated as a full 3D baroclinic ocean model.  It  advances     !
!  forward the primitive equations for all  nested  grids, if any,     !
!  for the specified time interval (seconds), RunInterval.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE mod_stepping
!
      USE bulk_flux_mod, ONLY : bulk_flux
      USE diag_mod, ONLY : diag
      USE ini_fields_mod, ONLY : ini_fields, ini_zeta
      USE lmd_vmix_mod, ONLY : lmd_vmix
      USE omega_mod, ONLY : omega
      USE rho_eos_mod, ONLY : rho_eos
      USE rhs3d_mod, ONLY : rhs3d
      USE set_avg_mod, ONLY : set_avg
      USE set_depth_mod, ONLY : set_depth
      USE set_massflux_mod, ONLY : set_massflux
      USE set_vbc_mod, ONLY : set_vbc
      USE set_zeta_mod, ONLY : set_zeta
      USE step2d_mod, ONLY : step2d
      USE step3d_t_mod, ONLY : step3d_t
      USE step3d_uv_mod, ONLY : step3d_uv
      USE wvelocity_mod, ONLY : wvelocity
!
      implicit none
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      integer :: ng, tile
      integer :: my_iif, next_indx1
      real(r8) :: MaxDT, my_StepTime
!
!=======================================================================
!  Time-step nonlinear 3D primitive equations by the specified time.
!=======================================================================
!
      my_StepTime=0.0_r8
      MaxDT=MAXVAL(dt)
      STEP_LOOP : DO WHILE (my_StepTime.le.(RunInterval+0.5_r8*MaxDT))
        my_StepTime=my_StepTime+MaxDT
!
!  Set time indices and time clock.
!
        DO ng=1,Ngrids
          iic(ng)=iic(ng)+1
          nstp(ng)=1+MOD(iic(ng)-ntstart(ng),2)
          nnew(ng)=3-nstp(ng)
          nrhs(ng)=nstp(ng)
          time(ng)=time(ng)+dt(ng)
          tdays(ng)=time(ng)*sec2day
          CALL time_string (time(ng), time_code(ng))
        END DO
!
!-----------------------------------------------------------------------
!  Read in required data, if any, from input NetCDF files.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          CALL get_data (ng)
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
!-----------------------------------------------------------------------
!  If applicable, process input data: time interpolate between data
!  snapshots.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL set_data (ng, tile)
          END DO
        END DO
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Initialize all time levels and compute other initial fields.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (iic(ng).eq.ntstart(ng)) THEN
!
!  Initialize free-surface and compute initial level thicknesses and
!  depths.
!
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ini_zeta (ng, tile, iNLM)
              CALL set_depth (ng, tile)
            END DO
!
!  Initialize other state variables.
!
            DO tile=last_tile(ng),first_tile(ng),-1
              CALL ini_fields (ng, tile, iNLM)
            END DO
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Compute horizontal mass fluxes (Hz*u/n and Hz*v/m), density related
!  quatities and report global diagnostics.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL set_massflux (ng, tile)
            CALL rho_eos (ng, tile)
            CALL diag (ng, tile)
          END DO
        END DO
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Set fields for vertical boundary conditions. Process tidal forcing,
!  if any.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL bulk_flux (ng, tile)
            CALL set_vbc (ng, tile)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Run ice model for one step
!-----------------------------------------------------------------------
!
          CALL seaice
!
!-----------------------------------------------------------------------
!  Compute time-dependent vertical/horizontal mixing coefficients for
!  momentum and tracers. Compute S-coordinate vertical velocity,
!  diagnostically from horizontal mass divergence.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=last_tile(ng),first_tile(ng),-1
            CALL lmd_vmix (ng, tile)
            CALL omega (ng, tile)
            CALL wvelocity (ng, tile, nstp(ng))
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Set free-surface to it time-averaged value.  If applicable,
!  accumulate time-averaged output data which needs a irreversible
!  loop in shared-memory jobs.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1         ! irreversible
            CALL set_zeta (ng, tile)
            CALL set_diags (ng, tile)
            CALL set_avg (ng, tile)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  Exit if last
!  time step.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          CALL output (ng)
          IF ((exit_flag.ne.NoError).or.                                &
     &        ((iic(ng).eq.(ntend(ng)+1)).and.(ng.eq.Ngrids))) RETURN
        END DO
!
!-----------------------------------------------------------------------
!  Compute right-hand-side terms for 3D equations.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=last_tile(ng),first_tile(ng),-1
            CALL rhs3d (ng, tile)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Solve the vertically integrated primitive equations for the
!  free-surface and barotropic momentum components.
!-----------------------------------------------------------------------
!
        LOOP_2D : DO my_iif=1,MAXVAL(nfast)+1
!
!  Set time indices for predictor step. The PREDICTOR_2D_STEP switch
!  it is assumed to be false before the first time-step.
!
          DO ng=1,Ngrids
            next_indx1=3-indx1(ng)
            IF (.not.PREDICTOR_2D_STEP(ng).and.                         &
     &          my_iif.le.(nfast(ng)+1)) THEN
              PREDICTOR_2D_STEP(ng)=.TRUE.
              iif(ng)=my_iif
              IF (iif(ng).eq.1) THEN
                kstp(ng)=indx1(ng)
              ELSE
                kstp(ng)=3-indx1(ng)
              END IF
              knew(ng)=3
              krhs(ng)=indx1(ng)
            END IF
!
!  Predictor step - Advance barotropic equations using 2D time-step
!  ==============   predictor scheme.  No actual time-stepping is
!  performed during the auxiliary (nfast+1) time-step. It is needed
!  to finalize the fast-time averaging of 2D fields, if any, and
!  compute the new time-evolving depths.
!
            IF (my_iif.le.(nfast(ng)+1)) THEN
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL step2d (ng, tile)
              END DO
            END IF
          END DO
!
!  Set time indices for corrector step.
!
          DO ng=1,Ngrids
            IF (PREDICTOR_2D_STEP(ng)) THEN
              PREDICTOR_2D_STEP(ng)=.FALSE.
              knew(ng)=next_indx1
              kstp(ng)=3-knew(ng)
              krhs(ng)=3
              IF (iif(ng).lt.(nfast(ng)+1)) indx1(ng)=next_indx1
            END IF
!
!  Corrector step - Apply 2D time-step corrector scheme.  Notice that
!  ==============   there is not need for a corrector step during the
!  auxiliary (nfast+1) time-step.
!
            IF (iif(ng).lt.(nfast(ng)+1)) THEN
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL step2d (ng, tile)
              END DO
            END IF
          END DO
        END DO LOOP_2D
!
!-----------------------------------------------------------------------
!  Recompute depths and thicknesses using the new time filtered
!  free-surface.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=last_tile(ng),first_tile(ng),-1
            CALL set_depth (ng, tile)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Time-step 3D momentum equations.
!-----------------------------------------------------------------------
!
!  Time-step 3D momentum equations and couple with vertically
!  integrated equations.
!
        DO ng=1,Ngrids
          DO tile=last_tile(ng),first_tile(ng),-1
            CALL step3d_uv (ng, tile)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Time-step vertical mixing turbulent equations and passive tracer
!  source and sink terms, if applicable.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL omega (ng, tile)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Time-step tracer equations.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=last_tile(ng),first_tile(ng),-1
            CALL step3d_t (ng, tile)
          END DO
        END DO
      END DO STEP_LOOP
      RETURN
      END SUBROUTINE main3d
