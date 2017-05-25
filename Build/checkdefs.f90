      SUBROUTINE checkdefs
!
!svn $Id: checkdefs.F 1484 2012-06-13 19:25:03Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine checks activated C-preprocessing options for        !
!  consistency.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE mod_strings
!
      USE strings_mod, ONLY : uppercase
!
      implicit none
!
!  Local variable declarations.
!
      integer :: ibbl = 0
      integer :: ibiology = 0
      integer :: idriver = 0
      integer :: itrcHadv = 0
      integer :: itrcVadv = 0
      integer :: itrcHadvtl = 0
      integer :: itrcVadvtl = 0
      integer :: ivelHadv = 0
      integer :: ivelVadv = 0
      integer :: ivmix = 0
      integer :: nearshore = 0
      integer :: is, lstr, ng
!
!-----------------------------------------------------------------------
!  Report activated C-preprocessing options.
!-----------------------------------------------------------------------
!
      Coptions=' '
      IF (Master) WRITE (stdout,10)
  10  FORMAT (/,' Activated C-preprocessing Options:',/)
  20  FORMAT (1x,a,t22,a)
!
      IF (Master) THEN
        WRITE (stdout,20) TRIM(ADJUSTL(MyAppCPP)), TRIM(ADJUSTL(title))
      END IF
      is=LEN_TRIM(Coptions)+1
      lstr=LEN_TRIM(MyAppCPP)
      Coptions(is:is+lstr)=TRIM(ADJUSTL(MyAppCPP))
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is)=','
!
      IF (Master) WRITE (stdout,20) 'ANA_BPFLUX',                       &
     &   'Analytical bottom passive tracers fluxes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BPFLUX,'
      IF (Master) WRITE (stdout,20) 'ANA_BSFLUX',                       &
     &   'Analytical kinematic bottom salinity flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BSFLUX,'
      IF (Master) WRITE (stdout,20) 'ANA_BTFLUX',                       &
     &   'Analytical kinematic bottom temperature flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BTFLUX,'
      IF (Master) WRITE (stdout,20) 'ANA_SPFLUX',                       &
     &   'Analytical surface passive tracer fluxes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_SPFLUX,'
      IF (Master) WRITE (stdout,20) 'ANA_SRFLUX',                       &
     &   'Analytical kinematic shortwave radiation flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_SRFLUX,'
      IF (Master) WRITE (stdout,20) 'ASSUMED_SHAPE',                    &
     &   'Using assumed-shape arrays.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' ASSUMED_SHAPE,'
      IF (Master) WRITE (stdout,20) 'AVERAGES',                         &
     &   'Writing out time-averaged nonlinear model fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' AVERAGES,'
      IF (Master) WRITE (stdout,20) 'BULK_FLUXES',                      &
     &   'Surface bulk fluxes parameterization.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' BULK_FLUXES,'
      IF (Master) WRITE (stdout,20) 'DIAGNOSTICS_TS',                   &
     &   'Computing and writing tracer diagnostic terms.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' DIAGNOSTICS_TS,'
      IF (Master) WRITE (stdout,20) 'DJ_GRADPS',                        &
     &   'Parabolic Splines density Jacobian (Shchepetkin, 2002).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' DJ_GRADPS,'
      IF (Master) WRITE (stdout,20) 'DOUBLE_PRECISION',                 &
     &   'Double precision arithmetic.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' DOUBLE_PRECISION,'
      IF (Master) WRITE (stdout,20) 'EMINUSP',                          &
     &   'Compute Salt Flux using E-P.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' EMINUSP,'
      IF (Master) WRITE (stdout,20) 'FORWARD_WRITE',                    &
     &   'Write out Forward solution for Tangent/Adjoint.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' FORWARD_WRITE,'
      IF (Master) WRITE (stdout,20) 'ICE_ADVECT',                       &
     &   'Advection of ice tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ICE_ADVECT,'
      IF (Master) WRITE (stdout,20) 'ICE_ALB_EC92',                     &
     &   'Ebert and Curry albedo formula.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ICE_ALB_EC92,'
      IF (Master) WRITE (stdout,20) 'ICE_BULK_FLUXES',                  &
     &   'Ice bulk fluxes from the atmosphere.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' ICE_BULK_FLUXES,'
      IF (Master) WRITE (stdout,20) 'ICE_CONVSNOW',                     &
     &   'Conversion of flooded snow to ice.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ICE_CONVSNOW,'
      IF (Master) WRITE (stdout,20) 'ICE_EVP',                          &
     &   'Elastic-viscous-plastic ice rheology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' ICE_EVP,'
      IF (Master) WRITE (stdout,20) 'ICE_MK',                           &
     &   'Mellor-Kantha ice thermodynamics.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' ICE_MK,'
      IF (Master) WRITE (stdout,20) 'ICE_MODEL',                        &
     &   'Include sea ice model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ICE_MODEL,'
      IF (Master) WRITE (stdout,20) 'ICE_MOMENTUM',                     &
     &   'Compute ice momentum equations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ICE_MOMENTUM,'
      IF (Master) WRITE (stdout,20) 'ICE_SMOLAR',                       &
     &   'Advect ice tracers with MPDATA scheme.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ICE_SMOLAR,'
      IF (Master) WRITE (stdout,20) 'ICE_THERMO',                       &
     &   'Include ice thermodynamics.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ICE_THERMO,'
      IF (Master) WRITE (stdout,20) 'ICESHELF',                         &
     &   'Include Ice Shelf Cavities.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ICESHELF,'
      IF (Master) WRITE (stdout,20) 'ICESHELF_3EQ',                     &
     &   'Include 3eq Ice Shelf Thermodynamics.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ICESHELF_3EQ,'
      IF (Master) WRITE (stdout,20) 'LMD_CONVEC',                       &
     &   'LMD convective mixing due to shear instability.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' LMD_CONVEC,'
      IF (Master) WRITE (stdout,20) 'LMD_MIXING',                       &
     &   'Large/McWilliams/Doney interior mixing.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' LMD_MIXING,'
      ivmix=ivmix+1
      IF (Master) WRITE (stdout,20) 'LMD_NONLOCAL',                     &
     &   'LMD convective nonlocal transport.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' LMD_NONLOCAL,'
      IF (Master) WRITE (stdout,20) 'LMD_RIMIX',                        &
     &   'LMD diffusivity due to shear instability.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' LMD_RIMIX,'
      IF (Master) WRITE (stdout,20) 'LMD_SKPP',                         &
     &   'KPP surface boundary layer mixing.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' LMD_SKPP,'
      IF (Master) WRITE (stdout,20) 'LONGWAVE',                         &
     &   'Compute net longwave radiation internally.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' LONGWAVE,'
      IF (Master) WRITE (stdout,20) 'MASKING',                          &
     &   'Land/Sea masking.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' MASKING,'
      IF (Master) WRITE (stdout,20) 'MIX_GEO_TS',                       &
     &   'Mixing of tracers along geopotential surfaces.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' MIX_GEO_TS,'
      IF (Master) WRITE (stdout,20) 'MIX_S_UV',                         &
     &   'Mixing of momentum along constant S-surfaces.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' MIX_S_UV,'
      IF (Master) WRITE (stdout,20) 'MPI',                              &
     &   'MPI distributed-memory configuration.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+4)=' MPI,'
      IF (Master) WRITE (stdout,20) 'NONLINEAR',                        &
     &   'Nonlinear Model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' NONLINEAR,'
      IF (Master) WRITE (stdout,20) 'NONLIN_EOS',                       &
     &   'Nonlinear Equation of State for seawater.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' NONLIN_EOS,'
      IF (Master) WRITE (stdout,20) 'POWER_LAW',                        &
     &   'Power-law shape time-averaging barotropic filter.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' POWER_LAW,'
      IF (Master) WRITE (stdout,20) 'PROFILE',                          &
     &   'Time profiling activated .'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' PROFILE,'
      IF (Master) WRITE (stdout,20) '!RST_SINGLE',                      &
     &   'Double precision fields in restart NetCDF file.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' !RST_SINGLE,'
      IF (Master) WRITE (stdout,20) 'SALINITY',                         &
     &   'Using salinity.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' SALINITY,'
      IF (Master) WRITE (stdout,20) 'SOLAR_SOURCE',                     &
     &   'Solar Radiation Source Term.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' SOLAR_SOURCE,'
      IF (Master) WRITE (stdout,20) 'SOLVE3D',                          &
     &   'Solving 3D Primitive Equations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' SOLVE3D,'
      IF (Master) WRITE (stdout,20) 'SPLINES',                          &
     &   'Conservative parabolic spline reconstruction.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' SPLINES,'
      IF (Master) WRITE (stdout,20) 'STATIONS',                         &
     &   'Writing out station data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' STATIONS,'
      IF (Master) WRITE (stdout,20) 'TCLIMATOLOGY',                     &
     &   'Processing tracer climatology data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' TCLIMATOLOGY,'
      IF (Master) WRITE (stdout,20) 'TCLM_NUDGING',                     &
     &   'Nudging toward tracer climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' TCLM_NUDGING,'
      IF (Master) WRITE (stdout,20) 'T_PASSIVE',                        &
     &   'Advecting and diffusing inert passive tracer.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' T_PASSIVE,'
      IF (Master) WRITE (stdout,20) 'TS_U3HADVECTION',                  &
     &   'Third-order upstream horizontal advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TS_U3HADVECTION,'
      itrcHadv=itrcHadv+1
      IF (Master) WRITE (stdout,20) 'TS_C4VADVECTION',                  &
     &   'Fourth-order centered vertical advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TS_C4VADVECTION,'
      itrcVadv=itrcVadv+1
      IF (Master) WRITE (stdout,20) 'TS_DIF2',                          &
     &   'Harmonic mixing of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' TS_DIF2,'
      IF (Master) WRITE (stdout,20) 'UV_ADV',                           &
     &   'Advection of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' UV_ADV,'
      IF (Master) WRITE (stdout,20) 'UV_COR',                           &
     &   'Coriolis term.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' UV_COR,'
      IF (Master) WRITE (stdout,20) 'UV_U3HADVECTION',                  &
     &   'Third-order upstream horizontal advection of 3D momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' UV_U3HADVECTION,'
      ivelHadv=ivelHadv+1
      IF (Master) WRITE (stdout,20) 'UV_C4VADVECTION',                  &
     &   'Fourth-order centered vertical advection of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' UV_C4VADVECTION,'
      ivelVadv=ivelVadv+1
      IF (Master) WRITE (stdout,20) 'UV_QDRAG',                         &
     &   'Quadratic bottom stress.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' UV_QDRAG,'
      ibbl=ibbl+1
      IF (Master) WRITE (stdout,20) 'UV_VIS2',                          &
     &   'Harmonic mixing of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' UV_VIS2,'
      IF (Master) WRITE (stdout,20) 'VAR_RHO_2D',                       &
     &   'Variable density barotropic mode.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' VAR_RHO_2D,'
!
!-----------------------------------------------------------------------
!  Stop if unsupported C-preprocessing options or report issues with
!  particular options.
!-----------------------------------------------------------------------
!
      CALL checkadj
!
!-----------------------------------------------------------------------
!  Check C-preprocessing options.
!-----------------------------------------------------------------------
!
!  Stop if more than one vertical closure scheme is selected.
!
      IF (Master.and.(ivmix.gt.1)) THEN
        WRITE (stdout,30)
  30    FORMAT (/,' CHECKDEFS - only one vertical closure scheme',      &
     &            ' is allowed.')
        exit_flag=5
      END IF
!
!  Stop if more that one bottom stress formulation is selected.
!
      IF (Master.and.(ibbl.gt.1)) THEN
        WRITE (stdout,40)
  40    FORMAT (/,' CHECKDEFS - only one bottom stress formulation is', &
     &            ' allowed.')
        exit_flag=5
      END IF
!
!  Stop if no bottom stress formulation is selected.
!
      IF (Master.and.(ibbl.eq.0)) THEN
        WRITE (stdout,50)
  50    FORMAT (/,' CHECKDEFS - no bottom stress formulation is',       &
     &            ' selected.')
        exit_flag=5
      END IF
!
!  Stop if more than one biological module is selected.
!
      IF (Master.and.(ibiology.gt.1)) THEN
        WRITE (stdout,60)
  60    FORMAT (/,' CHECKDEFS - only one biology MODULE is allowed.')
        exit_flag=5
      END IF
!
!  Stop if more that one model driver is selected.
!
      IF (Master.and.(idriver.gt.1)) THEN
        WRITE (stdout,70)
  70    FORMAT (/,' CHECKDEFS - only one model example is allowed.')
        exit_flag=5
      END IF
!
!  Stop if more than one advection scheme has been activated.
!
      IF (Master.and.(ivelHadv.gt.1)) THEN
        WRITE (stdout,140) 'horizontal','momentum','ivelHadv =',ivelHadv
        exit_flag=5
      END IF
      IF (Master.and.(ivelVadv.gt.1)) THEN
        WRITE (stdout,140) 'vertical','momentum','ivelVadv =',ivelVadv
        exit_flag=5
      END IF
      IF (Master.and.(itrcHadv.gt.1)) THEN
        WRITE (stdout,140) 'horizontal','tracers','itrcHadv =',itrcHadv
        exit_flag=5
      END IF
      IF (Master.and.(itrcVadv.gt.1)) THEN
        WRITE (stdout,140) 'vertical','tracers','itrcVadv =',itrcVadv
        exit_flag=5
      END IF
 140  FORMAT (/,' CHECKDEFS - only one ',a,' advection scheme',         &
     &        /,13x,'is allowed for ',a,', ',a,1x,i1)
      RETURN
      END SUBROUTINE checkdefs
