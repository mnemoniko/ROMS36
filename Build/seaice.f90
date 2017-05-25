      SUBROUTINE seaice
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This is the main driver routine for the sea ice portion of the
!  model.
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
      USE mod_stepping
      USE mod_ice
      USE mod_forces
      USE ice_spdiw_mod, ONLY : ice_spdiw
      USE ice_vbc_mod, ONLY : ice_vbc
      USE ice_thermo_mod, ONLY : ice_thermo
      USE ice_evp_mod, ONLY : ice_evp
      USE ice_evp_sig_mod, ONLY : ice_evp_sig
      USE ice_elastic_mod, ONLY : ice_elastic
      USE ice_advect_mod, ONLY : ice_advect
      USE ice_enthalpi_mod, ONLY : ice_enthalpi
      USE ice_limit_mod, ONLY : ice_limit
      implicit none
      integer :: thread, subs, tile
      integer :: i, nforc, ng, my_ievp, nelas, iter
      real(r8), parameter :: dt_large = 1.0E+23_r8
      real(r8) :: dtice_sav
      DO ng=1,Ngrids
        liold(ng) = linew(ng)
        linew(ng) = 3-liold(ng)
! ----------------------------------------------------------------------
!  Compute the ice-ocean shear.
! ----------------------------------------------------------------------
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ice_spdiw(ng, TILE)
        END DO
      END DO
! ----------------------------------------------------------------------
!  Compute the stresses on the ice from the air and water.
! ----------------------------------------------------------------------
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ice_vbc(ng, TILE)
        END DO
      END DO
! ----------------------------------------------------------------------
!  Compute the internal ice stresses according to the 
!  Elastic-Viscous-Plastic rheology (EVP).
! ----------------------------------------------------------------------
      DO ng=1,Ngrids
        nelas = nevp(ng)
        liuol(ng) = liunw(ng)
        liunw(ng) = 3-liuol(ng)
        dte(ng) = dtice(ng)/FLOAT(nevp(ng))
        DO my_ievp=1,nelas
          lieol(ng) = lienw(ng)
          lienw(ng) = 3-lieol(ng)
          ievp(ng)=my_ievp
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ice_evp(ng, TILE)
          END DO
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ice_evp_sig(ng, TILE)
          END DO
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ice_elastic(ng, TILE)
          END DO
        END DO
      END DO
! ----------------------------------------------------------------------
!  Compute the ice enthalpi before advection.
! ----------------------------------------------------------------------
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ice_enthalpi(ng, TILE)
        END DO
      END DO
! ----------------------------------------------------------------------
!  Compute the advection of the ice tracer fields.
! ----------------------------------------------------------------------
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ice_advect(ng, TILE)
          CALL ice_limit(ng, TILE)
        END DO
      END DO
! ----------------------------------------------------------------------
!  Compute the ice thermodynamics.
! ----------------------------------------------------------------------
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ice_thermo (ng, TILE)
!         CALL ice_limit(ng, TILE)
        END DO
      END DO
      RETURN
      END SUBROUTINE seaice
