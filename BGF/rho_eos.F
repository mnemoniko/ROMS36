#include "cppdefs.h"
      MODULE rho_eos_mod
#ifdef SOLVE3D
!
!svn $Id: rho_eos.F 443 2010-02-04 21:31:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes  "in situ" density and other associated       !
!  quantitites as a function of potential temperature,  salinity,      !
!  and pressure using two different polynomial expressions:            !
!                                                                      !  
!  The default algorithms are based on Jackett & McDougall (1992)      !
!  The polynomial expression was found from fitting to 248 values      !
!  in the oceanographic ranges of salinity, potential temperature      !
!  and  pressure. It  assumes  no  pressure  variation  along          !
!  geopotential surfaces, that is, depth (meters; negative) and        !
!  pressure (dbar; assumed negative here) are interchangeable.         !
!                                                                      !
!  Check Values: (T=3 C, S=35.5 PSU, Z=-5000 m)                        !
!                                                                      !
!     alpha = 2.1014611551470d-04 (1/Celsius)                          !
!     beta  = 7.2575037309946d-04 (1/PSU)                              !
!     gamma = 3.9684764511766d-06 (1/Pa)                               !
!     den   = 1050.3639165364     (kg/m3)                              !
!     den1  = 1028.2845117925     (kg/m3)                              !
!     sound = 1548.8815240223     (m/s)                                !
!     bulk  = 23786.056026320     (Pa)                                 !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Jackett, D. R. and T. J. McDougall, 1995, Minimal Adjustment of     !
!    Hydrostatic Profiles to Achieve Static Stability, J. of Atmos.    !
!    and Oceanic Techn., vol. 12, pp. 381-389.                         !
!                                                                      !
!  For JM2006 EOS:                                                     !
!  Algorithms have been updated to the latest polynomials supplied     !
!  in Jackett et al. (2003) and McDougall et al. (2003), by            !
!  Benjamin K. Galton-Fenzi (Ben.Galton-Fenzi@utas.edu.au). Pressure   !
!  is converted to be positive meters.                                 !
!                                                                      !
!  Check Values: (t1=20 C, s1=20, p1=1000 m)                           !
!     den   =  1017.728868019642  (kg/m3)                              !   
!     rho_s =  0.7510471164699279 (kg/m3/PSU)                          !
!     rho_t = -0.2570255211349140 (kg/m3/C)                            !
!                                                                      !
!  Jackett, D. R., T. J. McDougall, R. Feistel, D. G. Wright and S. M. !
!    Griffies, 2006, Algorithms for Density, Potential Temperature,    !
!    Conservative Temperature and Freezing Temperature of Seawater,    !
!    J. of Atmos. and Oceanic Techn., vol. 23, pp. 1709-1728.          !
!    and Oceanic Techn., vol. 12, pp. 381-389.                         !
!                                                                      !
!  McDougall, T. J., D. R. Jackett, D. G. Wright and R. Feistel, 2003, !
!    Accurate and Computationally Efficient Algorithms for Potential   !
!    Temperature and Density of Seawater, vol. 20, pp. 730-741.        !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: rho_eos

      CONTAINS
!
!***********************************************************************
      SUBROUTINE rho_eos (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 14)
# endif
      CALL rho_eos_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng),                                      &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
# ifdef VAR_RHO_2D
     &                   GRID(ng) % Hz,                                 &
# endif
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   OCEAN(ng) % t,                                 &
# ifdef VAR_RHO_2D
     &                   COUPLING(ng) % rhoA,                           &
     &                   COUPLING(ng) % rhoS,                           &
# endif
# ifdef BV_FREQUENCY
     &                   MIXING(ng) % bvf,                              &
# endif
# if defined LMD_SKPP    || defined LMD_BKPP         || \
     defined BULK_FLUXES || defined BALANCE_OPERATOR
     &                   MIXING(ng) % alpha,                            &
     &                   MIXING(ng) % beta,                             &
#  ifdef LMD_DDMIX
     &                   MIXING(ng) % alfaobeta,                        &
#  endif
# endif
# if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
     &                   OCEAN(ng) % pden,                              &
# endif
     &                   OCEAN(ng) % rho)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 14)
# endif
      RETURN
      END SUBROUTINE rho_eos
# ifdef NONLIN_EOS
!
!***********************************************************************
      SUBROUTINE rho_eos_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,                                    &
#  ifdef MASKING
     &                         rmask,                                   &
#  endif
#  ifdef VAR_RHO_2D
     &                         Hz,                                      &
#  endif
     &                         z_r, z_w, t,                             &
#  ifdef VAR_RHO_2D
     &                         rhoA, rhoS,                              &
#  endif
#  ifdef BV_FREQUENCY
     &                         bvf,                                     &
#  endif
#  if defined LMD_SKPP    || defined LMD_BKPP         || \
      defined BULK_FLUXES || defined BALANCE_OPERATOR
     &                         alpha, beta,                             &
#   ifdef LMD_DDMIX
     &                         alfaobeta,                               &
#   endif
#  endif
#  if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
     &                         pden,                                    &
#  endif
     &                         rho)
!***********************************************************************
!
      USE mod_param
      USE mod_eoscoef
      USE mod_scalars
#  ifdef SEDIMENT
      USE mod_sediment
#  endif
!
#  if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
      USE exchange_3d_mod
#  endif
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#   endif
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
#   ifdef VAR_RHO_2D
      real(r8), intent(out) :: rhoA(LBi:,LBj:)
      real(r8), intent(out) :: rhoS(LBi:,LBj:)
#   endif
#   ifdef BV_FREQUENCY
      real(r8), intent(out) :: bvf(LBi:,LBj:,0:)
#   endif
#   if defined LMD_SKPP    || defined LMD_BKPP         || \
       defined BULK_FLUXES || defined BALANCE_OPERATOR
      real(r8), intent(out) :: alpha(LBi:,LBj:)
      real(r8), intent(out) :: beta(LBi:,LBj:)
#    ifdef LMD_DDMIX
      real(r8), intent(out) :: alfaobeta(LBi:,LBj:,0:)
#    endif
#   endif
#   if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
      real(r8), intent(out) :: pden(LBi:,LBj:,:)
#   endif
      real(r8), intent(out) :: rho(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#   ifdef VAR_RHO_2D
      real(r8), intent(out) :: rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rhoS(LBi:UBi,LBj:UBj)
#   endif
#   ifdef BV_FREQUENCY
      real(r8), intent(out) :: bvf(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#   if defined LMD_SKPP    || defined LMD_BKPP         || \
       defined BULK_FLUXES || defined BALANCE_OPERATOR
      real(r8), intent(out) :: alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: beta(LBi:UBi,LBj:UBj)
#    ifdef LMD_DDMIX
      real(r8), intent(out) :: alfaobeta(LBi:UBi,LBj:UBj,0:N(ng))
#    endif
#   endif
#   if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
      real(r8), intent(out) :: pden(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(out) :: rho(LBi:UBi,LBj:UBj,N(ng))
#  endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: i, ised, itrc, j, k
      real(r8) :: Tp
      real(r8), dimension(IminS:ImaxS,N(ng)) :: SedDen 
# ifdef JM2006_EOS
! standard atmospheric pressure (dbar).
! Should be set to 0.0 if assume zero pressure for the
! overlying atmosphere.  But if have a realistic atmospheric
! pressure loading, then set press_standard=10.1325.
      real(r8) :: press_standard = 0.0_r8 !10.1325_r8
      real(r8) :: eps = 1.0E-16_r8
      real(r8) :: num1b, num2b, num1, num2
      real(r8) :: p1, p1t1
# else
      real(r8) :: Tpr10, Ts, Tt, sqrtTs
# endif
#  ifdef BV_FREQUENCY
#   ifdef JM2006_EOS
      real(r8) :: den_dn, den_up, denr_up,denr_dn
#   else
      real(r8) :: bulk_dn, bulk_up, den_dn, den_up
#   endif
#  endif
#  ifdef  JM2006_EOS
#   ifdef EOS_TDERIVATIVE
      real(r8) :: anum_s, aden_s, anum_t, aden_t
      real(r8) :: rho_s, rho_t, den1000
#   endif
      real(r8), dimension(IminS:ImaxS,N(ng)) :: t1
      real(r8), dimension(IminS:ImaxS,N(ng)) :: t2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: s1
      real(r8), dimension(IminS:ImaxS,N(ng)) :: sp5
      real(r8), dimension(IminS:ImaxS,N(ng)) :: num1a
      real(r8), dimension(IminS:ImaxS,N(ng)) :: num2a
      real(r8), dimension(IminS:ImaxS,N(ng)) :: denr
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Scof
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Tcof
#  else       
      real(r8) :: cff
      real(r8), dimension(0:9) :: C
#   ifdef EOS_TDERIVATIVE
      real(r8), dimension(0:9) :: dCdT(0:9)
      real(r8), dimension(IminS:ImaxS,N(ng)) :: DbulkDS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: DbulkDT
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Dden1DS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Dden1DT
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Scof
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Tcof
      real(r8), dimension(IminS:ImaxS,N(ng)) :: wrk
#   endif
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk0
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk1
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: den1
#  endif
      real(r8), dimension(IminS:ImaxS,N(ng)) :: den
      real(r8) :: cff1, cff2

#  include "set_bounds.h"
!
!=======================================================================
!  Nonlinear equation of state.Notice that the default equation of state
!  is only valid for potential temperature range of -2C to 40C and
!  a salinity range of 0 PSU to 42 PSU.
!  The JM2006 equation of state is valid for potential temperature range
!  of -10C to 50C,  a salinity range of 0 PSU to 50 PSU  and  a pressure
!  range of 0 dbar to 10E4 dbar.
!=======================================================================
!
      DO j=JstrR,JendR
        DO k=1,N(ng)
          DO i=IstrR,IendR
!
!  Check temperature and salinity lower values. Assign depth to the
!  pressure.
! 
#  ifdef JM2006_EOS
            t1(i,k)=MAX(-10.0_r8,t(i,j,k,nrhs,itemp))
#  else
            Tt=MAX(-2.0_r8,t(i,j,k,nrhs,itemp))
#  endif

#  ifdef SALINITY
#   ifdef JM2006_EOS
            s1(i,k)=MAX(0.0_r8,t(i,j,k,nrhs,isalt))
#   else
            Ts=MAX(0.0_r8,t(i,j,k,nrhs,isalt))
            sqrtTs=SQRT(Ts)
#   endif
#  else
#   ifdef JM2006_EOS
           s1(i,k)=0.0_r8
#   else
            Ts=0.0_r8
            sqrtTs=0.0_r8
#   endif
#  endif
#  ifdef JM2006_EOS
            Tp=-z_r(i,j,k)
#  else
            Tp=z_r(i,j,k)
            Tpr10=0.1_r8*Tp
#  endif

!
!-----------------------------------------------------------------------
!  Compute density (kg/m3) at standard one atmosphere pressure.
!-----------------------------------------------------------------------
!
#  ifdef JM2006_EOS
            t2(i,k)  = t1(i,k)*t1(i,k)
            sp5(i,k) = sqrt(s1(i,k))
                p1   = Tp - press_standard
                p1t1 = p1*t1(i,k)
          num1a(i,k) = a0 + t1(i,k)*(a1 + t1(i,k)*(a2+a3*t1(i,k)) )     &
     &               + s1(i,k)*(a4 + a5*t1(i,k)  + a6*s1(i,k))          

               num1b = p1*(a7 + a8*t2(i,k) + a9*s1(i,k)                 &
     &               + p1*(a10+a11*t2(i,k)))

          num2a(i,k) = b0 + t1(i,k)*(b1 + t1(i,k)*(b2                   &
     &               + t1(i,k)*(b3 + t1(i,k)*b4)))                      &
     &               + s1(i,k)*(b5 + t1(i,k)*(b6 + b7*t2(i,k))          &
     &               + sp5(i,k)*(b8 + b9*t2(i,k)))                     
               num2b = p1*(b10 + p1t1*(b11*t2(i,k) + b12*p1))
           denr(i,k) = 1.0/(eps+(num2a(i,k)+num2b))
            den(i,k) = (num1a(i,k)+num1b)*denr(i,k)

#  else
            C(0)=Q00+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))
            C(1)=U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))
            C(2)=V00+Tt*(V01+Tt*V02)
#   ifdef EOS_TDERIVATIVE
!
            dCdT(0)=Q01+Tt*(2.0_r8*Q02+Tt*(3.0_r8*Q03+Tt*(4.0_r8*Q04+   &
     &                      Tt*5.0_r8*Q05)))
            dCdT(1)=U01+Tt*(2.0_r8*U02+Tt*(3.0_r8*U03+Tt*4.0_r8*U04))
            dCdT(2)=V01+Tt*2.0_r8*V02
#   endif
!
            den1(i,k)=C(0)+Ts*(C(1)+sqrtTs*C(2)+Ts*W00)

#   ifdef EOS_TDERIVATIVE
!
!   Compute d(den1)/d(S) and d(den1)/d(T) derivatives used in the
!   computation of thermal expansion and saline contraction
!   coefficients.
!
            Dden1DS(i,k)=C(1)+1.5_r8*C(2)*sqrtTs+2.0_r8*W00*Ts
            Dden1DT(i,k)=dCdT(0)+Ts*(dCdT(1)+sqrtTs*dCdT(2))
#   endif
!
!-----------------------------------------------------------------------
!   Compute secant bulk modulus.
!-----------------------------------------------------------------------
!
            C(3)=A00+Tt*(A01+Tt*(A02+Tt*(A03+Tt*A04)))
            C(4)=B00+Tt*(B01+Tt*(B02+Tt*B03))
            C(5)=D00+Tt*(D01+Tt*D02)
            C(6)=E00+Tt*(E01+Tt*(E02+Tt*E03))
            C(7)=F00+Tt*(F01+Tt*F02)
            C(8)=G01+Tt*(G02+Tt*G03)
            C(9)=H00+Tt*(H01+Tt*H02)
#   ifdef EOS_TDERIVATIVE
!
            dCdT(3)=A01+Tt*(2.0_r8*A02+Tt*(3.0_r8*A03+Tt*4.0_r8*A04))
            dCdT(4)=B01+Tt*(2.0_r8*B02+Tt*3.0_r8*B03)
            dCdT(5)=D01+Tt*2.0_r8*D02
            dCdT(6)=E01+Tt*(2.0_r8*E02+Tt*3.0_r8*E03)
            dCdT(7)=F01+Tt*2.0_r8*F02
            dCdT(8)=G02+Tt*2.0_r8*G03
            dCdT(9)=H01+Tt*2.0_r8*H02
#   endif
!
            bulk0(i,k)=C(3)+Ts*(C(4)+sqrtTs*C(5))
            bulk1(i,k)=C(6)+Ts*(C(7)+sqrtTs*G00)
            bulk2(i,k)=C(8)+Ts*C(9)
            bulk (i,k)=bulk0(i,k)-Tp*(bulk1(i,k)-Tp*bulk2(i,k))

#   if defined LMD_SKPP    || defined LMD_BKPP         || \
       defined BULK_FLUXES || defined BALANCE_OPERATOR
!
!  Compute d(bulk)/d(S) and d(bulk)/d(T) derivatives used
!  in the computation of thermal expansion and saline contraction
!  coefficients.
!
            DbulkDS(i,k)=C(4)+sqrtTs*1.5_r8*C(5)-                       &
     &                   Tp*(C(7)+sqrtTs*1.5_r8*G00-Tp*C(9))
            DbulkDT(i,k)=dCdT(3)+Ts*(dCdT(4)+sqrtTs*dCdT(5))-           &
     &                   Tp*(dCdT(6)+Ts*dCdT(7)-                        &
     &                       Tp*(dCdT(8)+Ts*dCdT(9)))
#   endif
!
!-----------------------------------------------------------------------
!  Compute local "in situ" density anomaly (kg/m3 - 1000).
!-----------------------------------------------------------------------
!
            cff=1.0_r8/(bulk(i,k)+Tpr10)
            den(i,k)=den1(i,k)*bulk(i,k)*cff
#  endif

#  if defined SEDIMENT && defined SED_DENS
            SedDen(i,k)=0.0_r8
            DO ised=1,NST
              itrc=idsed(ised)
              cff1=1.0_r8/Srho(ised,ng)
              SedDen(i,k)=SedDen(i,k)+                                  &
     &               t(i,j,k,nrhs,itrc)*                                &
     &               (Srho(ised,ng)-den(i,k))*cff1
            END DO
            den(i,k)=den(i,k)+SedDen(i,k)
#  endif
            den(i,k)=den(i,k)-1000.0_r8
#  ifdef MASKING
            den(i,k)=den(i,k)*rmask(i,j)
#  endif
          END DO
        END DO

#  ifdef VAR_RHO_2D
!
!-----------------------------------------------------------------------
!  Compute vertical averaged density (rhoA) and density perturbation
!  (rhoS) used in barotropic pressure gradient.
!-----------------------------------------------------------------------
!
        DO i=IstrR,IendR
          cff1=den(i,N(ng))*Hz(i,j,N(ng))
          rhoS(i,j)=0.5_r8*cff1*Hz(i,j,N(ng))
          rhoA(i,j)=cff1
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrR,IendR
            cff1=den(i,k)*Hz(i,j,k)
            rhoS(i,j)=rhoS(i,j)+Hz(i,j,k)*(rhoA(i,j)+0.5_r8*cff1)
            rhoA(i,j)=rhoA(i,j)+cff1
          END DO
        END DO
        cff2=1.0_r8/rho0
        DO i=IstrR,IendR
          cff1=1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,0))
          rhoA(i,j)=cff2*cff1*rhoA(i,j)
          rhoS(i,j)=2.0_r8*cff1*cff1*cff2*rhoS(i,j)
        END DO
#  endif

#  if defined BV_FREQUENCY
!
!-----------------------------------------------------------------------
!  Compute Brunt-Vaisala frequency (1/s2) at horizontal RHO-points
!  and vertical W-points:
!
!                  bvf = - g/rho d(rho)/d(z).
!
!  The density anomaly difference is computed by lowering/rising the
!  water parcel above/below adiabatically at W-point depth "z_w".
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)-1
          DO i=IstrR,IendR
# ifdef JM2006_EOS
                p1   = -z_w(i,j,k) - press_standard
! Compute den_up
                p1t1 = p1*t1(i,k+1)

                num1 = num1a(i,k+1)                                     &
     &               + p1*(a7 + a8*t2(i,k+1) + a9*s1(i,k+1)             &
     &               + p1*(a10+a11*t2(i,k+1)))

                num2 = num2a(i,k+1)                                     &
     &               + p1*(b10 + p1t1*(b11*t2(i,k+1) + b12*p1))
               denr_up = 1.0/(eps+num2)
#  if defined SEDIMENT && defined SED_DENS
                den_up = num1*denr_up+SedDen(i,k+1)
#  else
                den_up = num1*denr_up
#  endif
! Compute den_dn
                p1t1 = p1*t1(i,k)

                num1 = num1a(i,k)                                       &
     &               + p1*(a7 + a8*t2(i,k) + a9*s1(i,k)                 &
     &               + p1*(a10+a11*t2(i,k)))

                num2 = num2a(i,k)                                       &
     &               + p1*(b10 + p1t1*(b11*t2(i,k) + b12*p1))
             denr_dn = 1.0/(eps+num2)
#  if defined SEDIMENT && defined SED_DENS
                den_dn = num1*denr_dn+SedDen(i,k)
#  else 
                den_dn = num1*denr_dn
#  endif
# else
            bulk_up=bulk0(i,k+1)-                                       &
     &              z_w(i,j,k)*(bulk1(i,k+1)-                           &
     &                          bulk2(i,k+1)*z_w(i,j,k))
            bulk_dn=bulk0(i,k  )-                                       &
     &              z_w(i,j,k)*(bulk1(i,k  )-                           &
     &                          bulk2(i,k  )*z_w(i,j,k))
            cff1=1.0_r8/(bulk_up+0.1_r8*z_w(i,j,k))
            cff2=1.0_r8/(bulk_dn+0.1_r8*z_w(i,j,k))
#  if defined SEDIMENT && defined SED_DENS
            den_up=cff1*(den1(i,k+1)*bulk_up)+SedDen(i,k+1)
            den_dn=cff2*(den1(i,k  )*bulk_dn)+SedDen(i,k)
#  else 
            den_up=cff1*(den1(i,k+1)*bulk_up)
            den_dn=cff2*(den1(i,k  )*bulk_dn)
#  endif
# endif
            bvf(i,j,k)=-g*(den_up-den_dn)/                              &
     &                 (0.5_r8*(den_up+den_dn)*                         &
     &                  (z_r(i,j,k+1)-z_r(i,j,k)))
          END DO
        END DO
        DO i=IstrR,IendR
!!        bvf(i,j,0)=bvf(i,j,1)
!!        bvf(i,j,N(ng))=bvf(i,j,N(ng)-1)
          bvf(i,j,0)=0.0_r8
          bvf(i,j,N(ng))=0.0_r8
        END DO
#  endif

#  if defined LMD_SKPP    || defined LMD_BKPP         || \
      defined BULK_FLUXES || defined BALANCE_OPERATOR
!
!-----------------------------------------------------------------------
!  Compute thermal expansion (1/Celsius) and saline contraction
!  (1/PSU) coefficients.
!-----------------------------------------------------------------------
!
#   ifdef LMD_DDMIX
        DO k=1,N(ng)
#   else
        DO k=N(ng),N(ng)
#   endif
          DO i=IstrR,IendR
!
!  Compute thermal expansion and saline contraction coefficients.
! 
#    ifdef JM2006_EOS
            p1 = -z_r(i,j,k) - press_standard

            den1000 = den(i,k)+1000.0_r8

            anum_s = a4 + t1(i,k)*a5 + s1(i,k)*two_a6 + p1*a9

            anum_t = a1 + t1(i,k)*(two_a2 + t1(i,k)*three_a3)           &
     &             + s1(i,k)*a5 + p1*t1(i,k)*(two_a8+p1*two_a11)

            aden_s = b5 + t1(i,k)*(b6 + t2(i,k)*b7)                     &
     &             + sp5(i,k)*(onep5_b8 + t2(i,k)*onep5_b9)

            aden_t = b1 + t1(i,k)*(two_b2 + t1(i,k)*(three_b3           &
     &             + t1(i,k)*four_b4))                                  &
     &             + s1(i,k)*(b6 + t2(i,k)*three_b7                     &
     &             + t1(i,k)*sp5(i,k)*two_b9)                           &
     &             - p1*p1*(-t2(i,k)*three_b11 - p1*b12)

            rho_s = (anum_s-aden_s*den1000)*denr(i,k)
            rho_t = (anum_t-aden_t*den1000)*denr(i,k)

            Scof(i,k)=rho_s/den1000
            Tcof(i,k)=-rho_t/den1000
#    else
            Tpr10=0.1_r8*z_r(i,j,k)
            cff=bulk(i,k)+Tpr10
            cff1=Tpr10*den1(i,k)
            cff2=bulk(i,k)*cff
            wrk(i,k)=(den(i,k)+1000.0_r8)*cff*cff
            Tcof(i,k)=-(DbulkDT(i,k)*cff1+                              &
     &                  Dden1DT(i,k)*cff2)
            Scof(i,k)= (DbulkDS(i,k)*cff1+                              &
     &                  Dden1DS(i,k)*cff2)
#    endif
#   ifdef LMD_DDMIX
            alfaobeta(i,j,k)=Tcof(i,k)/Scof(i,k)
#   endif
          END DO
          IF (k.eq.N(ng)) THEN
            DO i=IstrR,IendR
#   ifdef  JM2006_EOS
              alpha(i,j)=Tcof(i,N(ng))
              beta (i,j)=Scof(i,N(ng))
#   else
              cff=1.0_r8/wrk(i,N(ng))
              alpha(i,j)=cff*Tcof(i,N(ng))
              beta (i,j)=cff*Scof(i,N(ng))
#   endif
            END DO
          END IF
        END DO
#  endif
!
!-----------------------------------------------------------------------
!  Load "in situ" density anomaly (kg/m3 - 1000) and potential
!  density anomaly (kg/m3 - 1000) referenced to the surface into global
!  arrays. Notice that this is done in a separate (i,k) DO-loops to
!  facilitate the adjoint.
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)
          DO i=IstrR,IendR
            rho(i,j,k)=den(i,k)
#  if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
#   ifdef JM2006_EOS   
            pden(i,j,k)=(num1a(i,k)/(eps+num2a(i,k)))
#    if defined SEDIMENT && defined SED_DENS
            pden(i,j,k)=pden(i,j,k)+SedDen(i,k)
#    endif
            pden(i,j,k)=pden(i,j,k)-1000.0_r8
#   else
#    if defined SEDIMENT && defined SED_DENS
            pden(i,j,k)=pden(i,j,k)+SedDen(i,k)
#    endif
            pden(i,j,k)=den1(i,k)-1000.0_r8
#   endif
#   ifdef MASKING
            pden(i,j,k)=pden(i,j,k)*rmask(i,j)
#   endif
#  endif
          END DO
        END DO
      END DO
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
#   if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        rho)
#    if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        pden)
#    endif
#    if defined LMD_SKPP    || defined LMD_BKPP         || \
        defined BULK_FLUXES || defined BALANCE_OPERATOR

#     ifdef LMD_DDMIX
      CALL exchange_w3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        alfaobeta)
#     endif
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        alpha)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        beta)
#    endif
#    ifdef VAR_RHO_2D
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        rhoA)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        rhoS)
#    endif
#    ifdef BV_FREQUENCY
      CALL exchange_w3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        bvf)
#    endif
#   endif

#   ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    rho)
#    if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    pden)
#    endif
#    if defined LMD_SKPP    || defined LMD_BKPP         || \
        defined BULK_FLUXES || defined BALANCE_OPERATOR
#     ifdef LMD_DDMIX
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    alfaobeta)
#     endif
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    alpha, beta)
#    endif
#    ifdef VAR_RHO_2D
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    rhoA, rhoS)
#    endif
#    ifdef BV_FREQUENCY
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bvf)
#    endif
#   endif
#  endif

      RETURN
      END SUBROUTINE rho_eos_tile
# endif

# ifndef NONLIN_EOS
!
!***********************************************************************
      SUBROUTINE rho_eos_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,                                    &
#  ifdef MASKING
     &                         rmask,                                   &
#  endif
#  ifdef VAR_RHO_2D
     &                         Hz,                                      &
#  endif
     &                         z_r, z_w, t,                             &
#  ifdef VAR_RHO_2D
     &                         rhoA, rhoS,                              &
#  endif
#  ifdef BV_FREQUENCY
     &                         bvf,                                     &
#  endif
#  if defined LMD_SKPP    || defined LMD_BKPP         || \
      defined BULK_FLUXES || defined BALANCE_OPERATOR
     &                         alpha, beta,                             &
#   ifdef LMD_DDMIX
     &                         alfaobeta,                               &
#   endif
#  endif
#  if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
     &                         pden,                                    &
#  endif
     &                         rho)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
#  ifdef SEDIMENT
      USE mod_sediment
#  endif
!
#  if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
      USE exchange_3d_mod
#  endif
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#   endif
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)

#   ifdef VAR_RHO_2D
      real(r8), intent(out) :: rhoA(LBi:,LBj:)
      real(r8), intent(out) :: rhoS(LBi:,LBj:)
#   endif
#   ifdef BV_FREQUENCY
      real(r8), intent(out) :: bvf(LBi:,LBj:,0:)
#   endif
#   if defined LMD_SKPP    || defined LMD_BKPP         || \
       defined BULK_FLUXES || defined BALANCE_OPERATOR
      real(r8), intent(out) :: alpha(LBi:,LBj:)
      real(r8), intent(out) :: beta(LBi:,LBj:)
#    ifdef LMD_DDMIX
      real(r8), intent(out) :: alfaobeta(LBi:,LBj:,0:)
#    endif
#   endif
#   if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
      real(r8), intent(out) :: pden(LBi:,LBj:,:)
#   endif
      real(r8), intent(out) :: rho(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#   ifdef VAR_RHO_2D
      real(r8), intent(out) :: rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rhoS(LBi:UBi,LBj:UBj)
#   endif
#   ifdef BV_FREQUENCY
      real(r8), intent(out) :: bvf(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#   if defined LMD_SKPP    || defined LMD_BKPP         || \
       defined BULK_FLUXES || defined BALANCE_OPERATOR
      real(r8), intent(out) :: alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: beta(LBi:UBi,LBj:UBj)
#    ifdef LMD_DDMIX
      real(r8), intent(out) :: alfaobeta(LBi:UBi,LBj:UBj,0:N(ng))
#    endif
#   endif
#   if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
      real(r8), intent(out) :: pden(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(out) :: rho(LBi:UBi,LBj:UBj,N(ng))
#  endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: i, ised, itrc, j, k
      real(r8) :: SedDen, cff, cff1, cff2

#  include "set_bounds.h"
!
!=======================================================================
!  Linear equation of state.
!=======================================================================
!
!-----------------------------------------------------------------------
!  Compute "in situ" density anomaly (kg/m3 - 1000) using the linear
!  equation of state.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO k=1,N(ng)
          DO i=IstrR,IendR
            rho(i,j,k)=R0(ng)-                                          &
     &                 R0(ng)*Tcoef(ng)*(t(i,j,k,nrhs,itemp)-T0(ng))
#  ifdef SALINITY
            rho(i,j,k)=rho(i,j,k)+                                      &
     &                 R0(ng)*Scoef(ng)*(t(i,j,k,nrhs,isalt)-S0(ng))
#  endif
#  if defined SEDIMENT && defined SED_DENS
            SedDen=0.0_r8
            DO ised=1,NST
              itrc=idsed(ised)
              cff1=1.0_r8/Srho(ised,ng)
              SedDen=SedDen+                                            &
     &               t(i,j,k,nrhs,itrc)*                                &
     &               (Srho(ised,ng)-rho(i,j,k))*cff1
            END DO
            rho(i,j,k)=rho(i,j,k)+SedDen
#  endif
            rho(i,j,k)=rho(i,j,k)-1000.0_r8
#  ifdef MASKING
            rho(i,j,k)=rho(i,j,k)*rmask(i,j)
#  endif
#  if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
            pden(i,j,k)=rho(i,j,k)
#  endif
          END DO
        END DO

#  ifdef VAR_RHO_2D
!
!-----------------------------------------------------------------------
!  Compute vertical averaged density (rhoA) and density perturbation
!  used (rhoS) in barotropic pressure gradient.
!-----------------------------------------------------------------------
!
        DO i=IstrR,IendR
          cff1=rho(i,j,N(ng))*Hz(i,j,N(ng))
          rhoS(i,j)=0.5_r8*cff1*Hz(i,j,N(ng))
          rhoA(i,j)=cff1
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrR,IendR
            cff1=rho(i,j,k)*Hz(i,j,k)
            rhoS(i,j)=rhoS(i,j)+Hz(i,j,k)*(rhoA(i,j)+0.5_r8*cff1)
            rhoA(i,j)=rhoA(i,j)+cff1
          END DO
        END DO
        cff2=1.0_r8/rho0
        DO i=IstrR,IendR
          cff1=1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,0))
          rhoA(i,j)=cff2*cff1*rhoA(i,j)
          rhoS(i,j)=2.0_r8*cff1*cff1*cff2*rhoS(i,j)
        END DO
#  endif

#  ifdef BV_FREQUENCY
!
!-----------------------------------------------------------------------
!  Compute Brunt-Vaisala frequency (1/s2) at horizontal RHO-points
!  and vertical W-points.
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)-1
          DO i=IstrR,IendR
            bvf(i,j,k)=-gorho0*(rho(i,j,k+1)-rho(i,j,k))/               &
     &                         (z_r(i,j,k+1)-z_r(i,j,k))
          END DO
        END DO
#  endif

#  if defined LMD_SKPP    || defined LMD_BKPP         || \
      defined BULK_FLUXES || defined BALANCE_OPERATOR
!
!-----------------------------------------------------------------------
!  Compute thermal expansion (1/Celsius) and saline contraction
!  (1/PSU) coefficients.
!-----------------------------------------------------------------------
!
        DO i=IstrR,IendR
          alpha(i,j)=ABS(Tcoef(ng))
#   ifdef SALINITY
          beta(i,j)=ABS(Scoef(ng))
#   else
          beta(i,j)=0.0_r8
#   endif
        END DO
#   ifdef LMD_DDMIX
!
!  Compute ratio of thermal expansion and saline contraction
!  coefficients.
!
        IF (Scoef(ng).eq.0.0_r8) THEN
          cff=1.0_r8
        ELSE
          cff=1.0_r8/Scoef(ng)
        END IF
        DO k=1,N(ng)
          DO i=IstrR,IendR
            alfaobeta(i,j,k)=cff*Tcoef(ng)
          END DO
        END DO
#   endif
#  endif
      END DO
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
#   if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        rho)
#    if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        pden)
#    endif
#    if defined LMD_SKPP    || defined LMD_BKPP         || \
        defined BULK_FLUXES || defined BALANCE_OPERATOR
#     ifdef LMD_DDMIX
      CALL exchange_w3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        alfaobeta)
#     endif
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        alpha)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        beta)
#    endif
#    ifdef VAR_RHO_2D
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        rhoA)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        rhoS)
#    endif
#    ifdef BV_FREQUENCY
      CALL exchange_w3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        bvf)
#    endif
#   endif

#   ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    rho)
#    if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    pden)
#    endif
#    if defined LMD_SKPP    || defined LMD_BKPP         || \
        defined BULK_FLUXES || defined BALANCE_OPERATOR
#     ifdef LMD_DDMIX
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    alfaobeta)
#     endif
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    alpha, beta)
#    endif
#    ifdef VAR_RHO_2D
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    rhoA, rhoS)
#    endif
#    ifdef BV_FREQUENCY
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bvf)
#    endif
#   endif
#  endif
      RETURN
      END SUBROUTINE rho_eos_tile
# endif
#endif
      END MODULE rho_eos_mod

