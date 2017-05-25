      MODULE ncep_flux_mod
      END MODULE ncep_flux_mod
! -----------------------------------------------------------------
      SUBROUTINE bulktf(du,zu,ta,zt,qa,zq,ts,qs,cice,cd,ch,ce,wg2,      &
     &                  cd_i,cd_w)
!
! --- this routine computes momentum, sensible heat, and latent heat
! --- transfer coefficients by one interation of a bulk flux algorithm
! --- (Fairall et al.,1996)
!
      USE mod_kinds
      implicit none
! --- input variables:
! ---   du     - velocity difference (wind speed - current) [m/s]
! ---   zu     - measurement height of wind speed [m]
! ---   ta     - air temperature [K]
! ---   zt     - measurement height of air temperature [m]
! ---   qa     - specific humidity of air [kg/kg]
! ---   zq     - measurement height of specific humidity [m]
! ---   ts     - sea surface temperature [K]
! ---   qs     - specific humidity of air at the surface [kg/kg]
! ---   cice   - ice concentration (0 open water, 1 sea ice cover)
! ---   cd     - momentum transfer coeffisient
! ---   ch     - sensible heat transfer coeffisient
! ---   ce     - latent heat transfer coeffisient
! ---   wg2    - gustiness squared [m^2/s^2]
!
      real(r8), intent(in) :: du
      real(r8), intent(in) :: zu
      real(r8), intent(in) :: ta
      real(r8), intent(in) :: zt
      real(r8), intent(in) :: qa
      real(r8), intent(in) :: zq
      real(r8), intent(in) :: ts
      real(r8), intent(in) :: qs
      real(r8), intent(in) :: cice
      real(r8), intent(inout) :: cd
      real(r8), intent(inout) :: ch
      real(r8), intent(inout) :: ce
      real(r8), intent(inout) :: wg2
      real(r8), intent(out) :: cd_i
      real(r8), intent(out) :: cd_w
!
! --- parameters:
! ---   eps    - molecular weight ratio of dry air and water vapour
! ---   t0     - freezing point of water [K]
! ---   zi     - inversion height [m]
! ---   g      - acceleration due to gravity [m/s^2]
! ---   beta   - empirical constant relating gustiness to convective
! ---            scaling velocity
! ---   alpha  - Charnock constant
! ---   k      - von Karman constant
!
      real(r8), parameter :: eps=.62197_r8, cv=1._r8/eps-1._r8
      real(r8), parameter :: t0=273.16_r8
      real(r8), parameter :: zi=600._r8, g=9.8_r8, beta=1.2_r8
      real(r8), parameter :: athird=1._r8/3._r8
      real(r8), parameter :: alpha=.011_r8, gi=1._r8/g
      real(r8), parameter :: k=.4_r8, ki=2.5_r8
!
      real(r8) :: tv, tac, visca, dt, dq, du1, du2, s, ustar2, ustar
      real(r8) :: fac, tstar, qstar
      real(r8) :: tvstar, li, w3, wg, zetau, zetat, zetaq, z0, cd2
      real(r8) :: psiu, reu, ret, req
      real(r8) :: z0t, z0q, ct2, psitq, cq2, zice, zwat
      real(r8) :: cd2i, cd2w
!
! --- virtual temperature
      tv=ta*(1.+cv*qa)
!
! --- kinematic viscosity of dry air - Andreas (1989) CRREL Rep. 89-11
! --- [m^2/s]
      tac=ta-t0
      visca=1.326E-5_r8*(1._r8+tac*(6.542E-3_r8+tac*                    &
     &                            (8.301E-6_r8-tac*4.84E-9_r8)))
!
! --- temperature difference
      dt=ta-ts+.0098_r8*zt
!
! --- humidity difference
      dq=qa-qs
!
! --- initial values
      du1=max(du,1.E-2_r8)
      du2=du1*du1
      s=SQRT(du2+wg2)
      ustar2=cd*s*du1
      ustar=SQRT(ustar2)
      fac=ustar/(cd*du1)
      tstar=fac*ch*dt
      qstar=fac*ce*dq
!
! --- inverse Monin-Obukov lenght
      tvstar=tstar*(1+cv*qa)+cv*ta*qstar
      li=min(3._r8/zu,g*k*tvstar/(ustar2*tv))
!
! --- gustiness
      w3=-zi*g*ustar*tvstar/ta
      wg=max(.1_r8,beta*max(0._r8,w3)**athird)
!
! --- wind speed included gustiness
      s=sqrt(du2+wg*wg)
!
! --- nondimensional heights
      zetau=zu*li
      zetat=zt*li
      zetaq=zq*li
!
! --- roughness length - choose roughness lenght for sea ice
! --- corresponding to smooth multiyear ice, Guest and Davidson (1991)
! --- JGR 96, 4709-4721
      z0=cice*2.E-3_r8+(1._r8-cice)*                                    &
     &   (0.11_r8*visca/ustar+alpha*ustar2*gi)
!
! --- square root of the momentum transfer coefficient
      cd2=k/max(7._r8,LOG(zu/z0)-psiu(zetau))
!
! --- update friction velocity
      ustar=cd2*SQRT(s*du1)
!
! --- roughness Reynolds numbers
      reu=ustar*z0/visca
      CALL lkb(reu,ret,req)
!
! --- temperature and humidity roughness scales
      fac=visca/ustar
      z0t=fac*ret
      z0q=fac*req
!
! --- transfer coeffisient components for temperature and humidity
      ct2=k/max(7._r8,LOG(zt/z0t)-psitq(zetat))
      cq2=k/max(7._r8,LOG(zq/z0q)-psitq(zetaq))
!
! --- momentum transfer coeffisient
      cd=cd2*cd2
!
! --- heat transfer coeffisients
      ch=cd2*ct2
      ce=cd2*cq2
!
! --- gustiness squared
      wg2=wg*wg
! Get drag coefficients for open water, sea ice cover
      zwat=(0.11_r8*visca/ustar+alpha*ustar2*gi)
      zice=2.E-3_r8
      cd2i=k/MAX(7._r8,LOG(zu/zice)-psiu(zetau))
      cd2w=k/MAX(7._r8,LOG(zu/zwat)-psiu(zetau))
      cd_i = cd2i*cd2i
      cd_w = cd2w*cd2w
      RETURN
      END SUBROUTINE bulktf
!
! --- ------------------------------------------------------------------
!
      SUBROUTINE lkb(reu, ret, req)
!
! --- computes the roughness Reynolds for temperature and humidity
! --- following Liu, Katsaros and Businger (1979), J. Atmos.  Sci., 36,
! --- 1722-1735.
!
      USE mod_kinds
      real(r8), intent(in) :: reu
      real(r8), intent(out) :: ret
      real(r8), intent(out) :: req
!
      integer :: i
!
      real(r8), dimension(8,2) ::  a
      real(r8), dimension(8) :: re
      real(r8), dimension(8,2) ::  b
      data a /                                                          &
     &       0.177_r8, 1.376_r8, 1.026_r8, 1.625_r8,                    &
     &       4.661_r8, 34.904_r8, 1667.19_r8, 5.88E5_r8,                &
     &       0.292_r8, 1.808_r8, 1.393_r8, 1.956_r8,                    &
     &       4.994_r8, 30.709_r8, 1448.68_r8, 2.98E5_r8 /
      data b                                                            &
     &      /0._r8, 0.929_r8, -0.599_r8, -1.018_r8,                     &
     &      -1.475_r8, -2.067_r8, -2.907_r8, -3.935_r8,                 &
     &       0._r8, 0.826_r8, -0.528_r8, -0.870_r8,                     &
     &      -1.297_r8, -1.845_r8, -2.682_r8, -3.616_r8/
      data re                                                           &
     &      /0.11_r8, 0.825_r8, 3.0_r8, 10.0_r8,                        &
     &       30.0_r8, 100._r8, 300._r8, 1000._r8/
!
      i=1
 10   IF (reu.gt.re(i).and.i.lt.8) THEN
        i=i+1
        GOTO 10
      ENDIF
!
      ret=a(i,1)*reu**b(i,1)
      req=a(i,2)*reu**b(i,2)
!
      RETURN
      END SUBROUTINE lkb
!
! --- ------------------------------------------------------------------
!
      FUNCTION psiu(zeta)
!
! --- Monin-Obukhov similarity velocity profile function
!
      USE mod_kinds
      real(r8) :: psiu, zeta
!
      real(r8), parameter :: athird=1._r8/3._r8
      real(r8), parameter :: pi=3.141592653589793_r8
      real(r8), parameter :: sqrt3=1.732050807568877_r8
      real(r8), parameter :: sqrt3i=.5773502691896258_r8
!
      real(r8) :: x, psik, y, psic, f
!
      IF (zeta.eq.0._r8) THEN
        psiu=0.
      ELSEIF (zeta.gt.0._r8) THEN
        psiu=-4.7_r8*zeta
      ELSE
        x=(1._r8-16._r8*zeta)**.25_r8
        psik=2._r8*LOG((1._r8+x)*.5_r8)+LOG((1._r8+x*x)*.5_r8)          &
     &       -2._r8*ATAN(x)+pi*.5_r8
        y=(1.-12.87*zeta)**athird
        psic=1.5_r8*LOG((y*y+y+1._r8)*athird)                           &
     &      -sqrt3*ATAN((2._r8*y+1._r8)*sqrt3i)                         &
     &      +pi*sqrt3i
        f=1._r8/(1._r8+zeta*zeta)
        psiu=f*psik+(1._r8-f)*psic
      ENDIF
!
      RETURN
      END FUNCTION psiu
!
! --- ------------------------------------------------------------------
!
      FUNCTION psitq(zeta)
!
! --- Monin-Obukhov similarity temperature and humidity profile function
!
      USE mod_kinds
      real(r8) :: psitq, zeta
!
      real(r8), parameter :: athird=1._r8/3._r8
      real(r8), parameter :: pi=3.141592653589793_r8
      real(r8), parameter :: sqrt3=1.732050807568877_r8
      real(r8), parameter :: sqrt3i=.5773502691896258_r8
!
      real(r8) :: x, psik, y, psic, f
!
      IF (zeta.eq.0._r8) THEN
        psitq=0._r8
      ELSEIF (zeta.gt.0._r8) THEN
        psitq=-4.7_r8*zeta
      ELSE
        x=(1._r8-16._r8*zeta)**.25_r8
        psik=2._r8*LOG((1._r8+x*x)*.5_r8)
        y=(1._r8-12.87_r8*zeta)**athird
        psic=1.5_r8*LOG((y*y+y+1._r8)*athird)                           &
     &             -sqrt3*ATAN((2._r8*y+1._r8)*sqrt3i)                  &
     &      +pi*sqrt3i
        f=1./(1._r8+zeta*zeta)
        psitq=f*psik+(1._r8-f)*psic
      ENDIF
!
      RETURN
      END FUNCTION psitq
      FUNCTION qsatw(t,p)
!
! --- computes saturation specific humidity [kg/kg] over open water,
! --- Buck (1981) JAM 20, 1527-1532
!
! --- input variables:
! ---   t      - air temperature [K]
! ---   p      - atmospheric pressure [Pa]
!
      USE mod_kinds
      real(r8) :: qsatw, t, p
!
! --- parameters:
! ---   eps    - molecular weight ratio of dry air and water vapour
!
      real(r8), parameter :: eps=0.62197_r8
!
      real(r8) :: e
!
      e=611.21_r8*(1.0007_r8+3.46E-8_r8*p)*EXP(17.502_r8*(t-273.16_r8)/ &
     &            (t-32.19_r8))
      qsatw=eps*e/(p-(1._r8-eps)*e)
!
      RETURN
      END FUNCTION qsatw
!
! --- ------------------------------------------------------------------
!
      FUNCTION qsati(t,p)
!
! --- computes saturation specific humidity [kg/kg] over sea ice,
! --- Parkinson and Washington (1979) JGR 84, 311-337
!
! --- input variables:
! ---   t      - air temperature [K]
! ---   p      - atmospheric pressure [Pa]
!
      USE mod_kinds
      real(r8) :: qsati, t, p
!
! --- parameters:
! ---   eps    - molecular weight ratio of dry air and water vapour
!
      real(r8), parameter :: eps=0.62197_r8
!
      real(r8) :: e
!
      e=611._r8*10._r8**(9.5_r8*(t-273.16_r8)/(t-7.66_r8))
      qsati=eps*e/(p-(1._r8-eps)*e)
!
      RETURN
      END FUNCTION qsati
!
! --- ------------------------------------------------------------------
!
      FUNCTION rhoair(t,q,p)
!
! --- computes air density [kg/m^3]
!
! --- input variables:
! ---   t      - air temperature [K]
! ---   q      - specific humidity of air [kg/kg]
! ---   p      - atmospheric pressure [Pa]
!
      USE mod_kinds
      real(r8) :: rhoair, t, q, p
!
! --- parameters:
! ---   eps    - molecular weight ratio of dry air and water vapour
! ---   rgas   - gas constant for dry air [J/(kg*K)]
!
      real(r8), parameter :: eps=0.62197_r8
      real(r8), parameter :: cv=1.0_r8/eps-1.0_r8
      real(r8), parameter :: rgas=287.04_r8
!
      real(r8) :: tv
!
! --- virtual temperature
      tv=t*(1.0_r8+cv*q)
!
! --- moist air density [kg/m^3]
      rhoair=p/(rgas*tv)
!
      RETURN
      END FUNCTION rhoair
