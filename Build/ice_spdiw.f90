      MODULE ice_spdiw_mod
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This module computes the magnitude of the shear between the ice
!  and the surface water. In this case, the surface water is defined
!  as the water in a surface mixed layer, so that velocity must be
!  computed first.
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC ice_spdiw
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_spdiw (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean
      USE mod_ice
      USE mod_coupling
      USE mod_mixing
      USE mod_stepping
!
      implicit none
!
      integer, intent(in) :: ng, tile
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
      CALL wclock_on (ng, iNLM, 6)
      CALL ice_spdiw_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nrhs(ng),                                    &
     &                     liuol(ng),                                   &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w,                              &
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % v,                               &
     &                     MIXING(ng) % hsbl,                           &
     &                     ICE(ng) % ui,                                &
     &                     ICE(ng) % vi,                                &
     &                     ICE(ng) % uwater,                            &
     &                     ICE(ng) % vwater,                            &
     &                     ICE(ng) % spd_iw                             &
     &                     )
      CALL wclock_off (ng, iNLM, 6)
      RETURN
      END SUBROUTINE ice_spdiw
!
!***********************************************************************
      SUBROUTINE ice_spdiw_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs,                                  &
     &                           liuol,                                 &
     &                           Hz, z_r, z_w,                          &
     &                           u, v,                                  &
     &                           hsbl,                                  &
     &                           ui, vi,                                &
     &                           uwater, vwater,                        &
     &                           spd_iw)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod
      USE mp_exchange_mod, ONLY : mp_exchange2d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
      integer, intent(in) :: liuol
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: hsbl(LBi:,LBj:)
      real(r8), intent(in) :: ui(LBi:,LBj:,:)
      real(r8), intent(in) :: vi(LBi:,LBj:,:)
      real(r8), intent(out) :: uwater(LBi:,LBj:)
      real(r8), intent(out) :: vwater(LBi:,LBj:)
      real(r8), intent(out) :: spd_iw(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j
      integer :: nlio, nbotu, nbotv, k
      integer, dimension(IminS:ImaxS,JminS:JmaxS) :: nbot
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: uw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vw
      real(r8) :: mlio
      real(r8) :: dml
      real(r8) :: totml
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
        do j=MAX(Jstr-2,0),MIN(Jend+2,Mm(ng)+1)
        do i=MAX(Istr-2,0),MIN(Iend+2,Lm(ng)+1)
!             sl_dpth = lmd_epsilon*(z_w(i,j,N(ng))-hsbl(i,j))
              mlio = min(hsbl(i,j),-10._r8)
              nbot(i,j) = 1
              do k=N(ng),1,-1
                  if(z_r(i,j,k).lt.mlio) then
                     nbot(i,j) = min(k,N(ng))
                     nbot(i,j) = max(nbot(i,j),1)
                     goto 1111
                  endif
              enddo
 1111         continue
        enddo
        enddo
        do j=Jstr,Jend
        do i=MAX(Istr-1,1),Iend+1
             nlio = 0
             nbotu = NINT(0.5_r8*(nbot(i-1,j)+nbot(i,j)))
             nbotu = max(min(nbotu,N(ng)),1)
             uw(i,j) = 0._r8
             totml = 0._r8
             do k=N(ng),nbotu,-1
                nlio = nlio + 1
                dml = 0.5_r8*(z_w(i-1,j,k)-z_w(i-1,j,k-1)               &
     &                      + z_w(i,j,k)-z_w(i,j,k-1))
                uw(i,j) = uw(i,j) + u(i,j,k,nrhs)*dml
                totml = totml + dml
             enddo
             uw(i,j) = uw(i,j)/totml
!            uw(i,j) =  u(i,j,N,nrhs)
          enddo
        enddo
        do j=MAX(Jstr-1,1),Jend+1
          do i=Istr,Iend
             nlio = 0
             nbotv = NINT(0.5_r8*(nbot(i,j-1)+nbot(i,j)))
             nbotv = max(min(nbotv,N(ng)),1)
             vw(i,j) = 0._r8
             totml = 0._r8
             do k=N(ng),nbotv,-1
                nlio = nlio + 1
                dml = 0.5_r8*(z_w(i,j-1,k)-z_w(i,j-1,k-1)               &
     &                      + z_w(i,j,k)-z_w(i,j,k-1))
                vw(i,j) = vw(i,j) + v(i,j,k,nrhs)*dml
                totml = totml + dml
             enddo
             vw(i,j) = vw(i,j)/totml
!            vw(i,j) =  v(i,j,N,nrhs)
          enddo
        enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          spd_iw(i,j) = 0.5*sqrt((uw(i,j)-ui(i,j,liuol)                 &
     &                 +  uw(i+1,j)-ui(i+1,j,liuol))**2                 &
     &                  +(vw(i,j)-vi(i,j,liuol)                         &
     &                 +  vw(i,j+1)-vi(i,j+1,liuol))**2)
        enddo
      enddo
      do j=Jstr,Jend
        do i=IstrU,Iend
           uwater(i,j) = uw(i,j)
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
           vwater(i,j) = vw(i,j)
        enddo
      enddo
!
!  Apply boundary conditions.
!
        CALL bc_r2d_tile (ng, tile,                                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          spd_iw)
        CALL bc_u2d_tile (ng, tile,                                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          uwater)
        CALL bc_v2d_tile (ng, tile,                                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vwater)
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    spd_iw, uwater, vwater)
      RETURN
      END SUBROUTINE ice_spdiw_tile
      END MODULE ice_spdiw_mod
