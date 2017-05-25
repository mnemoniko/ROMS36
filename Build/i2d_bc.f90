      MODULE i2d_bc_mod
! 
!***********************************************************************
!  Compute lateral boundary conditions for any 2D ice variable.
!***********************************************************************
      implicit none
      PRIVATE
      PUBLIC i2d_bc_tile
      CONTAINS
!
!***********************************************************************
      SUBROUTINE i2d_bc_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        liold, linew,                             &
     &                        i_west, i_east, i_north, i_south,         &
     &                        ui, vi, ai, S, switchvar)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liold, linew
      TYPE(T_LBC), intent(in) :: S(4)
      integer, intent(in) :: switchvar   !SM 1/17/2013
      real(r8), intent(in)    :: ui(LBi:,LBj:,:)
      real(r8), intent(in)    :: vi(LBi:,LBj:,:)
      real(r8), intent(inout) :: ai(LBi:,LBj:,:)
      real(r8), intent(in)    :: i_west(LBj:)
      real(r8), intent(in)    :: i_east(LBj:)
      real(r8), intent(in)    :: i_north(LBi:)
      real(r8), intent(in)    :: i_south(LBi:)
!
!  Local variable declarations.
!
      integer :: i, j, know
      real(r8), parameter :: eps =1.0E-20_r8
      real(r8) :: Ce, Cx, cff, dTde, dTdt, dTdx, tau 
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad
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
!  Set time-indices
!-----------------------------------------------------------------------
!
        know=liold
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, implicit upstream radiation condition.
!
        IF (S(iwest)%radiation) THEN
          DO j=Jstr,Jend+1
            grad(Istr-1,j)=ai(Istr-1,j,know)-ai(Istr-1,j-1,know)
            grad(Istr-1,j)=grad(Istr-1,j)*GRID(ng)%vmask(Istr-1,j)
            grad(Istr,j)=ai(Istr,j,know)-ai(Istr,j-1,know)
            grad(Istr,j)=grad(Istr,j)*GRID(ng)%vmask(Istr,j)
          END DO
          DO j=Jstr,Jend
            dTdt=ai(Istr,j,know)-ai(Istr  ,j,linew)
            dTdx=ai(Istr,j,linew)-ai(Istr+1,j,linew)
            IF (S(iwest)%nudging) THEN
              tau=Tobc_out(1,ng,iwest)
              IF ((dTdt*dTdx).lt.0.0_r8) tau=Tobc_in(1,ng,iwest)
              tau=tau*dt(ng)
            END IF
            IF ((dTdt*dTdx).lt.0.0_r8) dTdt=0.0_r8
            IF ((dTdt*(grad(Istr,j)+grad(Istr,j+1))).gt.0.0_r8) THEN
              dTde=grad(Istr,j  )
            ELSE
              dTde=grad(Istr,j+1)
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
            Cx=dTdt*dTdx
            Ce=0.0_r8
            ai(Istr-1,j,linew)=(cff*ai(Istr-1,j,know)+                  &
     &                          Cx *ai(Istr  ,j,linew)-                 &
     &                          MAX(Ce,0.0_r8)*grad(Istr-1,j  )-        &
     &                          MIN(Ce,0.0_r8)*grad(Istr-1,j+1))/       &
     &                              (cff+Cx)
            IF (S(iwest)%nudging) THEN
              ai(Istr-1,j,linew)=ai(Istr-1,j,linew)+                    &
     &                       tau*(i_west(j)-ai(Istr-1,j,know))
            END IF
            ai(Istr-1,j,linew)=ai(Istr-1,j,linew)*                      &
     &                              GRID(ng)%rmask(Istr-1,j)
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (S(iwest)%clamped) THEN
          DO j=Jstr,Jend
             IF(switchvar==1) THEN !hi & hsn
                ai(0,j,linew)=i_west(j)*                                &

     &                        BOUNDARY(ng)%ai_west(j)   !SM 1/17/2013
             ELSE
                ai(0,j,linew)=i_west(j)
             END IF
            ai(0,j,linew)=ai(0,j,linew)*                                &
     &                 GRID(ng)%rmask(0,j)
          END DO
!
!  Western edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (S(iwest)%mixed) THEN
          DO j=Jstr,Jend
            IF (ui(1,j,linew).ge.0._r8) THEN
               IF(switchvar==1) THEN
                  ai(0,j,linew)=i_west(j)*                              &
     &                   BOUNDARY(ng)%ai_west(j) !SM 2013-06-06
               ELSE
                  ai(0,j,linew)=i_west(j)
               END IF
              ai(0,j,linew)=ai(0,j,linew)*                              &
     &                   GRID(ng)%rmask(0,j)
            ELSE
              ai(0,j,linew)=ai(1,j,liold)
              ai(0,j,linew)=ai(0,j,linew)*                              &
     &                   GRID(ng)%rmask(0,j)
            END IF
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (S(iwest)%closed) THEN
          DO j=Jstr,Jend
            ai(0,j,linew)=ai(1,j,linew)
            ai(0,j,linew)=ai(0,j,linew)*                                &
     &                   GRID(ng)%rmask(0,j)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!
        IF (S(ieast)%radiation) THEN
          DO j=Jstr,Jend+1
            grad(Iend,j)=ai(Iend,j,know)-ai(Iend,j-1,know)
            grad(Iend,j)=grad(Iend,j)*GRID(ng)%vmask(Iend  ,j)
            grad(Iend+1,j)=ai(Iend+1,j,know)-ai(Iend+1,j-1,know)
            grad(Iend+1,j)=grad(Iend+1,j)*GRID(ng)%vmask(Iend+1,j)
          END DO
          DO j=Jstr,Jend
            dTdt=ai(Iend,j,know)-ai(Iend  ,j,linew)
            dTdx=ai(Iend,j,linew)-ai(Iend-1,j,linew)
            IF (S(ieast)%nudging) THEN
              tau=Tobc_out(1,ng,ieast)
              IF ((dTdt*dTdx).lt.0.0_r8) tau=Tobc_in(1,ng,ieast)
              tau=tau*dt(ng)
            END IF
            IF ((dTdt*dTdx).lt.0.0_r8) dTdt=0.0_r8
              IF ((dTdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
              dTde=grad(Iend,j  )
            ELSE
              dTde=grad(Iend,j+1)
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
            Cx=dTdt*dTdx
            Ce=0.0_r8
            ai(Iend+1,j,linew)=(cff*ai(Iend+1,j,know)+                  &
     &                          Cx *ai(Iend  ,j,linew)-                 &
     &                          MAX(Ce,0.0_r8)*grad(Iend+1,j  )-        &
     &                          MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/       &
     &                              (cff+Cx)
            IF (S(ieast)%nudging) THEN
              ai(Iend+1,j,linew)=ai(Iend+1,j,linew)+                    &
     &             tau*(i_east(j)-ai(Iend+1,j,know))
            END IF
            ai(Iend+1,j,linew)=ai(Iend+1,j,linew)*                      &
     &                              GRID(ng)%rmask(Iend+1,j)
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (S(ieast)%clamped) THEN
          DO j=Jstr,Jend
             IF(switchvar==1) THEN !hi & hsn
                ai(Lm(ng)+1,j,linew)=i_east(j)*                         &

     &                        BOUNDARY(ng)%ai_east(j)   !SM 1/17/2013
             ELSE
                ai(Lm(ng)+1,j,linew)=i_east(j)
             END IF
            ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                  &
     &                        GRID(ng)%rmask(Lm(ng)+1,j)
          END DO
!
!  Eastern edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (S(iwest)%mixed) THEN
          DO j=Jstr,Jend
            IF (ui(Lm(ng)+1,j,linew).le.0._r8) THEN
               IF(switchvar==1) THEN
                  ai(Lm(ng)+1,j,linew)=i_east(j)*                       &
     &                          BOUNDARY(ng)%ai_east(j) !SM 2013-05-20
               ELSE
                  ai(Lm(ng)+1,j,linew)=i_east(j)
               END IF
              ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
            ELSE
              ai(Lm(ng)+1,j,linew)=ai(Lm(ng),j,liold)
              ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
            END IF
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (S(ieast)%closed) THEN
          DO j=Jstr,Jend
            ai(Lm(ng)+1,j,linew)=ai(Lm(ng),j,linew)
            ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                  &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
        IF (S(isouth)%radiation) THEN
          DO i=Istr,Iend+1
            grad(i,Jstr)=ai(i,Jstr,know)-ai(i-1,Jstr,know)
            grad(i,Jstr)=grad(i,Jstr)*GRID(ng)%umask(i,Jstr)
            grad(i,Jstr-1)=ai(i,Jstr-1,know)-ai(i-1,Jstr-1,know)
            grad(i,Jstr-1)=grad(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
          END DO
          DO i=Istr,Iend
            dTdt=ai(i,Jstr,know)-ai(i,Jstr  ,linew)
            dTde=ai(i,Jstr,linew)-ai(i,Jstr+1,linew)
            IF (S(isouth)%nudging) THEN
              tau=Tobc_out(1,ng,isouth)
              IF ((dTdt*dTde).lt.0.0_r8) tau=Tobc_in(1,ng,isouth)
              tau=tau*dt(ng)
            END IF
            IF ((dTdt*dTde).lt.0.0_r8) dTdt=0.0_r8
            IF ((dTdt*(grad(i,Jstr)+grad(i+1,Jstr))).gt.0.0_r8) THEN
              dTdx=grad(i  ,Jstr)
            ELSE
              dTdx=grad(i+1,Jstr)
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
            Cx=0.0_r8
            Ce=dTdt*dTde
            ai(i,Jstr-1,linew)=(cff*ai(i,Jstr-1,know)+                  &
     &                          Ce *ai(i,Jstr  ,linew)-                 &
     &                          MAX(Cx,0.0_r8)*grad(i  ,Jstr-1)-        &
     &                          MIN(Cx,0.0_r8)*grad(i+1,Jstr-1))/       &
     &                              (cff+Ce)
            IF (S(isouth)%nudging) THEN
              ai(i,Jstr-1,linew)=ai(i,Jstr-1,linew)+                    &
     &           tau*(i_south(i)-ai(i,Jstr-1,know))
            END IF
            ai(i,Jstr-1,linew)=ai(i,Jstr-1,linew)*                      &
     &                              GRID(ng)%rmask(i,Jstr-1)
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (S(isouth)%clamped) THEN
          DO i=Istr,Iend
             IF(switchvar==1) THEN !hi & hsn
                ai(i,0,linew)=i_south(i)*                               &
     &                        BOUNDARY(ng)%ai_south(i)
             ELSE
                ai(i,0,linew)=i_south(i)
             END IF
            ai(i,0,linew)=ai(i,0,linew)*                                &
     &                   GRID(ng)%rmask(i,0)
          END DO
!
!  Southern edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (S(isouth)%mixed) THEN
          DO i=Istr,Iend
            IF (vi(i,1,linew).ge.0._r8) THEN
               IF(switchvar==1) THEN !hi & hsn
                  ai(i,0,linew)=i_south(i)*                             &
     &                          BOUNDARY(ng)%ai_south(i)
               ELSE
                  ai(i,0,linew)=i_south(i)
               END IF
              ai(i,0,linew)=ai(i,0,linew)*                              &
     &                   GRID(ng)%rmask(i,0)
            ELSE
              ai(i,0,linew)=ai(i,1,liold)
              ai(i,0,linew)=ai(i,0,linew)*                              &
     &                   GRID(ng)%rmask(i,0)
            END IF
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (S(isouth)%closed) THEN
          DO i=Istr,Iend
            ai(i,0,linew)=ai(i,1,linew)
            ai(i,0,linew)=ai(i,0,linew)*                                &
     &                   GRID(ng)%rmask(i,0)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
        IF (S(inorth)%radiation) THEN
          DO i=Istr,Iend+1
            grad(i,Jend)=ai(i,Jend,know)-ai(i-1,Jend,know)
            grad(i,Jend)=grad(i,Jend)*GRID(ng)%umask(i,Jend)
            grad(i,Jend+1)=ai(i,Jend+1,know)-ai(i-1,Jend+1,know)
            grad(i,Jend+1)=grad(i,Jend+1)*GRID(ng)%umask(i,Jend+1)
          END DO
          DO i=Istr,Iend
            dTdt=ai(i,Jend,know)-ai(i,Jend  ,linew)
            dTde=ai(i,Jend,linew)-ai(i,Jend-1,linew)
            IF (S(inorth)%nudging) THEN
              tau=Tobc_out(1,ng,inorth)
              IF ((dTdt*dTde).lt.0.0_r8) tau=Tobc_in(1,ng,inorth)
              tau=tau*dt(ng)
            END IF
            IF ((dTdt*dTde).lt.0.0_r8) dTdt=0.0_r8
            IF ((dTdt*(grad(i,Jend)+grad(i+1,Jend))).gt.0.0_r8) THEN
              dTdx=grad(i  ,Jend)
            ELSE
              dTdx=grad(i+1,Jend)
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
            Cx=0.0_r8
            Ce=dTdt*dTde
            ai(i,Jend+1,linew)=(cff*ai(i,Jend+1,know)+                  &
     &                          Ce *ai(i,Jend  ,linew)-                 &
     &                          MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-        &
     &                          MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/       &
     &                              (cff+Ce)
            IF (S(inorth)%nudging) THEN
              ai(i,Jend+1,linew)=ai(i,Jend+1,linew)+                    &
     &              tau*(i_north(i)-ai(i,Jend+1,know))
            END IF
            ai(i,Jend+1,linew)=ai(i,Jend+1,linew)*                      &
     &                              GRID(ng)%rmask(i,Jend+1)
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (S(inorth)%clamped) THEN
          DO i=Istr,Iend
             IF(switchvar==1) THEN !hi & hsn
                ai(i,Mm(ng)+1,linew)=i_north(i)*                        &

     &                        BOUNDARY(ng)%ai_north(i)   !SM 1/17/2013
             ELSE
                ai(i,Mm(ng)+1,linew)=i_north(i)
             END IF
            ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                  &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
          END DO
!
!  Northern edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (S(inorth)%mixed) THEN
          DO i=Istr,Iend
            IF (vi(i,Mm(ng)+1,linew).le.0._r8) THEN
               IF(switchvar==1) THEN
                  ai(i,Mm(ng)+1,linew)=i_north(i)*                      &
     &                          BOUNDARY(ng)%ai_north(i) !SM 2013-06-06
               ELSE
                  ai(i,Mm(ng)+1,linew)=i_north(i)
               END IF
              ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
            ELSE
              ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng),liold)
              ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
            END IF
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (S(inorth)%closed) THEN
          DO i=Istr,Iend
            ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng),linew)
            ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                  &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng).and. .not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          ai(0,0,linew)=0.5_r8*(ai(1,0,linew)+                          &
     &                         ai(0,1,linew))
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          ai(Lm(ng)+1,0,linew)=0.5_r8*(ai(Lm(ng)+1,1,linew)+            &
     &                                ai(Lm(ng)  ,0,linew))
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          ai(0,Mm(ng)+1,linew)=0.5_r8*(ai(0,Mm(ng)  ,linew)+            &
     &                                ai(1,Mm(ng)+1,linew))
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          ai(Lm(ng)+1,Mm(ng)+1,linew)=0.5_r8*                           &
     &             (ai(Lm(ng)+1,Mm(ng)  ,linew)+                        &
     &              ai(Lm(ng)  ,Mm(ng)+1,linew))
        END IF
      END IF
      RETURN
      END SUBROUTINE i2d_bc_tile
      END MODULE i2d_bc_mod
