      MODULE uibc_mod
!***********************************************************************
!  Compute the lateral boundary conditions on the ice U-velocity.
!***********************************************************************
      implicit none
      PRIVATE
      PUBLIC uibc_tile
      CONTAINS
!
!***********************************************************************
      SUBROUTINE uibc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
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
      CALL  uibc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 liuol(ng), liunw(ng),                            &
     &                 ICE(ng) % ui)
      RETURN
      END SUBROUTINE uibc
!
!***********************************************************************
      SUBROUTINE uibc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liuol, liunw,                               &
     &                      ui)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
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
      integer, intent(in) :: liuol, liunw
      real(r8), intent(inout) :: ui(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, Imax, Imin, j, know
      real(r8), parameter :: eps =1.0E-20_r8
      real(r8) :: Ce, Cx
      real(r8) :: cff, dUde, dUdt, dUdx, tau
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
        know=liuol
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        IF (LBC(iwest,isUice,ng)%radiation) THEN
!
!  Western edge, implicit upstream radiation condition.
!
          DO j=Jstr,Jend+1
            grad(Istr  ,j)=ui(Istr  ,j  ,know)-                         &
     &                     ui(Istr  ,j-1,know)
            grad(Istr+1,j)=ui(Istr+1,j  ,know)-                         &
     &                     ui(Istr+1,j-1,know)
          END DO
          DO j=Jstr,Jend
            dUdt=ui(Istr+1,j,know)-ui(Istr+1,j,liunw)
            dUdx=ui(Istr+1,j,liunw)-ui(Istr+2,j,liunw)
            IF (LBC(iwest,isUice,ng)%nudging) THEN
              IF ((dUdt*dUdx).lt.0.0_r8) THEN
                tau=M2obc_in(ng,iwest)
              ELSE
                tau=M2obc_out(ng,iwest)
              END IF
              tau=tau*dt(ng)
            END IF
            IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
            IF ((dUdt*(grad(Istr+1,j)+grad(Istr+1,j+1))).gt.0.0_r8) THEN
              dUde=grad(Istr+1,j  )
            ELSE
              dUde=grad(Istr+1,j+1)
            END IF
            cff=MAX(dUdx*dUdx+dUde*dUde,eps)
            Cx=dUdt*dUdx
            Ce=0.0_r8
            ui(Istr,j,liunw)=(cff*ui(Istr  ,j,know)+                    &
     &                        Cx *ui(Istr+1,j,liunw)-                   &
     &                        MAX(Ce,0.0_r8)*grad(Istr,j  )-            &
     &                        MIN(Ce,0.0_r8)*grad(Istr,j+1))/           &
     &                       (cff+Cx)
            IF (LBC(iwest,isUice,ng)%nudging) THEN
              ui(Istr,j,liunw)=ui(Istr,j,liunw)+                        &
     &                       tau*(BOUNDARY(ng)%ui_west(j)-              &
     &                       ui(Istr,j,know))
            END IF
            ui(Istr,j,liunw)=ui(Istr,j,liunw)*                          &
     &                       GRID(ng)%umask(Istr,j)
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (LBC(iwest,isUice,ng)%clamped) THEN
          DO j=Jstr,Jend
            ui(1,j,liunw)=BOUNDARY(ng)%ui_west(j)
            ui(1,j,liunw)=ui(1,j,liunw)*                                &
     &                   GRID(ng)%umask(1,j)
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isUice,ng)%gradient) THEN
          DO j=Jstr,Jend
            ui(1,j,liunw)=ui(2,j,liunw)
            ui(1,j,liunw)=ui(1,j,liunw)*                                &
     &                   GRID(ng)%umask(1,j)
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(iwest,isUice,ng)%closed) THEN
          DO j=Jstr,Jend
            ui(1,j,liunw)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        IF (LBC(ieast,isUice,ng)%radiation) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!
          DO j=Jstr,Jend+1
            grad(Iend  ,j)=ui(Iend  ,j  ,know)-                         &
     &                     ui(Iend  ,j-1,know)
            grad(Iend+1,j)=ui(Iend+1,j  ,know)-                         &
     &                     ui(Iend+1,j-1,know)
          END DO
          DO j=Jstr,Jend
            dUdt=ui(Iend,j,know)-ui(Iend  ,j,liunw)
            dUdx=ui(Iend,j,liunw)-ui(Iend-1,j,liunw)
            IF (LBC(ieast,isUice,ng)%nudging) THEN
              IF ((dUdt*dUdx).lt.0.0_r8) THEN
                tau=M2obc_in(ng,ieast)
              ELSE
              tau=M2obc_out(ng,ieast)
              END IF
              tau=tau*dt(ng)
            END IF
            IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
            IF ((dUdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
              dUde=grad(Iend,j)
            ELSE
              dUde=grad(Iend,j+1)
            END IF
            cff=MAX(dUdx*dUdx+dUde*dUde,eps)
            Cx=dUdt*dUdx
            Ce=0.0_r8
            ui(Iend+1,j,liunw)=(cff*ui(Iend+1,j,know)+                  &
     &                         Cx *ui(Iend  ,j,liunw)-                  &
     &                         MAX(Ce,0.0_r8)*grad(Iend+1,j  )-         &
     &                         MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/        &
     &                        (cff+Cx)
            IF (LBC(ieast,isUice,ng)%nudging) THEN
              ui(Iend+1,j,liunw)=ui(Iend+1,j,liunw)+                    &
     &                           tau*(BOUNDARY(ng)%ui_east(j)-          &
     &                           ui(Iend+1,j,know))
            END IF
            ui(Iend+1,j,liunw)=ui(Iend+1,j,liunw)*                      &
     &                        GRID(ng)%umask(Iend+1,j)
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (LBC(ieast,isUice,ng)%clamped) THEN
          DO j=Jstr,Jend
            ui(Lm(ng)+1,j,liunw)=BOUNDARY(ng)%ui_east(j)
            ui(Lm(ng)+1,j,liunw)=ui(Lm(ng)+1,j,liunw)*                  &
     &                          GRID(ng)%umask(Lm(ng)+1,j)
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isUice,ng)%gradient) THEN
          DO j=Jstr,Jend
            ui(Lm(ng)+1,j,liunw)=ui(Lm(ng),j,liunw)
            ui(Lm(ng)+1,j,liunw)=ui(Lm(ng)+1,j,liunw)*                  &
     &                          GRID(ng)%umask(Lm(ng)+1,j)
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (LBC(ieast,isUice,ng)%closed) THEN
          DO j=Jstr,Jend
            ui(Lm(ng)+1,j,liunw)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        IF (LBC(isouth,isUice,ng)%radiation) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
          DO i=IstrU-1,Iend
            grad(i,Jstr-1)=ui(i+1,Jstr-1,know)-                         &
     &                     ui(i  ,Jstr-1,know)
            grad(i,Jstr  )=ui(i+1,Jstr  ,know)-                         &
     &                     ui(i  ,Jstr  ,know)
          END DO
          DO i=IstrU,Iend
            dUdt=ui(i,Jstr,know)-ui(i,Jstr  ,liunw)
            dUde=ui(i,Jstr,liunw)-ui(i,Jstr+1,liunw)
            IF (LBC(isouth,isUice,ng)%nudging) THEN
              IF ((dUdt*dUde).lt.0.0_r8) THEN
                tau=M2obc_in(ng,isouth)
              ELSE
                tau=M2obc_out(ng,isouth)
              END IF
              tau=tau*dt(ng)
            END IF
            IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
            IF ((dUdt*(grad(i-1,Jstr)+grad(i,Jstr))).gt.0.0_r8) THEN
              dUdx=grad(i-1,Jstr)
            ELSE
                dUdx=grad(i  ,Jstr)
            END IF
            cff=MAX(dUdx*dUdx+dUde*dUde,eps)
            Cx=0.0_r8
            Ce=dUdt*dUde
            ui(i,Jstr-1,liunw)=(cff*ui(i,Jstr-1,know)+                  &
     &                         Ce *ui(i,Jstr  ,liunw)-                  &
     &                         MAX(Cx,0.0_r8)*grad(i-1,Jstr-1)-         &
     &                         MIN(Cx,0.0_r8)*grad(i  ,Jstr-1))/        &
     &                        (cff+Ce)
            IF (LBC(isouth,isUice,ng)%nudging) THEN
              ui(i,Jstr-1,liunw)=ui(i,Jstr-1,liunw)+                    &
     &                          tau*(BOUNDARY(ng)%ui_south(i)-          &
     &                             ui(i,Jstr-1,know))
            END IF
            ui(i,Jstr-1,liunw)=ui(i,Jstr-1,liunw)*                      &
     &                        GRID(ng)%umask(i,Jstr-1)
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,isUice,ng)%clamped) THEN
          DO i=IstrU,Iend
            ui(i,0,liunw)=BOUNDARY(ng)%ui_south(i)
            ui(i,0,liunw)=ui(i,0,liunw)*                                &
     &                   GRID(ng)%umask(i,0)
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isUice,ng)%gradient) THEN
          DO i=IstrU,Iend
            ui(i,0,liunw)=ui(i,1,liunw)
            ui(i,0,liunw)=ui(i,0,liunw)*                                &
     &                   GRID(ng)%umask(i,0)
          END DO
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (LBC(isouth,isUice,ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrU
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendR
          END IF
          DO i=Imin,Imax
            ui(i,0,liunw)=gamma2(ng)*ui(i,1,liunw)
            ui(i,0,liunw)=ui(i,0,liunw)*                                &
     &                   GRID(ng)%umask(i,0)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        IF (LBC(inorth,isUice,ng)%radiation) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
          DO i=IstrU-1,Iend
            grad(i,Jend  )=ui(i+1,Jend  ,know)-                         &
     &                     ui(i  ,Jend  ,know)
            grad(i,Jend+1)=ui(i+1,Jend+1,know)-                         &
     &                     ui(i  ,Jend+1,know)
          END DO
          DO i=IstrU,Iend
            dUdt=ui(i,Jend,know)-ui(i,Jend  ,liunw)
            dUde=ui(i,Jend,liunw)-ui(i,Jend-1,liunw)
            IF (LBC(inorth,isUice,ng)%nudging) THEN
              IF ((dUdt*dUde).lt.0.0_r8) THEN
                tau=M2obc_in(ng,inorth)
              ELSE
                tau=M2obc_out(ng,inorth)
              END IF
              tau=tau*dt(ng)
            END IF
            IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
            IF ((dUdt*(grad(i-1,Jend)+grad(i,Jend))).gt.0.0_r8) THEN
              dUdx=grad(i-1,Jend)
            ELSE
              dUdx=grad(i  ,Jend)
            END IF
            cff=MAX(dUdx*dUdx+dUde*dUde,eps)
            Cx=0.0_r8
            Ce=dUdt*dUde
            ui(i,Jend+1,liunw)=(cff*ui(i,Jend+1,know)+                  &
     &                         Ce *ui(i,Jend  ,liunw)-                  &
     &                         MAX(Cx,0.0_r8)*grad(i-1,Jend+1)-         &
     &                         MIN(Cx,0.0_r8)*grad(i  ,Jend+1))/        &
     &                        (cff+Ce)
            ui(i,Jend+1,liunw)=ui(i,Jend+1,liunw)*                      &
     &                        GRID(ng)%umask(i,Jend+1)
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (LBC(inorth,isUice,ng)%clamped) THEN
          DO i=IstrU,Iend
            ui(i,Mm(ng)+1,liunw)=BOUNDARY(ng)%ui_north(i)
            ui(i,Mm(ng)+1,liunw)=ui(i,Mm(ng)+1,liunw)*                  &
     &                          GRID(ng)%umask(i,Mm(ng)+1)
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isUice,ng)%gradient) THEN
          DO i=IstrU,Iend
            ui(i,Mm(ng)+1,liunw)=ui(i,Mm(ng),liunw)
            ui(i,Mm(ng)+1,liunw)=ui(i,Mm(ng)+1,liunw)*                  &
     &                          GRID(ng)%umask(i,Mm(ng)+1)
          END DO
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (LBC(inorth,isUice,ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrU
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendR
          END IF
          DO i=Imin,Imax
            ui(i,Mm(ng)+1,liunw)=gamma2(ng)*ui(i,Mm(ng),liunw)
            ui(i,Mm(ng)+1,liunw)=ui(i,Mm(ng)+1,liunw)*                  &
     &                          GRID(ng)%umask(i,Mm(ng)+1)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          ui(1,0,liunw)=0.5_r8*(ui(2,0,liunw)+                          &
     &                         ui(1,1,liunw))
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          ui(Lm(ng)+1,0,liunw)=0.5_r8*(ui(Lm(ng)  ,0,liunw)+            &
     &                                ui(Lm(ng)+1,1,liunw))
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          ui(1,Mm(ng)+1,liunw)=0.5_r8*(ui(2,Mm(ng)+1,liunw)+            &
     &                                ui(1,Mm(ng)  ,liunw))
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          ui(Lm(ng)+1,Mm(ng)+1,liunw)=0.5_r8*                           &
     &                               (ui(Lm(ng)  ,Mm(ng)+1,liunw)+      &
     &                                ui(Lm(ng)+1,Mm(ng)  ,liunw))
        END IF
      END IF
      RETURN
      END SUBROUTINE uibc_tile
      END MODULE uibc_mod
