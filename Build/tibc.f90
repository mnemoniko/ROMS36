      MODULE tibc_mod
!***********************************************************************
!  Compute the lateral boundary conditions on the internal ice
!  temperature.
!***********************************************************************
      implicit none
      PRIVATE
      PUBLIC tibc_tile
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tibc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
      USE mod_scalars
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
      CALL tibc_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                liold(ng), linew(ng),                             &
     &                ICE(ng) % ui,                                     &
     &                ICE(ng) % vi,                                     &
     &                ICE(ng) % hi,                                     &
     &                ICE(ng) % ti,                                     &
     &                ICE(ng) % enthalpi)
      RETURN
      END SUBROUTINE tibc
!
!***********************************************************************
      SUBROUTINE tibc_tile (ng, tile,                                   &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           liold, linew,                          &
     &                           ui, vi, hi, ti, enthalpi)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_boundary
      USE mod_grid
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: liold, linew
      real(r8), intent(in)    :: ui(LBi:,LBj:,:)
      real(r8), intent(in)    :: vi(LBi:,LBj:,:)
      real(r8), intent(in)    :: hi(LBi:,LBj:,:)
      real(r8), intent(inout) :: ti(LBi:,LBj:,:)
      real(r8), intent(inout) :: enthalpi(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, j, know
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
!  Western edge, clamped boundary condition.
!
        IF (LBC(iwest,isTice,ng)%clamped) THEN
          DO j=Jstr,Jend
            enthalpi(0,j,linew)=BOUNDARY(ng)%hi_west(j)*                &
     &                           BOUNDARY(ng)%ti_west(j)
            enthalpi(0,j,linew)=enthalpi(0,j,linew)*                    &
     &                   GRID(ng)%rmask(0,j)
          END DO
!
!  Western edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (LBC(iwest,isTice,ng)%mixed) THEN
          DO j=Jstr,Jend
            IF (ui(1,j,linew).ge.0._r8) THEN
              enthalpi(0,j,linew)=BOUNDARY(ng)%hi_west(j)*              &
     &                           BOUNDARY(ng)%ti_west(j)
              enthalpi(0,j,linew)=enthalpi(0,j,linew)*                  &
     &                   GRID(ng)%rmask(0,j)
            ELSE
              enthalpi(0,j,linew)=enthalpi(1,j,liold)
              enthalpi(0,j,linew)=enthalpi(0,j,linew)*                  &
     &                   GRID(ng)%rmask(0,j)
            END IF
            ti(0,j,linew) = enthalpi(0,j,linew)/                        &
     &                       MAX(hi(0,j,linew),1.0E-6_r8)
            IF (hi(0,j,linew).LE.min_h(ng)) THEN
              enthalpi(0,j,linew) = 0.0_r8
              ti(0,j,linew) = 0.0_r8
            END IF
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isTice,ng)%gradient) THEN
          DO j=Jstr,Jend
            enthalpi(0,j,linew)=hi(1,j,linew)*ti(1,j,linew)
            enthalpi(0,j,linew)=enthalpi(0,j,linew)*                    &
     &                   GRID(ng)%rmask(0,j)
            ti(0,j,linew) = enthalpi(0,j,linew)/                        &
     &                       MAX(hi(1,j,linew),1.0E-6_r8)
            IF (hi(0,j,linew).LE.min_h(ng)) THEN
              enthalpi(0,j,linew) = 0.0_r8
              ti(0,j,linew) = 0.0_r8
            END IF
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(iwest,isTice,ng)%closed) THEN
          DO j=Jstr,Jend
            enthalpi(0,j,linew)=hi(1,j,linew)*ti(1,j,linew)
            enthalpi(0,j,linew)=enthalpi(0,j,linew)*                    &
     &                   GRID(ng)%rmask(0,j)
            ti(0,j,linew) = enthalpi(0,j,linew)/                        &
     &                       MAX(hi(0,j,linew),1.0E-6_r8)
            IF (hi(0,j,linew).LE.min_h(ng)) THEN
              enthalpi(0,j,linew) = 0.0_r8
              ti(0,j,linew) = 0.0_r8
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        IF (LBC(ieast,isTice,ng)%clamped) THEN
!
!  Eastern edge, clamped boundary condition.
!
          DO j=Jstr,Jend
            enthalpi(Lm(ng)+1,j,linew)=BOUNDARY(ng)%hi_east(j)*         &
     &                                  BOUNDARY(ng)%ti_east(j)
            enthalpi(Lm(ng)+1,j,linew)=enthalpi(Lm(ng)+1,j,linew)*      &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
          END DO
!
!  Eastern edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (LBC(ieast,isTice,ng)%mixed) THEN
          DO j=Jstr,Jend
            IF (ui(Lm(ng)+1,j,linew).le.0._r8) THEN
              enthalpi(Lm(ng)+1,j,linew)=BOUNDARY(ng)%hi_east(j)*       &
     &                                  BOUNDARY(ng)%ti_east(j)
              enthalpi(Lm(ng)+1,j,linew)=enthalpi(Lm(ng)+1,j,linew)*    &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
            ELSE
              enthalpi(Lm(ng)+1,j,linew)=hi(Lm(ng),j,liold)*            &
     &                                  ti(Lm(ng),j,liold)
              enthalpi(Lm(ng)+1,j,linew)=enthalpi(Lm(ng)+1,j,linew)*    &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
              ti(Lm(ng)+1,j,linew) = enthalpi(Lm(ng)+1,j,linew)/        &
     &                       MAX(hi(Lm(ng)+1,j,linew),1.0E-6_r8)
              IF (hi(Lm(ng)+1,j,linew).LE.min_h(ng)) THEN
                enthalpi(Lm(ng)+1,j,linew) = 0.0_r8
                ti(Lm(ng)+1,j,linew) = 0.0_r8
              END IF
            END IF
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isTice,ng)%gradient) THEN
          DO j=Jstr,Jend
            enthalpi(Lm(ng)+1,j,linew)=hi(Lm(ng),j,linew)*              &
     &                               ti(Lm(ng),j,linew)
            enthalpi(Lm(ng)+1,j,linew)=enthalpi(Lm(ng)+1,j,linew)*      &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
            ti(Lm(ng)+1,j,linew) = enthalpi(Lm(ng)+1,j,linew)/          &
     &                       MAX(hi(Lm(ng),j,linew),1.0E-6_r8)
            IF (hi(Lm(ng)+1,j,linew).LE.min_h(ng)) THEN
              enthalpi(Lm(ng)+1,j,linew) = 0.0_r8
              ti(Lm(ng)+1,j,linew) = 0.0_r8
            END IF
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (LBC(ieast,isTice,ng)%closed) THEN
          DO j=Jstr,Jend
            enthalpi(Lm(ng)+1,j,linew)=hi(Lm(ng),j,linew)*              &
     &                               ti(Lm(ng),j,linew)
            enthalpi(Lm(ng)+1,j,linew)=enthalpi(Lm(ng)+1,j,linew)*      &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
            ti(Lm(ng)+1,j,linew) = enthalpi(Lm(ng)+1,j,linew)/          &
     &                       MAX(hi(Lm(ng)+1,j,linew),1.0E-6_r8)
            IF (hi(Lm(ng)+1,j,linew).LE.min_h(ng)) THEN
              enthalpi(Lm(ng)+1,j,linew) = 0.0_r8
              ti(Lm(ng)+1,j,linew) = 0.0_r8
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        IF (LBC(isouth,isTice,ng)%clamped) THEN
!
!  Southern edge, clamped boundary condition.
!
          DO i=Istr,Iend
            enthalpi(i,0,linew)=BOUNDARY(ng)%hi_south(i)*               &
     &                          BOUNDARY(ng)%ti_south(i)            
            enthalpi(i,0,linew)=enthalpi(i,0,linew)*                    &
     &                   GRID(ng)%rmask(i,0)
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,isTice,ng)%mixed) THEN
          DO i=Istr,Iend
            IF (vi(i,1,linew).ge.0._r8) THEN
              enthalpi(i,0,linew)=BOUNDARY(ng)%hi_south(i)*             &
     &                          BOUNDARY(ng)%ti_south(i)            
              enthalpi(i,0,linew)=enthalpi(i,0,linew)*                  &
     &                   GRID(ng)%rmask(i,0)
            ELSE
              enthalpi(i,0,linew)=enthalpi(i,1,liold)
              enthalpi(i,0,linew)=enthalpi(i,0,linew)*                  &
     &                   GRID(ng)%rmask(i,0)
              ti(i,0,linew) = enthalpi(i,0,linew)/                      &
     &                       MAX(hi(i,0,linew),1.0E-6_r8)
              IF (hi(i,0,linew).LE.min_h(ng)) THEN
                enthalpi(i,0,linew) = 0.0_r8
                ti(i,0,linew) = 0.0_r8
              END IF
            ENDIF
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isTice,ng)%gradient) THEN
          DO i=Istr,Iend
            enthalpi(i,0,linew)=hi(i,1,linew)*ti(i,1,linew)
            enthalpi(i,0,linew)=enthalpi(i,0,linew)*                    &
     &                   GRID(ng)%rmask(i,0)
            ti(i,0,linew) = enthalpi(i,0,linew)/                        &
     &                       MAX(hi(i,1,linew),1.0E-6_r8)
            IF (hi(i,0,linew).LE.min_h(ng)) THEN
              enthalpi(i,0,linew) = 0.0_r8
              ti(i,0,linew) = 0.0_r8
            END IF
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (LBC(isouth,isTice,ng)%closed) THEN
          DO i=Istr,Iend
            enthalpi(i,0,linew)=enthalpi(i,1,linew)
            enthalpi(i,0,linew)=enthalpi(i,0,linew)*                    &
     &                   GRID(ng)%rmask(i,0)
            ti(i,0,linew) = enthalpi(i,0,linew)/                        &
     &                    MAX(hi(i,0,linew),1.0E-6_r8)
            IF (hi(i,0,linew).LE.min_h(ng)) THEN
              enthalpi(i,0,linew) = 0.0_r8
              ti(i,0,linew) = 0.0_r8
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        IF (LBC(inorth,isTice,ng)%clamped) THEN
!
!  Northern edge, clamped boundary condition.
!
          DO i=Istr,Iend
            enthalpi(i,Mm(ng)+1,linew)=BOUNDARY(ng)%hi_north(i)*        &
     &                                 BOUNDARY(ng)%ti_north(i)
            enthalpi(i,Mm(ng)+1,linew)=enthalpi(i,Mm(ng)+1,linew)*      &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
          END DO
!
!  Northern edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (LBC(inorth,isTice,ng)%mixed) THEN
          DO i=Istr,Iend
            IF (vi(i,Mm(ng)+1,linew).le.0._r8) THEN
              enthalpi(i,Mm(ng)+1,linew)=BOUNDARY(ng)%hi_north(i)*      &
     &                                 BOUNDARY(ng)%ti_north(i)
              enthalpi(i,Mm(ng)+1,linew)=enthalpi(i,Mm(ng)+1,linew)*    &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
          ELSE
              enthalpi(i,Mm(ng)+1,linew)=enthalpi(i,Mm(ng),liold)
              enthalpi(i,Mm(ng)+1,linew)=enthalpi(i,Mm(ng)+1,linew)*    &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
          ENDIF
            ti(i,Mm(ng)+1,linew) = enthalpi(i,Mm(ng)+1,linew)/          &
     &                       MAX(hi(i,Mm(ng)+1,linew),1.0E-6_r8)
            IF (hi(i,Mm(ng)+1,linew).LE.min_h(ng)) THEN
              enthalpi(i,Mm(ng)+1,linew) = 0.0_r8
              ti(i,Mm(ng)+1,linew) = 0.0_r8
            END IF
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isTice,ng)%gradient) THEN
          DO i=Istr,Iend
            enthalpi(i,Mm(ng)+1,linew)=hi(i,Mm(ng),linew)*              &
     &                               ti(i,Mm(ng),linew)
            enthalpi(i,Mm(ng)+1,linew)=enthalpi(i,Mm(ng)+1,linew)*      &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
            ti(i,Mm(ng)+1,linew) = enthalpi(i,Mm(ng)+1,linew)/          &
     &                       MAX(hi(i,Mm(ng),linew),1.0E-6_r8)
            IF (hi(i,Mm(ng)+1,linew).LE.min_h(ng)) THEN
              enthalpi(i,Mm(ng)+1,linew) = 0.0_r8
              ti(i,Mm(ng)+1,linew) = 0.0_r8
            END IF
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (LBC(inorth,isTice,ng)%closed) THEN
          DO i=Istr,Iend
            enthalpi(i,Mm(ng)+1,linew)=hi(i,Mm(ng),linew)*              &
     &                               ti(i,Mm(ng),linew)
            enthalpi(i,Mm(ng)+1,linew)=enthalpi(i,Mm(ng)+1,linew)*      &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
            ti(i,Mm(ng)+1,linew) = enthalpi(i,Mm(ng)+1,linew)/          &
     &                       MAX(hi(i,Mm(ng)+1,linew),1.0E-6_r8)
            IF (hi(i,Mm(ng)+1,linew).LE.min_h(ng)) THEN
              enthalpi(i,Mm(ng)+1,linew) = 0.0_r8
              ti(i,Mm(ng)+1,linew) = 0.0_r8
            END IF
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
          enthalpi(0,0,linew)=0.5_r8*(enthalpi(1,0,linew)+              &
     &                         enthalpi(0,1,linew))
          enthalpi(0,0,linew)=enthalpi(0,0,linew)*                      &
     &                   GRID(ng)%rmask(0,0)
          ti(0,0,linew) = enthalpi(0,0,linew)/                          &
     &                       MAX(hi(0,0,linew),1.0E-6_r8)
          IF (hi(0,0,linew).LE.min_h(ng)) THEN
            enthalpi(0,0,linew) = 0.0_r8
              ti(0,0,linew) = 0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          enthalpi(Lm(ng)+1,0,linew)=0.5_r8*(enthalpi(Lm(ng)+1,1,linew)+&
     &                                enthalpi(Lm(ng)  ,0,linew))
          enthalpi(Lm(ng)+1,0,linew)=enthalpi(Lm(ng)+1,0,linew)*        &
     &                   GRID(ng)%rmask(Lm(ng)+1,0)
          ti(Lm(ng)+1,0,linew) = enthalpi(Lm(ng)+1,0,linew)/            &
     &                       MAX(hi(Lm(ng)+1,0,linew),1.0E-6_r8)
          IF (hi(Lm(ng)+1,0,linew).LE.min_h(ng)) THEN
            enthalpi(Lm(ng)+1,0,linew) = 0.0_r8
            ti(Lm(ng)+1,0,linew) = 0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          enthalpi(0,Mm(ng)+1,linew)=0.5_r8*(enthalpi(0,Mm(ng)  ,linew)+&
     &                                enthalpi(1,Mm(ng)+1,linew))
          enthalpi(0,Mm(ng)+1,linew)=enthalpi(0,Mm(ng)+1,linew)*        &
     &                   GRID(ng)%rmask(0,Mm(ng)+1)
          ti(0,Mm(ng)+1,linew) = enthalpi(0,Mm(ng)+1,linew)/            &
     &                       MAX(hi(0,Mm(ng)+1,linew),1.0E-6_r8)
          IF (hi(0,Mm(ng)+1,linew).LE.min_h(ng)) THEN
            enthalpi(0,Mm(ng)+1,linew) = 0.0_r8
            ti(0,Mm(ng)+1,linew) = 0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          enthalpi(Lm(ng)+1,Mm(ng)+1,linew)=0.5_r8*                     &
     &             (enthalpi(Lm(ng)+1,Mm(ng)  ,linew)+                  &
     &              enthalpi(Lm(ng)  ,Mm(ng)+1,linew))
          enthalpi(Lm(ng)+1,Mm(ng)+1,linew)=                            &
     &             enthalpi(Lm(ng)+1,Mm(ng)+1,linew)*                   &
     &             GRID(ng)%rmask(Lm(ng)+1,Mm(ng)+1)
          ti(Lm(ng)+1,Mm(ng)+1,linew) =                                 &
     &         enthalpi(Lm(ng)+1,Mm(ng)+1,linew)/                       &
     &         MAX(hi(Lm(ng)+1,Mm(ng)+1,linew),1.0E-6_r8)
          IF (hi(Lm(ng)+1,Mm(ng)+1,linew).LE.min_h(ng)) THEN
            enthalpi(Lm(ng)+1,Mm(ng)+1,linew) = 0.0_r8
            ti(Lm(ng)+1,Mm(ng)+1,linew) = 0.0_r8
          END IF
        END IF
      END IF
      RETURN
      END SUBROUTINE tibc_tile
      END MODULE tibc_mod
