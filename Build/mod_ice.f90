      MODULE mod_ice
        USE mod_kinds
        implicit none
        TYPE T_ICE
          real(r8), pointer :: ui(:,:,:)
          real(r8), pointer :: vi(:,:,:)
          real(r8), pointer :: uie(:,:,:)
          real(r8), pointer :: vie(:,:,:)
          real(r8), pointer :: ai(:,:,:)
          real(r8), pointer :: hi(:,:,:)
          real(r8), pointer :: hsn(:,:,:)
          real(r8), pointer :: sfwat(:,:,:)
          real(r8), pointer :: ageice(:,:,:)
          real(r8), pointer :: io_mflux(:,:)
          real(r8), pointer :: wfr(:,:)
          real(r8), pointer :: wai(:,:)
          real(r8), pointer :: wao(:,:)
          real(r8), pointer :: wio(:,:)
          real(r8), pointer :: wro(:,:)
          real(r8), pointer :: tis(:,:)
          real(r8), pointer :: ti(:,:,:)
          real(r8), pointer :: enthalpi(:,:,:)
          real(r8), pointer :: hage(:,:,:)
          real(r8), pointer :: utau_iw(:,:)
          real(r8), pointer :: chu_iw(:,:)
          real(r8), pointer :: spd_iw(:,:)
          real(r8), pointer :: coef_ice_heat(:,:)
          real(r8), pointer :: rhs_ice_heat(:,:)
          real(r8), pointer :: pice(:,:)
          real(r8), pointer :: zetai(:,:)
          real(r8), pointer :: eta(:,:)
          real(r8), pointer :: tauaiu(:,:)
          real(r8), pointer :: tauaiv(:,:)
          real(r8), pointer :: uwater(:,:)
          real(r8), pointer :: vwater(:,:)
          real(r8), pointer :: sealev(:,:)
          real(r8), pointer :: s0mk(:,:)
          real(r8), pointer :: t0mk(:,:)
          real(r8), pointer :: stflx_wk(:,:,:)
          real(r8), pointer :: sig11(:,:,:)
          real(r8), pointer :: sig22(:,:,:)
          real(r8), pointer :: sig12(:,:,:)
        END TYPE T_ICE
        TYPE (T_ICE), allocatable :: ICE(:)
      CONTAINS
      SUBROUTINE allocate_ice (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( ICE(Ngrids) )
!
      allocate ( ICE(ng) % ui(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % vi(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % uie(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % vie(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % ai(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % hi(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % hsn(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % sfwat(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % ageice(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % io_mflux(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % wfr(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % wai(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % wao(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % wio(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % wro(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % tis(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % ti(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % enthalpi(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % hage(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % utau_iw(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % chu_iw(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % spd_iw(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % coef_ice_heat(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % rhs_ice_heat(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % pice(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % zetai(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % eta(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % tauaiu(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % tauaiv(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % uwater(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % vwater(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % sealev(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % s0mk(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % t0mk(LBi:UBi,LBj:UBj) )
      allocate ( ICE(ng) % stflx_wk(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % sig11(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % sig22(LBi:UBi,LBj:UBj,2) )
      allocate ( ICE(ng) % sig12(LBi:UBi,LBj:UBj,2) )
      RETURN
      END SUBROUTINE allocate_ice
      SUBROUTINE initialize_ice (ng, tile)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: i, j
      integer :: Imin, Imax, Jmin, Jmax
      real(r8), parameter :: IniVal = 0.0_r8
!
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
!  Set array initialization range.
!
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          ICE(ng) % ui(i,j,1) = IniVal
          ICE(ng) % ui(i,j,2) = IniVal
          ICE(ng) % vi(i,j,1) = IniVal
          ICE(ng) % vi(i,j,2) = IniVal
          ICE(ng) % uie(i,j,1) = IniVal
          ICE(ng) % uie(i,j,2) = IniVal
          ICE(ng) % vie(i,j,1) = IniVal
          ICE(ng) % vie(i,j,2) = IniVal
          ICE(ng) % ai(i,j,1) = IniVal
          ICE(ng) % ai(i,j,2) = IniVal
          ICE(ng) % hi(i,j,1) = 0.5_r8
          ICE(ng) % hi(i,j,2) = 0.5_r8
          ICE(ng) % hsn(i,j,1) = IniVal
          ICE(ng) % hsn(i,j,2) = IniVal
          ICE(ng) % sfwat(i,j,1) = IniVal
          ICE(ng) % sfwat(i,j,2) = IniVal
          ICE(ng) % ageice(i,j,1) = IniVal
          ICE(ng) % ageice(i,j,2) = IniVal
          ICE(ng) % io_mflux(i,j) = IniVal
          ICE(ng) % wfr(i,j) = IniVal
          ICE(ng) % wai(i,j) = IniVal
          ICE(ng) % wao(i,j) = IniVal
          ICE(ng) % wio(i,j) = IniVal
          ICE(ng) % wro(i,j) = IniVal
          ICE(ng) % tis(i,j) = IniVal
          ICE(ng) % ti(i,j,1) = IniVal
          ICE(ng) % ti(i,j,2) = IniVal
          ICE(ng) % enthalpi(i,j,1) = IniVal
          ICE(ng) % enthalpi(i,j,2) = IniVal
          ICE(ng) % hage(i,j,1) = IniVal
          ICE(ng) % hage(i,j,2) = IniVal
          ICE(ng) % utau_iw(i,j) = IniVal
          ICE(ng) % chu_iw(i,j) = IniVal
          ICE(ng) % spd_iw(i,j) = IniVal
          ICE(ng) % coef_ice_heat(i,j) = IniVal
          ICE(ng) % rhs_ice_heat(i,j) = IniVal
          ICE(ng) % pice(i,j) = IniVal
          ICE(ng) % zetai(i,j) = IniVal
          ICE(ng) % eta(i,j) = IniVal
          ICE(ng) % tauaiu(i,j) = IniVal
          ICE(ng) % tauaiv(i,j) = IniVal
          ICE(ng) % uwater(i,j) = IniVal
          ICE(ng) % vwater(i,j) = IniVal
          ICE(ng) % sealev(i,j) = IniVal
          ICE(ng) % s0mk(i,j) = IniVal
          ICE(ng) % t0mk(i,j) = IniVal
          ICE(ng) % stflx_wk(i,j,1) = IniVal
          ICE(ng) % stflx_wk(i,j,2) = IniVal
          ICE(ng) % sig11(i,j,1) = IniVal
          ICE(ng) % sig11(i,j,2) = IniVal
          ICE(ng) % sig22(i,j,1) = IniVal
          ICE(ng) % sig22(i,j,2) = IniVal
          ICE(ng) % sig12(i,j,1) = IniVal
          ICE(ng) % sig12(i,j,2) = IniVal
        END DO
      END DO
      RETURN
      END SUBROUTINE initialize_ice
      END MODULE mod_ice
