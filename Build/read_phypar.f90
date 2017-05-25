      SUBROUTINE read_PhyPar (model, inp, out, Lwrite)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads and reports physical model input parameters.     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
      USE mod_strings
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      logical :: find_file, got_Ngrids, obc_data
      logical :: Lvalue(1)
      logical, allocatable :: Linert(:,:)
      logical, allocatable :: Ltracer(:,:)
      integer :: Lstr, Npts, Nval, i, itrc, ivar, k, ng, status
      integer :: ifield, ifile, igrid, itracer, nline, mFfiles
      integer :: ibcfile, iclmfile
      integer :: Ivalue(1)
      integer :: decode_line, load_i, load_l, load_lbc, load_r
      integer :: load_s1d, load_s2d
      integer, allocatable :: Nfiles(:)
      integer, allocatable :: Ncount(:,:)
      integer, allocatable :: NBCcount(:,:)
      integer, allocatable :: NCLMcount(:,:)
      real(r8), allocatable :: Rtracer(:,:)
      real(r8), allocatable :: tracer(:,:)
      real(r8) :: Rvalue(1)
      real(r8), dimension(200) :: Rval
      character (len=1  ), parameter :: blank = ' '
      character (len=19 ) :: ref_att
      character (len=40 ) :: KeyWord
      character (len=50 ) :: label
      character (len=256) :: fname, line
      character (len=256), dimension(200) :: Cval
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      ifile=1                            ! multiple file counter
      ibcfile=1                          ! multiple BC file counter
      iclmfile=1                          ! multiple BC file counter
      igrid=1                            ! nested grid counter
      itracer=0                          ! LBC tracer counter
      nline=0                            ! LBC multi-line counter
      DO i=1,LEN(label)
        label(i:i)=blank
      END DO
      got_Ngrids=.FALSE.
!
!-----------------------------------------------------------------------
!  Read in physical model parameters. Then, load input data into module.
!  Take into account nested grid configurations.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('TITLE')
              IF (Nval.eq.1) THEN
                title=TRIM(ADJUSTL(Cval(Nval)))
              ELSE
                WRITE(title,'(a,1x,a)') TRIM(ADJUSTL(title)),           &
     &                                  TRIM(ADJUSTL(Cval(Nval)))
              END IF
            CASE ('MyAppCPP')
              DO i=1,LEN(MyAppCPP)
                MyAppCPP(i:i)=blank
              END DO
              MyAppCPP=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('VARNAME')
              DO i=1,LEN(varname)
                varname(i:i)=blank
              END DO
              varname=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('Ngrids')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Ngrids=Ivalue(1)
              IF (Ngrids.le.0) THEN
                IF (Master) WRITE (out,290) 'Ngrids', Ngrids,           &
     &            'must be greater than zero.'
                exit_flag=5
                RETURN
              END IF
              got_Ngrids=.TRUE.
              CALL allocate_param       ! Start allocating variables in
              CALL allocate_parallel    ! modules that solely depend on
              CALL allocate_iounits     ! the number of nested grids
              CALL allocate_stepping
              IF (.not.allocated(Nfiles)) THEN
                allocate ( Nfiles(Ngrids) )
                Nfiles(1:Ngrids)=0
              END IF
            CASE ('Lm')
              IF (.not.got_Ngrids) THEN
                IF (Master) WRITE (out,320) 'Ngrids',                   &
     &            'Add "Ngrids" keyword before grid dimension (Lm, Mm).'
                exit_flag=5
                RETURN
              END IF
              Npts=load_i(Nval, Rval, Ngrids, Lm)
              DO ng=1,Ngrids
                IF (Lm(ng).le.0) THEN
                  IF (Master) WRITE (out,300) 'Lm', ng,                 &
     &              'must be greater than zero.'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            CASE ('Mm')
              Npts=load_i(Nval, Rval, Ngrids, Mm)
              DO ng=1,Ngrids
                IF (Mm(ng).le.0) THEN
                  IF (Master) WRITE (out,300) 'Mm', ng,                 &
     &              'must be greater than zero.'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            CASE ('N')
              Npts=load_i(Nval, Rval, Ngrids, N)
              DO ng=1,Ngrids
                IF (N(ng).lt.0) THEN
                  IF (Master) WRITE (out,300) 'N', ng,                  &
     &              'must be greater than zero.'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            CASE ('NAT')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              NAT=Ivalue(1)
              IF ((NAT.lt.1).or.(NAT.gt.2)) THEN
                IF (Master) WRITE (out,290) 'NAT = ', NAT,              &
     &            'make sure that NAT is either 1 or 2.'
                exit_flag=5
                RETURN
              END IF
              IF (NAT.ne.2) THEN
                IF (Master) WRITE (out,290) 'NAT = ', NAT,              &
     &            'make sure that NAT is equal to 2.'
                exit_flag=5
                RETURN
              END IF
            CASE ('NPT')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              NPT=Ivalue(1)
              IF (NPT.le.0) THEN
                IF (Master) WRITE (out,290) 'NPT = ', NPT,              &
     &            'must be greater than zero.'
                exit_flag=5
                RETURN
              END IF
            CASE ('NtileI')
              Npts=load_i(Nval, Rval, Ngrids, NtileI)
              NtileX(1:Ngrids)=NtileI(1:Ngrids)
            CASE ('NtileJ')
              Npts=load_i(Nval, Rval, Ngrids, NtileJ)
              NtileE(1:Ngrids)=NtileJ(1:Ngrids)
              CALL initialize_param    ! Continue allocating/initalizing
              CALL allocate_scalars    ! variables since the application
              CALL initialize_scalars  ! number of nested grids and
              CALL allocate_ncparam    ! domain parameters are known
              CALL initialize_ncparam
              IF (.not.allocated(Ltracer)) THEN
                allocate (Ltracer(NAT+NPT,Ngrids))
              END IF
              IF (.not.allocated(Linert)) THEN
                allocate (Linert(NPT,Ngrids))
              END IF
              IF (MAXVAL(inert).eq.0) THEN
                IF (Master) WRITE (out,280) 'inert'
                exit_flag=5
                RETURN
              END IF
              IF (.not.allocated(Rtracer)) THEN
                allocate (Rtracer(NAT+NPT,Ngrids))
              END IF
              IF (.not.allocated(tracer)) THEN
                allocate (tracer(MT,Ngrids))
              END IF
            CASE ('LBC(isFsur)')
              Npts=load_lbc(Nval, Cval, line, nline, isFsur, igrid,     &
     &                      0, 0, LBC)
            CASE ('LBC(isUbar)')
              Npts=load_lbc(Nval, Cval, line, nline, isUbar, igrid,     &
     &                      0, 0, LBC)
            CASE ('LBC(isVbar)')
              Npts=load_lbc(Nval, Cval, line, nline, isVbar, igrid,     &
     &                      0, 0, LBC)
            CASE ('LBC(isUvel)')
              Npts=load_lbc(Nval, Cval, line, nline, isUvel, igrid,     &
     &                      0, 0, LBC)
            CASE ('LBC(isVvel)')
              Npts=load_lbc(Nval, Cval, line, nline, isVvel, igrid,     &
     &                      0, 0, LBC)
            CASE ('LBC(isUice)')
              Npts=load_lbc(Nval, Cval, line, nline, isUice, igrid,     &
     &                      0, 0, LBC)
            CASE ('LBC(isVice)')
              Npts=load_lbc(Nval, Cval, line, nline, isVice, igrid,     &
     &                      0, 0, LBC)
            CASE ('LBC(isAice)')
              Npts=load_lbc(Nval, Cval, line, nline, isAice, igrid,     &
     &                      0, 0, LBC)
            CASE ('LBC(isHice)')
              Npts=load_lbc(Nval, Cval, line, nline, isHice, igrid,     &
     &                      0, 0, LBC)
            CASE ('LBC(isTice)')
              Npts=load_lbc(Nval, Cval, line, nline, isTice, igrid,     &
     &                      0, 0, LBC)
            CASE ('LBC(isHsno)')
              Npts=load_lbc(Nval, Cval, line, nline, isHsno, igrid,     &
     &                      0, 0, LBC)
            CASE ('LBC(isSfwat)')
              Npts=load_lbc(Nval, Cval, line, nline, isSfwat, igrid,    &
     &                      0, 0, LBC)
            CASE ('LBC(isSig11)')
              Npts=load_lbc(Nval, Cval, line, nline, isSig11, igrid,    &
     &                      0, 0, LBC)
            CASE ('LBC(isSig12)')
              Npts=load_lbc(Nval, Cval, line, nline, isSig12, igrid,    &
     &                      0, 0, LBC)
            CASE ('LBC(isSig22)')
              Npts=load_lbc(Nval, Cval, line, nline, isSig22, igrid,    &
     &                      0, 0, LBC)
            CASE ('LBC(isTvar)')
              IF (itracer.lt.(NAT+NPT)) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              ifield=isTvar(itracer)
              Npts=load_lbc(Nval, Cval, line, nline, ifield, igrid,     &
     &                      1, NAT+NPT, LBC)
            CASE ('VolCons(west)')
              Npts=load_l(Nval, Cval, Ngrids, VolCons(iwest,:))
            CASE ('VolCons(east)')
              Npts=load_l(Nval, Cval, Ngrids, VolCons(ieast,:))
            CASE ('VolCons(south)')
              Npts=load_l(Nval, Cval, Ngrids, VolCons(isouth,:))
            CASE ('VolCons(north)')
              Npts=load_l(Nval, Cval, Ngrids, VolCons(inorth,:))
            CASE ('NTIMES')
              Npts=load_i(Nval, Rval, Ngrids, ntimes)
            CASE ('DT')
              Npts=load_r(Nval, Rval, Ngrids, dt)
            CASE ('NDTFAST')
              Npts=load_i(Nval, Rval, Ngrids, ndtfast)
            CASE ('ERstr')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              ERstr=Ivalue(1)
            CASE ('ERend')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              ERend=Ivalue(1)
            CASE ('Nouter')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Nouter=Ivalue(1)
            CASE ('Ninner')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Ninner=Ivalue(1)
            CASE ('Nintervals')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Nintervals=Ivalue(1)
            CASE ('NRREC')
              Npts=load_i(Nval, Rval, Ngrids, nrrec)
              DO ng=1,Ngrids
                IF (nrrec(ng).lt.0) THEN
                  LastRec(ng)=.TRUE.
                ELSE
                  LastRec(ng)=.FALSE.
                END IF
              END DO
            CASE ('LcycleRST')
              Npts=load_l(Nval, Cval, Ngrids, LcycleRST)
            CASE ('NRST')
              Npts=load_i(Nval, Rval, Ngrids, nRST)
            CASE ('NSTA')
              Npts=load_i(Nval, Rval, Ngrids, nSTA)
            CASE ('NFLT')
              Npts=load_i(Nval, Rval, Ngrids, nFLT)
            CASE ('NINFO')
              Npts=load_i(Nval, Rval, Ngrids, ninfo)
            CASE ('LDEFOUT')
              Npts=load_l(Nval, Cval, Ngrids, ldefout)
            CASE ('NHIS')
              Npts=load_i(Nval, Rval, Ngrids, nHIS)
            CASE ('NDEFHIS')
              Npts=load_i(Nval, Rval, Ngrids, ndefHIS)
            CASE ('NTSAVG')
              Npts=load_i(Nval, Rval, Ngrids, ntsAVG)
            CASE ('NAVG')
              Npts=load_i(Nval, Rval, Ngrids, nAVG)
            CASE ('NDEFAVG')
              Npts=load_i(Nval, Rval, Ngrids, ndefAVG)
            CASE ('NTSAVG2')
              Npts=load_i(Nval, Rval, Ngrids, ntsAVG2)
            CASE ('NAVG2')
              Npts=load_i(Nval, Rval, Ngrids, nAVG2)
            CASE ('NDEFAVG2')
              Npts=load_i(Nval, Rval, Ngrids, ndefAVG2)
            CASE ('NTSDIA')
              Npts=load_i(Nval, Rval, Ngrids, ntsDIA)
            CASE ('NDIA')
              Npts=load_i(Nval, Rval, Ngrids, nDIA)
            CASE ('NDEFDIA')
              Npts=load_i(Nval, Rval, Ngrids, ndefDIA)
            CASE ('LcycleTLM')
              Npts=load_l(Nval, Cval, Ngrids, LcycleTLM)
            CASE ('NTLM')
              Npts=load_i(Nval, Rval, Ngrids, nTLM)
            CASE ('NDEFTLM')
              Npts=load_i(Nval, Rval, Ngrids, ndefTLM)
            CASE ('LcycleADJ')
              Npts=load_l(Nval, Cval, Ngrids, LcycleADJ)
            CASE ('NADJ')
              Npts=load_i(Nval, Rval, Ngrids, nADJ)
            CASE ('NDEFADJ')
              Npts=load_i(Nval, Rval, Ngrids, ndefADJ)
            CASE ('NOBC')
              Npts=load_i(Nval, Rval, Ngrids, nOBC)
            CASE ('NSFF')
              Npts=load_i(Nval, Rval, Ngrids, nSFF)
            CASE ('LmultiGST')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              LmultiGST=Lvalue(1)
            CASE ('LrstGST')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              LrstGST=Lvalue(1)
            CASE ('MaxIterGST')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              MaxIterGST=Ivalue(1)
            CASE ('NGST')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              nGST=Ivalue(1)
            CASE ('TNU2')
              Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  nl_tnu2(itrc,ng)=Rtracer(itrc,ng)
                END DO
                DO i=1,NPT
                  itrc=inert(i)
                  nl_tnu2(itrc,ng)=Rtracer(NAT+i,ng)
                END DO
              END DO
            CASE ('TNU4')
              Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  nl_tnu4(itrc,ng)=Rtracer(itrc,ng)
                END DO
                DO i=1,NPT
                  itrc=inert(i)
                  nl_tnu4(itrc,ng)=Rtracer(NAT+i,ng)
                END DO
              END DO
            CASE ('ad_TNU2')
              Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  ad_tnu2(itrc,ng)=Rtracer(itrc,ng)
                  tl_tnu2(itrc,ng)=Rtracer(itrc,ng)
                END DO
                DO i=1,NPT
                  itrc=inert(i)
                  ad_tnu2(itrc,ng)=Rtracer(NAT+i,ng)
                  tl_tnu2(itrc,ng)=Rtracer(NAT+i,ng)
                END DO
              END DO
            CASE ('ad_TNU4')
              Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  ad_tnu4(itrc,ng)=Rtracer(itrc,ng)
                  tl_tnu4(itrc,ng)=Rtracer(itrc,ng)
                END DO
                DO i=1,NPT
                  itrc=inert(i)
                  ad_tnu4(itrc,ng)=Rtracer(NAT+i,ng)
                  tl_tnu4(itrc,ng)=Rtracer(NAT+i,ng)
                END DO
              END DO
            CASE ('VISC2')
              Npts=load_r(Nval, Rval, Ngrids, nl_visc2)
            CASE ('VISC4')
              Npts=load_r(Nval, Rval, Ngrids, nl_visc4)
            CASE ('ad_VISC2')
              Npts=load_r(Nval, Rval, Ngrids, ad_visc2)
              DO ng=1,Ngrids
                tl_visc2(ng)=ad_visc2(ng)
              END DO
            CASE ('ad_VISC4')
              Npts=load_r(Nval, Rval, Ngrids, ad_visc4)
              DO ng=1,Ngrids
                tl_visc4(ng)=ad_visc4(ng)
              END DO
            CASE ('AKT_BAK')
              Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Akt_bak(itrc,ng)=Rtracer(itrc,ng)
                END DO
                DO i=1,NPT
                  itrc=inert(i)
                  Akt_bak(itrc,ng)=Rtracer(NAT+i,ng)
                END DO
              END DO
            CASE ('ad_AKT_fac')
              Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  ad_Akt_fac(itrc,ng)=Rtracer(itrc,ng)
                  tl_Akt_fac(itrc,ng)=Rtracer(itrc,ng)
                END DO
                DO i=1,NPT
                  itrc=inert(i)
                  ad_Akt_fac(itrc,ng)=Rtracer(NAT+i,ng)
                  tl_Akt_fac(itrc,ng)=Rtracer(NAT+i,ng)
                END DO
              END DO
            CASE ('AKV_BAK')
              Npts=load_r(Nval, Rval, Ngrids, Akv_bak)
            CASE ('ad_AKV_fac')
              Npts=load_r(Nval, Rval, Ngrids, ad_Akv_fac)
              DO ng=1,Ngrids
                tl_Akv_fac(ng)=ad_AKv_fac(ng)
              END DO
            CASE ('AKK_BAK')
              Npts=load_r(Nval, Rval, Ngrids, Akk_bak)
            CASE ('AKP_BAK')
              Npts=load_r(Nval, Rval, Ngrids, Akp_bak)
            CASE ('TKENU2')
              Npts=load_r(Nval, Rval, Ngrids, tkenu2)
            CASE ('TKENU4')
              Npts=load_r(Nval, Rval, Ngrids, tkenu4)
            CASE ('GLS_P')
              Npts=load_r(Nval, Rval, Ngrids, gls_p)
            CASE ('GLS_M')
              Npts=load_r(Nval, Rval, Ngrids, gls_m)
            CASE ('GLS_N')
              Npts=load_r(Nval, Rval, Ngrids, gls_n)
            CASE ('GLS_Kmin')
              Npts=load_r(Nval, Rval, Ngrids, gls_Kmin)
            CASE ('GLS_Pmin')
              Npts=load_r(Nval, Rval, Ngrids, gls_Pmin)
            CASE ('GLS_CMU0')
              Npts=load_r(Nval, Rval, Ngrids, gls_cmu0)
            CASE ('GLS_C1')
              Npts=load_r(Nval, Rval, Ngrids, gls_c1)
            CASE ('GLS_C2')
              Npts=load_r(Nval, Rval, Ngrids, gls_c2)
            CASE ('GLS_C3M')
              Npts=load_r(Nval, Rval, Ngrids, gls_c3m)
            CASE ('GLS_C3P')
              Npts=load_r(Nval, Rval, Ngrids, gls_c3p)
            CASE ('GLS_SIGK')
              Npts=load_r(Nval, Rval, Ngrids, gls_sigk)
            CASE ('GLS_SIGP')
              Npts=load_r(Nval, Rval, Ngrids, gls_sigp)
            CASE ('CHARNOK_ALPHA')
              Npts=load_r(Nval, Rval, Ngrids, charnok_alpha)
            CASE ('ZOS_HSIG_ALPHA')
              Npts=load_r(Nval, Rval, Ngrids, zos_hsig_alpha)
            CASE ('SZ_ALPHA')
              Npts=load_r(Nval, Rval, Ngrids, sz_alpha)
            CASE ('CRGBAN_CW')
              Npts=load_r(Nval, Rval, Ngrids, crgban_cw)
            CASE ('RDRG')
              Npts=load_r(Nval, Rval, Ngrids, rdrg)
            CASE ('RDRG2')
              Npts=load_r(Nval, Rval, Ngrids, rdrg2)
            CASE ('Zob')
              Npts=load_r(Nval, Rval, Ngrids, Zob)
            CASE ('Zos')
              Npts=load_r(Nval, Rval, Ngrids, Zos)
            CASE ('BLK_ZQ')
              Npts=load_r(Nval, Rval, Ngrids, blk_ZQ)
            CASE ('BLK_ZT')
              Npts=load_r(Nval, Rval, Ngrids, blk_ZT)
            CASE ('BLK_ZW')
              Npts=load_r(Nval, Rval, Ngrids, blk_ZW)
            CASE ('DCRIT')
              Npts=load_r(Nval, Rval, Ngrids, Dcrit)
            CASE ('WTYPE')
              Npts=load_i(Nval, Rval, Ngrids, lmd_Jwt)
            CASE ('LEVSFRC')
              Npts=load_i(Nval, Rval, Ngrids, levsfrc)
            CASE ('LEVBFRC')
              Npts=load_i(Nval, Rval, Ngrids, levbfrc)
            CASE ('Vtransform')
              Npts=load_i(Nval, Rval, Ngrids, Vtransform)
              DO ng=1,Ngrids
                IF ((Vtransform(ng).lt.0).or.                           &
     &              (Vtransform(ng).gt.2)) THEN
                  IF (Master) WRITE (out,260) 'Vtransform = ',          &
     &                                        Vtransform(ng),           &
     &                                        'Must be either 1 or 2'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            CASE ('Vstretching')
              Npts=load_i(Nval, Rval, Ngrids, Vstretching)
              DO ng=1,Ngrids
                IF ((Vstretching(ng).lt.0).or.                          &
     &              (Vstretching(ng).gt.4)) THEN
                  IF (Master) WRITE (out,260) 'Vstretching = ',         &
     &                                        Vstretching(ng),          &
     &                                        'Must between 1 and 4'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            CASE ('THETA_S')
              Npts=load_r(Nval, Rval, Ngrids, theta_s)
            CASE ('THETA_B')
              Npts=load_r(Nval, Rval, Ngrids, theta_b)
            CASE ('TCLINE')
              Npts=load_r(Nval, Rval, Ngrids, Tcline)
              DO ng=1,Ngrids
                hc(ng)=Tcline(ng)
              END DO
            CASE ('RHO0')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              rho0=Rvalue(1)
            CASE ('BVF_BAK')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              bvf_bak=Rvalue(1)
            CASE ('DSTART')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              dstart=Rvalue(1)
            CASE ('TIDE_START')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              tide_start=Rvalue(1)
            CASE ('TIME_REF')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              time_ref=Rvalue(1)
              r_text=ref_att(time_ref,r_date)
            CASE ('TNUDG')
              Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Tnudg(itrc,ng)=Rtracer(itrc,ng)
                END DO
                DO i=1,NPT
                  itrc=inert(i)
                  Tnudg(itrc,ng)=Rtracer(NAT+i,ng)
                END DO
              END DO
            CASE ('ZNUDG')
              Npts=load_r(Nval, Rval, Ngrids, Znudg)
            CASE ('M2NUDG')
              Npts=load_r(Nval, Rval, Ngrids, M2nudg)
            CASE ('M3NUDG')
              Npts=load_r(Nval, Rval, Ngrids, M3nudg)
            CASE ('OBCFAC')
              Npts=load_r(Nval, Rval, Ngrids, obcfac)
            CASE ('R0')
              Npts=load_r(Nval, Rval, Ngrids, R0)
              DO ng=1,Ngrids
                IF (R0(ng).lt.100.0_r8) R0(ng)=R0(ng)+1000.0_r8
              END DO
            CASE ('T0')
              Npts=load_r(Nval, Rval, Ngrids, T0)
            CASE ('S0')
              Npts=load_r(Nval, Rval, Ngrids, S0)
            CASE ('TCOEF')
              Npts=load_r(Nval, Rval, Ngrids, Tcoef)
              DO ng=1,Ngrids
                Tcoef(ng)=ABS(Tcoef(ng))
              END DO
            CASE ('SCOEF')
              Npts=load_r(Nval, Rval, Ngrids, Scoef)
              DO ng=1,Ngrids
                Scoef(ng)=ABS(Scoef(ng))
              END DO
            CASE ('GAMMA2')
              Npts=load_r(Nval, Rval, Ngrids, gamma2)
            CASE ('Hout(idUvel)')
              IF (idUvel.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUvel,:))
            CASE ('Hout(idVvel)')
              IF (idVvel.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVvel,:))
            CASE ('Hout(idWvel)')
              IF (idWvel.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWvel,:))
            CASE ('Hout(idOvel)')
              IF (idOvel.eq.0) THEN
                IF (Master) WRITE (out,280) 'idOvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idOvel,:))
            CASE ('Hout(idUbar)')
              IF (idUbar.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUbar,:))
            CASE ('Hout(idVbar)')
              IF (idVbar.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVbar,:))
            CASE ('Hout(idFsur)')
              IF (idFsur.eq.0) THEN
                IF (Master) WRITE (out,280) 'idFsur'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idFsur,:))
            CASE ('Hout(idu2dE)')
              IF (idu2dE.eq.0) THEN
                IF (Master) WRITE (out,280) 'idu2dE'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idu2dE,:))
            CASE ('Hout(idv2dN)')
              IF (idv2dN.eq.0) THEN
                IF (Master) WRITE (out,280) 'idv2dN'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idv2dN,:))
            CASE ('Hout(idu3dE)')
              IF (idu3dE.eq.0) THEN
                IF (Master) WRITE (out,280) 'idu3dE'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idu3dE,:))
            CASE ('Hout(idv3dN)')
              IF (idv3dN.eq.0) THEN
                IF (Master) WRITE (out,280) 'idv3dN'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idv3dN,:))
            CASE ('Hout(idTvar)')
              IF (MAXVAL(idTvar).eq.0) THEN
                IF (Master) WRITE (out,280) 'idTvar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTvar(itrc)
                  Hout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idUsms)')
              IF (idUsms.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUsms'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUsms,:))
            CASE ('Hout(idVsms)')
              IF (idVsms.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVsms'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVsms,:))
            CASE ('Hout(idUbms)')
              IF (idUbms.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbms'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUbms,:))
            CASE ('Hout(idVbms)')
              IF (idVbms.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbms'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVbms,:))
            CASE ('Hout(idUbrs)')
              IF (idUbrs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbrs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUbrs,:))
            CASE ('Hout(idVbrs)')
              IF (idVbrs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbrs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVbrs,:))
            CASE ('Hout(idUbws)')
              IF (idUbws.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbws'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUbws,:))
            CASE ('Hout(idVbws)')
              IF (idVbws.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbws'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVbws,:))
            CASE ('Hout(idUbcs)')
              IF (idUbcs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbcs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUbcs,:))
            CASE ('Hout(idVbcs)')
              IF (idVbcs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbcs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVbcs,:))
            CASE ('Hout(idUbot)')
              IF (idUbot.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbot'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUbot,:))
            CASE ('Hout(idVbot)')
              IF (idVbot.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbot'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVbot,:))
            CASE ('Hout(idUbur)')
              IF (idUbur.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbur'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUbur,:))
            CASE ('Hout(idVbvr)')
              IF (idVbvr.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbvr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVbvr,:))
            CASE ('Hout(idW2xx)')
              IF (idW2xx.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW2xx'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW2xx,:))
            CASE ('Hout(idW2xy)')
              IF (idW2xy.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW2xy'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW2xy,:))
            CASE ('Hout(idW2yy)')
              IF (idW2yy.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW2yy'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW2yy,:))
            CASE ('Hout(idU2rs)')
              IF (idU2rs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idU2rs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idU2rs,:))
            CASE ('Hout(idV2rs)')
              IF (idV2rs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idV2rs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idV2rs,:))
            CASE ('Hout(idU2Sd)')
              IF (idU2Sd.eq.0) THEN
                IF (Master) WRITE (out,280) 'idU2Sd'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idU2Sd,:))
            CASE ('Hout(idV2Sd)')
              IF (idV2Sd.eq.0) THEN
                IF (Master) WRITE (out,280) 'idV2Sd'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idV2Sd,:))
            CASE ('Hout(idW3xx)')
              IF (idW3xx.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW3xx'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW3xx,:))
            CASE ('Hout(idW3xy)')
              IF (idW3xy.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW3xy'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW3xy,:))
            CASE ('Hout(idW3yy)')
              IF (idW3yy.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW3yy'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW3yy,:))
            CASE ('Hout(idW3zx)')
              IF (idW3zx.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW3zx'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW3zx,:))
            CASE ('Hout(idW3zy)')
              IF (idW3zy.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW3zy'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW3zy,:))
            CASE ('Hout(idU3rs)')
              IF (idU3rs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idU3rs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idU3rs,:))
            CASE ('Hout(idV3rs)')
              IF (idV3rs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idV3rs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idV3rs,:))
            CASE ('Hout(idU3Sd)')
              IF (idU3Sd.eq.0) THEN
                IF (Master) WRITE (out,280) 'idU3Sd'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idU3Sd,:))
            CASE ('Hout(idV3Sd)')
              IF (idV3Sd.eq.0) THEN
                IF (Master) WRITE (out,280) 'idV3Sd'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idV3Sd,:))
            CASE ('Hout(idWamp)')
              IF (idWamp.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWamp'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWamp,:))
            CASE ('Hout(idWlen)')
              IF (idWlen.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWlen'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWlen,:))
            CASE ('Hout(idWdir)')
              IF (idWdir.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWdir'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWdir,:))
            CASE ('Hout(idPair)')
              IF (idPair.eq.0) THEN
                IF (Master) WRITE (out,280) 'idPair'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idPair,:))
            CASE ('Hout(idUair)')
              IF (idUair.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUair'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUair,:))
            CASE ('Hout(idVair)')
              IF (idVair.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVair'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVair,:))
            CASE ('Hout(idUairE)')
              IF (idUairE.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUairE'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUairE,:))
            CASE ('Hout(idVairN)')
              IF (idVairN.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVairN'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVairN,:))
            CASE ('Hout(idTsur)')
              IF (idTsur(itemp).eq.0) THEN
                IF (Master) WRITE (out,280) 'idTsur'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTsur(itrc)
                  Hout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idLhea)')
              IF (idLhea.eq.0) THEN
                IF (Master) WRITE (out,280) 'idLhea'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idLhea,:))
            CASE ('Hout(idShea)')
              IF (idShea.eq.0) THEN
                IF (Master) WRITE (out,280) 'idShea'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idShea,:))
            CASE ('Hout(idLrad)')
              IF (idLrad.eq.0) THEN
                IF (Master) WRITE (out,280) 'idLrad'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idLrad,:))
            CASE ('Hout(idSrad)')
              IF (idSrad.eq.0) THEN
                IF (Master) WRITE (out,280) 'idSrad'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idSrad,:))
            CASE ('Hout(idEmPf)')
              IF (idEmPf.eq.0) THEN
                IF (Master) WRITE (out,280) 'idEmPf'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idEmPf,:))
            CASE ('Hout(idevap)')
              IF (idevap.eq.0) THEN
                IF (Master) WRITE (out,280) 'idevap'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idevap,:))
            CASE ('Hout(idrain)')
              IF (idrain.eq.0) THEN
                IF (Master) WRITE (out,280) 'idrain'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idrain,:))
            CASE ('Hout(idDano)')
              IF (idDano.eq.0) THEN
                IF (Master) WRITE (out,280) 'idDano'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idDano,:))
            CASE ('Hout(idVvis)')
              IF (idVvis.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVvis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVvis,:))
            CASE ('Hout(idTdif)')
              IF (idTdif.eq.0) THEN
                IF (Master) WRITE (out,280) 'idTdif'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTdif,:))
            CASE ('Hout(idSdif)')
              IF (idSdif.eq.0) THEN
                IF (Master) WRITE (out,280) 'idSdif'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idSdif,:))
            CASE ('Hout(idHsbl)')
              IF (idHsbl.eq.0) THEN
                IF (Master) WRITE (out,280) 'idHsbl'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idHsbl,:))
            CASE ('Hout(idHbbl)')
              IF (idHbbl.eq.0) THEN
                IF (Master) WRITE (out,280) 'idHbbl'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idHbbl,:))
            CASE ('Hout(idMtke)')
              IF (idMtke.eq.0) THEN
                IF (Master) WRITE (out,280) 'idMtke'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idMtke,:))
            CASE ('Hout(idMtls)')
              IF (idMtls.eq.0) THEN
                IF (Master) WRITE (out,280) 'idMtls'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idMtls,:))
            CASE ('Hout(inert)')
              Npts=load_l(Nval, Cval, NPT*Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=idTvar(inert(i))
                  Hout(itrc,ng)=Linert(i,ng)
                END DO
              END DO
            CASE ('Aout(idUvel)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUvel,:))
            CASE ('Aout(idVvel)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVvel,:))
            CASE ('Aout(idWvel)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWvel,:))
            CASE ('Aout(idOvel)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idOvel,:))
            CASE ('Aout(idUbar)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUbar,:))
            CASE ('Aout(idVbar)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVbar,:))
            CASE ('Aout(idFsur)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idFsur,:))
            CASE ('Aout(idu2dE)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idu2dE,:))
            CASE ('Aout(idv2dN)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idv2dN,:))
            CASE ('Aout(idu3dE)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idu3dE,:))
            CASE ('Aout(idv3dN)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idv3dN,:))
            CASE ('Aout(idTvar)')
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTvar(itrc)
                  Aout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Aout(idUsms)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUsms,:))
            CASE ('Aout(idVsms)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVsms,:))
            CASE ('Aout(idUbms)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUbms,:))
            CASE ('Aout(idVbms)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVbms,:))
            CASE ('Aout(idW2xx)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW2xx,:))
            CASE ('Aout(idW2xy)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW2xy,:))
            CASE ('Aout(idW2yy)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW2yy,:))
            CASE ('Aout(idU2rs)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idU2rs,:))
            CASE ('Aout(idV2rs)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idV2rs,:))
            CASE ('Aout(idU2Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idU2Sd,:))
            CASE ('Aout(idV2Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idV2Sd,:))
            CASE ('Aout(idW3xx)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW3xx,:))
            CASE ('Aout(idW3xy)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW3xy,:))
            CASE ('Aout(idW3yy)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW3yy,:))
            CASE ('Aout(idW3zx)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW3zx,:))
            CASE ('Aout(idW3zy)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW3zy,:))
            CASE ('Aout(idU3rs)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idU3rs,:))
            CASE ('Aout(idV3rs)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idV3rs,:))
            CASE ('Aout(idU3Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idU3Sd,:))
            CASE ('Aout(idV3Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idV3Sd,:))
            CASE ('Aout(idPair)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idPair,:))
            CASE ('Aout(idUair)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUair,:))
            CASE ('Aout(idVair)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVair,:))
            CASE ('Aout(idUairE)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUairE,:))
            CASE ('Aout(idVairN)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVairN,:))
            CASE ('Aout(idTsur)')
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTsur(itrc)
                  Aout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Aout(idLhea)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idLhea,:))
            CASE ('Aout(idShea)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idShea,:))
            CASE ('Aout(idLrad)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idLrad,:))
            CASE ('Aout(idSrad)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idSrad,:))
            CASE ('Aout(idevap)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idevap,:))
            CASE ('Aout(idrain)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idrain,:))
            CASE ('Aout(idDano)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idDano,:))
            CASE ('Aout(idVvis)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVvis,:))
            CASE ('Aout(idTdif)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idTdif,:))
            CASE ('Aout(idSdif)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idSdif,:))
            CASE ('Aout(idHsbl)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idHsbl,:))
            CASE ('Aout(idHbbl)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idHbbl,:))
            CASE ('Aout(id2dRV)')
              IF (id2dRV.eq.0) THEN
                IF (Master) WRITE (out,280) 'id2dRV'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(id2dRV,:))
            CASE ('Aout(id3dRV)')
              IF (id3dRV.eq.0) THEN
                IF (Master) WRITE (out,280) 'id3dRV'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(id3dRV,:))
            CASE ('Aout(id2dPV)')
              IF (id2dPV.eq.0) THEN
                IF (Master) WRITE (out,280) 'id2dPV'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(id2dPV,:))
            CASE ('Aout(id3dPV)')
              IF (id3dPV.eq.0) THEN
                IF (Master) WRITE (out,280) 'id3dPV'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(id3dPV,:))
            CASE ('Aout(idHUav)')
              IF (idHUav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idHUav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idHUav,:))
            CASE ('Aout(idHVav)')
              IF (idHVav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idHVav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idHVav,:))
            CASE ('Aout(idUUav)')
              IF (idUUav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUUav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUUav,:))
            CASE ('Aout(idUVav)')
              IF (idUVav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUVav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUVav,:))
            CASE ('Aout(idVVav)')
              IF (idVVav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVVav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVVav,:))
            CASE ('Aout(idU2av)')
              IF (idU2av.eq.0) THEN
                IF (Master) WRITE (out,280) 'idU2av'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idU2av,:))
            CASE ('Aout(idV2av)')
              IF (idV2av.eq.0) THEN
                IF (Master) WRITE (out,280) 'idV2av'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idV2av,:))
            CASE ('Aout(idZZav)')
              IF (idZZav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idZZav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idZZav,:))
            CASE ('Aout(idTTav)')
              IF (MAXVAL(idTTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'idTTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTTav(itrc)
                  Aout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Aout(idUTav)')
              IF (MAXVAL(idUTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'idUTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idUTav(itrc)
                  Aout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Aout(idVTav)')
              IF (MAXVAL(idVTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'idVTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idVTav(itrc)
                  Aout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iHUTav)')
              IF (MAXVAL(iHUTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'iHUTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=iHUTav(itrc)
                  Aout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iHVTav)')
              IF (MAXVAL(iHVTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'iHVTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=iHVTav(itrc)
                  Aout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Aout(inert)')
              Npts=load_l(Nval, Cval, NPT*Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=idTvar(inert(i))
                  Aout(itrc,ng)=Linert(i,ng)
                END DO
              END DO
            CASE ('Dout(iTrate)')
              IF (iTrate.le.0) THEN
                IF (Master) WRITE (out,280) 'iTrate'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iTrate),ng)=Ltracer(itrc,ng)
                END DO
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iTrate),ng)=Ltracer(NAT+i,ng)
                END DO
              END DO
            CASE ('Dout(iThadv)')
              IF (iThadv.le.0) THEN
                IF (Master) WRITE (out,280) 'iThadv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iThadv),ng)=Ltracer(itrc,ng)
                END DO
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iThadv),ng)=Ltracer(NAT+i,ng)
                END DO
              END DO
            CASE ('Dout(iTxadv)')
              IF (iTxadv.le.0) THEN
                IF (Master) WRITE (out,280) 'iTxadv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iTxadv),ng)=Ltracer(itrc,ng)
                END DO
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iTxadv),ng)=Ltracer(NAT+i,ng)
                END DO
              END DO
            CASE ('Dout(iTyadv)')
              IF (iTyadv.le.0) THEN
                IF (Master) WRITE (out,280) 'iTyadv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iTyadv),ng)=Ltracer(itrc,ng)
                END DO
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iTyadv),ng)=Ltracer(NAT+i,ng)
                END DO
              END DO
            CASE ('Dout(iTvadv)')
              IF (iTvadv.le.0) THEN
                IF (Master) WRITE (out,280) 'iTvadv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iTvadv),ng)=Ltracer(itrc,ng)
                END DO
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iTvadv),ng)=Ltracer(NAT+i,ng)
                END DO
              END DO
            CASE ('Dout(iThdif)')
              IF (iThdif.le.0) THEN
                IF (Master) WRITE (out,280) 'iThdif'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iThdif),ng)=Ltracer(itrc,ng)
                END DO
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iThdif),ng)=Ltracer(NAT+i,ng)
                END DO
              END DO
            CASE ('Dout(iTxdif)')
              IF (iTxdif.le.0) THEN
                IF (Master) WRITE (out,280) 'iTxdif'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iTxdif),ng)=Ltracer(itrc,ng)
                END DO
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iTxdif),ng)=Ltracer(NAT+i,ng)
                END DO
              END DO
            CASE ('Dout(iTydif)')
              IF (iTydif.le.0) THEN
                IF (Master) WRITE (out,280) 'iTydif'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iTydif),ng)=Ltracer(itrc,ng)
                END DO
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iTydif),ng)=Ltracer(NAT+i,ng)
                END DO
              END DO
            CASE ('Dout(iTsdif)')
              IF (iTsdif.le.0) THEN
                IF (Master) WRITE (out,280) 'iTsdif'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iTsdif),ng)=Ltracer(itrc,ng)
                END DO
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iTsdif),ng)=Ltracer(NAT+i,ng)
                END DO
              END DO
            CASE ('Dout(iTvdif)')
              IF (iTvdif.le.0) THEN
                IF (Master) WRITE (out,280) 'iTvdif'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iTvdif),ng)=Ltracer(itrc,ng)
                END DO
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iTvdif),ng)=Ltracer(NAT+i,ng)
                END DO
              END DO
            CASE ('NUSER')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Nuser=Ivalue(1)
            CASE ('USER')
              Npts=load_r(Nval, Rval, MAX(1,Nuser), user)
            CASE ('NC_SHUFFLE')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              shuffle=Ivalue(1)
            CASE ('NC_DEFLATE')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              deflate=Ivalue(1)
            CASE ('NC_DLEVEL')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              deflate_level=Ivalue(1)
            CASE ('GSTNAME')
              label='GST - generalized stability theory analysis'
              Npts=load_s1d(Nval, Cval, line, label, igrid, Nfiles, GST)
            CASE ('RSTNAME')
              label='RST - restart fields'
              Npts=load_s1d(Nval, Cval, line, label, igrid, Nfiles, RST)
            CASE ('HISNAME')
              label='HIS - nonlinear model history fields'
              Npts=load_s1d(Nval, Cval, line, label, igrid, Nfiles, HIS)
            CASE ('TLMNAME')
              label='TLM - tangent linear model history fields'
              Npts=load_s1d(Nval, Cval, line, label, igrid, Nfiles, TLM)
            CASE ('TLFNAME')
              label='TLF - tangent linear model impulse forcing'
              Npts=load_s1d(Nval, Cval, line, label, igrid, Nfiles, TLF)
            CASE ('ADJNAME')
              label='ADM - adjoint model history fields'
              Npts=load_s1d(Nval, Cval, line, label, igrid, Nfiles, ADM)
            CASE ('AVGNAME')
              label='AVG - time-averaged history fields'
              Npts=load_s1d(Nval, Cval, line, label, igrid, Nfiles, AVG)
            CASE ('AVG2NAME')
              label='AVG2 - time-averaged history fields'
              Npts=load_s1d(Nval, Cval, line, label, igrid, Nfiles,     &
     &                     AVG2)
            CASE ('DIANAME')
              label='DIA - time-averaged diagnostics fields'
              Npts=load_s1d(Nval, Cval, line, label, igrid, Nfiles, DIA)
            CASE ('STANAME')
              label='STA - stations time-series'
              Npts=load_s1d(Nval, Cval, line, label, igrid, Nfiles, STA)
            CASE ('FLTNAME')
              label='FLT - Lagragian particles trajectories'
              Npts=load_s1d(Nval, Cval, line, label, igrid, Nfiles, FLT)
            CASE ('GRDNAME')
              label='GRD - application grid'
              Npts=load_s1d(Nval, Cval, line, label, igrid, Nfiles, GRD)
            CASE ('ININAME')
              label='INI - nonlinear model initial conditions'
              Npts=load_s1d(Nval, Cval, line, label, igrid, Nfiles, INI)
            CASE ('IRPNAME')
              label='IRP - representer model initial conditions'
              Npts=load_s1d(Nval, Cval, line, label, igrid, Nfiles, IRP)
            CASE ('ITLNAME')
              label='ITL - tangent linear model initial conditions'
              Npts=load_s1d(Nval, Cval, line, label, igrid, Nfiles, ITL)
            CASE ('IADNAME')
              label='IAD - adjoint model initial conditions'
              Npts=load_s1d(Nval, Cval, line, label, igrid, Nfiles, IAD)
            CASE ('NFFILES')
              Npts=load_i(Nval, Rval, Ngrids, nFfiles)
              DO ng=1,Ngrids
                IF (nFfiles(ng).le.0) THEN
                  IF (Master) WRITE (out,260) 'NFFILES', nFfiles(ng),   &
     &              'Must be equal or greater than one.'
                  exit_flag=4
                  RETURN
                END IF
              END DO
              mFfiles=MAXVAL(nFfiles)
              allocate ( FRC(mFfiles,Ngrids) )
              allocate ( FRCids(mFfiles,Ngrids) )
              allocate ( Ncount(mFfiles,Ngrids) )
              FRCids(1:mFfiles,1:Ngrids)=-1
              Ncount(1:mFfiles,1:Ngrids)=0
            CASE ('FRCNAME')
              label='FRC - forcing fields'
              Npts=load_s2d(Nval, Cval, line, label, ifile, igrid,      &
     &                      Ncount, mFfiles, nFfiles, FRC)
            CASE ('NCLMFILES')
              Npts=load_i(Nval, Rval, Ngrids, nCLMfiles)
              DO ng=1,Ngrids
                IF (nCLMfiles(ng).le.0) THEN
                  IF (Master) WRITE (out,260) 'NCLMFILES',              &
     &               nCLMfiles(ng),                                     &
     &              'Must be equal or greater than one.'
                  exit_flag=4
                  RETURN
                END IF
              END DO
              mFfiles=MAXVAL(nCLMfiles)
              allocate ( CLM(mFfiles,Ngrids) )
              allocate ( CLMids(mFfiles,Ngrids) )
              allocate ( NCLMcount(mFfiles,Ngrids) )
              CLMids(1:mFfiles,1:Ngrids)=-1
              NCLMcount(1:mFfiles,1:Ngrids)=0
            CASE ('CLMNAME')
              label='CLM - climatology fields'
              Npts=load_s2d(Nval, Cval, line, label, iclmfile, igrid,   &
     &                      NCLMcount, mFfiles, nCLMfiles, CLM)
            CASE ('NBCFILES')
              Npts=load_i(Nval, Rval, Ngrids, nBCfiles)
              DO ng=1,Ngrids
                IF (nBCfiles(ng).le.0) THEN 
                  IF (Master) WRITE (out,260) 'NBCFILES', nBCfiles(ng), &
     &                            'Must be equal or greater than one.'
                  exit_flag=4
                  RETURN
                END IF
              END DO
              mFfiles=MAXVAL(nBCfiles)
              allocate ( BRY(mFfiles,Ngrids) )
              allocate ( BRYids(mFfiles,Ngrids) )
              allocate ( NBCcount(mFfiles,Ngrids) )
              BRYids(1:mFfiles,1:Ngrids)=-1
              NBCcount(1:mFfiles,1:Ngrids)=0
            CASE ('BRYNAME')
              label='BRY - lateral open boundary conditions'
              Npts=load_s2d(Nval, Cval, line, label, ibcfile, igrid,    &
     &                      NBCcount, mFfiles, nBCfiles, BRY)
            CASE ('FWDNAME')
              label='FWD - basic state forward trajectory'
              Npts=load_s1d(Nval, Cval, line, label, igrid, Nfiles, FWD)
            CASE ('ADSNAME')
              label='ADS - adjoint sensitivity functional'
              Npts=load_s1d(Nval, Cval, line, label, igrid, Nfiles, ADS)
            CASE ('APARNAM')
              DO i=1,LEN(aparnam)
                aparnam(i:i)=blank
              END DO
              aparnam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('SPOSNAM')
              DO i=1,LEN(sposnam)
                sposnam(i:i)=blank
              END DO
              sposnam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('FPOSNAM')
              DO i=1,LEN(fposnam)
                fposnam(i:i)=blank
              END DO
              fposnam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('BPARNAM')
              DO i=1,LEN(bparnam)
                bparnam(i:i)=blank
              END DO
              bparnam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('IPARNAM')
              DO i=1,LEN(iparnam)
                iparnam(i:i)=blank
              END DO
              iparnam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('SPARNAM')
              DO i=1,LEN(sparnam)
                sparnam(i:i)=blank
              END DO
              sparnam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('USRNAME')
              DO i=1,LEN(USRname)
                USRname(i:i)=blank
              END DO
              USRname=TRIM(ADJUSTL(Cval(Nval)))
          END SELECT
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
  10  IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
  20  CLOSE (inp)
!
!-----------------------------------------------------------------------
!  Process input parameters.
!-----------------------------------------------------------------------
!
!  Make sure that both component switches are activated when processing
!  (Eastward,Northward) momentum components at RHO-points.
!
      DO ng=1,Ngrids
        IF (.not.Hout(idu2dE,ng).and.Hout(idv2dN,ng)) THEN
          Hout(idu2dE,ng)=.TRUE.
        END IF
        IF (.not.Hout(idv2dN,ng).and.Hout(idu2dE,ng)) THEN
          Hout(idv2dN,ng)=.TRUE.
        END IF
        IF (.not.Hout(idu3dE,ng).and.Hout(idv3dN,ng)) THEN
          Hout(idu3dE,ng)=.TRUE.
        END IF
        IF (.not.Hout(idv3dN,ng).and.Hout(idu3dE,ng)) THEN
          Hout(idv3dN,ng)=.TRUE.
        END IF
        IF (.not.Aout(idu2dE,ng).and.Aout(idv2dN,ng)) THEN
          Aout(idu2dE,ng)=.TRUE.
        END IF
        IF (.not.Aout(idv2dN,ng).and.Aout(idu2dE,ng)) THEN
          Aout(idv2dN,ng)=.TRUE.
        END IF
        IF (.not.Aout(idu3dE,ng).and.Aout(idv3dN,ng)) THEN
          Aout(idu3dE,ng)=.TRUE.
        END IF
        IF (.not.Aout(idv3dN,ng).and.Aout(idu3dE,ng)) THEN
          Aout(idv3dN,ng)=.TRUE.
        END IF
      END DO
!
!  Set switch to create NetCDF file.
!
      DO ng=1,Ngrids
        DO i=1,NV
          IF (Hout(i,ng)) LdefHIS(ng)=.TRUE.
        END DO
        IF (((nrrec(ng).eq.0).and.(nAVG(ng).gt.ntimes(ng))).or.         &
     &      (nAVG(ng).eq.0)) THEN
          LdefAVG(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nAVG2(ng).gt.ntimes(ng))).or.        &
     &      (nAVG2(ng).eq.0)) THEN
          LdefAVG2(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nDIA(ng).gt.ntimes(ng))).or.         &
     &      (nDIA(ng).eq.0)) THEN
          LdefDIA(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nFLT(ng).gt.ntimes(ng))).or.         &
     &      (nFLT(ng).eq.0)) THEN
          LdefFLT(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nHIS(ng).gt.ntimes(ng))).or.         &
     &      (nHIS(ng).eq.0)) THEN
          LdefHIS(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nRST(ng).gt.ntimes(ng))).or.         &
     &      (nRST(ng).eq.0)) THEN
          LdefRST(ng)=.FALSE.
        END  IF
        IF (((nrrec(ng).eq.0).and.(nSTA(ng).gt.ntimes(ng))).or.         &
     &      (nSTA(ng).eq.0)) THEN
          LdefSTA(ng)=.FALSE.
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        WRITE (out,60) TRIM(title), TRIM(my_os), TRIM(my_cpu),          &
     &                 TRIM(my_fort), TRIM(my_fc), TRIM(my_fflags),     &
     &                 TRIM(Iname), TRIM(svn_url), TRIM(svn_rev),       &
     &                 TRIM(Rdir), TRIM(Hdir), TRIM(Hfile), TRIM(Adir)
        DO ng=1,Ngrids
!
!  Report grid size and domain decomposition.  Check for correct tile
!  decomposition.
!
          WRITE (out,70) ng, Lm(ng), Mm(ng), N(ng), numthreads,         &
     &                   NtileI(ng), NtileJ(ng)
          IF ((NtileI(ng)*NtileJ(ng)).ne.numthreads) THEN
            WRITE (out,80) ng
            exit_flag=6
            RETURN
          END IF
!
!  Report physical parameters.
!
          WRITE (out,110) ng
          WRITE (out,120) ntimes(ng), 'ntimes',                         &
     &          'Number of timesteps for 3-D equations.'
          WRITE (out,140) dt(ng), 'dt',                                 &
     &          'Timestep size (s) for 3-D equations.'
          WRITE (out,130) ndtfast(ng), 'ndtfast',                       &
     &          'Number of timesteps for 2-D equations between',        &
     &          'each 3D timestep.'
          WRITE (out,120) ERstr, 'ERstr',                               &
     &          'Starting ensemble/perturbation run number.'
          WRITE (out,120) ERend, 'ERend',                               &
     &          'Ending ensemble/perturbation run number.'
          WRITE (out,120) nrrec(ng), 'nrrec',                           &
     &          'Number of restart records to read from disk.'
          WRITE (out,170) LcycleRST(ng), 'LcycleRST',                   &
     &          'Switch to recycle time-records in restart file.'
          WRITE (out,130) nRST(ng), 'nRST',                             &
     &          'Number of timesteps between the writing of data',      &
     &          'into restart fields.'
          WRITE (out,130) ninfo(ng), 'ninfo',                           &
     &          'Number of timesteps between print of information',     &
     &          'to standard output.'
          WRITE (out,130) nSTA(ng), 'nSTA',                             &
     &          'Number of timesteps between the writing of data',      &
     &          'the stations file.'
          WRITE (out,170) ldefout(ng), 'ldefout',                       &
     &          'Switch to create a new output NetCDF file(s).'
          WRITE (out,130) nHIS(ng), 'nHIS',                             &
     &          'Number of timesteps between the writing fields',       &
     &          'into history file.'
          IF (ndefHIS(ng).gt.0) THEN
            WRITE (out,130) ndefHIS(ng), 'ndefHIS',                     &
     &            'Number of timesteps between creation of new',        &
     &            'history files.'
          END IF
          WRITE (out,130) ntsAVG(ng), 'ntsAVG',                         &
     &          'Starting timestep for the accumulation of output',     &
     &          'time-averaged data.'
          WRITE (out,130) nAVG(ng), 'nAVG',                             &
     &          'Number of timesteps between the writing of',           &
     &          'time-averaged data into averages file.'
          IF (ndefAVG(ng).gt.0) THEN
            WRITE (out,130) ndefAVG(ng), 'ndefAVG',                     &
     &            'Number of timesteps between creation of new',        &
     &            'time-averaged file.'
          END IF
          WRITE (out,130) ntsDIA(ng), 'ntsDIA',                         &
     &          'Starting timestep for the accumulation of output',     &
     &          'time-averaged diagnostics data.'
          WRITE (out,130) nDIA(ng), 'nDIA',                             &
     &          'Number of timesteps between the writing of',           &
     &          'time-averaged data into diagnostics file.'
          IF (ndefDIA(ng).gt.0) THEN
            WRITE (out,130) ndefDIA(ng), 'ndefDIA',                     &
     &            'Number of timesteps between creation of new',        &
     &            'diagnostic file.'
          END IF
          DO i=1,NAT+NPT
            itrc=i
            IF (i.gt.NAT) itrc=inert(i-NAT)
            WRITE (out,190) nl_tnu2(itrc,ng), 'nl_tnu2', itrc,          &
     &            'NLM Horizontal, harmonic mixing coefficient',        &
     &            '(m2/s) for tracer ', itrc,                           &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          WRITE (out,210) nl_visc2(ng), 'nl_visc2',                     &
     &          'NLM Horizontal, harmonic mixing coefficient',          &
     &          '(m2/s) for momentum.'
          DO i=1,NAT+NPT
            itrc=i
            IF (i.gt.NAT) itrc=inert(i-NAT)
            WRITE (out,190) Akt_bak(itrc,ng), 'Akt_bak', itrc,          &
     &            'Background vertical mixing coefficient (m2/s)',      &
     &            'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
          WRITE (out,210) Akv_bak(ng), 'Akv_bak',                       &
     &          'Background vertical mixing coefficient (m2/s)',        &
     &          'for momentum.'
          WRITE (out,200) rdrg(ng), 'rdrg',                             &
     &          'Linear bottom drag coefficient (m/s).'
          WRITE (out,200) rdrg2(ng), 'rdrg2',                           &
     &          'Quadratic bottom drag coefficient.'
          WRITE (out,200) Zob(ng), 'Zob',                               &
     &          'Bottom roughness (m).'
          WRITE (out,200) blk_ZQ(ng), 'blk_ZQ',                         &
     &          'Height (m) of surface air humidity measurement.'
          IF (blk_ZQ(ng).le.0.0_r8) THEN
            WRITE (out,265) 'blk_ZQ = ', blk_ZQ(ng),                    &
     &            'It must be greater than zero.'
            exit_flag=5
            RETURN
          END IF
          WRITE (out,200) blk_ZT(ng), 'blk_ZT',                         &
     &          'Height (m) of surface air temperature measurement.'
          IF (blk_ZT(ng).le.0.0_r8) THEN
            WRITE (out,265) 'blk_ZT = ', blk_ZT(ng),                    &
     &            'It must be greater than zero.'
            exit_flag=5
            RETURN
          END IF
          WRITE (out,200) blk_ZW(ng), 'blk_ZW',                         &
     &          'Height (m) of surface winds measurement.'
          IF (blk_ZW(ng).le.0.0_r8) THEN
            WRITE (out,265) 'blk_ZW = ', blk_ZW(ng),                    &
     &            'It must be greater than zero.'
            exit_flag=5
            RETURN
          END IF
          WRITE (out,120) lmd_Jwt(ng), 'lmd_Jwt',                       &
     &          'Jerlov water type.'
          IF ((lmd_Jwt(ng).lt.1).or.(lmd_Jwt(ng).gt.5)) THEN
            WRITE (out,260) 'lmd_Jwt = ', lmd_Jwt(ng),                  &
     &            'It must between one and five.'
            exit_flag=5
            RETURN
          END IF
          WRITE (out,120) Vtransform(ng), 'Vtransform',                 &
     &          'S-coordinate transformation equation.'
          WRITE (out,120) Vstretching(ng), 'Vstretching',               &
     &          'S-coordinate stretching function.'
          WRITE (out,200) theta_s(ng), 'theta_s',                       &
     &          'S-coordinate surface control parameter.'
          WRITE (out,200) theta_b(ng), 'theta_b',                       &
     &          'S-coordinate bottom  control parameter.'
          IF (Tcline(ng).gt.1.0E+5_r8) THEN
            WRITE (out,210) Tcline(ng), 'Tcline',                       &
     &            'S-coordinate surface/bottom layer width (m) used',   &
     &            'in vertical coordinate stretching.'
          ELSE
            WRITE (out,160) Tcline(ng), 'Tcline',                       &
     &            'S-coordinate surface/bottom layer width (m) used',   &
     &            'in vertical coordinate stretching.'
          END IF
          WRITE (out,140) rho0, 'rho0',                                 &
     &          'Mean density (kg/m3) for Boussinesq approximation.'
          WRITE (out,140) dstart, 'dstart',                             &
     &          'Time-stamp assigned to model initialization (days).'
          WRITE (out,150) time_ref, 'time_ref',                         &
     &          'Reference time for units attribute (yyyymmdd.dd)'
          DO i=1,NAT+NPT
            itrc=i
            IF (i.gt.NAT) itrc=inert(i-NAT)
            WRITE (out,190) Tnudg(itrc,ng), 'Tnudg', itrc,              &
     &            'Nudging/relaxation time scale (days)',               &
     &            'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
          WRITE (out,210) Znudg(ng), 'Znudg',                           &
     &          'Nudging/relaxation time scale (days)',                 &
     &          'for free-surface.'
          WRITE (out,210) M2nudg(ng), 'M2nudg',                         &
     &          'Nudging/relaxation time scale (days)',                 &
     &          'for 2D momentum.'
          WRITE (out,210) M3nudg(ng), 'M3nudg',                         &
     &          'Nudging/relaxation time scale (days)',                 &
     &          'for 3D momentum.'
          WRITE (out,210) obcfac(ng), 'obcfac',                         &
     &          'Factor between passive and active',                    &
     &          'open boundary conditions.'
          WRITE (out,170) VolCons(1,ng), 'VolCons(1)',                  &
     &          'NLM western  edge boundary volume conservation.'
          WRITE (out,170) VolCons(2,ng), 'VolCons(2)',                  &
     &          'NLM southern edge boundary volume conservation.'
          WRITE (out,170) VolCons(3,ng), 'VolCons(3)',                  &
     &          'NLM eastern  edge boundary volume conservation.'
          WRITE (out,170) VolCons(4,ng), 'VolCons(4)',                  &
     &          'NLM northern edge boundary volume conservation.'
          WRITE (out,140) T0(ng), 'T0',                                 &
     &          'Background potential temperature (C) constant.'
          WRITE (out,140) S0(ng), 'S0',                                 &
     &          'Background salinity (PSU) constant.'
          WRITE (out,160) gamma2(ng), 'gamma2',                         &
     &          'Slipperiness variable: free-slip (1.0) or ',           &
     &          '                     no-slip (-1.0).'
          IF (Hout(idFsur,ng)) WRITE (out,170) Hout(idFsur,ng),         &
     &       'Hout(idFsur)',                                            &
     &       'Write out free-surface.'
          IF (Hout(idUbar,ng)) WRITE (out,170) Hout(idUbar,ng),         &
     &       'Hout(idUbar)',                                            &
     &       'Write out 2D U-momentum component.'
          IF (Hout(idVbar,ng)) WRITE (out,170) Hout(idVbar,ng),         &
     &       'Hout(idVbar)',                                            &
     &       'Write out 2D V-momentum component.'
          IF (Hout(idu2dE,ng)) WRITE (out,170) Hout(idu2dE,ng),         &
     &       'Hout(idu2dE)',                                            &
     &       'Write out 2D U-eastward  at RHO-points.'
          IF (Hout(idv2dN,ng)) WRITE (out,170) Hout(idv2dN,ng),         &
     &       'Hout(idv2dN)',                                            &
     &       'Write out 2D V-northward at RHO-points.'
          IF (Hout(idUvel,ng)) WRITE (out,170) Hout(idUvel,ng),         &
     &       'Hout(idUvel)',                                            &
     &       'Write out 3D U-momentum component.'
          IF (Hout(idVvel,ng)) WRITE (out,170) Hout(idVvel,ng),         &
     &       'Hout(idVvel)',                                            &
     &       'Write out 3D V-momentum component.'
          IF (Hout(idu3dE,ng)) WRITE (out,170) Hout(idu3dE,ng),         &
     &       'Hout(idu3dE)',                                            &
     &       'Write out 3D U-wastward  component at RHO-points.'
          IF (Hout(idv3dN,ng)) WRITE (out,170) Hout(idv3dN,ng),         &
     &       'Hout(idv3dN)',                                            &
     &       'Write out 3D V-northward component at RHO-points.'
          IF (Hout(idWvel,ng)) WRITE (out,170) Hout(idWvel,ng),         &
     &       'Hout(idWvel)',                                            &
     &       'Write out W-momentum component.'
          IF (Hout(idOvel,ng)) WRITE (out,170) Hout(idOvel,ng),         &
     &       'Hout(idOvel)',                                            &
     &       'Write out omega vertical velocity.'
          DO itrc=1,NAT
            IF (Hout(idTvar(itrc),ng)) WRITE (out,180)                  &
     &          Hout(idTvar(itrc),ng), 'Hout(idTvar)',                  &
     &          'Write out tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
          IF (Hout(idUsms,ng)) WRITE (out,170) Hout(idUsms,ng),         &
     &       'Hout(idUsms)',                                            &
     &       'Write out surface U-momentum stress.'
          IF (Hout(idVsms,ng)) WRITE (out,170) Hout(idVsms,ng),         &
     &       'Hout(idVsms)',                                            &
     &       'Write out surface V-momentum stress.'
          IF (Hout(idUbms,ng)) WRITE (out,170) Hout(idUbms,ng),         &
     &       'Hout(idUbms)',                                            &
     &       'Write out bottom U-momentum stress.'
          IF (Hout(idVbms,ng)) WRITE (out,170) Hout(idVbms,ng),         &
     &       'Hout(idVbms)',                                            &
     &       'Write out bottom V-momentum stress.'
          DO itrc=1,NPT
            IF (Hout(idTvar(inert(itrc)),ng)) WRITE (out,180)           &
     &          Hout(idTvar(inert(itrc)),ng), 'Hout(inert)',            &
     &          'Write out inert passive tracer ', itrc,                &
     &          TRIM(Vname(1,idTvar(inert(itrc))))
          END DO
          IF (Hout(idPair,ng)) WRITE (out,170) Hout(idPair,ng),         &
     &       'Hout(idPair)',                                            &
     &       'Write out surface air pressure.'
          IF (Hout(idUair,ng)) WRITE (out,170) Hout(idUair,ng),         &
     &       'Hout(idUair)',                                            &
     &       'Write out surface U-wind component.'
          IF (Hout(idVair,ng)) WRITE (out,170) Hout(idVair,ng),         &
     &       'Hout(idVair)',                                            &
     &       'Write out surface V-wind component.'
          IF (Hout(idUairE,ng)) WRITE (out,170) Hout(idUairE,ng),       &
     &       'Hout(idUairE)',                                           &
     &       'Write out surface Eastward U-wind component.'
          IF (Hout(idVairN,ng)) WRITE (out,170) Hout(idVairN,ng),       &
     &       'Hout(idVairN)',                                           &
     &       'Write out surface Northward V-wind component.'
          IF (Hout(idTsur(itemp),ng)) WRITE (out,170)                   &
     &        Hout(idTsur(itemp),ng), 'Hout(idTsur)',                   &
     &       'Write out surface net heat flux.'
          IF (Hout(idTsur(isalt),ng)) WRITE (out,170)                   &
     &        Hout(idTsur(isalt),ng), 'Hout(idTsur)',                   &
     &       'Write out surface net salt flux.'
          IF (Hout(idSrad,ng)) WRITE (out,170) Hout(idSrad,ng),         &
     &       'Hout(idSrad)',                                            &
     &       'Write out shortwave radiation flux.'
          IF (Hout(idLrad,ng)) WRITE (out,170) Hout(idLrad,ng),         &
     &       'Hout(idLrad)',                                            &
     &       'Write out longwave radiation flux.'
          IF (Hout(idLhea,ng)) WRITE (out,170) Hout(idLhea,ng),         &
     &       'Hout(idLhea)',                                            &
     &       'Write out latent heat flux.'
          IF (Hout(idShea,ng)) WRITE (out,170) Hout(idShea,ng),         &
     &       'Hout(idShea)',                                            &
     &       'Write out sensible heat flux.'
          IF (Hout(idEmPf,ng)) WRITE (out,170) Hout(idEmPf,ng),         &
     &       'Hout(idEmPf)',                                            &
     &       'Write out E-P flux.'
          IF (Hout(idevap,ng)) WRITE (out,170) Hout(idevap,ng),         &
     &       'Hout(idevap)',                                            &
     &       'Write out evaporation rate.'
          IF (Hout(idrain,ng)) WRITE (out,170) Hout(idrain,ng),         &
     &       'Hout(idrain)',                                            &
     &       'Write out rain rate.'
          IF (Hout(idDano,ng)) WRITE (out,170) Hout(idDano,ng),         &
     &       'Hout(idDano)',                                            &
     &       'Write out density anomaly.'
          IF (Hout(idVvis,ng)) WRITE (out,170) Hout(idVvis,ng),         &
     &       'Hout(idVvis)',                                            &
     &       'Write out vertical viscosity: AKv.'
          IF (Hout(idTdif,ng)) WRITE (out,170) Hout(idTdif,ng),         &
     &       'Hout(idTdif)',                                            &
     &       'Write out vertical diffusion: AKt(itemp).'
          IF (Hout(idSdif,ng)) WRITE (out,170) Hout(idSdif,ng),         &
     &       'Hout(idSdif)',                                            &
     &       'Write out vertical diffusion: AKt(isalt).'
          IF (Hout(idHsbl,ng)) WRITE (out,170) Hout(idHsbl,ng),         &
     &       'Hout(idHsbl)',                                            &
     &       'Write out depth of surface boundary layer.'
          WRITE (out,'(1x)')
          IF (Aout(idFsur,ng)) WRITE (out,170) Aout(idFsur,ng),         &
     &       'Aout(idFsur)',                                            &
     &       'Write out averaged free-surface.'
          IF (Aout(idUbar,ng)) WRITE (out,170) Aout(idUbar,ng),         &
     &       'Aout(idUbar)',                                            &
     &       'Write out averaged 2D U-momentum component.'
          IF (Aout(idVbar,ng)) WRITE (out,170) Aout(idVbar,ng),         &
     &       'Aout(idVbar)',                                            &
     &       'Write out averaged 2D V-momentum component.'
          IF (Aout(idu2dE,ng)) WRITE (out,170) Aout(idu2dE,ng),         &
     &       'Aout(idu2dE)',                                            &
     &       'Write out averaged 2D U-eastward  at RHO-points.'
          IF (Aout(idv2dN,ng)) WRITE (out,170) Aout(idv2dN,ng),         &
     &       'Aout(idv2dN)',                                            &
     &       'Write out averaged 2D V-northward at RHO-points.'
          IF (Aout(idUvel,ng)) WRITE (out,170) Aout(idUvel,ng),         &
     &       'Aout(idUvel)',                                            &
     &       'Write out averaged 3D U-momentum component.'
          IF (Aout(idVvel,ng)) WRITE (out,170) Aout(idVvel,ng),         &
     &       'Aout(idVvel)',                                            &
     &       'Write out averaged 3D V-momentum component.'
          IF (Aout(idu3dE,ng)) WRITE (out,170) Aout(idu3dE,ng),         &
     &       'Aout(idu3dE)',                                            &
     &       'Write out averaged 3D U-eastward  at RHO-points.'
          IF (Aout(idv3dN,ng)) WRITE (out,170) Aout(idv3dN,ng),         &
     &       'Aout(idv3dN)',                                            &
     &       'Write out averaged 3D V-northward at RHO-points.'
          IF (Aout(idWvel,ng)) WRITE (out,170) Aout(idWvel,ng),         &
     &       'Aout(idWvel)',                                            &
     &       'Write out averaged W-momentum component.'
          IF (Aout(idOvel,ng)) WRITE (out,170) Aout(idOvel,ng),         &
     &       'Aout(idOvel)',                                            &
     &       'Write out averaged omega vertical velocity.'
          DO itrc=1,NAT
            IF (Aout(idTvar(itrc),ng)) WRITE (out,180)                  &
     &          Aout(idTvar(itrc),ng), 'Aout(idTvar)',                  &
     &          'Write out averaged tracer ', itrc,                     &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
          IF (Aout(idUsms,ng)) WRITE (out,170) Aout(idUsms,ng),         &
     &       'Aout(idUsms)',                                            &
     &       'Write out averaged surface U-momentum stress.'
          IF (Aout(idVsms,ng)) WRITE (out,170) Aout(idVsms,ng),         &
     &       'Aout(idVsms)',                                            &
     &       'Write out averaged surface V-momentum stress.'
          IF (Aout(idUbms,ng)) WRITE (out,170) Aout(idUbms,ng),         &
     &       'Aout(idUbms)',                                            &
     &       'Write out averaged bottom U-momentum stress.'
          IF (Aout(idVbms,ng)) WRITE (out,170) Aout(idVbms,ng),         &
     &       'Aout(idVbms)',                                            &
     &       'Write out averaged bottom V-momentum stress.'
          DO itrc=1,NPT
            IF (Aout(idTvar(inert(itrc)),ng)) WRITE (out,180)           &
     &          Aout(idTvar(inert(itrc)),ng), 'Aout(inert)',            &
     &          'Write out  averaged inert passive tracer ', itrc,      &
     &          TRIM(Vname(1,idTvar(inert(itrc))))
          END DO
          IF (Aout(idPair,ng)) WRITE (out,170) Aout(idPair,ng),         &
     &       'Aout(idPair)',                                            &
     &       'Write out averaged surface air pressure.'
          IF (Aout(idUair,ng)) WRITE (out,170) Aout(idUair,ng),         &
     &       'Aout(idUair)',                                            &
     &       'Write out averaged surface U-wind component.'
          IF (Aout(idVair,ng)) WRITE (out,170) Aout(idVair,ng),         &
     &       'Aout(idVair)',                                            &
     &       'Write out averaged surface V-wind component.'
          IF (Aout(idUairE,ng)) WRITE (out,170) Aout(idUairE,ng),       &
     &       'Aout(idUairE)',                                           &
     &       'Write out averaged Eastward surface U-wind component.'
          IF (Aout(idVairN,ng)) WRITE (out,170) Aout(idVairN,ng),       &
     &       'Aout(idVairN)',                                           &
     &       'Write out averaged Northward surface V-wind component.'
          IF (Aout(idTsur(itemp),ng)) WRITE (out,170)                   &
     &        Aout(idTsur(itemp),ng), 'Aout(idTsur)',                   &
     &       'Write out averaged surface net heat flux.'
          IF (Aout(idTsur(isalt),ng)) WRITE (out,170)                   &
     &        Aout(idTsur(isalt),ng), 'Aout(idTsur)',                   &
     &       'Write out averaged surface net salt flux.'
          IF (Aout(idSrad,ng)) WRITE (out,170) Aout(idSrad,ng),         &
     &       'Aout(idSrad)',                                            &
     &       'Write out averaged shortwave radiation flux.'
          IF (Aout(idLrad,ng)) WRITE (out,170) Aout(idLrad,ng),         &
     &       'Aout(idLrad)',                                            &
     &       'Write out averaged longwave radiation flux.'
          IF (Aout(idLhea,ng)) WRITE (out,170) Aout(idLhea,ng),         &
     &       'Aout(idLhea)',                                            &
     &       'Write out averaged latent heat flux.'
          IF (Aout(idShea,ng)) WRITE (out,170) Aout(idShea,ng),         &
     &       'Aout(idShea)',                                            &
     &       'Write out averaged sensible heat flux.'
          IF (Aout(idevap,ng)) WRITE (out,170) Aout(idevap,ng),         &
     &       'Aout(idevap)',                                            &
     &       'Write out averaged evaporation rate.'
          IF (Aout(idrain,ng)) WRITE (out,170) Aout(idrain,ng),         &
     &       'Aout(idrain)',                                            &
     &       'Write out averaged rain rate.'
          IF (Aout(idDano,ng)) WRITE (out,170) Aout(idDano,ng),         &
     &       'Aout(idDano)',                                            &
     &       'Write out averaged density anomaly.'
          IF (Aout(idVvis,ng)) WRITE (out,170) Aout(idVvis,ng),         &
     &       'Aout(idVvis)',                                            &
     &       'Write out averaged vertical viscosity: AKv.'
          IF (Aout(idTdif,ng)) WRITE (out,170) Aout(idTdif,ng),         &
     &       'Aout(idTdif)',                                            &
     &       'Write out averaged vertical diffusion: AKt(itemp).'
          IF (Aout(idSdif,ng)) WRITE (out,170) Aout(idSdif,ng),         &
     &       'Aout(idSdif)',                                            &
     &       'Write out averaged vertical diffusion: AKt(isalt).'
          IF (Aout(idHsbl,ng)) WRITE (out,170) Aout(idHsbl,ng),         &
     &       'Aout(idHsbl)',                                            &
     &       'Write out averaged depth of surface boundary layer.'
          IF (Aout(id2dRV,ng)) WRITE (out,170) Aout(id2dRV,ng),         &
     &       'Aout(id2dRV)',                                            &
     &       'Write out averaged 2D relative vorticity.'
          IF (Aout(id2dPV,ng)) WRITE (out,170) Aout(id2dPV,ng),         &
     &       'Aout(id2dPV)',                                            &
     &       'Write out averaged 2D potential vorticity.'
          IF (Aout(id3dRV,ng)) WRITE (out,170) Aout(id3dRV,ng),         &
     &       'Aout(id3dRV)',                                            &
     &       'Write out averaged 3D relative vorticity.'
          IF (Aout(id3dPV,ng)) WRITE (out,170) Aout(id3dPV,ng),         &
     &       'Aout(id3dPV)',                                            &
     &       'Write out averaged 3D potential vorticity.'
          IF (Aout(idZZav,ng)) WRITE (out,170) Aout(idZZav,ng),         &
     &       'Aout(idZZav)',                                            &
     &       'Write out averaged quadratic <zeta*zeta> term.'
          IF (Aout(idU2av,ng)) WRITE (out,170) Aout(idU2av,ng),         &
     &       'Aout(idU2av)',                                            &
     &       'Write out averaged quadratic <ubar*ubar> term.'
          IF (Aout(idV2av,ng)) WRITE (out,170) Aout(idV2av,ng),         &
     &       'Aout(idV2av)',                                            &
     &       'Write out averaged quadratic <vbar*vbar> term.'
          IF (Aout(idHUav,ng)) WRITE (out,170) Aout(idHUav,ng),         &
     &       'Aout(idHUav)',                                            &
     &       'Write out averaged u-volume flux, Huon.'
          IF (Aout(idHVav,ng)) WRITE (out,170) Aout(idHVav,ng),         &
     &       'Aout(idHVav)',                                            &
     &       'Write out averaged v-volume flux, Hvom.'
          IF (Aout(idUUav,ng)) WRITE (out,170) Aout(idUUav,ng),         &
     &       'Aout(idUUav)',                                            &
     &       'Write out averaged quadratic <u*u> term.'
          IF (Aout(idUVav,ng)) WRITE (out,170) Aout(idUVav,ng),         &
     &       'Aout(idUVav)',                                            &
     &       'Write out averaged quadratic <u*v> term.'
          IF (Aout(idVVav,ng)) WRITE (out,170) Aout(idVVav,ng),         &
     &       'Aout(idVVav)',                                            &
     &       'Write out averaged quadratic <v*v> term.'
          DO itrc=1,NAT
            IF (Aout(idTTav(itrc),ng)) WRITE (out,180)                  &
     &          Aout(idTTav(itrc),ng), 'Aout(idTTav)',                  &
     &          'Write out averaged <t*t> for tracer ', itrc,           &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Aout(idUTav(itrc),ng)) WRITE (out,180)                  &
     &          Aout(idUTav(itrc),ng), 'Aout(idUTav)',                  &
     &          'Write out averaged <u*t> for tracer ', itrc,           &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Aout(idVTav(itrc),ng)) WRITE (out,180)                  &
     &          Aout(idVTav(itrc),ng), 'Aout(idVTav)',                  &
     &          'Write out averaged <v*t> for tracer ', itrc,           &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Aout(iHUTav(itrc),ng)) WRITE (out,180)                  &
     &          Aout(iHUTav(itrc),ng), 'Aout(iHUTav)',                  &
     &          'Write out averaged <Huon*t> for tracer ', itrc,        &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Aout(iHVTav(itrc),ng)) WRITE (out,180)                  &
     &          Aout(iHVTav(itrc),ng), 'Aout(iHVTav)',                  &
     &          'Write out averaged <Hvom*t> for tracer ', itrc,        &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
          WRITE (out,'(1x)')
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iTrate),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTrate)',                   &
     &            'Write out rate of change of tracer ', itrc,          &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(itrc,iTrate),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTrate)',                   &
     &            'Write out rate of change of tracer ', itrc,          &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iThadv),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iThadv)',                   &
     &            'Write out horizontal advection, tracer ', itrc,      &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(itrc,iThadv),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iThadv)',                   &
     &            'Write out horizontal advection, tracer ', itrc,      &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iTxadv),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTxadv)',                   &
     &            'Write out horizontal X-advection, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(itrc,iTxadv),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTxadv)',                   &
     &            'Write out horizontal X-advection, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iTyadv),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTyadv)',                   &
     &            'Write out horizontal Y-advection, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(itrc,iTyadv),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTyadv)',                   &
     &            'Write out horizontal Y-advection, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iTvadv),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTvadv)',                   &
     &            'Write out vertical advection, tracer ', itrc,        &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(itrc,iTvadv),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTvadv)',                   &
     &            'Write out vertical advection, tracer ', itrc,        &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iThdif),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iThdif)',                   &
     &            'Write out horizontal diffusion, tracer ', itrc,      &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(itrc,iThdif),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iThdif)',                   &
     &            'Write out horizontal diffusion, tracer ', itrc,      &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iTxdif),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTxdif)',                   &
     &            'Write out horizontal X-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(i,iTxdif),ng))                              &
     &        WRITE (out,180) .TRUE., 'Dout(iTxdif)',                   &
     &            'Write out horizontal X-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iTydif),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTydif)',                   &
     &            'Write out horizontal Y-diffusion , tracer ', itrc,   &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(itrc,iTydif),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTydif)',                   &
     &            'Write out horizontal Y-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iTsdif),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTsdif)',                   &
     &            'Write out horizontal S-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(itrc,iTsdif),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTsdif)',                   &
     &            'Write out horizontal S-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iTvdif),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTvdif)',                   &
     &            'Write out vertical diffusion, tracer ', itrc,        &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(itrc,iTvdif),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTvdif)',                   &
     &            'Write out vertical diffusion, tracer ', itrc,        &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
!
!-----------------------------------------------------------------------
!  Report output/input files and check availability of input files.
!-----------------------------------------------------------------------
!
          WRITE (out,220)
          WRITE (out,230) '           Output Restart File:  ',          &
     &                    TRIM(RST(ng)%name)
          IF (LdefHIS(ng)) THEN
            IF (ndefHIS(ng).eq.0) THEN
              WRITE (out,230) '           Output History File:  ',      &
     &                        TRIM(HIS(ng)%name)
            ELSE
              WRITE (out,230) '      Prefix for History Files:  ',      &
     &                        TRIM(HIS(ng)%base)
            END IF
          END IF
          WRITE (out,230) '     Output Forward State File:  ',          &
     &                    TRIM(FWD(ng)%name)
          IF (ndefAVG(ng).eq.0) THEN
            WRITE (out,230) '          Output Averages File:  ',        &
     &                      TRIM(AVG(ng)%name)
          ELSE
            WRITE (out,230) '     Prefix for Averages Files:  ',        &
     &                      TRIM(AVG(ng)%base)
          END IF
          IF (ndefDIA(ng).eq.0) THEN
            WRITE (out,230) '       Output Diagnostics File:  ',        &
     &                      TRIM(DIA(ng)%name)
          ELSE
            WRITE (out,230) '  Prefix for Diagnostics Files:  ',        &
     &                      TRIM(DIA(ng)%base)
          END IF
          WRITE (out,230) '          Output Stations File:  ',          &
     &                    TRIM(STA(ng)%name)
          fname=GRD(ng)%name
          IF (.not.find_file(ng, fname)) GO TO 30
          WRITE (out,230) '               Input Grid File:  ',          &
     &                    TRIM(fname)
          fname=INI(ng)%name
          IF (.not.find_file(ng, fname)) GO TO 30
          WRITE (out,230) '  Input Nonlinear Initial File:  ',          &
     &                    TRIM(fname)
          DO i=1,nFfiles(ng)
            DO ifile=1,FRC(i,ng)%Nfiles
              fname=FRC(i,ng)%files(ifile)
              IF (.not.find_file(ng, fname)) GO TO 30
              IF (ifile.eq.1) THEN
                WRITE (out,310) '         Input Forcing File ', i,      &
     &                          ':  ', TRIM(fname)
              ELSE
                WRITE (out,'(35x,a)') TRIM(fname)
              END IF
            END DO
          END DO
          DO i=1,nCLMfiles(ng)
            DO ifile=1,CLM(i,ng)%Nfiles
              fname=CLM(i,ng)%files(ifile)
              IF (.not.find_file(ng, fname)) GO TO 30
              IF (ifile.eq.1) THEN
                WRITE (out,230) '        Input Climatology File:  ',    &
     &                        TRIM(fname)
              ELSE
                WRITE (out,'(35x,a)') TRIM(fname)
              END IF
            END DO
          END DO
          ObcData(ng)=ObcData(ng).or.ANY(LBC(:,isFsur,ng)%acquire)
          ObcData(ng)=ObcData(ng).or.ANY(LBC(:,isUbar,ng)%acquire)      &
     &                           .or.ANY(LBC(:,isVbar,ng)%acquire)
          ObcData(ng)=ObcData(ng).or.ANY(LBC(:,isUvel,ng)%acquire)      &
     &                           .or.ANY(LBC(:,isVvel,ng)%acquire)
          ObcData(ng)=ObcData(ng).or.ANY(LBC(:,isTvar(:),ng)%acquire)
!
          IF (ObcData(ng)) THEN
            DO i=1,nBCfiles(ng)
              DO ifile=1,BRY(i,ng)%Nfiles
                fname=BRY(i,ng)%files(ifile)
                IF (.not.find_file(ng, fname)) GO TO 30
                IF (ifile.eq.1) THEN
                  WRITE (out,310) '       Input Boundary File ', i,     &
     &                            ':  ', TRIM(fname)
                ELSE
                  WRITE (out,'(35x,a)') TRIM(fname)
                END IF
              END DO
            END DO
          END IF
          fname=sposnam
          IF (.not.find_file(ng, fname)) GO TO 30
          WRITE (out,230) '        Station positions File:  ',          &
     &                    TRIM(fname)
          fname=iparnam
          IF (.not.find_file(ng, fname)) GO TO 30
          WRITE (out,230) '     Ice Model Parameters File:  ',          &
     &                    TRIM(fname)
          fname=varname
          IF (.not.find_file(ng, fname)) GO TO 30
          GO TO 40
  30      IF (Master) WRITE (out,270) TRIM(fname)
          exit_flag=4
          RETURN
  40      CONTINUE
        END DO
        IF (Nuser.gt.0) THEN
          WRITE (out,230) '        Input/Output USER File:  ',          &
     &                    TRIM(USRname)
        END IF
!
!-----------------------------------------------------------------------
!  Report generic USER parameters.
!-----------------------------------------------------------------------
!
        IF (Nuser.gt.0) THEN
          WRITE (out,240)
          DO i=1,Nuser
            WRITE (out,250) user(i), i, i
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Rescale active tracer parameters
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO i=1,NAT+NPT
          itrc=i
          IF (i.gt.NAT) itrc=inert(i-NAT)
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
          nl_tnu4(itrc,ng)=SQRT(ABS(nl_tnu4(itrc,ng)))
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
          IF (Tnudg(itrc,ng).gt.0.0_r8) THEN
            Tnudg(itrc,ng)=1.0_r8/(Tnudg(itrc,ng)*86400.0_r8)
          ELSE
            Tnudg(itrc,ng)=0.0_r8
          END IF
        END DO
      END DO
  50  FORMAT (/,' READ_PhyPar - Error while processing line: ',/,a)
  60  FORMAT (/,1x,a,/,                                                 &
     &        /,1x,'Operating system : ',a,                             &
     &        /,1x,'CPU/hardware     : ',a,                             &
     &        /,1x,'Compiler system  : ',a,                             &
     &        /,1x,'Compiler command : ',a,                             &
     &        /,1x,'Compiler flags   : ',a,/,                           &
     &        /,1x,'Input Script  : ',a,/,                              &
     &        /,1x,'SVN Root URL  : ',a,                                &
     &        /,1x,'SVN Revision  : ',a,/,                              &
     &        /,1x,'Local Root    : ',a,                                &
     &        /,1x,'Header Dir    : ',a,                                &
     &        /,1x,'Header file   : ',a,                                &
     &        /,1x,'Analytical Dir: ',a)
  70  FORMAT (/,' Resolution, Grid ',i2.2,': ',i4.4,'x',i4.4,'x',i3.3,  &
     &        ',',2x,'Parallel Nodes: ',i3,',',2x,'Tiling: ',i3.3,      &
     &        'x',i3.3)
  80  FORMAT (/,' ROMS/TOMS: Wrong choice of domain ',i2.2,1x,          &
     &        'partition or number of parallel threads.',               &
     &        /,12x,'NtileI * NtileJ  must be equal to the number of ', &
     &        'parallel nodes.',                                        &
     &        /,12x,'Change -np value to mpirun or',                    &
     &        /,12x,'change domain partition in input script.')
  90  FORMAT (/,' Resolution, Grid ',i2.2,': ',i4.4,'x',i4.4,'x',i3.3,  &
     &        ',',2x,'Parallel Threads: ',i2,',',2x,'Tiling: ',i3.3,    &
     &        'x',i3.3)
 100  FORMAT (/,' ROMS/TOMS: Wrong choice of domain ',i3.3,1x,          &
     &        'partition or number of parallel threads.',               &
     &        /,12x,'NtileI*NtileJ must be a positive multiple of the', &
     &        ' number of threads.',                                    &
     &        /,12x,'Change number of threads (environment variable) ', &
     &        'or',/,12x,'change domain partition in input script.')
 110  FORMAT (/,/,' Physical Parameters, Grid: ',i2.2,                  &
     &        /,  ' =============================',/)
 120  FORMAT (1x,i10,2x,a,t30,a)
 130  FORMAT (1x,i10,2x,a,t30,a,/,t32,a)
 140  FORMAT (f11.3,2x,a,t30,a)
 150  FORMAT (f11.2,2x,a,t30,a)
 160  FORMAT (f11.3,2x,a,t30,a,/,t32,a)
 170  FORMAT (10x,l1,2x,a,t30,a)
 180  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
 185  FORMAT (10x,l1,2x,a,'(',i2.2,')',t30,a,i2.2,':',1x,a)
 190  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,/,t32,a,i2.2,':',1x,a)
 195  FORMAT (1p,e11.4,2x,a,t30,a,i2.2,':',1x,a)
 200  FORMAT (1p,e11.4,2x,a,t30,a)
 210  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)
 220  FORMAT (/,' Output/Input Files:',/)
 230  FORMAT (2x,a,a)
 240  FORMAT (/,' Generic User Parameters:',/)
 250  FORMAT (1p,e11.4,2x,'user(',i2.2,')',t30,                         &
     &        'User parameter ',i2.2,'.')
 260  FORMAT (/,' READ_PHYPAR - Invalid input parameter, ',a,           &
     &        i4,/,15x,a)
 265  FORMAT (/,' READ_PHYPAR - Invalid input parameter, ',a,           &
     &        1p,e11.4,/,15x,a)
 270  FORMAT (/,' READ_PHYPAR - could not find input file:  ',a)
 280  FORMAT (/,' READ_PHYPAR - variable info not yet loaded, ', a)
 290  FORMAT (/,' READ_PHYPAR - Invalid dimension parameter, ',a,i4,    &
     &        /,15x,a)
 300  FORMAT (/,' READ_PHYPAR - Invalid dimension parameter, ',a,'(',   &
     &        i2.2,')',/,15x,a)
 310  FORMAT (2x,a,i2.2,a,a)
 320  FORMAT (/,' READ_PHYPAR - could not find input parameter: ', a,   &
     &        /,15x,'in ROMS standard input script.',/,15x,a)
 330  FORMAT (/,' READ_PHYPAR - Invalid input parameter, ',a,i4,/,15x,a)
      RETURN
      END SUBROUTINE read_PhyPar
