      MODULE mod_kinds
!
!svn $Id: mod_kinds.F 1451 2012-02-02 20:56:14Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!
        implicit none
!
        integer, parameter :: i1b= selected_int_kind(2)        !  1-byte
        integer, parameter :: i2b= selected_int_kind(4)        !  2-byte
        integer, parameter :: i4b= selected_int_kind(9)        !  4-byte
        integer, parameter :: c8 = selected_real_kind(6,30)    ! 32-bit
        integer, parameter :: r4 = selected_real_kind(6,30)    ! 32-bit
        integer, parameter :: r8 = selected_real_kind(12,300)  ! 64-bit
        integer, parameter :: r16 = selected_real_kind(15,300) !128-bit
        integer, parameter :: char_len = 256
      END MODULE mod_kinds
