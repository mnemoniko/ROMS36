/*
** svn $Id: flt_test.h 1451 2012-02-02 20:56:14Z kate $
*******************************************************************************
** Copyright (c) 2002-2012 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Floats Tracking Test.
**
** Application flag:   FLT_TEST
** Input script:       ocean_flt_test2d.in,  ocean_flt_test3d.in
**                     floats_flt_test2d.in, floats_flt_test3d.in
*/

#define UV_ADV
#define UV_QDRAG
#define UV_VIS2
#define MIX_S_UV
#define FLOATS
#define MASKING
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX

#ifdef SOLVE3D
# define DJ_GRADPS
# define TS_A4HADVECTION
# define TS_A4VADVECTION
# define BODYFORCE
# define SPLINES
# define ANA_BTFLUX
# define ANA_STFLUX
#endif

