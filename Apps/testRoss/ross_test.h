/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2009 The ROMS/TOMS Group
**
**   Licensed under a MIT/X style license
**
**   See License_ROMS.txt
**
*******************************************************************************
**
**  Options for Ross Sea Model
*/

#define UV_ADV
#define DJ_GRADPS
#define UV_QDRAG
#define UV_VIS2
#define UV_COR
#define TS_DIF2
#define TS_U3HADVECTION
#define TS_C4VADVECTION  /* Subject to change */

#define SALINITY
#define NONLIN_EOS
#define SOLVE3D
#define ANA_SRFLUX
#define ANA_BTFLUX
#define ANA_BSFLUX
#define MIX_S_UV         /* Note: This is different from old Ross */
#define MIX_GEO_TS
#define LMD_MIXING
#define LMD_RIMIX
#define LMD_CONVEC
#define LMD_SKPP
#define LMD_NONLOCAL

#define SPLINES
#define MASKING
#define ICESHELF

#undef BULK_FLUXES
#ifdef BULK_FLUXES
# ifdef ICE_MODEL
#  define ICE_BULK_FLUXES
# endif
# undef LONGWAVE
# undef EMINUSP
# undef SOLAR_SOURCE
#endif

#define AVERAGES
#undef  CURVGRID        /* not needed here...be careful with wind rotation */

#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SSFLUX
#define ANA_WINDS
#define EASTERN_WALL
#define WESTERN_WALL
#define SOUTHERN_WALL
#define NORTHERN_WALL

