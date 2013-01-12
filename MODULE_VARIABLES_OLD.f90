MODULE MODULE_VARIABLES_OLD

! ====================================================================
! SNOW.h
! ====================================================================
      integer SNOW_RUN
      parameter (SNOW_RUN=1)
! ====================================================================
! SNOW_RUN	0 if you want to shut the snow model off.
!		1 if you want the snow mode to run.
! ====================================================================

! ====================================================================
! wgtpar.h
! ====================================================================
! ====================================================================
! In editing this file : please make sure there are spaces between
! each word and character in the lines defining the variables.
! This is in order to allow the program calc_memory.! to run.
! Thus keep the format:
! parameter ( MAX_VAR = x )
! Also make sure parameter is spelled in lower case.
! Do not change the line with SNW_FLG.
! ====================================================================

      integer MAX_PIX,MAX_SPP,MAX_STA,MAX_SOI,MAX_VEG
      integer MAX_CAT,MAX_ROW,MAX_COL,MAX_FIL,MAX_SER
      integer MAX_ATB,MAX_TST,MAX_VST,MAX_TOP,MAX_LAN
      integer MAX_PP1,MAX_PP2,MOS_FLG,UST_FLG,SNW_FLG
      integer LAK_FLG

      parameter ( MAX_PIX = 3960)
      parameter ( MAX_PP1 = 3960)
      parameter ( MAX_PP2 = 1 )
      parameter ( MAX_SPP = 1 )
      parameter ( MAX_STA = 1 )
      parameter ( MAX_SOI = 3960)
      parameter ( MAX_VEG = 3960)
      parameter ( MAX_VST = 1 )
      parameter ( MAX_LAN = 1 )
      parameter ( MAX_CAT = 1 )
      parameter ( MAX_ROW = 66)
      parameter ( MAX_COL = 60)
      parameter ( MAX_FIL = 400 )
      parameter ( MAX_SER = 1 )
      parameter ( MAX_ATB = 30 )
      parameter ( MAX_TOP = 30 )
      parameter ( MAX_TST = 1 )
      parameter ( MOS_FLG = 0 )
      parameter ( UST_FLG = 0 )
      parameter ( LAK_FLG = 0 )

      parameter ( SNW_FLG = SNOW_RUN*(MOS_FLG+UST_FLG-MOS_FLG*UST_FLG) )

! ====================================================================
! Parameters for array dimensions for weighting and station files.
!
! MAX_PIX - Maximum number of pixels.
! MAX_PP1 - Maximum number of pixels in statistical mode.  If running
!	    in statistical mode this number should equal MAX_PIX.
!	    In distributed mode this number should equal 1.
! MAX_PP2 - Maximum number of pixels in statistical mode for fractional
!	    coverage.  If running in distributed mode or in statistical
!	    mode without representation of rainfall fractional coverage
!	    this number should equal 1.  If rainfall fractional coverage
!	    is represented this number should equal MAX_PIX.
! MAX_SPP - Maximum number of stations used per pixel.
! MAX_STA - Maximum number of stations for each forcing variable.
! MAX_SOI - Maximum number of soil classes.
! MAX_VEG - Maximum number of land cover classes.
! MAX_VST - Maximum number of vegetation classes per pixel in
!	    statistical mode.
!	    If running the model in statistical mode, this number
!	    should equal MAX_VEG or less.
!	    If running the model in distributed mode, this number
!	    should equal 1.
! MAX_LAN - Maximum number of vegetation classes per pixel.
!	    If running the model in distributed mode, this number
!	    should equal 1.
!	    If considering fractional coverage of rainfall than this
!	    number should equal MAX_VST.
!	    If assuming rainfall to be uniformly distributed over
!	    each grid than this paramters should equal 1.
! MAX_CAT - Maximum number of catchments.
! MAX_ROW - Maximum number of rows in the image.
! MAX_COL - Maximum number of columns in the image.
! MAX_FIL - Highest input or output file number.
! MAX_SER - Maximum number of series of printing the output
!           images (e.g. one serie is from time step 2 through 7, a
!           second serie from time step 25 through 33, ...).
! MAX_ATB - Maximum number of intervals in the atb distributions.
! MAX_TOP - Maximum number of topindex intervals per pixel.
!	    If considering fractional coverage of rainfall than this
!	    number should equal MAX_ATB.
!	    If assuming rainfall to be uniformly distributed over
!	    each grid than this paramters should equal 1.
! MAX_TST - Maximum number of time steps to be solved for.
! MOS_FLG - 1 if there is a moss layer in at least one land cover
!	    class.  zero if all land cover classes have no moss layer.
! UST_FLG - 1 if there is an understory layer in at least one land cover
!	    class.  zero if all land cover classes have no understory layer.
! LAK_FLG - 1 if there is at least one land cover class that is a lake.
!           0 if none of the land cover classes are lakes.
! ====================================================================

! ====================================================================
! LAKE.h
! ====================================================================

      integer MAX_NOD
      parameter ( MAX_NOD = 5 )

! ====================================================================
! GIS.h
! ====================================================================

! TOPMODEL parameters

      real*8 q0(MAX_CAT),ff(MAX_CAT),zw(MAX_PIX),atanb(MAX_PIX)
      real*8 xlamda(MAX_CAT),area(MAX_CAT),dd(MAX_CAT),xlength(MAX_CAT)
      real*8 psicav(MAX_CAT),basink(MAX_CAT),dtil(MAX_CAT)
      real*8 zbar(MAX_CAT),zbar1(MAX_CAT),zbarrg,zbar1rg,qb0(MAX_CAT)

! Soil moisture variables

      real*8 rzsmold
      real*8 r_mossm1(1+MOS_FLG*(MAX_PIX-1))
      real*8 tzsmold,smold

! Drainage parameters

      real*8 capflx

END MODULE MODULE_VARIABLES_OLD
