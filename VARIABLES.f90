MODULE VARIABLES

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
      real*8 dz_lake,delta,rhowat,rhosnow,rhoice,rhosurf,Le,Lei
      real*8 fusion,surf,fracmin,fraclim,qwtau,cpw_ice,beta,dm,pi
      real*8 snocrit
      parameter (dz_lake=.1)
      parameter ( MAX_NOD = 5 )

      logical iceflag
      parameter (iceflag = .true.)

      logical ar_flag
      parameter (ar_flag = .true.)
 
      parameter (delta=5.67e-8)   ! s-b constant
      parameter (rhowat = 1000.)
      parameter (rhosnow = 330.)  ! densities of water and snow 
      parameter (rhoice = 917.)   ! ice density
      parameter (rhosurf=1.275)   ! surface air density
      parameter (Le = 2.25e6)
      parameter (Lei = 2.5e6)     ! latent heats 
      parameter (fusion=3.34e5)   ! heat of fusion
      parameter (surf = 0.1)     ! surface thickness for E-B computation
      parameter (fracmin= 0.01)   ! min ice thickness in meters
      parameter (fraclim = 0.02)  ! lower limit on fractional ice cover
      parameter (qwtau=86400./2.) ! D. Pollard sub-ice time constant
      parameter (cpw_ice = 4200.) ! specific heat of ice
      parameter (beta = 0.4)      ! portion of s-w absorbed at the surface 
      real*8 kv 
      parameter (kv= 0.4)         ! vonkarman constant
      parameter (dm=1.38889E-07)  ! molecular diffusivity of water
      parameter (pi= 3.141592654)
      parameter (snocrit = 0.05)  ! for albedo, in m

! ====================================================================
! GIS.h
! ====================================================================

! Parameters describing the domain and the timing

      integer nsoil,i,ndata,m_lc,m_px,u_lc,u_px,s_lc,s_px,sw_lc,sw_px
      integer l_lc,l_px
      integer ixpix(MAX_PIX),iypix(MAX_PIX),icatch(MAX_PIX)
      integer nrow,ncol,npix,ipix,ncatch,ic
      integer ilandc(MAX_PIX),ivgtyp(MAX_VEG)
      integer i_und(1+UST_FLG*(MAX_VEG-1)),i_moss(1+MOS_FLG*(MAX_VEG-1)),i_2l
      integer inc_frozen,nlandc,ipixnum(MAX_ROW,MAX_COL),nlcs
      integer lakpix,MODE,FRCOV,img_opt,iyear,iday,ihour,dtveg
      integer iihour,iiday,iimonth,iiyear,jday
      integer ibeginhour,ibeginday,ibeginmonth,ibeginyear

      real*8 dt,pixsiz,r_lakearea(MAX_CAT),frcbeta

! Different options

      integer ioppet,iopwv,iopstab,irestype,iopbf
      integer iopthermc,iopgveg,iopthermc_v,ikopt,iopsmini
      integer iopveg,iopwt0

! Timing parameters

      real*8 r_minstep,day,hour,djday


      real*8 twet(MAX_PIX),twet_ic(MAX_PIX),rh_ic(MAX_PIX),rnetpn(MAX_PIX)
      real*8 gbspen(MAX_PIX),Tslope1(MAX_VEG),Tint1(MAX_VEG),Tincan(MAX_PIX)

! Parameters describing the storm/interstorm status

      integer intstp_moss(1+MOS_FLG*(MAX_PIX-1))
      integer istmst_moss(1+MOS_FLG*(MAX_PIX-1))
      integer intstm_moss(1+MOS_FLG*(MAX_PIX-1))
      integer istorm_moss(1+MOS_FLG*(MAX_PIX-1))
      integer intstp(MAX_PIX)
      integer istmst(MAX_PIX)
      integer intstm(MAX_PIX)
      integer istorm(MAX_PIX)

      real*8 xintst(MAX_PIX),endstm,xintst_moss(1+MOS_FLG*(MAX_PIX-1))

! TOPMODEL parameters

      real*8 q0(MAX_CAT),ff(MAX_CAT),zw(MAX_PIX),atanb(MAX_PIX)
      real*8 xlamda(MAX_CAT),area(MAX_CAT),dd(MAX_CAT),xlength(MAX_CAT)
      real*8 psicav(MAX_CAT),basink(MAX_CAT),dtil(MAX_CAT)
      real*8 zbar(MAX_CAT),zbar1(MAX_CAT),zbarrg,zbar1rg,qb0(MAX_CAT)
      real*8 wslp(MAX_PIX,2)
      integer iwel(MAX_PIX)

! Parameters describing the printing requirements of outputfiles

      integer iprn(MAX_FIL),ioutsp(MAX_FIL,MAX_SER)
      integer ioutst(MAX_FIL,MAX_SER),iouten(MAX_FIL,MAX_SER)
      integer nseries(MAX_FIL),icurser(MAX_FIL)

! Energy balance parameters

      integer maxnri

      real*8 toleb,rib(MAX_PIX),rib_us(1+UST_FLG*(MAX_PIX-1))
      real*8 rib_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 tolinf

! Soil temperature parameters

      real*8 zmid(MAX_SOI),tmid0(MAX_SOI),smpet0
      real*8 tmid0_moss(1+MOS_FLG*(MAX_VEG-1))
      real*8 tk0moss(1+MOS_FLG*(MAX_VEG-1))
      real*8 r_mossmpet0(1+MOS_FLG*(MAX_VEG-1))
      real*8 amp(MAX_SOI),phase(MAX_SOI),shift(MAX_SOI)
      real*8 zdeep(MAX_SOI),tdeep(MAX_SOI),Tdeepstep(MAX_SOI)

! Moss parameters

      real*8 bulk_dens_moss(1+MOS_FLG*(MAX_VEG-1))
      real*8 thetas_moss(1+MOS_FLG*(MAX_VEG-1))
      real*8 srespar1_moss(1+MOS_FLG*(MAX_VEG-1))
      real*8 srespar2_moss(1+MOS_FLG*(MAX_VEG-1))
      real*8 srespar3_moss(1+MOS_FLG*(MAX_VEG-1))
      real*8 alb_moss(1+MOS_FLG*(MAX_VEG-1))
      real*8 r_moss_depth(1+MOS_FLG*(MAX_VEG-1))

! Soil parameters

      integer ifcoarse(MAX_SOI),isoil(MAX_PIX)

      real*8 bcbeta(MAX_SOI),bcgamm(MAX_SOI),psic(MAX_SOI)
      real*8 thetas(MAX_SOI),thetar(MAX_SOI),quartz(MAX_SOI)
      real*8 srespar1(MAX_SOI),srespar2(MAX_SOI),srespar3(MAX_SOI)
      real*8 fbs(MAX_CAT),bulk_dens(MAX_SOI)
      real*8 fbsrg,rocpsoil(MAX_SOI),xk0(MAX_SOI)
      real*8 a_ice(MAX_SOI),b_ice(MAX_SOI),eps(1+MOS_FLG*(MAX_VEG-1))
      real*8 a_ice_moss(1+MOS_FLG*(MAX_VEG-1))
      real*8 b_ice_moss(1+MOS_FLG*(MAX_VEG-1))

! Vegetation paramters

      real*8 emiss(MAX_VEG),albd(MAX_VEG),albw(MAX_VEG)
      real*8 z0h(MAX_VEG),z0m(MAX_VEG),zpd(MAX_VEG),za(MAX_VEG)
      real*8 zww(MAX_VEG),rescan(MAX_VEG)
      real*8 emiss_us(1+UST_FLG*(MAX_VEG-1))
      real*8 albd_us(1+UST_FLG*(MAX_VEG-1))
      real*8 albw_us(1+UST_FLG*(MAX_VEG-1))
      real*8 z0h_us(1+UST_FLG*(MAX_VEG-1))
      real*8 z0m_us(1+UST_FLG*(MAX_VEG-1))
      real*8 zpd_us(1+UST_FLG*(MAX_VEG-1))
      real*8 rescan_us(1+UST_FLG*(MAX_VEG-1))
      real*8 emiss_moss(1+MOS_FLG*(MAX_VEG-1))
      real*8 z0h_moss(1+MOS_FLG*(MAX_VEG-1))
      real*8 z0m_moss(1+MOS_FLG*(MAX_VEG-1))
      real*8 zpd_moss(1+MOS_FLG*(MAX_VEG-1))
      real*8 rtres(MAX_VEG)
      real*8 tc(MAX_VEG),tw(MAX_VEG),rtact(MAX_VEG),rtdens(MAX_VEG)
      real*8 psicri(MAX_VEG),respla(MAX_VEG),xlai(MAX_VEG)
      real*8 xlai_wsc(MAX_VEG),rsmin(MAX_VEG),rsmax(MAX_VEG)
      real*8 Rpl(MAX_VEG),f3vpdpar(MAX_VEG),f4temppar(MAX_VEG)
      real*8 trefk(MAX_VEG),tcbeta(MAX_VEG)
      real*8 rsmin_us(1+UST_FLG*(MAX_VEG-1))
      real*8 xlai_us(1+UST_FLG*(MAX_VEG-1))
      real*8 xlai_wsc_us(1+UST_FLG*(MAX_VEG-1))
      real*8 rsmax_us(1+UST_FLG*(MAX_VEG-1))
      real*8 Rpl_us(1+UST_FLG*(MAX_VEG-1))
      real*8 f3vpdpar_us(1+UST_FLG*(MAX_VEG-1))
      real*8 trefk_us(1+UST_FLG*(MAX_VEG-1))
      real*8 f4temppar_us(1+UST_FLG*(MAX_VEG-1))
      real*8 tw_us(1+UST_FLG*(MAX_VEG-1))
      real*8 tcbeta_us(1+UST_FLG*(MAX_VEG-1))
      real*8 tc_us(1+UST_FLG*(MAX_VEG-1))
      real*8 rtact_us(1+UST_FLG*(MAX_VEG-1))
      real*8 rtdens_us(1+UST_FLG*(MAX_VEG-1))
      real*8 rtres_us(1+UST_FLG*(MAX_VEG-1))
      real*8 psicri_us(1+UST_FLG*(MAX_VEG-1))
      real*8 respla_us(1+UST_FLG*(MAX_VEG-1))
      real*8 canclos(MAX_VEG)

! Snow pack variables

      real*8 PackWater(1+SNOW_RUN*(MAX_PIX-1))
      real*8 SurfWater(1+SNOW_RUN*(MAX_PIX-1))
      real*8 Swq(1+SNOW_RUN*(MAX_PIX-1))
      real*8 VaporMassFlux(1+SNOW_RUN*(MAX_PIX-1))
      real*8 TPack(1+SNOW_RUN*(MAX_PIX-1))
      real*8 TSurf(1+SNOW_RUN*(MAX_PIX-1))
      real*8 r_MeltEnergy(1+SNOW_RUN*(MAX_PIX-1))
      real*8 Outflow(1+SNOW_RUN*(MAX_PIX-1))
      real*8 xleact_snow(1+SNOW_RUN*(MAX_PIX-1))
      real*8 hact_snow(1+SNOW_RUN*(MAX_PIX-1))
      real*8 rn_snow(1+SNOW_RUN*(MAX_PIX-1))
      real*8 PackWater_us(1+SNW_FLG*(MAX_PIX-1))
      real*8 SurfWater_us(1+SNW_FLG*(MAX_PIX-1))
      real*8 Swq_us(1+SNW_FLG*(MAX_PIX-1))
      real*8 VaporMassFlux_us(1+SNW_FLG*(MAX_PIX-1))
      real*8 TPack_us(1+SNW_FLG*(MAX_PIX-1))
      real*8 TSurf_us(1+SNW_FLG*(MAX_PIX-1))
      real*8 r_MeltEnergy_us(1+SNW_FLG*(MAX_PIX-1))
      real*8 Outflow_us(1+SNW_FLG*(MAX_PIX-1))
      real*8 xleact_snow_us(1+SNW_FLG*(MAX_PIX-1))
      real*8 hact_snow_us(1+SNW_FLG*(MAX_PIX-1))
      real*8 rn_snow_us(1+SNW_FLG*(MAX_PIX-1))
      real*8 alb_snow(MAX_PIX)
      real*8 precip_o(MAX_PIX)
      real*8 precip_u(MAX_PIX)
      real*8 dens(1+SNOW_RUN*(MAX_PIX-1))
      real*8 dens_us(1+SNW_FLG*(MAX_PIX-1))
      real*8 Swqsum,Swq_ussum
      real*8 Sdepth_us(1+SNW_FLG*(MAX_PIX-1))
      real*8 Sdepthsum,Sdepth_ussum
      real*8 dsty(1+SNOW_RUN*(MAX_PIX-1))
      real*8 Sdepth(1+SNOW_RUN*(MAX_PIX-1))
      real*8 dsty_us(1+SNW_FLG*(MAX_PIX-1))


! Parameters describing the understory

      real*8 extinct(MAX_VEG),f_und(1+UST_FLG*(MAX_VEG-1))
      real*8 f_moss(1+MOS_FLG*(MAX_VEG-1))

! Exfiltration parameters

      integer ievcon_us(1+UST_FLG*(MAX_PIX-1))
      integer ievcon_moss(1+MOS_FLG*(MAX_PIX-1))

      integer ievcon(MAX_PIX)
      integer idifind(MAX_SOI)

      real*8 par(MAX_SOI),corr(MAX_SOI)
      real*8 sesq(MAX_PIX),evtact(MAX_PIX),etpix(MAX_PIX),bsdew
      real*8 evtact_us(1+UST_FLG*(MAX_PIX-1))
      real*8 evtact_moss(1+MOS_FLG*(MAX_PIX-1))

! Infiltration parameters

      integer irntyp(MAX_PIX)

      real*8 cuminf(MAX_PIX),sorp(MAX_PIX),cc(MAX_PIX)
      real*8 xinact(MAX_PIX),runtot(MAX_PIX),xinfxr,satxr

! Actual energy fluxes

      real*8 rnact(MAX_PIX),xleact(MAX_PIX),hact(MAX_PIX)
      real*8 gact(MAX_PIX),tkact(MAX_PIX),tkmid(MAX_PIX)
      real*8 dshact(MAX_PIX)
      real*8 tkact_us(1+UST_FLG*(MAX_PIX-1))
      real*8 tkmid_us(1+UST_FLG*(MAX_PIX-1))
      real*8 dshact_us(1+UST_FLG*(MAX_PIX-1))
      real*8 rnact_us(1+UST_FLG*(MAX_PIX-1))
      real*8 xleact_us(1+UST_FLG*(MAX_PIX-1))
      real*8 hact_us(1+UST_FLG*(MAX_PIX-1))
      real*8 gact_us(1+UST_FLG*(MAX_PIX-1))
      real*8 tkact_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 tkmid_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 dshact_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 rnact_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 xleact_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 hact_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 gact_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 tskinact_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 rnsum,xlesum,hsum,gsum,dshsum
      real*8 tksum,tkmidsum,tkdeepsum

! Potential energy fluxes

      real*8 rnpet(MAX_PIX),xlepet(MAX_PIX),hpet(MAX_PIX)
      real*8 gpet(MAX_PIX),tkpet(MAX_PIX),tkmidpet(MAX_PIX)
      real*8 dspet(MAX_PIX),rnpet_us(1+UST_FLG*(MAX_PIX-1))
      real*8 xlepet_us(1+UST_FLG*(MAX_PIX-1))
      real*8 hpet_us(1+UST_FLG*(MAX_PIX-1))
      real*8 gpet_us(1+UST_FLG*(MAX_PIX-1))
      real*8 tkpet_us(1+UST_FLG*(MAX_PIX-1))
      real*8 tkmidpet_us(1+UST_FLG*(MAX_PIX-1))
      real*8 dspet_us(1+UST_FLG*(MAX_PIX-1))
      real*8 rnpet_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 xlepet_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 hpet_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 gpet_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 tkpet_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 tkmidpet_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 tskinpet_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 dspet_moss(1+MOS_FLG*(MAX_PIX-1))
      real*8 rnpetsum,xlepetsum,hpetsum,gpetsum,tkpetsum
      real*8 tkmidpetsum,dshpetsum

! Potential evapotranspiration variables

      real*8 ebspot(MAX_PIX)

! Canopy water balance variables

      real*8 wcip1(MAX_PIX),wsc(MAX_VEG),pnet(MAX_PIX),fwcat(MAX_CAT)
      real*8 wcip1sum,wcsum,epwms,dswc,wcrhs,dswcsum,wcrhssum
      real*8 wcip1_us(1+UST_FLG*(MAX_PIX-1)),wsc_us(1+UST_FLG*(MAX_VEG-1))
      real*8 epwms_us,dswc_us,wcrhs_us

! Catchment and regional average water balance variables

      real*8 ettot(MAX_CAT),etstsum(MAX_CAT),etwtsum(MAX_CAT)
      real*8 etbssum(MAX_CAT),etdcsum(MAX_CAT),etwcsum(MAX_CAT)
      real*8 contot(MAX_CAT),qsurf(MAX_CAT),sxrtot(MAX_CAT)
      real*8 xixtot(MAX_CAT),conrun(MAX_CAT),ranrun(MAX_CAT)
      real*8 gwtsum(MAX_CAT),capsum(MAX_CAT),pptsum(MAX_CAT)
      real*8 pnetsum(MAX_CAT),rzpsum(MAX_CAT),tzpsum(MAX_CAT)
      real*8 etlakesum(MAX_CAT)
      real*8 ettotrg,etstsumrg,etwtsumrg,etbssumrg,etdcsumrg
      real*8 etwcsumrg,contotrg,qsurfrg,sxrtotrg,xixtotrg
      real*8 conrunrg,ranrunrg,gwtsumrg,capsumrg,pptsumrg
      real*8 pnetsumrg,qbreg,dssum,svarhssum,grzsumrg,gtzsumrg
      real*8 difrzsumrg,etlakesumrg
      real*8 qb24sum

! Soil moisture variables

      real*8 rzsm(MAX_PIX),rzsm1(MAX_PIX),zrz,zrzmax,rzsmav
      real*8 dsrz,rzrhs,dsrzsum,rzrhssum,rzsmold,rzsm_u(MAX_PIX)
      real*8 rzsm1_u(MAX_PIX),rzsm_f(MAX_PIX),rzsm1_f(MAX_PIX)
      real*8 rzdthetaidt(MAX_PIX),rzdthetaudtemp(MAX_PIX)
      real*8 tzdthetaidt(MAX_PIX),zmoss(1+MOS_FLG*(MAX_PIX-1))
      real*8 r_mossm(1+MOS_FLG*(MAX_PIX-1))
      real*8 r_mossm1(1+MOS_FLG*(MAX_PIX-1))
      real*8 r_mossm_u(1+MOS_FLG*(MAX_PIX-1))
      real*8 r_mossm1_u(1+MOS_FLG*(MAX_PIX-1))
      real*8 r_mossm_f(1+MOS_FLG*(MAX_PIX-1))
      real*8 r_mossm1_f(1+MOS_FLG*(MAX_PIX-1))
      real*8 tzsm(MAX_PIX),tzsm1(MAX_PIX),ztz
      real*8 tzsm_u(MAX_PIX),tzsm1_u(MAX_PIX),tzsm_f(MAX_PIX)
      real*8 tzsm1_f(MAX_PIX),tzsmav,dstz,tzrhs,dstzsum
      real*8 tzrhssum,tzsmold,smold

! Lake parameters

      real*8 surface(1+LAK_FLG*(MAX_PIX-1),1+LAK_FLG*(MAX_NOD-1))
      real*8 temp_a(1+LAK_FLG*(MAX_PIX-1),1+LAK_FLG*(MAX_NOD-1))
      real*8 tp_in(1+LAK_FLG*(MAX_PIX-1))
      real*8 hice_in(1+LAK_FLG*(MAX_PIX-1))
      real*8 hsnw_in(1+LAK_FLG*(MAX_PIX-1))
      real*8 xlat_a(1+LAK_FLG*(MAX_PIX-1))
      real*8 xlong_a(1+LAK_FLG*(MAX_PIX-1))
      real*8 eta_a(1+LAK_FLG*(MAX_PIX-1))
      real*8 tempi_a(1+LAK_FLG*(MAX_PIX-1))
      real*8 hice_a(1+LAK_FLG*(MAX_PIX-1))
      real*8 hsnow_a(1+LAK_FLG*(MAX_PIX-1))
      real*8 preca_a(1+LAK_FLG*(MAX_PIX-1))
      real*8 fraci_a(1+LAK_FLG*(MAX_PIX-1))
      real*8 precacc(1+LAK_FLG*(MAX_PIX-1))
      real*8 Tcutoff,rhostp

      integer mixmax(1+LAK_FLG*(MAX_PIX-1))
      integer numnod(1+LAK_FLG*(MAX_PIX-1))

! Drainage parameters

      real*8 capflx,difrz,diftz,grz,gtz

! Parameters describing wetness of canopy

      real*8 dc,dc_us,fw,fw_us,fwreg

! Regional saturation parameters

      real*8 perrg1,perrg2,pr3sat
      real*8 pr2sat,pr2uns,pr1sat
      real*8 pr1rzs,pr1tzs,pr1uns
      real*8 persac,peruac,perusc,perixr,persxr

! Geographi! variables

      real*8 lat_deg,lat_min,lng_deg,lng_min,lng_mer
      real*8 rlatitude,rlongitude,rlng_merid

! File names

      character*100 fnimg(MAX_FIL)

! Miscilanuous constants

      real*8 zero,one,two,three,four,five,six
      real*8 etstore,etwt,gwt,dstore,svarhs

      integer icount(MAX_CAT)

! Water and energy balance constants

      real*8 cp,row,vkc,sbc,astab,cph2o,GRAV,roi

! Values of different constants

      data cp,row,vkc,sbc,astab,cph2o,GRAV,roi/1005.d0,&
            997.d0,0.4d0,5.6696d-8,1.d0,4186.d0,9.81d0,850.d0/
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
            3.d0,4.d0,5.d0,6.d0/
      data tolinf/1.0d-09/

! ====================================================================
! STAT.h
! ====================================================================

! Variables describing the topographi! configuration

      integer natb,int,landc,veg(MAX_PP1,MAX_VST),ivn,m_ivn,u_ivn

      real*8 atb(MAX_PP1,MAX_ATB),atb_pdf(MAX_PP1,MAX_ATB)
      real*8 veg_pdf(MAX_PP1,MAX_VST),f_lake(MAX_PP1)
      real*8 pixlamda(MAX_PP1)

! Rainfall fractional coverage paramters

      integer NEWSTORM(MAX_PP2,MAX_TST),NEWSTORM_VALUE

      real*8 FRC(MAX_PP2,MAX_TST),rainfall,FRC_VALUE,FRC_VALUE_OLD

! Variables describing the storm/interstorm situation

      integer s_istmst(MAX_ATB,MAX_VST,MAX_PP1)
      integer s_intstm(MAX_ATB,MAX_VST,MAX_PP1)
      integer s_istmst_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      integer s_intstm_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      integer s_intstp(MAX_ATB,MAX_VST,MAX_PP1)
      integer s_istorm(MAX_ATB,MAX_VST,MAX_PP1)
      integer s_intstp_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      integer s_istorm_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))

      real*8 s_xintst_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_xintst(MAX_ATB,MAX_VST,MAX_PP1)

! Snow pack variables

      real*8 s_PackWater(MAX_ATB,1+SNOW_RUN*(MAX_VST-1),1+SNOW_RUN*(MAX_PP1-1))
      real*8 s_SurfWater(MAX_ATB,1+SNOW_RUN*(MAX_VST-1),1+SNOW_RUN*(MAX_PP1-1))
      real*8 s_Swq(MAX_ATB,1+SNOW_RUN*(MAX_VST-1),1+SNOW_RUN*(MAX_PP1-1))
      real*8 s_VaporMassFlux(MAX_ATB,1+SNOW_RUN*(MAX_VST-1),1+SNOW_RUN*(MAX_PP1-1))
      real*8 s_TPack(MAX_ATB,1+SNOW_RUN*(MAX_VST-1),1+SNOW_RUN*(MAX_PP1-1))
      real*8 s_TSurf(MAX_ATB,1+SNOW_RUN*(MAX_VST-1),1+SNOW_RUN*(MAX_PP1-1))
      real*8 s_r_MeltEnergy(MAX_ATB,1+SNOW_RUN*(MAX_VST-1),1+SNOW_RUN*(MAX_PP1-1))
      real*8 s_Outflow(MAX_ATB,1+SNOW_RUN*(MAX_VST-1),1+SNOW_RUN*(MAX_PP1-1))
      real*8 s_xleact_snow(MAX_ATB,1+SNOW_RUN*(MAX_VST-1),1+SNOW_RUN*(MAX_PP1-1))
      real*8 s_hact_snow(MAX_ATB,1+SNOW_RUN*(MAX_VST-1),1+SNOW_RUN*(MAX_PP1-1))
      real*8 s_rn_snow(MAX_ATB,1+SNOW_RUN*(MAX_VST-1),1+SNOW_RUN*(MAX_PP1-1))
      real*8 s_PackWater_us(MAX_ATB,1+SNW_FLG*(MAX_VST-1),1+SNW_FLG*(MAX_PP1-1))
      real*8 s_SurfWater_us(MAX_ATB,1+SNW_FLG*(MAX_VST-1),1+SNW_FLG*(MAX_PP1-1))
      real*8 s_Swq_us(MAX_ATB,1+SNW_FLG*(MAX_VST-1),1+SNW_FLG*(MAX_PP1-1))
      real*8 s_VaporMassFlux_us(MAX_ATB,1+SNW_FLG*(MAX_VST-1),1+SNW_FLG*(MAX_PP1-1))
      real*8 s_TPack_us(MAX_ATB,1+SNW_FLG*(MAX_VST-1),1+SNW_FLG*(MAX_PP1-1))
      real*8 s_TSurf_us(MAX_ATB,1+SNW_FLG*(MAX_VST-1),1+SNW_FLG*(MAX_PP1-1))
      real*8 s_r_MeltEnergy_us(MAX_ATB,1+SNW_FLG*(MAX_VST-1),1+SNW_FLG*(MAX_PP1-1))
      real*8 s_Outflow_us(MAX_ATB,1+SNW_FLG*(MAX_VST-1),1+SNW_FLG*(MAX_PP1-1))
      real*8 s_xleact_snow_us(MAX_ATB,1+SNW_FLG*(MAX_VST-1),1+SNW_FLG*(MAX_PP1-1))
      real*8 s_hact_snow_us(MAX_ATB,1+SNW_FLG*(MAX_VST-1),1+SNW_FLG*(MAX_PP1-1))
      real*8 s_rn_snow_us(MAX_ATB,1+SNW_FLG*(MAX_VST-1),1+SNW_FLG*(MAX_PP1-1))
      real*8 s_alb_snow(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_precip_o(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_precip_u(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_dens(MAX_ATB,1+SNOW_RUN*(MAX_VST-1),1+SNOW_RUN*(MAX_PP1-1))
      real*8 s_dens_us(MAX_ATB,1+SNW_FLG*(MAX_VST-1),1+SNW_FLG*(MAX_PP1-1))
      real*8 s_dsty(MAX_ATB,1+SNOW_RUN*(MAX_VST-1),1+SNOW_RUN*(MAX_PP1-1))
      real*8 s_dsty_us(MAX_ATB,1+SNOW_RUN*(MAX_VST-1),1+SNOW_RUN*(MAX_PP1-1))
      real*8 s_Sdepth(MAX_ATB,1+SNOW_RUN*(MAX_VST-1),1+SNOW_RUN*(MAX_PP1-1))
      real*8 s_Sdepth_us(MAX_ATB,1+SNOW_RUN*(MAX_VST-1),1+SNOW_RUN*(MAX_PP1-1))

! Exfiltration parameters

      integer s_ievcon_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))
      integer s_ievcon_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      integer s_ievcon(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_sesq(MAX_ATB,MAX_VST,MAX_PP1)

! Infiltration parameters

      integer s_irntyp(MAX_ATB,MAX_VST,MAX_PP1)

      real*8 s_cuminf(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_sorp(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_cc(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_xinact(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_runtot(MAX_ATB,MAX_VST,MAX_PP1)

! Actual energy fluxes

      real*8 s_rnact(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_xleact(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_hact(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_gact(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_tkact(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_tkmid(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_dshact(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_tkact_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))
      real*8 s_tkmid_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))
      real*8 s_dshact_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))
      real*8 s_rnact_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))
      real*8 s_xleact_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))
      real*8 s_hact_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))
      real*8 s_gact_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))
      real*8 s_tkact_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_tkmid_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_dshact_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_rnact_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_xleact_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_hact_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_gact_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_tskinact_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_evtact(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_evtact_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))
      real*8 s_evtact_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))

! Potential energy fluxes

      real*8 s_rnpet(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_xlepet(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_hpet(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_gpet(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_tkpet(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_tkmidpet(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_dspet(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_rnpet_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))
      real*8 s_xlepet_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))
      real*8 s_hpet_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))
      real*8 s_gpet_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))
      real*8 s_tkpet_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))
      real*8 s_tkmidpet_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))
      real*8 s_dspet_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))
      real*8 s_rnpet_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_xlepet_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_hpet_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_gpet_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_tkpet_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_tkmidpet_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_tskinpet_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_dspet_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))

! Potential evapotranspiration variables

      real*8 s_ebspot(MAX_ATB,MAX_VST,MAX_PP1)

! Canopy water balance variables

      real*8 s_wcip1(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_wsc(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_pnet(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_wcip1_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))
      real*8 s_wsc_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))

! Soil moisture variables

      real*8 s_rzsm(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_rzsm1(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_rzsm_u(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_rzsm1_u(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_rzsm_f(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_rzsm1_f(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_rzdthetaidt(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_rzdthetaudtemp(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_tzdthetaidt(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_zmoss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_r_mossm(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_r_mossm1(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_r_mossm_u(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_r_mossm1_u(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_r_mossm_f(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_r_mossm1_f(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))
      real*8 s_tzsm(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_tzsm1(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_tzsm_u(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_tzsm1_u(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_tzsm_f(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_tzsm1_f(MAX_ATB,MAX_VST,MAX_PP1)

! Energy balance variables

      real*8 s_rib(MAX_ATB,MAX_VST,MAX_PP1)
      real*8 s_rib_us(MAX_ATB,1+UST_FLG*(MAX_VST-1),1+UST_FLG*(MAX_PP1-1))
      real*8 s_rib_moss(MAX_ATB,1+MOS_FLG*(MAX_VST-1),1+MOS_FLG*(MAX_PP1-1))

! Topmodel variables

      real*8 s_zw(MAX_ATB,MAX_VST,MAX_PP1)

! Lake variables

      real*8 s_preca_a(MAX_ATB,1+LAK_FLG*(MAX_VST-1),1+LAK_FLG*(MAX_PP1-1))
      real*8 s_tp_in(MAX_ATB,1+LAK_FLG*(MAX_VST-1),1+LAK_FLG*(MAX_PP1-1))
      real*8 s_hice_in(MAX_ATB,1+LAK_FLG*(MAX_VST-1),1+LAK_FLG*(MAX_PP1-1))
      real*8 s_hsnw_in(MAX_ATB,1+LAK_FLG*(MAX_VST-1),1+LAK_FLG*(MAX_PP1-1))
      real*8 s_temp_a(MAX_ATB,1+LAK_FLG*(MAX_VST-1),1+LAK_FLG*(MAX_PP1-1),1+LAK_FLG*(MAX_NOD-1))
      real*8 s_tempi_a(MAX_ATB,1+LAK_FLG*(MAX_VST-1),1+LAK_FLG*(MAX_PP1-1))
      real*8 s_hice_a(MAX_ATB,1+LAK_FLG*(MAX_VST-1),1+LAK_FLG*(MAX_PP1-1))
      real*8 s_hsnow_a(MAX_ATB,1+LAK_FLG*(MAX_VST-1),1+LAK_FLG*(MAX_PP1-1))
      real*8 s_xlat_a(MAX_ATB,1+LAK_FLG*(MAX_VST-1),1+LAK_FLG*(MAX_PP1-1))
      real*8 s_xlong_a(MAX_ATB,1+LAK_FLG*(MAX_VST-1),1+LAK_FLG*(MAX_PP1-1))
      real*8 s_eta_a(MAX_ATB,1+LAK_FLG*(MAX_VST-1),1+LAK_FLG*(MAX_PP1-1))
      real*8 s_surface(MAX_ATB,1+LAK_FLG*(MAX_VST-1),1+LAK_FLG*(MAX_PP1-1),1+LAK_FLG*(MAX_NOD-1))
      real*8 s_fraci_a(MAX_ATB,1+LAK_FLG*(MAX_VST-1),1+LAK_FLG*(MAX_PP1-1))
      real*8 s_precacc(MAX_ATB,1+LAK_FLG*(MAX_VST-1),1+LAK_FLG*(MAX_PP1-1))

      integer s_numnod(MAX_ATB,1+LAK_FLG*(MAX_VST-1),1+LAK_FLG*(MAX_PP1-1))
      integer s_mixmax(MAX_ATB,1+LAK_FLG*(MAX_VST-1),1+LAK_FLG*(MAX_PP1-1))

! ====================================================================
! STAT_NR.h
! ====================================================================

! Variables describing the storm/interstorm situation

      integer s_nr_istmst(MAX_TOP,MAX_LAN,MAX_PP2)
      integer s_nr_intstm(MAX_TOP,MAX_LAN,MAX_PP2)
      integer s_nr_istmst_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      integer s_nr_intstm_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      integer s_nr_intstp(MAX_TOP,MAX_LAN,MAX_PP2)
      integer s_nr_istorm(MAX_TOP,MAX_LAN,MAX_PP2)
      integer s_nr_intstp_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      integer s_nr_istorm_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))

      real*8 s_nr_xintst_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_xintst(MAX_TOP,MAX_LAN,MAX_PP2)

! Snow pack variables

      real*8 s_nr_PackWater(MAX_TOP,1+SNOW_RUN*(MAX_LAN-1),1+SNOW_RUN*(MAX_PP2-1))
      real*8 s_nr_SurfWater(MAX_TOP,1+SNOW_RUN*(MAX_LAN-1),1+SNOW_RUN*(MAX_PP2-1))
      real*8 s_nr_Swq(MAX_TOP,1+SNOW_RUN*(MAX_LAN-1),1+SNOW_RUN*(MAX_PP2-1))
      real*8 s_nr_VaporMassFlux(MAX_TOP,1+SNOW_RUN*(MAX_LAN-1),1+SNOW_RUN*(MAX_PP2-1))
      real*8 s_nr_TPack(MAX_TOP,1+SNOW_RUN*(MAX_LAN-1),1+SNOW_RUN*(MAX_PP2-1))
      real*8 s_nr_TSurf(MAX_TOP,1+SNOW_RUN*(MAX_LAN-1),1+SNOW_RUN*(MAX_PP2-1))
      real*8 s_nr_r_MeltEnergy(MAX_TOP,1+SNOW_RUN*(MAX_LAN-1),1+SNOW_RUN*(MAX_PP2-1))
      real*8 s_nr_Outflow(MAX_TOP,1+SNOW_RUN*(MAX_LAN-1),1+SNOW_RUN*(MAX_PP2-1))
      real*8 s_nr_xleact_snow(MAX_TOP,1+SNOW_RUN*(MAX_LAN-1),1+SNOW_RUN*(MAX_PP2-1))
      real*8 s_nr_hact_snow(MAX_TOP,1+SNOW_RUN*(MAX_LAN-1),1+SNOW_RUN*(MAX_PP2-1))
      real*8 s_nr_rn_snow(MAX_TOP,1+SNOW_RUN*(MAX_LAN-1),1+SNOW_RUN*(MAX_PP2-1))
      real*8 s_nr_PackWater_us(MAX_TOP,1+SNW_FLG*(MAX_LAN-1),1+SNW_FLG*(MAX_PP2-1))
      real*8 s_nr_SurfWater_us(MAX_TOP,1+SNW_FLG*(MAX_LAN-1),1+SNW_FLG*(MAX_PP2-1))
      real*8 s_nr_Swq_us(MAX_TOP,1+SNW_FLG*(MAX_LAN-1),1+SNW_FLG*(MAX_PP2-1))
      real*8 s_nr_VaporMassFlux_us(MAX_TOP,1+SNW_FLG*(MAX_LAN-1),1+SNW_FLG*(MAX_PP2-1))
      real*8 s_nr_TPack_us(MAX_TOP,1+SNW_FLG*(MAX_LAN-1),1+SNW_FLG*(MAX_PP2-1))
      real*8 s_nr_TSurf_us(MAX_TOP,1+SNW_FLG*(MAX_LAN-1),1+SNW_FLG*(MAX_PP2-1))
      real*8 s_nr_r_MeltEnergy_us(MAX_TOP,1+SNW_FLG*(MAX_LAN-1),1+SNW_FLG*(MAX_PP2-1))
      real*8 s_nr_Outflow_us(MAX_TOP,1+SNW_FLG*(MAX_LAN-1),1+SNW_FLG*(MAX_PP2-1))
      real*8 s_nr_xleact_snow_us(MAX_TOP,1+SNW_FLG*(MAX_LAN-1),1+SNW_FLG*(MAX_PP2-1))
      real*8 s_nr_hact_snow_us(MAX_TOP,1+SNW_FLG*(MAX_LAN-1),1+SNW_FLG*(MAX_PP2-1))
      real*8 s_nr_rn_snow_us(MAX_TOP,1+SNW_FLG*(MAX_LAN-1),1+SNW_FLG*(MAX_PP2-1))
      real*8 s_nr_alb_snow(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_precip_o(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_precip_u(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_dens(MAX_TOP,1+SNOW_RUN*(MAX_LAN-1),1+SNOW_RUN*(MAX_PP2-1))
      real*8 s_nr_dens_us(MAX_TOP,1+SNW_FLG*(MAX_LAN-1),1+SNW_FLG*(MAX_PP2-1))
      real*8 s_nr_dsty(MAX_TOP,1+SNOW_RUN*(MAX_LAN-1),1+SNOW_RUN*(MAX_PP2-1))
      real*8 s_nr_dsty_us(MAX_TOP,1+SNOW_RUN*(MAX_LAN-1),1+SNOW_RUN*(MAX_PP2-1))      
      real*8 s_nr_Sdepth(MAX_TOP,1+SNOW_RUN*(MAX_LAN-1),1+SNOW_RUN*(MAX_PP2-1))
      real*8 s_nr_Sdepth_us(MAX_TOP,1+SNOW_RUN*(MAX_LAN-1),1+SNOW_RUN*(MAX_PP2-1))

! Exfiltration parameters

      integer s_nr_ievcon_us(MAX_TOP,1+UST_FLG*(MAX_LAN-1),1+UST_FLG*(MAX_PP2-1))
      integer s_nr_ievcon_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      integer s_nr_ievcon(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_sesq(MAX_TOP,MAX_LAN,MAX_PP2)

! Infiltration parameters

      integer s_nr_irntyp(MAX_TOP,MAX_LAN,MAX_PP2)

      real*8 s_nr_cuminf(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_sorp(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_cc(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_xinact(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_runtot(MAX_TOP,MAX_LAN,MAX_PP2)

! Actual energy fluxes

      real*8 s_nr_rnact(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_xleact(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_hact(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_gact(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_tkact(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_tkmid(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_dshact(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_tkact_us(MAX_TOP,1+UST_FLG*(MAX_LAN-1),1+UST_FLG*(MAX_PP2-1))
      real*8 s_nr_tkmid_us(MAX_TOP,1+UST_FLG*(MAX_LAN-1),1+UST_FLG*(MAX_PP2-1))
      real*8 s_nr_dshact_us(MAX_TOP,1+UST_FLG*(MAX_LAN-1),1+UST_FLG*(MAX_PP2-1))
      real*8 s_nr_rnact_us(MAX_TOP,1+UST_FLG*(MAX_LAN-1),1+UST_FLG*(MAX_PP2-1))
      real*8 s_nr_xleact_us(MAX_TOP,1+UST_FLG*(MAX_LAN-1),1+UST_FLG*(MAX_PP2-1))
      real*8 s_nr_hact_us(MAX_TOP,1+UST_FLG*(MAX_LAN-1),1+UST_FLG*(MAX_PP2-1))
      real*8 s_nr_gact_us(MAX_TOP,1+UST_FLG*(MAX_LAN-1),1+UST_FLG*(MAX_PP2-1))
      real*8 s_nr_tkact_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_tkmid_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_dshact_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_rnact_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_xleact_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_hact_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_gact_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_tskinact_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_evtact(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_evtact_us(MAX_TOP,1+UST_FLG*(MAX_LAN-1),1+UST_FLG*(MAX_PP2-1))
      real*8 s_nr_evtact_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))

! Potential energy fluxes

      real*8 s_nr_rnpet(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_xlepet(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_hpet(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_gpet(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_tkpet(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_tkmidpet(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_dspet(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_rnpet_us(MAX_TOP,1+UST_FLG*(MAX_LAN-1),1+UST_FLG*(MAX_PP2-1))
      real*8 s_nr_xlepet_us(MAX_TOP,1+UST_FLG*(MAX_LAN-1),1+UST_FLG*(MAX_PP2-1))
      real*8 s_nr_hpet_us(MAX_TOP,1+UST_FLG*(MAX_LAN-1),1+UST_FLG*(MAX_PP2-1))
      real*8 s_nr_gpet_us(MAX_TOP,1+UST_FLG*(MAX_LAN-1),1+UST_FLG*(MAX_PP2-1))
      real*8 s_nr_tkpet_us(MAX_TOP,1+UST_FLG*(MAX_LAN-1),1+UST_FLG*(MAX_PP2-1))
      real*8 s_nr_tkmidpet_us(MAX_TOP,1+UST_FLG*(MAX_LAN-1),1+UST_FLG*(MAX_PP2-1))
      real*8 s_nr_dspet_us(MAX_TOP,1+UST_FLG*(MAX_LAN-1),1+UST_FLG*(MAX_PP2-1))
      real*8 s_nr_rnpet_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_xlepet_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_hpet_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_gpet_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_tkpet_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_tkmidpet_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_tskinpet_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_dspet_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))

! Potential evapotranspiration variables

      real*8 s_nr_ebspot(MAX_TOP,MAX_LAN,MAX_PP2)

! Canopy water balance variables

      real*8 s_nr_wcip1(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_wsc(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_pnet(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_wcip1_us(MAX_TOP,1+UST_FLG*(MAX_LAN-1),1+UST_FLG*(MAX_PP2-1))
      real*8 s_nr_wsc_us(MAX_TOP,1+UST_FLG*(MAX_LAN-1),1+UST_FLG*(MAX_PP2-1))

! Soil moisture variables

      real*8 s_nr_rzsm(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_rzsm1(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_rzsm_u(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_rzsm1_u(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_rzsm_f(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_rzsm1_f(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_rzdthetaidt(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_rzdthetaudtemp(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_tzdthetaidt(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_zmoss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_r_mossm(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_r_mossm1(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_r_mossm_u(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_r_mossm1_u(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_r_mossm_f(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_r_mossm1_f(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))
      real*8 s_nr_tzsm(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_tzsm1(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_tzsm_u(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_tzsm1_u(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_tzsm_f(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_tzsm1_f(MAX_TOP,MAX_LAN,MAX_PP2)

! Energy balance variables

      real*8 s_nr_rib(MAX_TOP,MAX_LAN,MAX_PP2)
      real*8 s_nr_rib_us(MAX_TOP,1+UST_FLG*(MAX_LAN-1),1+UST_FLG*(MAX_PP2-1))
      real*8 s_nr_rib_moss(MAX_TOP,1+MOS_FLG*(MAX_LAN-1),1+MOS_FLG*(MAX_PP2-1))

! Topmodel variables

      real*8 s_nr_zw(MAX_TOP,MAX_LAN,MAX_PP2)

! Lake variables

      real*8 s_nr_preca_a(MAX_TOP,1+LAK_FLG*(MAX_LAN-1),1+LAK_FLG*(MAX_PP2-1))
      real*8 s_nr_tp_in(MAX_TOP,1+LAK_FLG*(MAX_LAN-1),1+LAK_FLG*(MAX_PP2-1))
      real*8 s_nr_hice_in(MAX_TOP,1+LAK_FLG*(MAX_LAN-1),1+LAK_FLG*(MAX_PP2-1))
      real*8 s_nr_hsnw_in(MAX_TOP,1+LAK_FLG*(MAX_LAN-1),1+LAK_FLG*(MAX_PP2-1))
      real*8 s_nr_temp_a(MAX_TOP,1+LAK_FLG*(MAX_LAN-1),1+LAK_FLG*(MAX_PP2-1),1+LAK_FLG*(MAX_NOD-1))
      real*8 s_nr_tempi_a(MAX_TOP,1+LAK_FLG*(MAX_LAN-1),1+LAK_FLG*(MAX_PP2-1))
      real*8 s_nr_hice_a(MAX_TOP,1+LAK_FLG*(MAX_LAN-1),1+LAK_FLG*(MAX_PP2-1))
      real*8 s_nr_hsnow_a(MAX_TOP,1+LAK_FLG*(MAX_LAN-1),1+LAK_FLG*(MAX_PP2-1))
      real*8 s_nr_xlat_a(MAX_TOP,1+LAK_FLG*(MAX_LAN-1),1+LAK_FLG*(MAX_PP2-1))
      real*8 s_nr_xlong_a(MAX_TOP,1+LAK_FLG*(MAX_LAN-1),1+LAK_FLG*(MAX_PP2-1))
      real*8 s_nr_eta_a(MAX_TOP,1+LAK_FLG*(MAX_LAN-1),1+LAK_FLG*(MAX_PP2-1))
      real*8 s_nr_surface(MAX_TOP,1+LAK_FLG*(MAX_LAN-1),1+LAK_FLG*(MAX_PP2-1),1+LAK_FLG*(MAX_NOD-1))
      real*8 s_nr_fraci_a(MAX_TOP,1+LAK_FLG*(MAX_LAN-1),1+LAK_FLG*(MAX_PP2-1))
      real*8 s_nr_precacc(MAX_TOP,1+LAK_FLG*(MAX_LAN-1),1+LAK_FLG*(MAX_PP2-1))

      integer s_nr_mixmax(MAX_TOP,1+LAK_FLG*(MAX_LAN-1),1+LAK_FLG*(MAX_PP2-1))
      integer s_nr_numnod(MAX_TOP,1+LAK_FLG*(MAX_LAN-1),1+LAK_FLG*(MAX_PP2-1))

! Water table updating variables

      real*8 s_nr_zbar(MAX_CAT),s_nr_gwtsum(MAX_CAT)
      real*8 s_nr_capsum(MAX_CAT),s_nr_etwtsum(MAX_CAT)
      real*8 s_nr_rzpsum(MAX_CAT),s_nr_tzpsum(MAX_CAT)
      real*8 s_nr_zbar1(MAX_CAT),FRT
      integer numc

END MODULE VARIABLES
