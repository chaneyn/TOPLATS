MODULE MODULE_VARIABLES

!MAPPING STRUCTURES
type MAP_template
  integer :: ilat,ilon
endtype

!IO STRUCTURES
type FILE_template
  !Filename
  character(len=400) :: fname
  !pointer
  integer :: fp
  !spatial resolution of data
  real*8 :: spatial_res
  !undefined value
  real*8 :: undef
  !dimensions
  real*8 :: minlat,minlon
  integer :: nlat,nlon
end type FILE_template

!GRID DATA

type GRID_MET_template
        real*8 :: tdry,rh,press,pptms,rld,rsd,uzw
        type(MAP_template) :: tdry_MAP,rh_MAP,press_MAP,pptms_MAP,rld_MAP,rsd_MAP,uzw_MAP
end type GRID_MET_template

type GRID_SOIL_template
        real*8 :: bcbeta,psic,thetas,thetar,xk0,zdeep,tdeep,&
                zmid,tmid0,rocpsoil,quartz,srespar1,srespar2,srespar3,a_ice,&
                b_ice,bulk_dens,amp,phase,shift,bcgamm,par,corr,Tdeepstep
        integer :: isoil,ifcoarse,idifind
        type(MAP_template) :: MAP
        type(MAP_template) :: K0_MAP
end type GRID_SOIL_template

type GRID_VEG_template
        !Overstory
        real*8 :: xlai,xlai_wsc,albd,albw,emiss,za,zww,z0m,z0h,&
                zpd,rsmin,rsmax,Rpl,f3vpdpar,f4temppar,trefk,tcbeta,tc,tw,extinct,&
                canclos,Tslope1,Tint1,Tslope2,Tint2,Twslope1,Twint1,Twslope2,Twint2,&
                Tsep,Twsep,f_und,rtact,rtdens,rtres,psicri,respla,rescan,wsc
        !Understory
        real*8 :: xlai_us,xlai_wsc_us,albd_us,&
                albw_us,emiss_us,z0m_us,z0h_us,zpd_us,f3vpdpar_us,f4temppar_us,trefk_us,&
                rsmin_us,rsmax_us,Rpl_us,tcbeta_us,tc_us,tw_us,rtact_us,respla_us,&
                rtdens_us,rtres_us,rescan_us,wsc_us,psicri_us,wcip1_us
        !Moss layer
        real*8 :: zpd_moss,z0m_moss,z0h_moss,&
                emiss_moss,f_moss,srespar1_moss,srespar2_moss,alb_moss,&
                r_moss_depth,eps,srespar3_moss,thetas_moss,tk0moss,tmid0_moss,&
                r_mossmpet0,a_ice_moss,b_ice_moss,bulk_dens_moss
        !Miscellanous
        !real*8 :: wcip1
        !Temperature Variables
        real*8 :: tkd,tkmidd,tkw,tkmidw
        !Energy fluxes
        real*8 :: rnetd,xled,hd,gd,dshd,rnetw,xlew,hw,gw,dshw
        !Indices
        integer :: i_und,i_moss,ivgtyp,ilandc
        !Mapping NNI
        type(MAP_template) :: static_MAP
        type(MAP_template) :: dynamic_MAP

end type GRID_VEG_template

type GRID_VARS_template

        !Water Balance Variables
        real*8 :: rzsm,tzsm,rzsm1,tzsm1,rzsm_f,tzsm_f,rzsm1_f,tzsm1_f,&
                rzdthetaudtemp,rzdthetaidt,tzdthetaidt,zw,pnet,xinact,&
                runtot
        !Soil Moisture Variables
        real*8 :: tzsm1_u,rzsm1_u,zmoss,r_mossm1,r_mossm,&
                r_mossm1_u,r_mossm_u,r_mossm1_f,r_mossm_f
        !Meteorological Variables
        real*8 :: Tincan,rh_ic,precip_o,precip_u
        !Temperature Variables
        real*8 :: tkmid,tkact,tkmidpet,tkpet
        !Energy Fluxes
        real*8 :: dshact,rnetpn,gbspen,evtact,gact,rnact,xleact,&
                hact,ebspot,dspet,rnpet,xlepet,hpet,gpet
        integer :: ievcon,irntyp
        !Energy Balance
        real*8 :: rib
        !Evapotranspiration
        real*8 :: etpix
        real*8 :: epetw
        !Snow Variables
        real*8 :: PackWater_us,&
                SurfWater_us,Swq_us,VaporMassFlux_us,r_MeltEnergy_us,&
                Outflow_us,alb_snow
        real*8 :: TPack,TSurf,xleact_snow,hact_snow,rn_snow,&
                TPack_us,TSurf_us,xleact_snow_us,hact_snow_us,&
                rn_snow_us,dens,dens_us,dsty,dsty_us,Sdepth,Sdepth_us
        real*8 :: PackWater,&
                SurfWater,Swq,VaporMassFlux,r_MeltEnergy,Outflow
        !Old POINT_VARS!
        !Water balance variables
        real*8 :: zrz,ztz,smold,rzsmold,tzsmold,capflx,difrz,diftz,grz,&
                gtz,satxr,xinfxr,dc,fw,dsrz,rzrhs,dstz,tzrhs,dswc,wcrhs
        !Energy fluxes
        real*8 :: epwms
        !Constants
        real*8 :: row,cph2o,cp,roi
        !!!!!!!!!!!!!!!!!!!
        !Catchment grid
        integer :: icatch
        !STORM PARAM
        integer :: istorm,intstm,istmst,intstp
        integer :: istorm_moss,intstm_moss,istmst_moss,intstp_moss
        !TOPMODEL PARAM
        real*8 :: atanb,GSTI,TI,T0
        !INFILTRATION PARAM
        real*8 :: xintst,cuminf,sorp,cc,sesq,qb0,xintst_moss
        !Vegetation
        real*8 :: wcip1
        !NEW
        !Water balance
        real*8 :: rzsm_zrzmax !Volumetric soil water content up to zrzmax (zrz*zrzmax + thetas*(zrzmax-zrz))/zrzmax
        real*8 :: sm(2),sm1(2),sm_f(2),sm1_f(2),sm1_u(2),smdthetaidt(2),sm_old(2)
        real*8 :: z_layer(2)

end type GRID_VARS_template

type GRID_template
        type (GRID_MET_template) :: MET
        type (GRID_SOIL_template) :: SOIL
        type (GRID_VEG_template) :: VEG
        type (GRID_VARS_template) :: VARS
end type GRID_template 

!CATCHMENT DATA

type CATCHMENT_template
        !Vegetation
        real*8 :: fbs
        !Soil
        real*8 :: psicav
        !Evaporation and condensation
        real*8 :: ettot,etstsum,etwtsum,etbssum,&
                etdcsum,etwcsum,etlakesum,contot
        !Infiltration/runoff/precipitation
        real*8 :: pptsum,pnetsum,sxrtot,xixtot,&
                qsurf,ranrun,conrun
        !Vertical soil moisture fluxes and water table updating
        real*8 :: zbar,zbar1,capsum,gwtsum,&
                rzpsum,tzpsum,smpsum(2)
        !State variables
        real*8 :: fwcat
        !TOPMODEL PARAM
        real*8 :: q0,ff,dd,area,dtil,xlength,basink,xlamda,qb0
        real*8 :: zbar0,zbrpor
        !GENERALIZED TOPMODEL PARAMETERS
        real*8 :: lambda
        !Runoff
        real*8 :: qb
        !Generalized TOPMODEL paramters
        real*8 :: n
end type

!REGIONAL DATA

type REGIONAL_template
        !Vegetation
        real*8 :: fbsrg
        !Snow water equivalent sums
        real*8 :: Swqsum,Swq_ussum,Sdepthsum,Sdepth_ussum
        !Regional state variables
        real*8 :: fwreg,rzsmav,tzsmav,wcsum,wcip1sum
        !Regional evaporation sums
        real*8 :: ettotrg,etstsumrg,etwtsumrg,etbssumrg,&
                etdcsumrg,etwcsumrg,etlakesumrg
        !Regional condensation and precipitations sums
        real*8 :: pptsumrg,pnetsumrg,contotrg
        !Regional infiltration/runoff and baseflow sums
        real*8 :: sxrtotrg,xixtotrg,qsurfrg,ranrunrg,&
                conrunrg,qbreg
        !Regional vertical moisture flux sums and water table
        real*8 :: capsumrg,difrzsumrg,gwtsumrg,grzsumrg,&
                gtzsumrg,zbarrg,zbar1rg
        !Regional water balance variables
        real*8 :: dswcsum,dsrzsum,dstzsum,dssum,wcrhssum,&
                rzrhssum,tzrhssum,svarhssum
        !Regional regional actual and potential energy fluxes
        real*8 :: rnsum,xlesum,hsum,gsum,tksum,dshsum,tkmidsum,&
                tkdeepsum,rnpetsum,xlepetsum,hpetsum,gpetsum,&
                tkpetsum,tkmidpetsum,dshpetsum
        !Variables telling percent land coverage of various surface states
        real*8 :: perrg1,perrg2,pr3sat,pr2sat,pr2uns,pr1sat,pr1rzs,&
                pr1tzs,pr1uns,persac,peruac,perusc,persxr,perixr
end type REGIONAL_template

!GLOBAL PARAMETERS

type GLOBAL_template
        !Vegetation
        real*8 :: zrzmax
        !Soil
        real*8 :: smpet0
        integer :: nlayer
        !OPTIONS template
        integer :: ndata,nlandc,iopveg,inc_frozen,maxnri,iopbf,iopwt0
        integer :: ncatch,nrow,ncol,npix,i_2l,nsoil,irestype,ikopt
        integer :: ioppet,iopwv,iopstab,iopgveg,iopthermc,iopthermc_v
        integer :: iopsmini
        integer :: iopflg,istflg,iophd
        real*8 :: frcbeta
        !STORM PARAM
        real*8 :: endstm,toleb,pixsiz,dt
        integer :: newstorm
        !Vegetation
        integer :: ntdveg
        real*8 :: wc0
        integer :: dtveg
        integer :: dynamic_vegetation
        !Time
        integer :: time,itime,old_date(9),new_date(9)
        !Misc
        real*8 :: mul_fac
        !TOPMODEL
        integer :: KS_TYPE
        !OpenMP parameters
        integer :: nthreads
        !Files
        type(FILE_template) :: GENERAL_FILE
        type(FILE_template) :: TI_FILE
        type(FILE_template) :: Subbasin_FILE
        type(FILE_template) :: K0_FILE
        type(FILE_template) :: CL_table_FILE
        type(FILE_template) :: SOIL_FILE
        type(FILE_template) :: VEG_FILE
        type(FILE_template) :: DVEG_FILE
        type(FILE_template) :: FORCING_FILE
        type(FILE_template) :: OUTPUT_FILE
        type(FILE_template) :: REGIONAL_FILE
        type(FILE_template) :: CATCHMENT_FILE
        type(FILE_template) :: GSTI_FILE
        type(FILE_template) :: ET_FILE
        type(FILE_template) :: QS_FILE
  !spatial resolution of data
  real*8 :: spatial_res
  !dimensions
  real*8 :: minlat,minlon
  integer :: nlat,nlon
  !Water balance
  real*8 :: zmax_layer(2)

        
end type GLOBAL_template

type CELL_VARS_template

        real*8 tcel,vppa,psychr,xlhv,tkel,uzw,appa,vpsat,tcel_ic
        real*8 vppa_ic,psychr_ic,xlhv_ic,tkel_ic,vpsat_ic,twet_ic
        real*8 twet,qv,qv_ic,ra,ra_ic
        
        real*8 tkmid_us,tkact_us,tskinact_moss,tkact_moss,tkmid_moss
        real*8 epetd,epetd_us,dshact_moss,xle_act_moss
        real*8 tskinactd_moss,tkactd_moss,tkmidactd_moss,ds_p_moss
        real*8 dshact_us,rnetw_us,xlew_us,hw_us,gw_us
        real*8 dshw_us,tkw_us,tkmidw_us,epetw_us
        real*8 rnetd_us,xled_us,hd_us,gd_us,dshd_us,tkd_us
        real*8 tkmidd_us,rnet_pot_moss,xle_p_moss
        real*8 h_p_moss,g_p_moss,tk_p_moss,tkmid_p_moss,tskin_p_moss,eact_moss
        real*8 tsoilold,tkmidpet_us,tkmidpet_moss,dspet_us,dspet_moss
        real*8 rib_moss,epet_moss
        real*8 xleact_us,hact_us,gact_us,evtact_us,bsdew_moss
        real*8 evtact_moss,rnact_moss,xleact_moss,hact_moss,gact_moss
        real*8 rnact_us,bsdew
        real*8 f1par,f3vpd,f4temp,f1par_us,f3vpd_us
        real*8 f4temp_us,f1,f2,f3,f3vpdpar_us,f4temppar_us,roa,roa_ic
        
        real*8 ravd,rahd,ravd_us,rahd_us,rav_moss,rah_moss,RaSnow,rib_us
        real*8 ravw,ravw_us,rahw,rahw_us
        
        real*8 rzsm_u,tzsm_u,r_mossmold,deltrz,dc_us,fw_us,dewrun
        
        integer ievcon_us,ievcon_moss
end type CELL_VARS_template

!General variables
integer :: i,ipix

!General constants
real*8 zero,one,two,three,four,five,six
data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
     3.d0,4.d0,5.d0,6.d0/

END MODULE MODULE_VARIABLES
