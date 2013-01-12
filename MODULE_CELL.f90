MODULE MODULE_CELL

USE MODULE_VARIABLES_OLD

USE MODULE_VARIABLES

USE MODULE_LAND

USE MODULE_ATMOS

USE MODULE_CANOPY

USE MODULE_SNOW

contains

!#####################################################################
!
!                        subroutine Update_Cells
!
!#####################################################################
!
! Solve the water and energy budget for all land surface area
!
!#####################################################################

  subroutine Update_Cells(GRID,CAT,GLOBAL,i)
  
    implicit none
    type (GRID_template),dimension(:),intent(inout) :: GRID
    type (CATCHMENT_template),dimension(:),intent(inout) :: CAT
    type (GLOBAL_template),intent(in) :: GLOBAL
    integer,intent(in) :: i
    integer :: ipix,isoil,icatch,ilandc

!#####################################################################
! Update each grid cell
!#####################################################################

    call OMP_SET_NUM_THREADS(GLOBAL%nthreads)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ipix,isoil,ilandc,icatch) 
!$OMP DO SCHEDULE(DYNAMIC) ORDERED

    do ipix=1,GLOBAL%npix

      isoil = GRID(ipix)%SOIL%isoil
      ilandc = GRID(ipix)%VEG%ilandc
      icatch = GRID(ipix)%VARS%icatch

      call Update_Cell(ipix,i,GRID(ipix)%MET,GRID(isoil)%SOIL,&
         GRID(ilandc)%VEG,GRID(ipix)%VARS,GRID(ipix)%VARS%wcip1,&
         CAT(icatch),GLOBAL)

    enddo

!$OMP END DO
!$OMP END PARALLEL

  end subroutine Update_Cells


! ====================================================================
!
!                        subroutine Update_Cell
!
! ====================================================================
!
! Solve the water and energy budget for a land surface area
!
! ====================================================================


      subroutine Update_Cell(ipix,i,GRID_MET,GRID_SOIL,GRID_VEG,&
               GRID_VARS,wcip1,CAT,GLOBAL)

      implicit none
      include 'help/land_lake.h'
      type (GRID_VEG_template) :: GRID_VEG
      type (GRID_SOIL_template) :: GRID_SOIL
      type (GRID_MET_template) :: GRID_MET
      type (GRID_VARS_template) :: GRID_VARS
      type (CATCHMENT_template) :: CAT
      type (GLOBAL_template) :: GLOBAL
      GLOBAL%mul_fac = 1.0d0

! TEMPORARY LOCATION TO PASS STRUCTURE INFORMATION TO OLD FORMAT

!VEGETATION
xlai = GRID_VEG%xlai
emiss = GRID_VEG%emiss
zpd = GRID_VEG%zpd
z0m = GRID_VEG%z0m
z0h = GRID_VEG%z0h
rescan = GRID_VEG%rescan
tc = GRID_VEG%tc
tw = GRID_VEG%tw
rtact = GRID_VEG%rtact
rtdens = GRID_VEG%rtdens
psicri = GRID_VEG%psicri
respla = GRID_VEG%respla
rsmax = GRID_VEG%rsmax
Rpl = GRID_VEG%Rpl
f3vpdpar = GRID_VEG%f3vpdpar
rsmin = GRID_VEG%rsmin
trefk = GRID_VEG%trefk
f4temppar = GRID_VEG%f4temppar
canclos = GRID_VEG%canclos
extinct = GRID_VEG%extinct
albd = GRID_VEG%albd
albw = GRID_VEG%albw
zww = GRID_VEG%zww
za = GRID_VEG%za
Tslope1 = GRID_VEG%Tslope1
Tint1 = GRID_VEG%Tint1
Tslope2 = GRID_VEG%Tslope2
Tint2 = GRID_VEG%Tint2
Tsep = GRID_VEG%Tsep
Twslope1 = GRID_VEG%Twslope1
Twint1 = GRID_VEG%Twint1
Twslope2 = GRID_VEG%Twslope2
Twint2 = GRID_VEG%Twint2
Twsep = GRID_VEG%Twsep
wsc = GRID_VEG%wsc
tcbeta = GRID_VEG%tcbeta
ivgtyp = GRID_VEG%ivgtyp

!SOIL PROPERTIES
bcbeta = GRID_SOIL%bcbeta
psic = GRID_SOIL%psic
thetas = GRID_SOIL%thetas
thetar = GRID_SOIL%thetar
xk0 = GRID_SOIL%xk0
zdeep = GRID_SOIL%zdeep
tdeep = GRID_SOIL%tdeep
zmid = GRID_SOIL%zmid
tmid0 = GRID_SOIL%tmid0
rocpsoil = GRID_SOIL%rocpsoil
quartz = GRID_SOIL%quartz
srespar1 = GRID_SOIL%srespar1
srespar2 = GRID_SOIL%srespar2
srespar3 = GRID_SOIL%srespar3
a_ice = GRID_SOIL%a_ice
b_ice = GRID_SOIL%b_ice
bulk_dens = GRID_SOIL%bulk_dens
amp = GRID_SOIL%amp
phase = GRID_SOIL%phase
shift = GRID_SOIL%shift
bcgamm = GRID_SOIL%bcgamm
par = GRID_SOIL%par
corr = GRID_SOIL%corr
idifind = GRID_SOIL%idifind
par = GRID_SOIL%par
Tdeepstep = GRID_SOIL%Tdeepstep
ifcoarse = GRID_SOIL%ifcoarse

!METEOROLOGY
rsd = GRID_MET%rsd
rld = GRID_MET%rld
tdry = GRID_MET%tdry
rh = GRID_MET%rh
uzw = GRID_MET%uzw
press = GRID_MET%press
pptms = GRID_MET%pptms

!GRID VARIABLES
!Water Balance Variables
rzsm = GRID_VARS%rzsm
tzsm = GRID_VARS%tzsm
rzsm1 = GRID_VARS%rzsm1
tzsm1 = GRID_VARS%tzsm1
rzsm_f = GRID_VARS%rzsm_f
tzsm_f = GRID_VARS%tzsm_f
rzsm1_f = GRID_VARS%rzsm1_f
tzsm1_f = GRID_VARS%tzsm1_f
rzdthetaidt = GRID_VARS%rzdthetaidt
tzdthetaidt = GRID_VARS%tzdthetaidt
rzdthetaudtemp = GRID_VARS%rzdthetaudtemp
pnet = GRID_VARS%pnet
xinact = GRID_VARS%xinact
runtot = GRID_VARS%runtot
irntyp = GRID_VARS%irntyp
!Meteorological Variables
Tincan = GRID_VARS%Tincan
rh_ic = GRID_VARS%rh_ic
precip_o = GRID_VARS%precip_o
precip_u = GRID_VARS%precip_u
!Temperature Variables
tkmid = GRID_VARS%tkmid
tkact = GRID_VARS%tkact
tkmidpet = GRID_VARS%tkmidpet
tkpet = GRID_VARS%tkpet
!Energy Fluxes
dshact = GRID_VARS%dshact
rnetpn = GRID_VARS%rnetpn
gbspen = GRID_VARS%gbspen
evtact = GRID_VARS%evtact
ievcon = GRID_VARS%ievcon
gact = GRID_VARS%gact
rnact = GRID_VARS%rnact
xleact = GRID_VARS%xleact
hact = GRID_VARS%hact
ebspot = GRID_VARS%ebspot
dspet = GRID_VARS%dspet
rnpet = GRID_VARS%rnpet
xlepet = GRID_VARS%xlepet
hpet = GRID_VARS%hpet
gpet = GRID_VARS%gpet
rib = GRID_VARS%rib
!Snow
PackWater = GRID_VARS%Packwater
SurfWater = GRID_VARS%SurfWater
VaporMassFlux = GRID_VARS%VaporMassFlux
r_MeltEnergy = GRID_VARS%r_MeltEnergy
Outflow = GRID_VARS%Outflow
PackWater_us = GRID_VARS%PackWater_us
Swq = GRID_VARS%Swq
TPack = GRID_VARS%Tpack
TSurf = GRID_VARS%TSurf
xleact_snow = GRID_VARS%xleact_snow
hact_snow = GRID_VARS%hact_snow
rn_snow = GRID_VARS%rn_snow
Swq_us = GRID_VARS%Swq_us
TPack_us = GRID_VARS%TPack_us
TSurf_us = GRID_VARS%TSurf_us
xleact_snow_us = GRID_VARS%xleact_snow_us
hact_snow_us = GRID_VARS%hact_snow_us
rn_snow_us = GRID_VARS%rn_snow_us
dens = GRID_VARS%dens
dens_us = GRID_VARS%dens_us
dsty = GRID_VARS%dsty
dsty_us = GRID_VARS%dsty_us
Sdepth = GRID_VARS%Sdepth
Sdepth_us = GRID_VARS%Sdepth_us

!wcip1 = GRID_VARS%wcip1
cuminf = GRID_VARS%cuminf
sorp = GRID_VARS%sorp
cc = GRID_VARS%cc
sesq = GRID_VARS%sesq
alb_snow = GRID_VARS%alb_snow
istmst = GRID_VARS%istmst
intstm = GRID_VARS%intstm
intstp = GRID_VARS%intstp
istorm = GRID_VARS%istorm
xintst = GRID_VARS%xintst
atanb = GRID_VARS%atanb

!Global variables
dt = GLOBAL%dt
inc_frozen = GLOBAL%inc_frozen
zrzmax = GLOBAL%zrzmax
toleb = GLOBAL%toleb
maxnri = GLOBAL%maxnri
smpet0 = GLOBAL%smpet0
endstm = GLOBAL%endstm
iopthermc = GLOBAL%iopthermc
iopgveg = GLOBAL%iopgveg
iopthermc_v = GLOBAL%iopthermc_v
iopsmini = GLOBAL%iopsmini
ikopt = GLOBAL%ikopt
irestype = GLOBAL%irestype
ioppet = GLOBAL%ioppet
iopveg = GLOBAL%iopveg
iopstab = GLOBAL%iopstab
iopwv = GLOBAL%iopwv
i_2l = GLOBAL%i_2l
newstorm = GLOBAL%newstorm
mul_fac = GLOBAL%mul_fac

!Point Data
!Water Balance
zrz = 0.d0
ztz = 0.d0
smold = 0.d0
rzsmold = 0.d0
tzsmold = 0.d0
capflx = 0.d0
difrz = 0.d0
diftz = 0.d0
grz = 0.d0
gtz = 0.d0
satxr = 0.d0
xinfxr = 0.d0
dc = 0.d0!d
fw = 0.d0!fw
dsrz = 0.d0!dsrz
rzrhs = 0.d0!rzrhs
dstz = 0.d0!dstz
tzrhs = 0.d0!tzrhs
dswc = 0.d0!dswc
wcrhs = 0.d0!wcrhs
!Energy Fluxes
epwms = 0.d0!epwms
!Constants
row = GRID_VARS%row!row
cph2o = GRID_VARS%cph2o!cph2o
cp = GRID_VARS%cp!cp
roi = GRID_VARS%roi!roi

!Catchment
zbar = CAT%zbar
ff = CAT%ff
xlamda = CAT%xlamda

! ====================================================================
! If the vegetation type is greater than or equal to zero then
! solve the water and energy balance for a land area.
! ====================================================================

      if (ivgtyp.ge.0) then

! ....................................................................
! Calculate the local energy fluxes and set
! up the storm/interstorm event times and flags.
! ..................................................................

       GRID_VARS%PackWater = PackWater        
       GRID_VARS%SurfWater = SurfWater
       GRID_VARS%Swq = Swq
       GRID_VARS%VaporMassFlux = VaporMassFlux
       GRID_VARS%r_MeltEnergy = r_MeltEnergy
       GRID_VARS%Outflow = Outflow
       GRID_VARS%PackWater_us = PackWater_us
       GRID_VARS%SurfWater_us = SurfWater_us
       GRID_VARS%Swq_us = Swq_us
       GRID_VARS%VaporMassFlux_us = VaporMassFlux_us
       GRID_VARS%r_MeltEnergy_us = r_MeltEnergy_us
       GRID_VARS%Outflow_us = Outflow_us
       GRID_VARS%PackWater = PackWater        
       GRID_VARS%SurfWater = SurfWater
       GRID_VEG%i_und = i_und
       GLOBAL%zrzmax = zrzmax
       GLOBAL%toleb = toleb
       GLOBAL%maxnri = maxnri
       GLOBAL%smpet0 = smpet0
       GLOBAL%iopthermc = iopthermc
       GLOBAL%iopgveg = iopgveg
       GLOBAL%iopthermc_v = iopthermc_v
       GLOBAL%iopstab = iopstab
       GLOBAL%ioppet = ioppet
       GLOBAL%iopwv = iopwv
       GLOBAL%iopsmini = iopsmini 


         call atmos(ipix,i,dt,inc_frozen,i_2l,&

! General vegetation parameters

       GRID_VEG,&

! Snow pack variables

       TPack,TSurf,&
       xleact_snow,hact_snow,rn_snow,&
       TPack_us,TSurf_us,&
       xleact_snow_us,hact_snow_us,rn_snow_us,dens,dens_us,&

! Albedos of the over story, under story,&
! and moss layer

        

! Meteorological data
       GRID_MET,tcel,vppa,psychr,xlhv,tkel,zww,za,uzw,press,&
       appa,vpsat,tcel_ic,vppa_ic,psychr_ic,xlhv_ic,tkel_ic,vpsat_ic,&
       Tslope1,Tint1,Tslope2,Tint2,Tsep,Tincan,tdry,Twslope1,Twint1,&
       Twslope2,Twint2,Twsep,twet_ic,twet,rh,rh_ic,qv,qv_ic,ra,ra_ic,&

! Temperature variables

       GRID_VARS,tkmid,tkact,tkmid_us,tkact_us,tskinact_moss,tkact_moss,&
       tkmid_moss,Tdeepstep,amp,phase,shift,tdeep,tmid0,tmid0_moss,tk0moss,&

! Energy fluxes and states

       dshact,epetd,gact,epetd_us,dshact_moss,xle_act_moss,rnetd,xled,hd,&
       gd,dshd,tkd,tkmidd,rnetw,xlew,hw,gw,dshw,tkw,&
       tkmidw,tskinactd_moss,tkactd_moss,tkmidactd_moss,ds_p_moss,epetw,&
       dshact_us,rnetw_us,xlew_us,hw_us,gw_us,&
       dshw_us,tkw_us,tkmidw_us,epetw_us,&
       rnetd_us,xled_us,hd_us,gd_us,dshd_us,tkd_us,&
       tkmidd_us,rnet_pot_moss,xle_p_moss,&
       h_p_moss,g_p_moss,tk_p_moss,tkmid_p_moss,tskin_p_moss,eact_moss,ebspot,&
       tsoilold,tkmidpet,tkpet,tkmidpet_us,tkmidpet_moss,&
       dspet,dspet_us,dspet_moss,rnetpn,gbspen,&

! Soil parameters

       GRID_SOIL,thetar,thetas,psic,bcbeta,quartz,ifcoarse,rocpsoil,&
       tcbeta,tcbeta_us,zdeep,zmid,zrzmax,&

! Moss parameters

       r_moss_depth,eps,emiss_moss,zpd_moss,rib_moss,&
       z0m_moss,z0h_moss,epet_moss,&

! Vegetation parameters

       xlai,xlai_us,emiss,zpd,zpd_us,z0m,z0h,z0m_us,z0h_us,&
       f1par,f3vpd,f4temp,f1par_us,f3vpd_us,f4temp_us,rescan,&
       rescan_us,f1,f2,f3,emiss_us,rsmin,rsmax,rsmin_us,rsmax_us,Rpl,&
       Rpl_us,f3vpdpar,f3vpdpar_us,trefk,f4temppar,trefk_us,f4temppar_us,&

! Constants

       row,cph2o,roa,cp,roi,toleb,maxnri,roa_ic,&

! Energy balance variables

       ravd,rahd,ravd_us,rahd_us,rav_moss,rah_moss,rib,RaSnow,rib_us,&
       ravw,ravw_us,rahw,rahw_us,&

! Water balance variables

       rzsm,tzsm,rzsm1,tzsm1,r_mossm,zrz,smold,rzdthetaudtemp,smpet0,&

! Different option paramters

       GLOBAL,iopthermc,iopgveg,iopthermc_v,iopstab,ioppet,iopwv,iopsmini)

! ....................................................................
! Calculate local wet canopy water balance.
! ....................................................................

         GRID_VARS%epetw = epetw
         GRID_VEG%rnetd = rnetd
         GRID_VEG%xled = xled
         GRID_VEG%hd = hd
         GRID_VEG%gd = gd
         GRID_VEG%dshd = dshd
         GRID_VEG%rnetw = rnetw
         GRID_VEG%xlew = xlew
         GRID_VEG%hw = hw
         GRID_VEG%gw = gw
         GRID_VEG%dshw = dshw
         GRID_VEG%tkd = tkd
         GRID_VEG%tkmidd = tkmidd
         GRID_VEG%tkw = tkw
         GRID_VEG%tkmidw = tkmidw 
 
         call canopy(ipix,GRID_VARS,GRID_VEG,GRID_MET,GLOBAL)

         epetw = GRID_VARS%epetw
         rnetd = GRID_VEG%rnetd
         xled = GRID_VEG%xled
         hd = GRID_VEG%hd
         gd = GRID_VEG%gd
         dshd = GRID_VEG%dshd
         rnetw = GRID_VEG%rnetw
         xlew = GRID_VEG%xlew
         hw = GRID_VEG%hw
         gw = GRID_VEG%gw
         dshw = GRID_VEG%dshw
         tkd = GRID_VEG%tkd
         tkmidd = GRID_VEG%tkmidd
         tkw = GRID_VEG%tkw
         tkmidw = GRID_VEG%tkmidw
! ....................................................................
! Calculate the local land surface water/energy balance.
! ....................................................................

! ....................................................................
! Option 2 : the incoming long wave radiation for both under and over
! story is equal and is the atmospheri! incoming long wave radiation.
! The uncouples the radiation balances for both layers from each
! other.  This option is also used when under story is not represented.
! ....................................................................
       
       !GRID_VARS%dc = dc
       !GRID_VARS%fw = fw
       GRID_MET%uzw = uzw
       !GRID_VARS%precip_o = precip_o
       GRID_VARS%tkmid = tkmid

       pnet = GRID_VARS%pnet
       epwms = GRID_VARS%epwms

            call land(newstorm,ipix,i,dt,inc_frozen,i_2l,&

! Factor to multiply the regional parameters with

       mul_fac,&

! General vegetation parameters

       canclos,extinct,i_und,i_moss,ivgtyp,f_moss,f_und,&

! Snow pack variables

       PackWater,SurfWater,Swq,VaporMassFlux,TPack,TSurf,&
       r_MeltEnergy,Outflow,xleact_snow,hact_snow,rn_snow,PackWater_us,&
       SurfWater_us,Swq_us,VaporMassFlux_us,TPack_us,TSurf_us,r_MeltEnergy_us,&
       Outflow_us,xleact_snow_us,hact_snow_us,rn_snow_us,dens,dens_us,&

! Albedos of the over story, under story,&
! and moss layer

       albd_us,alb_moss,alb_snow,albd,&

! Meteorological data

       rsd,rld,tcel,vppa,psychr,xlhv,tkel,zww,za,uzw,press,pptms,appa,&
       vpsat,tcel_ic,vppa_ic,psychr_ic,xlhv_ic,tkel_ic,vpsat_ic,precip_o,&
       precip_u,&

! Temperature variables

       tkmid,tkact,tkmid_us,tkact_us,tskinact_moss,tkact_moss,&
       tkmid_moss,tkmidpet,tkmidpet_us,tkmidpet_moss,tsoilold,Tdeepstep,&

! Energy fluxes

       dshact,rnetpn,gbspen,epetd,evtact,ievcon,bsdew,gact,&
       rnact,xleact,hact,epetd_us,dshact_moss,xle_act_moss,rnetd,xled,hd,&
       gd,dshd,tkd,tkmidd,rnetw,xlew,hw,gw,dshw,tkw,tkmidw,ievcon_us,rnact_us,&
       xleact_us,hact_us,gact_us,dshact_us,rnetw_us,xlew_us,hw_us,gw_us,&
       dshw_us,tkw_us,tkmidw_us,evtact_us,rnetd_us,xled_us,hd_us,gd_us,dshd_us,&
       tkd_us,tkmidd_us,ievcon_moss,bsdew_moss,evtact_moss,rnet_pot_moss,&
       xle_p_moss,h_p_moss,g_p_moss,tk_p_moss,tkmid_p_moss,&
       tskin_p_moss,eact_moss,rnact_moss,xleact_moss,hact_moss,gact_moss,&
       ds_p_moss,&

! Soil parameters

       thetar,thetas,psic,bcbeta,quartz,ifcoarse,rocpsoil,tcbeta,&
       tcbeta_us,bulk_dens,a_ice,b_ice,xk0,bcgamm,&
       srespar1,srespar2,srespar3,zdeep,zmid,zrzmax,&

! Moss parameters

       r_moss_depth,thetas_moss,srespar1_moss,srespar2_moss,srespar3_moss,&
       eps,emiss_moss,zpd_moss,rib_moss,z0m_moss,z0h_moss,epet_moss,&
       a_ice_moss,b_ice_moss,bulk_dens_moss,&

! Vegetation parameters

       xlai,xlai_us,emiss,zpd,zpd_us,z0m,z0h,&
       f1par,f3vpd,f4temp,f1par_us,f3vpd_us,f4temp_us,rescan,&
       tc,tw,tc_us,tw_us,rescan_us,rtact,rtdens,psicri,&
       respla,f1,f2,f3,emiss_us,&

! Constants

       row,cph2o,roa,cp,roi,toleb,maxnri,roa_ic,&

! Energy balance variables

       ravd,rahd,ravd_us,rahd_us,rav_moss,rah_moss,rib,RaSnow,&

! Water balance variables

       rzsm,tzsm,rzsm1,tzsm1,rzsm_u,tzsm_u,rzsm1_u,tzsm1_u,rzsm_f,&
       tzsm_f,rzsm1_f,tzsm1_f,r_mossm,r_mossm1,r_mossm_f,r_mossm1_f,r_mossm_u,&
       r_mossm1_u,zrz,ztz,r_mossmold,smold,rzsmold,tzsmold,rzdthetaudtemp,&
       rzdthetaidt,tzdthetaidt,zw,zbar,zmoss,&
       capflx,difrz,diftz,grz,gtz,pnet,cuminf,sorp,cc,deltrz,&
       xinact,satxr,xinfxr,runtot,irntyp,sesq,corr,&
       idifind,dc,fw,dc_us,fw_us,wcip1,par,dewrun,dsrz,rzrhs,dstz,tzrhs,&

! Storm parameters

       istmst,intstm,istmst_moss,intstm_moss,intstp,istorm,&

! Topmodel parameters

       ff,atanb,xlamda,&

! Regional saturation parameters

       fwcat,fwreg,pr3sat,perrg2,pr2sat,pr2uns,perrg1,pr1sat,&
       pr1rzs,pr1tzs,pr1uns,persxr,perixr,persac,peruac,perusc,&

! Different option paramters

       iopthermc,iopgveg,iopthermc_v,iopsmini,ikopt,&
       irestype,ioppet,iopveg,GRID_MET,GRID_VEG,GRID_VARS,GRID_SOIL,CAT,GLOBAL)

! ====================================================================
! Calculate the density and depth of the snow layers.
! ====================================================================

         if ( (Swq.gt.0.d0).and.(SNOW_RUN.eq.1) ) then

            call calcrain (tcel,snow,rain,precip_o,dt)
            call snow_density(dsty,snow,tcel,Swq,Sdepth,TSurf,dt)

         else

            Sdepth=0.d0
            dsty=0.d0

         endif

      endif

! ====================================================================
! In the vegetation type is lower than zero then solve the open 
! water energy and water balance.
! ====================================================================

      if (ivgtyp.eq.(-1)) then

! ....................................................................
! Calculate the deep soil temperature.
! ....................................................................

         if ( (amp.eq.(0.d0)).and.&
              (phase.eq.(0.d0)).and.&
              (shift.eq.(0.d0)) ) then

            Tdeepstep=tdeep

         else

            rrr=real(i)

            Tdeepstep=tdeep + amp*cos ( rrr*phase - shift )

         endif

      endif

!TEMPORARY - Convert variables back to structure
!GRID VARIABLES
GRID_VARS%rzsm = rzsm
GRID_VARS%tzsm = tzsm
GRID_VARS%rzsm1 = rzsm1
GRID_VARS%tzsm1 = tzsm1
GRID_VARS%rzsm_f = rzsm_f
GRID_VARS%tzsm_f = tzsm_f
GRID_VARS%rzsm1_f = rzsm1_f
GRID_VARS%tzsm1_f = tzsm1_f
GRID_VARS%rzdthetaidt = rzdthetaidt
GRID_VARS%tzdthetaidt = tzdthetaidt
GRID_VARS%rzdthetaudtemp = rzdthetaudtemp
GRID_VARS%pnet = pnet
GRID_VARS%xinact = xinact
GRID_VARS%runtot = runtot
GRID_VARS%irntyp = irntyp
!Meteorological Variables
GRID_VARS%Tincan = Tincan
GRID_VARS%rh_ic = rh_ic
GRID_VARS%precip_o = precip_o
GRID_VARS%precip_u = precip_u
!Temperature Variables
GRID_VARS%tkmid = tkmid
GRID_VARS%tkact = tkact
GRID_VARS%tkmidpet = tkmidpet
GRID_VARS%tkpet = tkpet
!Energy Fluxes
GRID_VARS%dshact = dshact
GRID_VARS%rnetpn = rnetpn
GRID_VARS%gbspen = gbspen
GRID_VARS%evtact = evtact
GRID_VARS%ievcon = ievcon
GRID_VARS%gact = gact
GRID_VARS%rnact = rnact
GRID_VARS%xleact = xleact
GRID_VARS%hact = hact
GRID_VARS%ebspot = ebspot
GRID_VARS%dspet = dspet
GRID_VARS%rnpet = rnpet
GRID_VARS%xlepet = xlepet
GRID_VARS%hpet = hpet
GRID_VARS%gpet = gpet
!Point Data
GRID_VARS%zrz = zrz
GRID_VARS%ztz = ztz
GRID_VARS%smold = smold
GRID_VARS%rzsmold = rzsmold
GRID_VARS%tzsmold = tzsmold
GRID_VARS%capflx = capflx
GRID_VARS%difrz = difrz
GRID_VARS%diftz = diftz
GRID_VARS%grz = grz
GRID_VARS%gtz = gtz
GRID_VARS%satxr = satxr
GRID_VARS%xinfxr = xinfxr
GRID_VARS%dc = dc
GRID_VARS%fw = fw
GRID_VARS%dsrz = dsrz
GRID_VARS%rzrhs = rzrhs
GRID_VARS%dstz = dstz
GRID_VARS%tzrhs = tzrhs
!GRID_VARS%dswc = dswc
!GRID_VARS%wcrhs = wcrhs
!Energy Fluxes
GRID_VARS%epwms = epwms
!Constants
GRID_VARS%row = row
GRID_VARS%cph2o = cph2o
GRID_VARS%cp = cp
GRID_VARS%roi = roi
!GRID_VARS%wcip1 = wcip1

!SNOW
GRID_VARS%SurfWater = SurfWater
GRID_VARS%Outflow = Outflow
GRID_VARS%Swq = Swq
GRID_VARS%Tpack = Tpack
GRID_VARS%TSurf = Tsurf
GRID_VARS%xleact_snow = xleact_snow
GRID_VARS%hact_snow = hact_snow
GRID_VARS%rn_snow = rn_snow
GRID_VARS%Swq_us = Swq_us
GRID_VARS%TPack_us = TPack_us
GRID_VARS%TSurf_us = TSurf_us
GRID_VARS%xleact_snow_us = xleact_snow_us
GRID_VARS%hact_snow_us = hact_snow_us
GRID_VARS%rn_snow_us = rn_snow_us
GRID_VARS%dens = dens
GRID_VARS%dens_us = dens_us
GRID_VARS%dsty = dsty
GRID_VARS%dsty_us = dsty_us
GRID_VARS%Sdepth = Sdepth
GRID_VARS%Sdepth_us = Sdepth_us

!SOIL
GRID_SOIL%Tdeepstep = Tdeepstep


      return

      end subroutine Update_Cell

END MODULE MODULE_CELL
