      integer maxnri,i

      real*8 thermc2,heatcap_us,heatcap2,heatcapold,rs_under
      real*8 rain,snow,albd_us,ravd_us,rahd_us,f1par_us
      real*8 f3vpd_us,f4temp_us,rescan_us,tkd_us,tkmidd_us,rnetd_us
      real*8 xled_us,epetd_us,hd_us,gd_us,dshd_us,tcel_ic,vppa_ic
      real*8 roa_ic,psychr_ic,xlhv_ic,zdeep,zmid,Tdeepstep
      real*8 toleb,dt,PackWater_us,Swq_us,SurfWater_us,TPack_us
      real*8 VaporMassFlux_us,TSurf_us,r_MeltEnergy_us,Outflow_us
      real*8 xleact_snow_us,hact_snow_us,rn_snow_us,dens_us,RaSnow
      real*8 za,zpd_us,z0h_us,alb_snow,appa,vpsat_ic,uzw,gact,row
      real*8 thermc_us,rld,emiss_us
      real*8 zero,one,two,three,four,five,six
      real*8 dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,dum10
      real*8 dum11,dum12,tsnow,dum
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
