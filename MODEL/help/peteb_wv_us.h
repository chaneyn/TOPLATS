      integer maxnri,i,iopstab,iter

      real*8 thermc2,heatcap_us,heatcap2,heatcapold,rs_under
      real*8 rain,snow,Swq_us,albw_us,dt,emiss_us,thermc_us
      real*8 ravw_us,rahw_us,tkw_us,tkmidw_us,rnetw_us,xlew_us,epetw_us
      real*8 hw_us,gw_us,dshw_us,tcel_ic,vppa_ic,roa_ic,psychr_ic
      real*8 xlhv_ic,zdeep,Tdeepstep,zmid,rld,toleb,tkact_us
      real*8 zww,za,uzw,zpd_us,z0m_us,tkel_ic,press,rib_us,z0h_us
      real*8 PackWater_us,SurfWater_us,VaporMassFlux_us,TPack_us,TSurf_us
      real*8 r_MeltEnergy_us,Outflow_us,xleact_snow_us,hact_snow_us
      real*8 rn_snow_us,dens_us,row,RaSnow,alb_snow,appa,vpsat_ic,gact
      real*8 zero,one,two,three,four,five,six,tktmp
      real*8 dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,dum10
      real*8 dum11,dum12,tsnow,dum

      real*8 calcra
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
