      integer i,iopstab,maxnri,iter

      real*8 thermc1,thermc2,heatcap_moss,heatcap1,heatcap2
      real*8 heatcap_us,rs_under,rain,snow,tskin_p_moss
      real*8 tk_p_moss,tkmid_p_moss,thermc_moss,Tdeepstep,xle_act_moss
      real*8 h_p_moss,g_p_moss,xle_p_moss,ds_p_moss,tkel_ic,rav_moss
      real*8 rah_moss,r_moss_depth,zmid,zdeep,eps,dt,toleb
      real*8 rld,alb_moss,rnet_pot_moss,vppa_ic,emiss_moss,roa_ic
      real*8 psychr_ic,xlhv_ic,epet_moss,tskinact_moss,zww
      real*8 za,uzw,zpd_moss,press,z0m_moss,rib_moss,z0h_moss
      real*8 PackWater_us,SurfWater_us,VaporMassFlux_us,TPack_us
      real*8 TSurf_us,r_MeltEnergy_us,Outflow_us,xleact_snow_us
      real*8 hact_snow_us,rn_snow_us,dens_us,RaSnow,alb_snow
      real*8 appa,tcel_ic,gact,epetw_us,row,Swq_us
      real*8 zero,one,two,three,four,five,six
      real*8 vpsat_ic,told0,told1,told2,hold0,hold1,hold2
      real*8 tcold0,tcold1,tcold2,tktmp,dum1,dum2,dum3,dum4,dum5
      real*8 dum6,dum7,dum8,dum9,dum10,dum11,dum12,tsnow
      real*8 dum

      real*8 calcra
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
