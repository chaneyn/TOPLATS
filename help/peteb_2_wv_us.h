      integer maxnri,i_und,i,iopstab,i_moss,iter

      real*8 zero,one,two,three,four,five,six,rlwdn,rswdn
      real*8 rlwup,ccc,thermc1,thermc2,thermc2_us,heatcap1,heatcap2
      real*8 heatcapold,heatcap_us,heatcap_moss,rain,snow
      real*8 tkmidw_us_tmp,tkw_us_tmp,tskin_p_moss_tmp,tk_p_moss_tmp
      real*8 tkmid_p_moss_tmp,tkactd_moss,tskinactd_moss
      real*8 tkmidactd_moss,Swq_us,tkw_us,tkmidw_us,f3
      real*8 albw_us,emiss_us,thermc_us,ravw_us,rahw_us,rnetw_us
      real*8 xlew_us,epetw_us,hw_us,gw_us,dshw_us,tcel_ic,vppa_ic
      real*8 roa_ic,psychr_ic,xlhv_ic,dt,zdeep,Tdeepstep,zmid,toleb
      real*8 tkact_us,zww,za,uzw,zpd_us,z0m_us,tkel_ic,press,rib_us
      real*8 z0h_us,tskin_p_moss,tk_p_moss,tkmid_p_moss,xle_act_moss,h_p_moss
      real*8 g_p_moss,xle_p_moss,ds_p_moss,eps,alb_moss,rnet_pot_moss
      real*8 emiss_moss,epet_moss,row,tskinact_moss,z0m_moss,rib_moss
      real*8 z0h_moss,PackWater_us,SurfWater_us,VaporMassFlux_us,TPack_us
      real*8 TSurf_us,r_MeltEnergy_us,Outflow_us,xleact_snow_us
      real*8 hact_snow_us,rn_snow_us,dens_us,alb_snow,RaSnow
      real*8 appa,vpsat_ic,gact,r_moss_depth,zpd_moss,rav_moss
      real*8 rah_moss,thermc_moss,albsnow
      real*8 dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,dum10
      real*8 dum11,dum12,tsnow,told0,told1,told2,hold0,hold1,hold2
      real*8 tcold0,tcold1,tcold2,tktmp,dum

      real*8 calcra

      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/

