      integer maxnri,i_und,i,i_moss,iopstab,iter

      real*8 rlwdn,rswdn,rlwup,ccc,thermc1,thermc2,thermc2_us
      real*8 heatcap1,heatcap2,heatcapold,heatcap_us,rain,snow
      real*8 tkmidd_us_tmp,tskin_p_moss_tmp,tk_p_moss_tmp,tkmid_p_moss_tmp
      real*8 Swq_us,tkmidd_us,f3,albd_us,thermc_us,emiss_us
      real*8 f1par_us,f3vpd_us,f4temp_us,rescan_us,ravd_us,rahd_us
      real*8 tkd_us,rnetd_us,xled_us,epetd_us,hd_us,gd_us,dshd_us
      real*8 zdeep,zmid,Tdeepstep,toleb,dt,tskin_p_moss,tk_p_moss
      real*8 tkmid_p_moss,heatcap_moss,thermc_moss,xle_act_moss,h_p_moss
      real*8 g_p_moss,xle_p_moss,ds_p_moss,tkel_ic,rav_moss,rah_moss
      real*8 r_moss_depth,eps,alb_moss,rnet_pot_moss
      real*8 vppa_ic,roa_ic,psychr_ic,emiss_moss,xlhv_ic,epet_moss
      real*8 tskinact_moss,zww,za,uzw,zpd_moss,z0m_moss,press,rib_moss
      real*8 z0h_moss,PackWater_us,SurfWater_us,VaporMassFlux_us,TPack_us
      real*8 TSurf_us,r_MeltEnergy_us,Outflow_us,xleact_snow_us
      real*8 hact_snow_us,rn_snow_us,dens_us,alb_snow,zpd_us,z0h_us
      real*8 RaSnow,tcel_ic,vpsat_ic,gact,row,heatcap2_us
      real*8 zero,one,two,three,four,five,six
      real*8 told0,told1,told2,hold0,hold1,hold2
      real*8 tcold0,tcold1,tcold2,tktmp
      real*8 dum1,dum2,dum3,dum4,dum5,dum6,dum7
      real*8 dum8,dum9,dum10,dum11,dum12,albsnow,appa,tsnow,dum
      real*8 tskinactd_moss,tkactd_moss,tkmidactd_moss
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/

      real*8 calcra
