      integer ipix,i,inc_frozen,i_und,i_moss,ivgtyp
      integer iopthermc,iopgveg,iopthermc_v,iopstab,maxnri
      integer ifcoarse,istart,iffroz,istop,ii,iopsmini

      real*8 dt,canclos,extinct,PackWater,SurfWater,Swq
      real*8 VaporMassFlux,TPack,TSurf,r_MeltEnergy,Outflow
      real*8 xleact_snow,hact_snow,rn_snow,PackWater_us
      real*8 SurfWater_us,Swq_us,VaporMassFlux_us,TPack_us
      real*8 TSurf_us,r_MeltEnergy_us,Outflow_us,xleact_snow_us
      real*8 hact_snow_us,rn_snow_us,dens,dens_us,albd_us
      real*8 alb_moss,alb_snow,albd,albw,albw_us,rsd,rld,tcel
      real*8 vppa,psychr,xlhv,tkel,zww,za,uzw,press,appa,vpsat
      real*8 tcel_ic,vppa_ic,psychr_ic,xlhv_ic,tkel_ic,vpsat_ic
      real*8 tkmid,tkact,tkmid_us,tkact_us,tskinact_moss,tkact_moss
      real*8 tkmid_moss,Tdeepstep,dshact,epetd,gact,epetd_us
      real*8 dshact_moss,xle_act_moss,rnetd,xled,hd,gd,dshd,tkd
      real*8 tkmidd,rnetw,xlew,hw,gw,dshw,tkw,tkmidw
      real*8 tskinactd_moss,tkactd_moss,tkmidactd_moss,ds_p_moss
      real*8 epetw,dshact_us,rnetw_us,xlew_us,hw_us,gw_us,dshw_us
      real*8 tkw_us,tkmidw_us,epetw_us,rnetd_us,xled_us,hd_us,gd_us
      real*8 dshd_us,tkd_us,tkmidd_us,rnet_pot_moss,xle_p_moss,h_p_moss
      real*8 g_p_moss,tk_p_moss,tkmid_p_moss,tskin_p_moss,eact_moss
      real*8 thetar,thetas,psic,bcbeta,quartz,rocpsoil,tcbeta
      real*8 tcbeta_us,zdeep,zmid,zrzmax,r_moss_depth,eps,emiss_moss
      real*8 zpd_moss,rib_moss,z0m_moss,z0h_moss,epet_moss,xlai
      real*8 xlai_us,emiss,zpd,zpd_us,z0m,z0h,z0m_us,z0h_us,f1par,f3vpd
      real*8 f4temp,f1par_us,f3vpd_us,f4temp_us,rescan,rescan_us,f1
      real*8 f2,f3,emiss_us,rsmin,rsmax,rsmin_us,rsmax_us,Rpl,Rpl_us
      real*8 f3vpdpar,f3vpdpar_us,trefk,f4temppar,trefk_us
      real*8 f4temppar_us,row,cph2o,roa,cp,roi,toleb,roa_ic
      real*8 ravd,rahd,ravd_us,rahd_us,rav_moss,rah_moss
      real*8 rib,RaSnow,rib_us,ravw,ravw_us,rahw,rahw_us,rzsm
      real*8 tzsm,rzsm1,tzsm1,r_mossm,zrz,smold,rzdthetaudtemp,smpet0
      real*8 zero,one,two,three,four,five,six
      real*8 ccc,rain,snow,smtmp,thermc1,thermc2,heatcap1,heatcap2
      real*8 heatcapold,thermc,heatcap,heatcap2_us,thermc_moss
      real*8 thermc_us,heatcap_us,thermc2_us,r_mindiff,tt,rlwdn,rswdn
      real*8 rlwup,tfin,tkmidw_us_tmp,tkw_us_tmp,tkmidw_tmp,tkw_tmp,r_diff
      real*8 tkmidd_us_tmp,tskin_p_moss_tmp,tk_p_moss_tmp
      real*8 tkmid_p_moss_tmp,tkmidd_tmp,heatcap_moss

      real*8 clcf1par,clcf3vpd,clcf4temp

      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
