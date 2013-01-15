      integer ipix,i,inc_frozen,i_und,i_moss,ivgtyp
      integer iopthermc,iopgveg,iopthermc_v,iopstab
      integer ifcoarse,ioppet,i_2l,iopwv,maxnri,iopsmini

      real*8 dt,canclos,extinct,PackWater,SurfWater,Swq,VaporMassFlux
      real*8 TPack,TSurf,r_MeltEnergy,Outflow,xleact_snow,hact_snow
      real*8 rn_snow,PackWater_us,SurfWater_us,Swq_us,VaporMassFlux_us
      real*8 TPack_us,TSurf_us,r_MeltEnergy_us,Outflow_us,xleact_snow_us
      real*8 hact_snow_us,rn_snow_us,dens,dens_us,albd_us,alb_moss,alb_snow
      real*8 albd,albw,albw_us,rsd,rld,tcel,vppa,psychr,xlhv,tkel,zww
      real*8 za,uzw,press,appa,vpsat,tcel_ic,vppa_ic,psychr_ic,xlhv_ic
      real*8 tkel_ic,vpsat_ic,Tslope1,Tint1,Tslope2,Tint2,Tsep,Tincan
      real*8 tdry,Twslope1,Twint1,Twslope2,Twint2,Twsep
      real*8 rh,rh_ic,qv,qv_ic,ra,ra_ic,tkmid,tkact,tkmid_us,tkact_us
      real*8 tskinact_moss,tkact_moss,tkmid_moss,Tdeepstep,amp,phase
      real*8 shift,tdeep,tmid0,tmid0_moss,tk0moss,dshact,epetd,gact
      real*8 epetd_us,dshact_moss,xle_act_moss,rnetd,xled,hd,gd,dshd
      real*8 tkd,tkmidd,rnetw,xlew,hw,gw,dshw,tkw,tkmidw,tskinactd_moss
      real*8 tkactd_moss,tkmidactd_moss,ds_p_moss,epetw,dshact_us
      real*8 rnetw_us,xlew_us,hw_us,gw_us,dshw_us,tkw_us,tkmidw_us
      real*8 epetw_us,rnetd_us,xled_us,hd_us,gd_us,dshd_us,tkd_us,tkmidd_us
      real*8 rnet_pot_moss,xle_p_moss,h_p_moss,g_p_moss,tk_p_moss
      real*8 tkmid_p_moss,tskin_p_moss,eact_moss,ebspot,tsoilold
      real*8 tkmidpet,tkpet,tkmidpet_us,tkmidpet_moss,dspet,dspet_us
      real*8 dspet_moss,rnetpn,gbspen,thetar,thetas,psic,bcbeta,quartz
      real*8 rocpsoil,tcbeta,tcbeta_us,zdeep,zmid,zrzmax,r_moss_depth
      real*8 eps,emiss_moss,zpd_moss,rib_moss,z0m_moss,z0h_moss,epet_moss
      real*8 xlai,xlai_us,emiss,zpd,zpd_us,z0m,z0h,z0m_us,z0h_us,f1par
      real*8 f3vpd,f4temp,f1par_us,f3vpd_us,f4temp_us,rescan,rescan_us
      real*8 f1,f2,f3,emiss_us,rsmin,rsmax,rsmin_us,rsmax_us,Rpl,Rpl_us
      real*8 f3vpdpar,f3vpdpar_us,trefk,f4temppar,trefk_us,f4temppar_us
      real*8 row,cph2o,roa,cp,roi,toleb,roa_ic,ravd,rahd,ravd_us
      real*8 rahd_us,rav_moss,rah_moss,rib,RaSnow,rib_us,ravw,ravw_us
      real*8 rahw,rahw_us,rzsm,tzsm,rzsm1,tzsm1,r_mossm,zrz,smold
      real*8 rzdthetaudtemp,smpet0,twet,twet_ic
      real*8 zero,one,two,three,four,five,six,rrr,rrrr,vpdef

      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
            3.d0,4.d0,5.d0,6.d0/

