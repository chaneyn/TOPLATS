      integer ipix,i,inc_frozen,i_2l,lakpix,i_und,i_moss,ivgtyp
      integer iopthermc,iopgveg,iopthermc_v,iopsmini,ikopt
      integer irestype,ioppet,iopveg,iopstab,iopwv
      integer istmst,intstm,istmst_moss,intstm_moss,intstp
      integer istorm,intstp_moss,istorm_moss
      integer maxnri,ievcon,ievcon_us,ievcon_moss
      integer ifcoarse,irntyp,idifind,newstorm
      real*8 snow,rain,dsty,dsty_us,Sdepth,Sdepth_us

      integer mixmax,numnod

      real*8 mul_fac
      real*8 dt,f_moss,f_und,endstm,xintst_moss,xintst
      real*8 PackWater,SurfWater,Swq,VaporMassFlux
      real*8 TPack,TSurf,r_MeltEnergy,Outflow,xleact_snow
      real*8 hact_snow,rn_snow,PackWater_us,SurfWater_us
      real*8 Swq_us,VaporMassFlux_us,TPack_us,TSurf_us
      real*8 r_MeltEnergy_us,Outflow_us,xleact_snow_us,hact_snow_us
      real*8 rn_snow_us,dens,dens_us,albd_us,alb_moss,alb_snow
      real*8 albd,albw,albw_us,rsd,rld,tdry,Tincan,twet,twet_ic
      real*8 rh,rh_ic,zww,za,uzw,press,pptms,precip_o,precip_u
      real*8 Tslope1,Tint1,Tslope2,Tint2,Tsep,Twslope1,Twint1
      real*8 Twslope2,Twint2,Twsep,tkmid,tkact,tkmid_us,tkact_us
      real*8 tskinact_moss,tkact_moss,tkmid_moss,tkmidpet,tkmidpet_us
      real*8 tkmidpet_moss,Tdeepstep,tkpet,amp,phase,shift,tdeep
      real*8 tmid0,tmid0_moss,tk0moss,tkpet_us,dshact,rnetpn
      real*8 gbspen,evtact,bsdew,gact,rnact,xleact,hact
      real*8 dshact_moss,ebspot,rnact_us,xleact_us
      real*8 hact_us,gact_us,dshact_us,evtact_us,evtact_moss
      real*8 rnact_moss,xleact_moss,hact_moss,gact_moss,dspet
      real*8 dspet_us,dspet_moss,epwms,epwms_us,rnpet,xlepet
      real*8 hpet,gpet,xlepet_us,hpet_us,gpet_us,rnpet_moss
      real*8 gpet_moss,hpet_moss,xlepet_moss,tkpet_moss,tskinpet_moss
      real*8 rnpet_us,thetar,thetas,psic,bcbeta,quartz
      real*8 rocpsoil,tcbeta,tcbeta_us,bulk_dens,a_ice,b_ice
      real*8 xk0,bcgamm,srespar1,srespar2,srespar3,zdeep
      real*8 zmid,zrzmax,r_moss_depth,thetas_moss,srespar1_moss
      real*8 srespar2_moss,srespar3_moss,eps,emiss_moss,zpd_moss
      real*8 rib_moss,z0m_moss,z0h_moss,xlai,xlai_us,emiss,zpd
      real*8 zpd_us,z0m,z0h,z0m_us,z0h_us,rescan,tc,tw,tc_us,tw_us
      real*8 rescan_us,rtact,rtdens,psicri,respla,emiss_us,rsmax
      real*8 rsmin_us,rsmax_us,Rpl,Rpl_us,f3vpdpar,f3vpdpar_us
      real*8 rsmin,trefk,f4temppar,trefk_us,f4temppar_us,rib
      real*8 rib_us,rzsm,tzsm,rzsm1,tzsm1,rzsm_u,tzsm_u,rzsm1_u
      real*8 tzsm1_u,rzsm_f,tzsm_f,rzsm1_f,tzsm1_f,r_mossm
      real*8 r_mossm1,r_mossm_f,r_mossm1_f,r_mossm_u,r_mossm1_u
      real*8 zrz,ztz,smold,rzsmold,tzsmold,rzdthetaudtemp
      real*8 rzdthetaidt,tzdthetaidt,zw,zbar,zmoss,capflx
      real*8 difrz,diftz,grz,gtz,pnet,cuminf,sorp,cc,xinact
      real*8 satxr,xinfxr,runtot,sesq,corr
      real*8 dc,fw,dc_us,fw_us,wcip1,par,smpet0
      real*8 wsc,dsrz,rzrhs,dstz,tzrhs,wcip1_us,wsc_us,dswc,wcrhs
      real*8 dswc_us,wcrhs_us,ff,atanb,xlamda,fwcat,fwreg,pr3sat
      real*8 perrg2,pr2sat,pr2uns,perrg1,pr1sat,pr1rzs,pr1tzs,pr1uns
      real*8 persxr,perixr,persac,peruac,perusc
      real*8 canclos,extinct,cp,cph2o,row,roi,toleb

      real*8 tcel,vppa,psychr,xlhv,tkel,appa,vpsat,tcel_ic
      real*8 vppa_ic,psychr_ic,xlhv_ic,tkel_ic,vpsat_ic,qv
      real*8 qv_ic,ra,ra_ic,epetd,epetd_us,xle_act_moss
      real*8 rnetd,xled,hd,gd,dshd,tkd,tkmidd,rnetw,xlew,hw
      real*8 gw,dshw,tkw,tkmidw,tskinactd_moss,tkactd_moss,tkmidactd_moss
      real*8 ds_p_moss,epetw,rnetw_us,xlew_us,hw_us,gw_us,dshw_us
      real*8 tkw_us,tkmidw_us,epetw_us,rnetd_us,xled_us,hd_us
      real*8 gd_us,dshd_us,tkd_us,tkmidd_us,rnet_pot_moss,xle_p_moss
      real*8 h_p_moss,g_p_moss,tk_p_moss,tkmid_p_moss,tskin_p_moss,eact_moss
      real*8 tsoilold,epet_moss,f1par,f3vpd,f4temp,f1par_us,f3vpd_us
      real*8 f4temp_us,f1,f2,f3,roa,roa_ic,ravd,rahd,ravd_us,rahd_us
      real*8 rav_moss,rah_moss,rasnow,ravw,ravw_us,rahw,rahw_us,wc
      real*8 wc_us,dc_moss,bsdew_moss,r_mossmold,deltrz,dewrun
      real*8 a_ice_moss,b_ice_moss,bulk_dens_moss

      real*8 preca_a,tp_in,hice_in,hsnw_in,tempi_a,hice_a,hsnow_a
      real*8 xlat_a,xlong_a,eta_a,fraci_a,precacc,avpx,tax,rha
      real*8 psurfx,qax,rlwdx,swx,uax

      real*8 temp_a(1+LAK_FLG*(MAX_PIX-1),1+LAK_FLG*(4))
      real*8 surface(1+LAK_FLG*(MAX_PIX-1),1+LAK_FLG*(MAX_NOD-1))
      real*8 rrr,tprev
