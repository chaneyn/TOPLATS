      integer ievcon_moss,i_2l,ipix,initer,iffroz_us,iopthermc,ifcoarse
      integer inc_frozen,i,maxnri

      real*8 canclos,xleactd_moss,xlhv_ic,row,bsdew_moss
      real*8 evtact_moss,tkmid,zmid,zrzmax,smtmp,rzsm,tzsm
      real*8 smold,rzsmold,tzsmold,thermc1,thermc2,thetar,thetas
      real*8 psic,bcbeta,quartz,heatcap1,heatcap2,heatcapold,rocpsoil
      real*8 cph2o,roa,cp,roi,thermc,heatcap,rzdthetaudtemp,tcbeta_us
      real*8 xlai_us,thermc_us,heatcap_us,thermc_moss,heatcap_moss
      real*8 r_mossm,tskinactd_moss,tskinact_moss,tkactd_moss,tkact_moss
      real*8 tkmidactd_moss,tkmid_moss,Tdeepstep,hactd_moss,gactd_moss
      real*8 dshactd_moss,tkel_ic,rav_moss,rah_moss,r_moss_depth,zdeep
      real*8 eps,dt,toleb,r_sdn,r_ldn,f3,alb_moss,rnactd_moss,vppa_ic
      real*8 emiss_moss,roa_ic,psychr_ic,eact_moss,rld,rnet_pot_moss
      real*8 xle_p_moss,h_p_moss,g_p_moss,ds_p_moss,tk_p_moss,tkmid_p_moss
      real*8 tskin_p_moss,zmoss,thetas_moss,rnact_moss
      real*8 xleact_moss,hact_moss,gact_moss,dshact_moss,gact
      real*8 zero,one,two,three,four,five,six,trlup,ccc,tau_us
      real*8 told0,told1,told2,hold0,hold1,hold2,tcold0,tcold1
      real*8 tcold2,tdum1,tdum2,tdum3,dummy
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
