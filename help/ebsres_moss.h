      integer i,inc_frozen,irestype,ioppet,iopthermc,ifcoarse,i_2l
      integer ievcon_moss,maxiter,iffroz,iter,maxnri

      real*8 rsoil,r_mossm,srespar1_moss,zrzmax,srespar2_moss,smtmp
      real*8 srespar3_moss,rzsm,thetas_moss,tzsm,tkact_moss,r_mossm_u
      real*8 rav_moss,rah_moss,tkmid,smold,zmid,rzsmold,tzsmold
      real*8 thermc1,thermc2,thermc,thetar,psic,bcbeta,quartz
      real*8 heatcap1,heatcap2,heatcapold,heatcap,roi,rocpsoil,row
      real*8 cph2o,roa,cp,rzdthetaudtemp,thermc_us,tcbeta_us,thermc_moss
      real*8 xlai_us,heatcap_moss,tskinact_moss,tkmid_moss
      real*8 dshact_moss,Tdeepstep,xle_act_moss,tcel_ic,r_moss_depth,zdeep
      real*8 eps,dt,toleb,r_sdn,r_ldn,f3,canclos,vppa_ic,roa_ic
      real*8 alb_moss,psychr_ic,xlhv_ic,ebscap,emiss_moss,rld,zww
      real*8 tkel_ic,za,uzw,press,zpd_moss,rib_moss,z0m_moss
      real*8 epet_moss,z0h_moss,zmoss,evtact_moss,bsdew_moss,gbspen
      real*8 rnetpn,r_mossm_f,thetas,tolstab
      real*8 zero,one,two,three,four,five,six
      real*8 ccc,raveff,heatcap_us,dumtskin,dumtk,dumtkmid,dumds
      real*8 tskintemp,ttemp,tmidtemp,told0,told1,told2,hold0,hold1
      real*8 hold2,tcold0,tcold1,tcold2,dumh,dumg,dumxle,dum_rnet
      real*8 tacttemp,vpsat,vpdef,dvpsdt,pstar,ebsmf

      real*8 calcra
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
