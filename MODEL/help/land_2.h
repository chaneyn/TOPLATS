      integer ievcon,ivgtyp,ioppet,iffroz,iopthermc,ifcoarse
      integer iopgveg,iopthermc_v,i_2l,i,inc_frozen,ipix,initer
      integer maxnri

      real*8 ccc,rain,snow,tt,r_diff,canclos,xlhv,row,xleactd
      real*8 evtact,bsdew,tkmid,zmid,zrzmax,smtmp,rzsm,tzsm
      real*8 smold,rzsmold,tzsmold,thermc1,thermc2,thetar,thetas
      real*8 psic,bcbeta,quartz,heatcap1,heatcap2,heatcapold
      real*8 rocpsoil,cph2o,roa,cp,roi,thermc,heatcap
      real*8 rzdthetaudtemp,tcbeta,xlai,tkact,tkactd,tkmidactd
      real*8 f1,f2,f3,emiss,rescan,ravd,rahd,rnactd,hactd,gactd
      real*8 dshactd,tcel,vppa,psychr,zdeep,Tdeepstep,rsd,r_lup
      real*8 rld,toleb,dt,albd,r_sdn,rnetpn,gbspen,rnetd
      real*8 xled,hd,gd,dshd,tkd,tkmidd,rnact,xleact,hact,gact
      real*8 dshact,rnetw,xlew,hw,gw,dshw,tkw,tkmidw,dc,fw,tdiff
      real*8 Swq,precip_o,PackWater,SurfWater,TPack,TSurf,r_MeltEnergy
      real*8 Outflow,xleact_snow,dens,hact_snow,rn_snow,alb_snow,za
      real*8 zpd,z0h,RaSnow,appa,vpsat,uzw
      real*8 zero,one,two,three,four,five,six,dum1,dum2,dum3,dum4,dum5
      real*8 dum6,dum7,dum8,dum9,dum10,dum11,dum12,albsnow
      real*8 VaporMassFlux
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
