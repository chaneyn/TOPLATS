      integer ievcon,ivgtyp,ioppet,iffroz,iopthermc,ifcoarse,iopgveg,i_2l
      integer iopthermc_v,i,inc_frozen,ipix,initer,maxnri

      real*8 rain,snow,thermc2,heatcap,heatcap2,heatcapold
      real*8 tkactd,tkmidactd,canclos,xlhv,row,xleactd
      real*8 evtact,bsdew,tkmid,zmid,zrzmax,smtmp,rzsm,tzsm
      real*8 smold,rzsmold,tzsmold,thermc1,thetar,thetas
      real*8 psic,bcbeta,quartz,heatcap1
      real*8 rocpsoil,cph2o,roa,cp,roi,thermc,rzdthetaudtemp
      real*8 tcbeta,xlai,tkact,f1,f2,f3,emiss
      real*8 rescan,ravd,rahd,rnactd,hactd,gactd,dshactd,tcel
      real*8 vppa,psychr,zdeep,Tdeepstep,rsd,r_lup,rld,toleb,dt
      real*8 albd,r_sdn,rnetpn,gbspen,rnetd,xled,hd,gd,dshd,tkd
      real*8 tkmidd,rnact,xleact,hact,gact,dshact,rnetw,xlew,hw
      real*8 gw,dshw,tkw,tkmidw,dc,fw,tdiff,PackWater,SurfWater
      real*8 Swq,VaporMassFlux,TPack,TSurf,r_MeltEnergy,Outflow
      real*8 xleact_snow,hact_snow,dens,precip_o,za,zpd,z0h
      real*8 RaSnow,appa,vpsat,uzw,rn_snow,alb_snow
      real*8 zero,one,two,three,four,five,six,tsnow
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
