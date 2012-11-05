      integer ievcon,ivgtyp,ioppet,iffroz,iopthermc,ifcoarse
      integer iopgveg,iopthermc_v,i_2l,i,inc_frozen
      integer ipix,initer,maxnri

      real*8 ccc,rain,snow,thermc2,heatcap,heatcap2,heatcapold
      real*8 tkactd,tkmidactd,canclos,xlhv,row,xleactd,evtact
      real*8 bsdew,tkmid,zmid,zrzmax,smtmp,rzsm,tzsm,smold
      real*8 rzsmold,tzsmold,thermc1,thetar,thetas,psic,bcbeta
      real*8 quartz,heatcap1,rocpsoil,cph2o,roa,cp,roi,thermc
      real*8 rzdthetaudtemp,tcbeta,xlai,tkact,f1,f2,f3,emiss
      real*8 rescan,ravd,rahd,rnactd,hactd,gactd,dshactd,tcel
      real*8 vppa,psychr,zdeep,Tdeepstep,rsd,r_lup,rld,toleb
      real*8 dt,albd,r_sdn,rnetpn,gbspen,rnetd,xled,hd
      real*8 gd,dshd,tkd,tkmidd,rnact,xleact,hact,gact,dshact
      real*8 rnetw,xlew,hw,gw,dshw,tkw,tkmidw,dc,fw,tdiff,Swq
      real*8 precip_o,alb_snow,za,zpd,z0h,RaSnow,appa,vpsat,uzw
      real*8 PackWater,SurfWater,VaporMassFlux,TPack,TSurf
      real*8 r_MeltEnergy,Outflow,xleact_snow,hact_snow,rn_snow,dens
      real*8 zero,one,two,three,four,five,six,albsnow,tsnow
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/

