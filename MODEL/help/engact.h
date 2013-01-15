      integer ievcon,inc_frozen,ipix,initer,ivgtyp,ioppet,iffroz
      integer iopthermc,ifcoarse,iopgveg,iopthermc_v,i_2l,i
      integer maxnri

      real*8 canclos,xlhv,row,xleactd,evtact,bsdew,tkmid,zmid
      real*8 zrzmax,smtmp,rzsm,tzsm,smold,rzsmold,tzsmold
      real*8 thermc1,thermc2,thetar,thetas,psic,bcbeta,quartz
      real*8 heatcap1,heatcap2,heatcapold,rocpsoil,cph2o,roa,cp
      real*8 roi,thermc,heatcap,rzdthetaudtemp,tcbeta,xlai
      real*8 tkact,tkactd,tkmidactd,f1,f2,f3,emiss,rescan,ravd
      real*8 rahd,rnactd,hactd,gactd,dshactd,tcel,vppa,psychr
      real*8 zdeep,Tdeepstep,rsd,r_lup,rld,toleb,dt,albd
      real*8 r_sdn,rnetpn,gbspen,rnetd,xled,hd,gd,dshd,tkd,tkmidd
      real*8 rnact,xleact,hact,gact,dshact,rnetw,xlew,hw,gw,dshw
      real*8 tkw,tkmidw,dc,fw,tdiff
      real*8 zero,one,two,three,four,five,six,ccc,tdum1,tdum2
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
