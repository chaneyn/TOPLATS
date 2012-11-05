      integer inc_frozen,irestype,iffroz,iopthermc,ifcoarse,i
      integer ievcon,ioppet,maxiter,iter,maxnri

      real*8 srespar1,tkact,srespar2,rzsm_u,srespar3,ravd,rzsm
      real*8 thetas,tkmid,zmid,zrzmax,smtmp,tzsm,smold,rzsmold
      real*8 tzsmold,thermc1,thermc2,thetar,heatcapold,psic,bcbeta
      real*8 quartz,heatcap1,heatcap2,rocpsoil,row,cph2o,roa,cp
      real*8 roi,thermc,heatcap,rzdthetaudtemp,dshact,albd,emiss
      real*8 rahd,ebscap,tcel,vppa,psychr,xlhv,zdeep,Tdeepstep,rsd
      real*8 rld,toleb,dt,tkel,zww,za,uzw,zpd,z0m,press
      real*8 rib,rnetpn,gbspen,epetd,evtact,bsdew,z0h
      real*8 zero,one,two,three,four,five,six
      real*8 tolstab
      real*8 ebsmf,pstar,dvpsdt,vpdef,vpsat,dumg,dumh,dumxle
      real*8 dumrn,tmidtemp,ttemp,dumds,dumtkmid,dumtk
      real*8 raveff,rsoil,tacttemp

      real*8 calcra
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
