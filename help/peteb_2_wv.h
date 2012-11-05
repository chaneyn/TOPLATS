      integer maxnri,i,iopstab,iter

      real*8 rlwdn,rswdn,rlwup,ccc,thermc2,heatcap,heatcap2
      real*8 heatcapold,tkmidw_tmp,tkw_tmp,r_diff,tt,rain,snow
      real*8 Swq,tkw,tkmidw,f1,f2,f3,thermc,dt,emiss,ravw,rahw
      real*8 rnetw,xlew,epetw,hw,gw,dshw,tcel,vppa,roa
      real*8 psychr,xlhv,zdeep,rld,Tdeepstep,zmid,toleb,tkact,zww
      real*8 za,uzw,zpd,z0m,press,tkel,rib,z0h,PackWater,SurfWater
      real*8 VaporMassFlux,TPack,TSurf,r_MeltEnergy,Outflow
      real*8 xleact_snow,hact_snow,rn_snow,dens,alb_snow
      real*8 RaSnow,appa,vpsat,gact,row
      real*8 zero,one,two,three,four,five,six
      real*8 dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9
      real*8 dum10,dum11,dum12,tsnow,dum,albsnow,tktmp

      real*8 calcra

      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/

