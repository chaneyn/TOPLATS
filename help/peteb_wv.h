      integer maxnri,i,iopstab,iter

      real*8 thermc2,heatcap,heatcap2,heatcapold,rs_over,rain
      real*8 snow,Swq,albw,emiss,thermc,ravw,rahw,tkw,tkmidw
      real*8 rnetw,xlew,epetw,hw,gw,dshw,tcel,vppa,roa,psychr
      real*8 xlhv,zdeep,Tdeepstep,zmid,rld,toleb,dt
      real*8 tkact,zww,za,uzw,zpd,z0m,tkel,press,rib
      real*8 z0h,PackWater,SurfWater,VaporMassFlux,TPack,TSurf,r_MeltEnergy
      real*8 Outflow,xleact_snow,hact_snow,rn_snow,dens,RaSnow
      real*8 alb_snow,appa,vpsat,gact,row
      real*8 zero,one,two,three,four,five,six
      real*8 dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,dum10
      real*8 dum11,dum12,tsnow,dum,tktmp

      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
