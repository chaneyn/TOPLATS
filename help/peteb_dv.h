      integer i,maxnri

      real*8 thermc2,heatcap,heatcap2,heatcapold,rs_over,rain
      real*8 snow,Swq,albd,emiss,thermc,f1par,f3vpd,f4temp
      real*8 rescan,ravd,rahd,tkd,tkmidd,rnetd,xled,epetd
      real*8 hd,gd,dshd,tcel,vppa,roa,psychr,xlhv,zdeep
      real*8 Tdeepstep,zmid,rld,toleb,dt,PackWater,SurfWater,VaporMassFlux
      real*8 TPack,TSurf,r_MeltEnergy,Outflow,xleact_snow,hact_snow
      real*8 rn_snow,dens,za,zpd,z0h,RaSnow,appa,uzw,gact,alb_snow,row
      real*8 zero,one,two,three,four,five,six
      real*8 dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,dum10
      real*8 dum11,dum12,tsnow,dum,vpsat
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
