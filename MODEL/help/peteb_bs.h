      integer maxnri,i

      real*8 thermc1,thermc2,heatcap1,heatcap2,heatcapold,rain
      real*8 snow,Swq,albd,emiss,ravd,rahd,tkd,tkmidd,rnetd
      real*8 xled,epetd,hd,gd,dshd,tcel,vppa,roa,psychr,xlhv
      real*8 zdeep,Tdeepstep,zmid,rsd,rld,toleb,dt,tkw,tkmidw
      real*8 rnetw,xlew,epetw,hw,gw,dshw,ravw,rahw,PackWater
      real*8 SurfWater,VaporMassFlux,TPack,TSurf,r_MeltEnergy
      real*8 Outflow,xleact_snow,hact_snow,rn_snow,dens,za,zpd
      real*8 albw,z0h,RaSnow,alb_snow,appa,vpsat,uzw,gact,row
      real*8 zero,one,two,three,four,five,six
      real*8 dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9
      real*8 dum10,dum11,dum12,tsnow,dum
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
