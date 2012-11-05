      integer maxnri,i

      real*8 rlwdn,rswdn,rlwup,ccc,thermc2,heatcap,heatcap2
      real*8 heatcapold,tkmidd_tmp,r_diff,tt,rain,snow,Swq
      real*8 tkmidd,f1,f2,f3,emiss,thermc,f1par,f3vpd,f4temp
      real*8 rescan,ravd,rahd,tkd,rnetd,xled,epetd,hd,gd,dshd
      real*8 tcel,vppa,roa,psychr,xlhv,zdeep,Tdeepstep,zmid
      real*8 rld,toleb,dt,PackWater,SurfWater,VaporMassFlux
      real*8 TPack,TSurf,r_MeltEnergy,Outflow,xleact_snow,hact_snow
      real*8 rn_snow,dens,alb_snow,za,zpd,z0h,RaSnow,appa,vpsat
      real*8 uzw,gact,row
      real*8 zero,one,two,three,four,five,six
      real*8 dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9
      real*8 dum10,dum11,dum12,tsnow,dum,albsnow
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
