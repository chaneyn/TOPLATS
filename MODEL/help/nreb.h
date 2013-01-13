      integer ipetopt,itime,itmpflg,ihflg,ileflg,iconv,iter,maxnri

      real*8 alb,emiss,thermc1,thermc2,heatcap1,heatcap2,heatcapold
      real*8 rav,rah,tskink,tmidk,rn,xle,epot,h,g,ds
      real*8 tairc,vppa,roa,psychr,xlhv,zdeep
      real*8 tdeep,zmid,rsd,rld,toleb,dt
      real*8 cp,row,vkc,sbc,astab,cph2o,HMIN,DTMAX,GRAV
      real*8 zero,one,two,three,four,five,six
      real*8 deltnr,dftdt,dstmp,ft,tmidknew,gtmp,htmp
      real*8 dxledt,xletmp,drndt,rntmp,dftfac,vpsat,vpdef
      real*8 ddsdt,dgdt,gdenom,dzdeep,dhdt,dsold
      real*8 tskinc,toldmidk,toldk
      real*8 T0new,T1new,deltnew
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
