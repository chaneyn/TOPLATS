      integer maxiter,itime,ipetopt,iter

      real*8 T0,T1,T2,T3,Tn0,Tn1,Tn2,tc0,tc1,tc2,tcn0,tcn1,tcn2
      real*8 heatcap1,heatcap2,heatcap3,heatcapn1,heatcapn2,heatcapn3
      real*8 r_LE_act,H,G,r_LE,ds,Ta,rav,rah,zmoss,z1,z2,eps,delt,toleb
      real*8 rsd,rld,alb,Rn,vappres,emiss,rhoa,psychr,r_lhv,epot
      real*8 SIGMA,row,Cp,d1,d2,d3,dz0,dz1,dz2,dz01,dz12
      real*8 c1,c2,cn1,cn2,dsold,dold,Told,dtest,deltnr
      real*8 dftdt,ft,dhdt,dledt,dftfac,drndt
      real*8 r_letmp,gtmp,htmp,rntmp,dstmp,vpsat,vpdef
      real*8 dsnew,ddsdt,dgdt,deltnew

      real*8 T0new,T1new,T2new
