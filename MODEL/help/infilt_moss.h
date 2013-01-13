      integer istmst_moss,inc_frozen,intstm_moss,irntyp

      real*8 pnet,dt,PackWater_us,SurfWater_us,Swq_us,Outflow_us
      real*8 cuminf,rzsmst,rzsm,rzsm_u,thetas,thetar,tolinf,sorp,xk0
      real*8 psic,bcgamm,bcbeta,deltrz,cc,zw,xinact,satxr,xinfxr
      real*8 xinfcp,runtot,r_moss_depth,thetas_moss,zmoss
      real*8 r_mossm,r_mossm_u,r_mossm_f
      real*8 zero,one,two,three,four,five,six,precipitation
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
