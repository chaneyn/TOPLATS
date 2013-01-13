      integer i_moss,i_und,istmst,inc_frozen,intstm,irntyp

      real*8 pnet,PackWater_us,SurfWater_us,Swq_us,Outflow_us,dt
      real*8 PackWater,SurfWater,Swq,Outflow,cuminf
      real*8 rzsmst,rzsm,rzsm_u,thetas,thetar,tolinf,sorp,xk0
      real*8 psic,bcgamm,bcbeta,deltrz,cc,zw,xinact,satxr,xinfxr
      real*8 xinfcp,runtot,precipitation
      real*8 zero,one,two,three,four,five,six
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
