      integer ntimes

      real*8 sw,tp,tcutoff,hice,hsnow,ta,qa,ua
      real*8 psurf,rhosurf,delq,evap,hsen,rlwd
      real*8 Le,Lei,stefbol,emice
      real*8 x,tcutc,a,b,tposs,t,tlat,qsen,qmet
      real*8 qlat,t4,condbar,val,val2,q0,t0

      parameter ( Le = 2.25e6, Lei = 2.5e6 )
      parameter ( stefbol = 5.6696e-8 , emice = 0.97 )
