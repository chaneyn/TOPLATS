      real*8 radlwd,tempice,qsen,qlat,tcutoff,sw,hi,hs,twater,qbot
      real*8 qw,ds,evapi,qnetice,fracice,evaps,dt

      real*8 tmelt,emice,t4,q0t0,x,condqw,condbar,val,evapl,tprev
      real*8 qmet,q0,qmelts,qf,qmeltb,qmeltsx,disurf,dibot,hiprv,extradi
      real*8 df,xfrac,di,diextra,extraf,val2

      parameter ( tmelt = 0.0 )
      parameter ( emice = 0.97 )
