      integer lakpix,mixmax,numnod

      real*8 ta,ua,qa,psurf,sw,rlwd,lnet,rh,Qe,Qh
      real*8 xlat,xlong,eta_a,dt
      real*8 tempi_a,hice_a,hsnow_a,preca_a,fraci_a,precacc,rnet
      real*8 surface(MAX_PIX,MAX_NOD),temp_a(MAX_PIX,MAX_NOD)

      integer i_shuf,iwater,mixdep,k

      real*8 luw,lui,lnetw,lneti,lu,tin
      real*8 T(MAX_NOD,2), de(MAX_NOD), dnsty(MAX_NOD)
      real*8 Ti(MAX_NOD,2),sworig,Tcutoff,rhostp,taC,prec
      real*8 addprec,albs,albi,albw,swi,sww,Tcutk,hicedum,delq
      real*8 evapw,Qhw,Qew,evapi,Qhi,Qei,tkw,tki,fracprv,evap
      real*8 eflux,eadd,qbot,qw,snowmlt,qnetice,evaps
      real*8 u2,fracadd,tw1,tw2,r1,r2,rtu
