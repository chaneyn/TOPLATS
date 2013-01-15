      integer i,newstorm,inc_frozen,ikopt,ivgtyp,num
      integer iter,max_iter,non_linear_flag
      integer i_und,i_moss,num_exp_iter

      real*8 rzsm,tzsm,tzsm_u,rzsm_u,rzsm1,tzsm1
      real*8 zrz,ztz,zw0,zrzmax,dt
      real*8 evtact,evtact_us,bsdew,dewrun,grz,gtz,diftz,difrz
      real*8 satxr,runtot,xinact,cuminf
      real*8 ff,thetar,thetas,bcbeta,xk0,psic
      real*8 Swq,Swq_us,dc,fw,dc_us,fw_us,evrz_moss,f_und
      real*8 rzsm0,tzsm0,rzsm1old,tzsm1new,tzsm1old,rzsm1new
      real*8 evtran_rz,evtran_tz,cor_flx_tz,cor_flx_rz,ddifrzdth1
      real*8 ddifrzdth2,dgrzdth1,dgrzdth2,dgtzdth1,dgtzdth2,ddiftzdth1
      real*8 ddiftzdth2,explicit_weight,implicit_weight,tol
      real*8 srzflx,dewrz,stzflx,F11,F22,dF1dtheta1,dF1dtheta2
      real*8 dF2dtheta1,dF2dtheta2,del_rzsm,del_tzsm
      real*8 grz_sum,difrz_sum,gtz_sum,diftz_sum
      real*8 dummy,d2rzsmdt2,d2tzsmdt2
      real*8 dtsrzflx,dtstzflx,case_flag,xksrz,xkstz
      real*8 dstz,dsrz,rzrhs,tzrhs
      real*8 zero,one,two,three,four,five,six
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
