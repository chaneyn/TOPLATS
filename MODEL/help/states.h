      integer inc_frozen,i_moss,i,iopsmini

      real*8 tkmid_moss,r_mossm_u,r_mossm_f,r_mossm,zw,zbar,ff
      real*8 atanb,xlamda,psic,zrz,ztz,rzsm1,tzsm1,thetas
      real*8 zrzmax,thetar,bcbeta,rzsm1_u,tzsm1_u,rzsm1_f,tzsm1_f
      real*8 tsoilold,bulk_dens,a_ice,b_ice,row,rzsmold,tzsmold,r_mossmold
      real*8 rzsm,tzsm,r_mossm1,zmoss,r_moss_depth,thetas_moss,rzsm_u,rzsm_f
      real*8 tzsm_u,tzsm_f,r_mossm1_u,r_mossm1_f
      real*8 a_ice_moss,b_ice_moss,bulk_dens_moss
      real*8 zero,one,two,three,four,five,six
      real*8 zw0
      data zero,one,two,three,four,five,six/0.d0,1.d0,2.d0,&
             3.d0,4.d0,5.d0,6.d0/
