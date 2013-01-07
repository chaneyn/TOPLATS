MODULE MODULE_UNIT_TESTS

USE FRUIT

USE MODULE_ATMOS
USE MODULE_LAND

implicit none

contains

  subroutine clcf1par_test1()
    implicit none
    real*8 :: f1par_result,f1par_true
    call set_unit_name ('clcf1par_test1')
    f1par_result = clcf1par(0.3d0,5.0d0,500.0d0,100.0d0,5000.0d0,100.0d0)
    f1par_true = 2.2405063029073640d0
    call assert_equals (f1par_true,f1par_result)
  end subroutine

  subroutine clcf1par_test2()
    implicit none
    real*8 :: f1par_result,f1par_true
    call set_unit_name ('clcf1par_test2')
    f1par_result = clcf1par(0.1d0,1.0d0,100.0d0,200.0d0,3000.0d0,10.0d0)
    f1par_true = 1.0936454829336983d0
    call assert_equals (f1par_true,f1par_result)
  end subroutine

  subroutine clcf3vpd_test1()
    implicit none
    real*8 :: f3vpd_result,f3vpd_true
    call set_unit_name ('clcf3vpd_test1')
    f3vpd_result = clcf3vpd(275.0d0,1.0d0,2.0d0)
    f3vpd_true = 4.d0
    call assert_equals (f3vpd_true,f3vpd_result)
  end subroutine

  subroutine clcf3vpd_test2()
    implicit none
    real*8 :: f3vpd_result,f3vpd_true
    call set_unit_name ('clcf3vpd_test2')
    f3vpd_result = clcf3vpd(1.0d0,1000.0d0,250.0d0)
    f3vpd_true = 4.d0
    call assert_equals (f3vpd_true,f3vpd_result)
  end subroutine

  subroutine clcf4temp_test1()
    implicit none
    real*8 :: f4temp_result,f4temp_true
    call set_unit_name ('clcf4temp_test1')
    f4temp_result = clcf4temp(2.0d0,27.0d0,5.0d0)
    f4temp_true = -0.0097428548611568806d0
    call assert_equals (f4temp_true,f4temp_result)
  end subroutine

  subroutine clcf4temp_test2()
    implicit none
    real*8 :: f4temp_result,f4temp_true
    call set_unit_name ('clcf4temp_test2')
    f4temp_result = clcf4temp(200.0d0,1.0d0,20.0d0)
    f4temp_true = -0.0021964570949596605d0
    call assert_equals (f4temp_true,f4temp_result)
  end subroutine

  subroutine sm_cen_dif_test1()
    implicit none
    integer :: iffroz
    real*8 :: sm_cen_dif_result,sm_cen_dif_true,smtmp,smold
    call set_unit_name ('sm_cen_dif_test1')
    smold = 0.0d0
    smtmp = 0.0d0
    iffroz = 0
    call sm_cen_dif(iffroz,280.0d0,5.0d0,&
      10.0d0,smtmp,15.0d0,5.0d0,smold,25.0d0,35.0d0)
    sm_cen_dif_result = smtmp
    sm_cen_dif_true = 15.0d0
    call assert_equals (sm_cen_dif_true,sm_cen_dif_result)
  end subroutine

  subroutine sm_cen_dif_test2()
    implicit none
    integer :: iffroz
    real*8 :: sm_cen_dif_result,sm_cen_dif_true,smtmp,smold
    call set_unit_name ('sm_cen_dif_test2')
    call set_unit_name ('sm_cen_dif_test1')
    smold = 0.0d0
    smtmp = 0.0d0
    iffroz = 0
    call sm_cen_dif(iffroz,270.0d0,7.0d0,&
      2.0d0,smtmp,23.0d0,4.0d0,smold,100.0d0,7.0d0)
    sm_cen_dif_result = smold
    sm_cen_dif_true = 53.5d0
    call assert_equals (sm_cen_dif_true,sm_cen_dif_result)
  end subroutine

  subroutine soiladapt_test1()
    implicit none
    integer :: iopgveg,iopthermc_v
    real*8 :: tcbeta,xlai,thermc1,heatcap,heatcap1,zero,tau
    real*8 :: soiladapt_result,soiladapt_true,thermc
    iopgveg = 10
    thermc = 18.0d0
    iopthermc_v = 5
    tcbeta = 7.0d0
    xlai = 2.0d0
    thermc1 = 17.0d0
    heatcap = 4.0d0
    heatcap1 = 3.0d0
    zero = 0.0d0
    call set_unit_name ('soiladapt_test1')
    call soiladapt(iopgveg,thermc,iopthermc_v,tcbeta,xlai,&
         thermc1,heatcap,heatcap1,zero)
    soiladapt_result = thermc
    soiladapt_true = 0.000014135988224760654d0
    call assert_equals (soiladapt_true,soiladapt_result)
  end subroutine

  subroutine soiladapt_test2()
    implicit none
    integer :: iopgveg,iopthermc_v
    real*8 :: tcbeta,xlai,thermc1,heatcap,heatcap1,zero,tau
    real*8 :: soiladapt_result,soiladapt_true,thermc
    iopgveg = 5
    thermc = 18.0d0
    iopthermc_v = 5
    tcbeta = 7.0d0
    xlai = 1.0d0
    thermc1 = 17.0d0
    heatcap = 4.0d0
    heatcap1 = 3.0d0
    zero = 0.0d0
    call set_unit_name ('soiladapt_test2')
    call soiladapt(iopgveg,thermc,iopthermc_v,tcbeta,xlai,&
         thermc1,heatcap,heatcap1,zero)
    soiladapt_result = thermc
    soiladapt_true = 0.015501993414426776d0
    call assert_equals (soiladapt_true,soiladapt_result)
  end subroutine

  subroutine esat_test1()
    implicit none
    real*8 :: esat_result,esat_true
    call set_unit_name ('esat_test1')
    esat_result = esat(300.0d0)
    esat_true = 3535.2433230660672d0
    call assert_equals (esat_true, esat_result)
  end subroutine

  subroutine esat_test2()
    implicit none
    real*8 :: esat_result,esat_true
    call set_unit_name ('esat_test2')
    esat_result = esat(260.0d0)
    esat_true = 221.83511852146523d0
    call assert_equals (esat_true, esat_result)
  end subroutine

  subroutine nreb_snow_test1()
    implicit none
    real*8 :: nreb_snow_result,nreb_snow_true,dzdeep
    real*8 :: thermc1,thermc2,heatcap1,heatcap2,heatcapold
    real*8 :: tskink,tmidknew,tairc,zdeep,tdeep,zmid,dt,gtmp
    thermc1 = 5.0d0
    thermc2 = 18.0d0
    heatcap1 = 2.0d0
    heatcap2 = 0.0d0
    heatcapold = 14.0d0
    tskink = 200.0d0
    tmidknew = 50.0d0
    tairc = 3.0d0
    zdeep = 47.0d0
    tdeep = 25.0d0
    zmid = 13.0d0
    dt = 2.0d0 
    gtmp = 6.0d0
    call set_unit_name ('nreb_snow_test1')
    call nreb_snow(thermc1,thermc2,heatcap1,heatcap2,heatcapold,&
                   tskink,tmidknew,tairc,zdeep,tdeep,zmid,dt,gtmp)
    nreb_snow_result = gtmp
    nreb_snow_true = -4.9009900990099009d0
    call assert_equals (nreb_snow_true,nreb_snow_result)
  end subroutine

  subroutine nreb_snow_test2()
    implicit none
    real*8 :: nreb_snow_result,nreb_snow_true,dzdeep
    real*8 :: thermc1,thermc2,heatcap1,heatcap2,heatcapold
    real*8 :: tskink,tmidknew,tairc,zdeep,tdeep,zmid,dt,gtmp
    thermc1 = 5.0d0
    thermc2 = 18.0d0
    heatcap1 = 2.0d0
    heatcap2 = 0.0d0
    heatcapold = 14.0d0
    tskink = 200.0d0
    tmidknew = 50.0d0
    tairc = 3.0d0
    zdeep = 47.0d0
    tdeep = 25.0d0
    zmid = 13.0d0
    dt = 2.0d0
    gtmp = 6.0d0
    call set_unit_name ('nreb_snow_test2')
    call nreb_snow(thermc1,thermc2,heatcap1,heatcap2,heatcapold,&
                   tskink,tmidknew,tairc,zdeep,tdeep,zmid,dt,gtmp)
    nreb_snow_result = tmidknew
    nreb_snow_true = 20.371287128712872d0
    call assert_equals (nreb_snow_true,nreb_snow_result)
  end subroutine

  subroutine stabcor_test1()
    implicit none
    real*8 :: zwind,zhum,wind_s,zpdis,r_m,tair
    real*8 :: airp,t_skin,vappres,richn
    real*8 :: stabcor_result,stabcor_true,uza
    zwind = 10.0d0
    zhum = 4.0d0
    wind_s = 45.0d0
    zpdis = 2.5d0
    r_m = 5.0d0
    tair = 150.0d0
    airp = 2.0d0
    t_skin = 11.0d0
    vappres = 2.5d0
    richn = 1.5d0
    call set_unit_name ('stabcor_test1')
    call stabcor(zwind,zhum,wind_s,zpdis,r_m,tair,airp,&
                 t_skin,vappres,richn)
    stabcor_result = richn
    stabcor_true = 3636.2399999999993d0
    call assert_equals (stabcor_true,stabcor_result)
  end subroutine

  subroutine stabcor_test2()
    implicit none
    real*8 :: zwind,zhum,wind_s,zpdis,r_m,tair
    real*8 :: airp,t_skin,vappres,richn
    real*8 :: stabcor_result,stabcor_true,uza
    zwind = 12.0d0
    zhum = 5.0d0
    wind_s = 35.0d0
    zpdis = 4.0d0
    r_m = 0.3d0
    tair = 6.0d0
    airp = 18.0d0
    t_skin = 37.0d0
    vappres = 4.0d0
    richn = 2.2d0
    call set_unit_name ('stabcor_test2')
    call stabcor(zwind,zhum,wind_s,zpdis,r_m,tair,airp,&
                 t_skin,vappres,richn)
    stabcor_result = richn
    stabcor_true = -1.5386200769623049d0
    call assert_equals (stabcor_true,stabcor_result)
  end subroutine

  subroutine calcra_test1()
    implicit none
    real*8 :: calcra_result,calcra_true
    calcra_result = calcra(0.5d0,5.0d0,8.0d0,4.2d0,10.0d0,&
      15.0d0,2.0d0)
    calcra_true = 173.39747643603178d0
    call set_unit_name ('calcra_test1')
    call assert_equals (calcra_true,calcra_result)
  end subroutine

  subroutine calcra_test2()
    implicit none
    real*8 :: calcra_result,calcra_true
    calcra_result = calcra(3.0d0,7.0d0,11.0d0,2.2d0,9.0d0,&
      13.0d0,1.0d0)
    calcra_true = 1.3890511481178962d0
    call set_unit_name ('calcra_test2')
    call assert_equals (calcra_true,calcra_result)
  end subroutine

  subroutine calcrain_test1()
    implicit none
    real*8 :: tcel,snow,rain,precip_o,dt,snow_true,rain_true    
    tcel = 10.0d0
    snow = 0.0d0
    rain = 0.0d0
    precip_o = 0.1d0
    dt = 3600.0d0
    snow_true = 0.0d0
    rain_true = 360.0d0

    call calcrain(tcel,snow,rain,precip_o,dt)
        
    call set_unit_name ('calcrain_test1')
    call assert_equals (snow_true,snow)
    call assert_equals (rain_true,rain)
  end subroutine 

  subroutine calcrain_test2()
    implicit none
    real*8 :: tcel,snow,rain,precip_o,dt,snow_true,rain_true    
    tcel = -0.1d0
    snow = 1.0d0
    rain = 12.0d0
    precip_o = 0.01d0
    dt = 3600.0d0
    snow_true = 36.0d0
    rain_true = 0.0d0

    call calcrain(tcel,snow,rain,precip_o,dt)
        
    call set_unit_name ('calcrain_test2')
    call assert_equals (snow_true,snow)
    call assert_equals (rain_true,rain)
  end subroutine 

  subroutine calcrsoil_test1()

    integer :: irestype
    real*8 :: rsoil,srespar1,srespar2,srespar3
    real*8 :: thetas,rzsm,tkact,rsoil_true

    irestype = 2
    rsoil = 0.0d0
    srespar1 = 0.1d0
    srespar2 = 0.1d0
    srespar3 = 0.1d0
    thetas = 0.1d0
    rzsm = 0.1d0 
    tkact = 0.1d0
    rsoil_true = 0.2d0
    
    call calcrsoil(irestype,rsoil,srespar1,srespar2,srespar3,&
                   thetas,rzsm,tkact)
    call set_unit_name ('calcrsoil_test1')
    call assert_equals (rsoil_true,rsoil)

  end subroutine

  subroutine calcrsoil_test2()
  
    integer :: irestype
    real*8 :: rsoil,srespar1,srespar2,srespar3
    real*8 :: thetas,rzsm,tkact,rsoil_true


    irestype = 3
    rsoil = 0.0d0
    srespar1 = 0.1d0
    srespar2 = 1.0d0
    srespar3 = 0.1d0
    thetas = 1.0d0
    rzsm = 0.1d0 
    tkact = 10.0d0
    rsoil_true = 5470.7735443438487d0
    
    call calcrsoil(irestype,rsoil,srespar1,srespar2,srespar3,&
                   thetas,rzsm,tkact)
    call set_unit_name ('calcrsoil_test2')
    call assert_equals (rsoil_true,rsoil)
  end subroutine

  subroutine calcrsoil_test3()

    integer :: irestype
    real*8 :: rsoil,srespar1,srespar2,srespar3
    real*8 :: thetas,rzsm,tkact,rsoil_true

    irestype = 4
    rsoil = 0.0d0
    srespar1 = 0.1d0
    srespar2 = 0.1d0
    srespar3 = 0.1d0
    thetas = 2.0d0
    rzsm = 0.1d0 
    tkact = 10.0d0
    rsoil_true = 0.09d0
    
    call calcrsoil(irestype,rsoil,srespar1,srespar2,srespar3,&
                   thetas,rzsm,tkact)
    call set_unit_name ('calcrsoil_test3')
    call assert_equals (rsoil_true,rsoil)
  end subroutine

  subroutine calcrsoil_test4()
    integer :: irestype 
    real*8 :: rsoil,srespar1,srespar2,srespar3
    real*8 :: thetas,rzsm,tkact,rsoil_true

    irestype = 5
    rsoil = 0.0d0
    srespar1 = 0.1d0
    srespar2 = 0.1d0
    srespar3 = 0.1d0
    thetas = 0.1d0
    rzsm = 0.1d0 
    tkact = 0.1d0
    rsoil_true = 0.11051709180756478d0
    
    call calcrsoil(irestype,rsoil,srespar1,srespar2,srespar3,&
                   thetas,rzsm,tkact)
    call set_unit_name ('calcrsoil_test4')
    call assert_equals (rsoil_true,rsoil)
  end subroutine
  
  subroutine reset_inf_pars_test1()
  
    real*8 :: cuminf,zero,rzsmt,rzsm,thetas,tolinf
    real*8 :: sorp,two,xk0,psic,thetar,bcgamm,bcbeta,deltrz,cc,one
    real*8 :: cuminf_true,rzsmt_true,sorp_true,deltrz_true,cc_true

    cuminf = 1.0d0
    zero = 0.0d0
    one = 1.0d0
    two = 2.0d0
    rzsmt = 0.0d0
    rzsm  = 1.0d0
    thetas = 2.0d0
    tolinf = 0.0d0
    sorp = 0.0d0
    xk0 = 0.01d0
    psic = 2.0d0
    thetar = 0.5d0
    bcgamm = 17.0d0
    bcbeta = 5.0d0
    deltrz = 0.0d0
    cc = 0.0d0

    cuminf_true = 0.0d0
    rzsmt_true = 1.0d0
    sorp_true = 0.20357171080835726d0
    deltrz_true = 0.5d0
    cc_true = 0.51193322249957873d0

    call reset_inf_pars(cuminf,zero,rzsmt,rzsm,thetas,tolinf,sorp,two,&
                        xk0,psic,thetar,bcgamm,bcbeta,deltrz,cc,one)
    call set_unit_name ('reset_inf_pars_test1')
    call assertequals(cuminf_true,cuminf)
    call assertequals(rzsmt_true,rzsmt)
    call assertequals(sorp_true,sorp)
    call assertequals(deltrz_true,deltrz)
    call assertequals(cc_true,cc)
    
  end subroutine

  subroutine reset_inf_pars_test2()
  
    real*8 :: cuminf,zero,rzsmt,rzsm,thetas,tolinf
    real*8 :: sorp,two,xk0,psic,thetar,bcgamm,bcbeta,deltrz,cc,one
    real*8 :: cuminf_true,rzsmt_true,sorp_true,deltrz_true,cc_true

    cuminf = 1.0d0
    zero = 0.0d0
    one = 1.0d0
    two = 2.0d0
    rzsmt = 0.0d0
    rzsm  = 2.0d0
    thetas = 2.0d0
    tolinf = 0.1d0
    sorp = 0.0d0
    xk0 = 0.01d0
    psic = 2.0d0
    thetar = 0.5d0
    bcgamm = 17.0d0
    bcbeta = 5.0d0
    deltrz = 0.0d0
    cc = 0.0d0

    cuminf_true = 0.0d0
    rzsmt_true = 1.9d0
    sorp_true = 0.063359406676628677d0
    deltrz_true = 1.4d0
    cc_true = 0.89545315507131318d0

    call reset_inf_pars(cuminf,zero,rzsmt,rzsm,thetas,tolinf,sorp,two,&
                        xk0,psic,thetar,bcgamm,bcbeta,deltrz,cc,one)
    call set_unit_name ('reset_inf_pars_test2')
    call assertequals(cuminf_true,cuminf)
    call assertequals(rzsmt_true,rzsmt)
    call assertequals(sorp_true,sorp)
    call assertequals(deltrz_true,deltrz)
    call assertequals(cc_true,cc)
    
  end subroutine

  subroutine calcsmcond_test1()
    
    real*8 :: rzsm,tc,smcond,one,tw,zero
    real*8 :: smcond_true

    one = 1.0d0
    zero = 0.0d0
    rzsm = 0.5d0
    tc = 0.75d0
    tw = 0.1d0
    smcond = 0.0d0

    smcond_true = 0.61538461538461542d0

    call calcsmcond(rzsm,tc,smcond,one,tw,zero)
    call set_unit_name ('calcsmcond_test1')
    call assertequals(smcond_true,smcond)
  end subroutine

  subroutine calcsmcond_test2()
    
    real*8 :: rzsm,tc,smcond,one,tw,zero
    real*8 :: smcond_true

    one = 1.0d0
    zero = 0.0d0
    rzsm = 0.5d0
    tc = 0.4d0
    tw = 0.1d0
    smcond = 0.0d0

    smcond_true = 1.0d0

    call calcsmcond(rzsm,tc,smcond,one,tw,zero)
    call set_unit_name ('calcsmcond_test2')
    call assertequals(smcond_true,smcond)
  end subroutine
 
  subroutine new_dwnflx_test1()

    integer :: ikopt,inc_frozen
    real*8 :: zrz,rzsm,ff,thetar,thetas,rzsm_u,grz,bcbeta
    real*8 :: ztz,gtz,tzsm,tzsm_u,xksrz,relsrz,xkstz,relstz
    real*8 :: grz_true,gtz_true,relsrz_true,relstz_true
   
    ikopt = 1
    inc_frozen = 1
    zrz = 1.0d0
    rzsm = 0.5d0
    ff = 1.0d0
    thetar = 0.3d0
    thetas = 0.7d0
    rzsm_u = 0.0d0
    grz = 0.15d0
    bcbeta = 10.0d0
    ztz = 0.0d0
    gtz = 0.0d0
    tzsm = 0.2d0
    tzsm_u = 0.03d0
    xksrz = 0.25d0
    relsrz = 0.5d0
    xkstz = 0.1d0
    relstz = 0.05d0
    
    grz_true = 0.0d0
    relsrz_true = 0.5d0
    gtz_true = 0.0d0
    relstz_true = 0.05d0

    call new_dwnflx(zrz,ikopt,xksrz,xkstz,rzsm,ff,inc_frozen,&
                    thetar,thetas,rzsm_u,grz,bcbeta,ztz,gtz,&
                    tzsm,tzsm_u)
    call set_unit_name ('new_dwnflx_test1')
    call assertequals(grz_true,grz)
    call assertequals(relsrz_true,relsrz)
    call assertequals(gtz_true,gtz)
    call assertequals(relstz_true,relstz)
  end subroutine
        
  subroutine new_dwnflx_test2()

    integer :: ikopt,inc_frozen
    real*8 :: zrz,rzsm,ff,thetar,thetas,rzsm_u,grz,bcbeta
    real*8 :: ztz,gtz,tzsm,tzsm_u,xksrz,relsrz,xkstz,relstz
    real*8 :: grz_true,gtz_true,relsrz_true,relstz_true
   
    ikopt = 1
    inc_frozen = 0
    zrz = 1.0d0
    rzsm = 0.5d0
    ff = 1.0d0
    thetar = 0.1d0
    thetas = 0.7d0
    rzsm_u = 0.0d0
    grz = 0.15d0
    bcbeta = 7.0d0
    ztz = 0.0d0
    gtz = 0.0d0
    tzsm = 0.37d0
    tzsm_u = 0.03d0
    xksrz = 0.26d0
    relsrz = 0.5d0
    xkstz = 0.7d0
    relstz = 0.05d0
    
    grz_true = 0.068610057374659361d0
    relsrz_true = 0.5d0
    gtz_true = 0.0d0
    relstz_true = 0.05d0

    call new_dwnflx(zrz,ikopt,xksrz,xkstz,rzsm,ff,inc_frozen,&
                    thetar,thetas,rzsm_u,grz,bcbeta,ztz,gtz,&
                    tzsm,tzsm_u)
    call set_unit_name ('new_dwnflx_test2')
    call assertequals(grz_true,grz)
    call assertequals(relsrz_true,relsrz)
    call assertequals(gtz_true,gtz)
    call assertequals(relstz_true,relstz)
  end subroutine

  subroutine new_difflx_test1()
    
    integer :: ikopt,inc_frozen
    real*8 :: xksrz,xkstz,ff,zrz,ztz,rzsm,tzsm,thetas,thetar
    real*8 :: bcbeta,psic,difrz,rzsm_u,tzsm_u,diftz
    real*8 :: difrz_true,diftz_true
    ikopt = 1
    inc_frozen = 1
    xksrz = 0.26d0
    xkstz = 0.7d0
    ff = 1.0d0
    zrz = 0.5d0
    ztz = 1.0d0
    rzsm = 0.5d0
    tzsm = 0.3d0
    thetas = 0.7d0
    thetar = 0.1d0
    bcbeta = 7.0d0
    psic = 2.0d0 
    difrz = 1.0d0
    rzsm_u = 0.5d0
    tzsm_u = 0.4d0
    diftz = 0.5d0

    difrz_true = 0.0085829351856200782d0
    diftz_true = -0.035652192309498787d0

    call new_difflx(ikopt,xksrz,xkstz,ff,zrz,ztz,inc_frozen,&
                    rzsm,tzsm,thetas,thetar,bcbeta,psic,difrz,&
                    rzsm_u,tzsm_u,diftz)
    call set_unit_name('new_difflx_test1')
    call assertequals(difrz_true,difrz)
    call assertequals(diftz_true,diftz)
  end subroutine

  subroutine new_difflx_test2()
    
    integer :: ikopt,inc_frozen
    real*8 :: xksrz,xkstz,ff,zrz,ztz,rzsm,tzsm,thetas,thetar
    real*8 :: bcbeta,psic,difrz,rzsm_u,tzsm_u,diftz
    real*8 :: difrz_true,diftz_true
    ikopt = 1
    inc_frozen = 0
    xksrz = 0.26d0
    xkstz = 0.8d0
    ff = 1.0d0
    zrz = 0.3d0
    ztz = 0.9d0
    rzsm = 0.5d0
    tzsm = 0.3d0
    thetas = 0.7d0
    thetar = 0.1d0
    bcbeta = 7.0d0
    psic = 2.0d0 
    difrz = 1.0d0
    rzsm_u = 0.5d0
    tzsm_u = 0.4d0
    diftz = 0.5d0

    difrz_true = 0.011936146649400578d0
    diftz_true = -0.0081018518518518705d0

    call new_difflx(ikopt,xksrz,xkstz,ff,zrz,ztz,inc_frozen,&
                    rzsm,tzsm,thetas,thetar,bcbeta,psic,difrz,&
                    rzsm_u,tzsm_u,diftz)
    call set_unit_name('new_difflx_test2')
    call assertequals(difrz_true,difrz)
    call assertequals(diftz_true,diftz)
  end subroutine

  subroutine clc_evrz_test1()
    
    integer :: ivgtyp,i_und,i_moss
    real*8 :: evrz,Swq,Swq_us,evtact,dc,fw,evtact_us,dc_us
    real*8 :: fw_us,evrz_moss,dummy,f_und
    real*8 :: evrz_true

    ivgtyp = 0
    i_und = 0
    i_moss = 0
    evrz = 0.0d0
    Swq = -0.1d0
    Swq_us = 0.0d0
    evtact = 0.5d0
    evtact_us = 0.0d0
    dc = 0.5d0
    dc_us = 0.0d0
    fw = 0.5d0
    fw_us = 0.0d0
    evrz_moss = 0.0d0
    dummy = 0.0d0
    f_und = 0.0d0
    
    evrz_true = 0.25d0
        
    call clc_evrz(evrz,Swq,Swq_us,ivgtyp,evtact,dc,i_und,&
                  i_moss,fw,evtact_us,dc_us,fw_us,evrz_moss,&
                  dummy,f_und)
    call set_unit_name('clc_evrz_test1')
    call assertequals(evrz_true,evrz)     
  end subroutine   
 
  subroutine clc_evrz_test2()
    
    integer :: ivgtyp,i_und,i_moss
    real*8 :: evrz,Swq,Swq_us,evtact,dc,fw,evtact_us,dc_us
    real*8 :: fw_us,evrz_moss,dummy,f_und
    real*8 :: evrz_true

    ivgtyp = 1
    i_und = 0
    i_moss = 0
    evrz = 0.0d0
    Swq = -0.1d0
    Swq_us = 0.0d0
    evtact = 0.5d0
    evtact_us = 0.1d0
    dc = 0.5d0
    dc_us = 0.1d0
    fw = 0.5d0
    fw_us = 0.1d0
    evrz_moss = 0.0d0
    dummy = 0.0d0
    f_und = 0.0d0
    
    evrz_true = 0.125d0
        
    call clc_evrz(evrz,Swq,Swq_us,ivgtyp,evtact,dc,i_und,&
                  i_moss,fw,evtact_us,dc_us,fw_us,evrz_moss,&
                  dummy,f_und)
    call set_unit_name('clc_evrz_test2')
    call assertequals(evrz_true,evrz)     
  end subroutine   
      
  subroutine run_unit_tests()
  
    !Driver to run all the unit tests of the key subroutine in TOPLATS
   
    !Initialize the unit tests library 
    call init_fruit

    !Run each test
    call clcf1par_test1()
    call clcf1par_test2()
    call clcf3vpd_test1()
    call clcf3vpd_test2()
    call clcf4temp_test1()
    call clcf4temp_test2()
    call sm_cen_dif_test1()
    call sm_cen_dif_test2()
    call soiladapt_test1()
    call soiladapt_test2()
    call esat_test1()
    call esat_test2()
    call nreb_snow_test1()
    call nreb_snow_test2()
    call stabcor_test1()
    call stabcor_test2()
    call calcra_test1()
    call calcra_test2()
    call calcrain_test1()
    call calcrain_test2()
    call calcrsoil_test1()
    call calcrsoil_test2()
    call calcrsoil_test3()
    call calcrsoil_test4()
    call reset_inf_pars_test1()    
    call reset_inf_pars_test2()
    call calcsmcond_test1()
    call calcsmcond_test2()    
    call new_dwnflx_test1()
    call new_dwnflx_test2()
    call new_difflx_test1()
    call new_difflx_test2()
    call clc_evrz_test1()
    call clc_evrz_test2()   
    !Summarize and finalize the unit tests
    call fruit_summary
    call fruit_finalize

  end subroutine run_unit_tests

END MODULE MODULE_UNIT_TESTS

PROGRAM TESTS_DRIVER

USE MODULE_UNIT_TESTS

implicit none

call run_unit_tests()

END PROGRAM TESTS_DRIVER
