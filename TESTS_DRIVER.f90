MODULE MODULE_UNIT_TESTS

USE FRUIT

USE MODULE_ATMOS

USE MODULE_CANOPY

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

  subroutine calcfw_test1()
    implicit none
    real*8 :: calcfw_result,calcfw_true
    real*8 :: Swq,wc,zero,fw,wsc
    Swq = 0.0d0
    wc = 3.1d0
    zero = 0.0d0
    fw = 2.2d0
    wsc = 1.7d0
    call calcfw(Swq,wc,zero,fw,wsc)
    calcfw_result = fw
    calcfw_true = 1.000000000
    call set_unit_name ('calcfw_test1')
    call assert_equals (calcfw_result, calcfw_true)
  end subroutine
 
  subroutine calcdc_test1()
    implicit none
    real*8 :: calcdc_result, calcdc_true
    real*8 :: dc,one,epetw,zero
    dc = 2.38
    one = 1.0d0
    epetw = 2.0d0
    zero = 0.0d0
    call calcdc(dc,one,epetw,zero)
    calcdc_result = dc
    calcdc_true = 1.00000000 
    call set_unit_name ('calcdc_test1')
    call assert_equals (calcdc_result, calcdc_true)
  end subroutine

  subroutine calcepw_test1()
    implicit none
    real*8 :: calcepw_result, calcepw_true
    real*8 :: epwms,epetw,one,dc,fw,dt,wc
    epwms = 2.5d0
    epetw = 3.2d0
    one = 1.0d0
    dc = 1.5d0
    fw = 0.5d0
    dt = 1.0d0
    wc = 0.2d0
    call calcepw(epwms,epetw,one,dc,fw,dt,wc)
    calcepw_result = fw
    calcepw_true = 0.12500000000000000
    call set_unit_name ('calcepw_test1')
    call assert_equals (calcepw_result, calcepw_true)
  end subroutine
  
  subroutine calcepw_test2()
    implicit none
    real*8 :: calcepw_result, calcepw_true
    real*8 :: epwms,epetw,one,dc,fw,dt,wc
    epwms = 2.5d0
    epetw = 3.2d0
    one = 1.0d0
    dc = 1.5d0
    fw = 0.5d0
    dt = 1.0d0
    wc = 0.2d0
    call calcepw(epwms,epetw,one,dc,fw,dt,wc)
    calcepw_result = epwms
    calcepw_true = -1.000000000000
    call set_unit_name ('calcepw_test2')
    call assert_equals (calcepw_result, calcepw_true)
  end subroutine

  subroutine interstorm_test1()
    implicit none
    real*8 :: interstorm_result, interstorm_true
    integer :: ipix,intstp,istmst,istorm,intstm
    real*8 :: precipi,outf,snowp,xintst,dt,endstm,r_input
    ipix = 1
    intstp = 1
    istmst = 1
    istorm = 1
    intstm = 1
    precipi = 2.0d0
    outf = 0.5d0
    snowp = 0.1d0
    xintst = 1.0d0
    dt = 0.5d0
    endstm = 2.0d0
    r_input = 1.0d0
    call interstorm(ipix,precipi,outf,snowp,xintst,&
        dt,intstp,endstm,istmst,istorm,intstm)
    interstorm_result = istmst
    interstorm_true = 2.00000000000
    call set_unit_name ('interstorm_test1')
    call assert_equals (interstorm_result, interstorm_true)
  end subroutine

  subroutine interstorm_test2()
    implicit none
    real*8 :: interstorm_result, interstorm_true
    integer :: ipix,intstp,istmst,istorm,intstm
    real*8 :: precipi,outf,snowp,xintst,dt,endstm,r_input
    ipix = 1
    intstp = 1
    istmst = 1
    istorm = 1
    intstm = 1
    precipi = 2.0d0
    outf = 0.5d0
    snowp = 0.1d0
    xintst = 1.0d0
    dt = 0.5d0
    endstm = 2.0d0
    r_input = 1.0d0
    call interstorm(ipix,precipi,outf,snowp,xintst,&
        dt,intstp,endstm,istmst,istorm,intstm)
    interstorm_result = intstp
    interstorm_true = 1.00000000000
    call set_unit_name ('interstorm_test2')
    call assert_equals (interstorm_result, interstorm_true)
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
    call calcfw_test1()
    call calcdc_test1()
    call calcepw_test1()
    call calcepw_test2()
    call interstorm_test1()
    call interstorm_test2() 

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
