MODULE MODULE_UNIT_TESTS

USE FRUIT

USE MODULE_ATMOS

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
    f4temp_true = -0.00974285486
    call assert_equals (f4temp_true,f4temp_result)
  end subroutine

  subroutine clcf4temp_test2()
    implicit none
    real*8 :: f4temp_result,f4temp_true
    call set_unit_name ('clcf4temp_test2')
    f4temp_result = clcf4temp(200.0d0,1.0d0,20.0d0)
    f4temp_true = -0.00219645709
    call assert_equals (f4temp_true,f4temp_result)
  end subroutine

  subroutine sm_cen_dif_test1()
    implicit none
    real*8 :: sm_cen_dif_result,sm_cen_dif_true,smtmp
    call set_unit_name ('sm_cen_dif_test1')
    call sm_cen_dif(0,280.0d0,5.0d0,&
      10.0d0,0.0d0,15.0d0,5.0d0,0.0d0,25.0d0,35.0d0)
    sm_cen_dif_result = smtmp
    sm_cen_dif_true = 15.0d0
    call assert_equals (sm_cen_dif_true,sm_cen_dif_result)
  end subroutine

  subroutine sm_cen_dif_test2()
    implicit none
    real*8 :: sm_cen_dif_result,sm_cen_dif_true,smold
    call set_unit_name ('sm_cen_dif_test2')
    call sm_cen_dif(0,270.0d0,7.0d0,&
      2.0d0,0.0d0,23.0d0,4.0d0,0.0d0,100.0d0,7.0d0)
    sm_cen_dif_result = smold
    sm_cen_dif_true = 53.5d0
    call assert_equals (sm_cen_dif_true,sm_cen_dif_result)
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
