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

  subroutine run_unit_tests()
  
    !Driver to run all the unit tests of the key subroutine in TOPLATS
   
    !Initialize the unit tests library 
    call init_fruit

    !Run each test
    call clcf1par_test1()
    call clcf1par_test2()

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
