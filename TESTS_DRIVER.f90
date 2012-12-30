MODULE MODULE_UNIT_TESTS

USE FRUIT

IMPLICIT NONE

type unit_tests
  contains
  procedure :: one_plus_one
end type

contains

  subroutine one_plus_one()
   print*,'hello'
  end subroutine

  subroutine run_unit_tests()
  
    !Driver to run all the unit tests of the key subroutine in TOPLATS
    
    call init_fruit !Initialize the fruit library
    !
    call set_unit_name ('test')
    call assert_equals (1,1)
    call fruit_summary !Summarize the fruit output for this time step
    call fruit_finalize !Finalize the fruit library

  end subroutine run_unit_tests

END MODULE MODULE_UNIT_TESTS

PROGRAM TESTS_DRIVER

USE MODULE_UNIT_TESTS

IMPLICIT NONE


call run_unit_tests()

END PROGRAM TESTS_DRIVER
