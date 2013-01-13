#GFC = /home/ice/nchaney/UTILS/gcc-4.8/bin/gfortran
GFC = gfortran
IFC = /opt/intel/bin/ifort
FFLAGS = -g -fbounds-check -finit-local-zero
LIBS = -fopenmp
TOPLATS = TOPLATS_3.0
TESTS = TESTS_DRIVER
CATCHMENT_TESTS = CATCHMENT_TESTS_DRIVER
TOPLATS_OBJS = MODEL/MODULE_VARIABLES.f90 MODEL/MODULE_SNOW.f90 MODEL/MODULE_IO.f90 MODEL/MODULE_SHARED.f90 MODEL/MODULE_FROZEN_SOIL.f90 MODEL/MODULE_ATMOS.f90 MODEL/MODULE_LAND.f90 MODEL/MODULE_CANOPY.f90 MODEL/MODULE_CELL.f90 MODEL/MODULE_TOPMODEL.f90 MODEL/MODULE_REGIONAL.f90 TOPLATS_DRIVER.f90
TEST_OBJS = fruit/fruit_util.f90 fruit/fruit.f90 MODULE_VARIABLES.f90 MODULE_SNOW.f90 MODULE_IO.f90 MODULE_SHARED.f90 MODULE_FROZEN_SOIL.f90 MODULE_ATMOS.f90 MODULE_LAND.f90 MODULE_CANOPY.f90 MODULE_CELL.f90 MODULE_TOPMODEL.f90 MODULE_REGIONAL.f90 TESTS_DRIVER.f90
CATCHMENT_TEST_OBJS = fruit/fruit_util.f90 fruit/fruit.f90 MODULE_VARIABLES.f90 CATCHMENT_TESTS_DRIVER.f90

all: $(TOPLATS) $(TESTS) $(CATCHMENT_TESTS)

$(TOPLATS) : $(TOPLATS_OBJS)
	$(GFC) $(FFLAGS) -o $@ $^ $(LIBS)

$(TESTS) : $(TEST_OBJS)
	$(GFC) $(FFLAGS) -o $@ $^ $(LIBS)

$(CATCHMENT_TESTS) : $(CATCHMENT_TEST_OBJS)
	$(GFC) $(FFLAGS) -Wall -o $@ $^ $(LIBS)

clean:
	rm $(TOPLATS) $(TESTS) $(CATCHMENT_TESTS) *.mod

