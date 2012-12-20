FC = gfortran
OPTS = 
FFLAGS = -g -fbounds-check -finit-local-zero
LIBS = -fopenmp
PROG_SSP = TOPLATS_3.0
OBJS_SSP = fruit/fruit_util.f90 fruit/fruit.f90 VARIABLES.f90 MODULE_TESTS.f90 MODULE_IO.f90 MODULE_SHARED.f90 MODULE_ATMOS.f90 MODULE_LAND.f90 MODULE_CELL.f90 MODULE_TOPMODEL.f90 *.f90

all: $(PROG_SSP)

$(PROG_SSP) : $(OBJS_SSP)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

clean:
	rm $(PROG_SSP) *.mod

