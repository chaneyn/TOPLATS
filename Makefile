FC = gfortran
OPTS = 
FFLAGS = -g -fbounds-check -finit-local-zero
LIBS = -fopenmp
PROG_SSP = TOPLATS_3.0
OBJS_SSP = fruit/fruit_util.f90 fruit/fruit.f90 VARIABLES.f90 MODULE_IO.f90 MODULE_POINT.f90 MODULE_TOPLATS.f90 *.f90

all: $(PROG_SSP)

$(PROG_SSP) : $(OBJS_SSP)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

clean:
	rm $(PROG_SSP) *.mod

