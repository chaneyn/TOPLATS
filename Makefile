FC = gfortran
OPTS = 
FFLAGS = -g -fbounds-check -finit-local-zero
LIBS = -fopenmp
PROG_SSP = TOPLATS_3.0
OBJS_SSP = *.f90

all: $(PROG_SSP)

$(PROG_SSP) : $(OBJS_SSP)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

clean:
	rm $(PROG_SSP) *.mod

