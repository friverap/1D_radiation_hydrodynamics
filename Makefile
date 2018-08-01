#FC = mpif90
FC=gfortran
#FC=/opt/intel/bin/ifort

FFLAGS = -g -O3 -fopenmp
#FFLAGS = -O3
#-O -Wall -fcheck=all -g -fbacktrace

#LNK = mpif90
LNK=gfortran
#LNK=/opt/intel/bin/ifort

OBJS = evolve.o allocate.o grid.o initial.o rhs.o main.o check_parameters.o diagnostics.o save0Ddata.o save1Ddata.o saveLum.o recv_primitives.o flux_hlle.o sigma.o reconstructor.o lpm.o BC.o sources.o save1DRaddata.o

MODS = arrays.o global_numbers.o

$(OBJS):	$(MODS)

all:  $(OBJS) $(MODS)
		$(LNK) $(FFLAGS) -o xRHD $(OBJS) $(MODS) 
	@ mkdir -p xxx
	@ mv xRHD xxx
#	@ cp input.par xxx

.PHONY:	clean

clean:
	-rm -f *.o *.mod xxx/xRHD xxx/*xyl xxx/*xzl xxx/*yzl xxx/*tl xxx/*xl xxx/*yl xxx/*zl

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@
