
###############################################################################
#                                                                             #
# Creates the LAPACK library from separate Fortran source files obtained from #
# http://people.scs.fsu.edu/~burkardt/f_src/lapack/lapack.html                #
# and split by f90split.                                                      #
#                                                                             #
# A. Kuronen, 2007, 2009, 2013                                                #
#                                                                             #
###############################################################################


#########################
# GNU Fortran
F90=gfortran
FFLAGS=-O3 #-ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops -funit-at-a-time
AR=ar
#########################
# Intel Fortran
#F90=ifort
#FFLAGS=-O3
#AR=xiar
#########################
# Open64 Fortran 
#F90=openf95
#FFLAGS=-O
#AR=ar
#########################

LIB=liblapack.a
OBJ=$(patsubst %.f90,%.o,$(wildcard *.f90))

lib: $(OBJ)
	rm -f $(LIB)
	$(AR) qc $(LIB) $(OBJ)


# gfortran + optimization + DLAMC1 -> problems
DLAMC1.o: DLAMC1.f90
	$(F90) -O0 -c $<

%.o: %.f90
	$(F90) $(FFLAGS) -c $<


clean:
	rm -f *.o

cleanall:
	rm -f *.o $(LIB)
