# -*-makefile-*-

F90=gfortran

FFLAGS=-O0 -p -g

OBJ=simpleprofiling.f90

simpprof: $(OBJ)
	$(F90) $(FFLAGS) -o simpprof $(OBJ)

.PHONY: clean

clean:
	rm -f simpprof
	rm -rf *.dSYM
