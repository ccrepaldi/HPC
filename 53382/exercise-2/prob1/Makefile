# -*-makefile-*-

F90=gfortran

FFLAGS=-O0 -g

OBJ=array.f90

ex2p1: $(OBJ)
	$(F90) $(FFLAGS) -o ex2p1 $(OBJ) 

.PHONY: clean

clean:
	rm -f ex2p1
	rm -rf *.dSYM
