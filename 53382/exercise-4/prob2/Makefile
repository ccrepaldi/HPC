# -*-makefile-*-

F90=gfortran

.PHONY: clean

ex4p2: ex4p2.f90
	$(F90) -O0 -o ex4p2 $<
	./ex4p2 > ex4p2-O0.output
	$(F90) -O3 -o ex4p2 $<
	./ex4p2 > ex4p2-O3.output
	rm -f ex4p2
	cat *.output | less

clean:
	rm -f *.o *.output
