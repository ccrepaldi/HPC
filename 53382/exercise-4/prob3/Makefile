# -*-makefile-*-

F90=gfortran

.PHONY: clean

ex4p3: ex4p3.f90
	$(F90) -O0 -o ex4p3 $<
	./ex4p3 > ex4p3-O0.output
	$(F90) -O3 -o ex4p3 $<
	./ex4p3 > ex4p3-O3.output
	rm -f ex4p3
	cat *.output | less

clean:
	rm -f *.o *.output
