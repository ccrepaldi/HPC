# -*-makefile-*-

F90=gfortran

.PHONY: clean

ex4p4: ex4p4.f90
	$(F90) -O2 -ftree-vectorize -ftree-vectorizer-verbose=1 -o ex4p4 $<
	./ex4p4 > ex4p4-vectorized.output
	$(F90) -O2 -fno-tree-vectorize -o ex4p4 $<
	./ex4p4 > ex4p4-no-tree-vectorized.output
	rm -f ex4p4
	cat *.output | less

clean:
	rm -f *.o *.output
