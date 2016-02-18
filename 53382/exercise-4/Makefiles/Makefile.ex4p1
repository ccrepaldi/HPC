# -*-makefile-*-

F90=gfortran
PROBLEMS=ex4p1a ex4p1a-original ex4p1b ex4p1b-original

.PHONY: all clean print

all: $(PROBLEMS) print

ex4p1a: ex4p1a.f90
	$(F90) -O0 -o ex4p1a $<
	./ex4p1a > ex4p1a-O0.output
	$(F90) -O3 -o ex4p1a $<
	./ex4p1a > ex4p1a-O3.output
	rm -f ex4p1a

ex4p1a-original: ex4p1a-original.f90
	$(F90) -O0 -o ex4p1a-original $<
	./ex4p1a-original > ex4p1a-original-O0.output
	$(F90) -O3 -o ex4p1a-original $<
	./ex4p1a-original > ex4p1a-original-O3.output
	rm -f ex4p1a-original

ex4p1b: ex4p1b.f90
	$(F90) -O0 -o ex4p1b $<
	./ex4p1b > ex4p1b-O0.output
	$(F90) -O3 -o ex4p1b $<
	./ex4p1b > ex4p1b-O3.output
	rm -f ex4p1b

ex4p1b-original: ex4p1b-original.f90
	$(F90) -O0 -o ex4p1b-original $<
	./ex4p1b-original > ex4p1b-original-O0.output
	$(F90) -O3 -o ex4p1b-original $<
	./ex4p1b-original > ex4p1b-original-O3.output
	rm -f ex4p1b-original

print:
	cat *.output | less

clean:
	rm -f *.output *.o
