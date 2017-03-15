# Makefile for building the C language shared library for the CalcMatrixAvg demonstration package.
C = gcc
OPTS = -c -fPIC
LOADER = gcc

OBJECTS = Pmax2.o

 Pmax2.so: $(OBJECTS)
	R CMD SHLIB -o Pmax2.so $(OBJECTS)

.c.o: ; $(C) $(OPTS) -c $<

clean:
	-rm *.o *.so