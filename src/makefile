# Makefile for building the C language shared library
C = gcc
OPTS = -c -fPIC
LOADER = gcc

OBJECTS = graphclass.o

 graphclass.so: $(OBJECTS)
	R CMD SHLIB -o graphclass.so $(OBJECTS)

.c.o: ; $(C) $(OPTS) -c $<

clean:
	-rm *.o *.so