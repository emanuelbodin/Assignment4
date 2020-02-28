CFLAGS= -g -pg -Wall -lm -Ofast

galsim: galsim.o
	gcc -o galsim galsim.o

galsim.o: galsim.c 
	gcc $(CFLAGS) -c galsim.c

clean:
	rm -f ./galsim *.o
