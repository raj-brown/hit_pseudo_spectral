FC=gfortran 
FFLAGS=-O3 -Wall -Wextra -fdefault-real-8 -fimplicit-none -std=f2003 
INCLUDE=-I/usr/local/include


FFT_LIB:=/usr/local/lib
LDFLAGS+=-L$(FFT_LIB) -lfftw3 -lm

all: main

main.o: main.f90
	$(FC) $(INCLUDE) $(FFLAGS) -c main.f90

module.o: fft_module.f90
	$(FC) $(INCLUDE) $(FFLAGS) -c main.f90


main: main.o
	$(FC) -o main main.o $(LDFLAGS)

clean:
	rm -f *.o *.mod main *.tec *.dat
