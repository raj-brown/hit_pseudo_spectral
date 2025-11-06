FC=gfortran-11 -Wall -fimplicit-none -J.
FFLAGS=-O3 #-Wall -Wextra -fdefault-real-8 -fimplicit-none -std=f2003 
INCLUDE=-I/mnt/projects/gits/hit_pseudo_spectral/tpl/fftw-3.3.10/include


FFT_LIB=/mnt/projects/gits/hit_pseudo_spectral/tpl/fftw-3.3.10/lib
LDFLAGS=-L$(FFT_LIB) -Wl,-rpath=$(FFT_LIB) -lfftw3 -lm
$(info LDFLAGS is $(LDFLAGS))

# SOURCE File
MODULE_SRC = fft_module.f90
MAIN_SRC = main.f90

#OBJECT File
MODULE_OBJ = $(MODULE_SRC:.f90=.o)
MAIN_OBJ = $(MAIN_SRC:.f90=.o)

# Executable
EXE = hit_cpu

all: $(EXE)

# Link Executable
$(EXE): $(MODULE_OBJ) $(MAIN_OBJ)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Compile Module
$(MODULE_OBJ): $(MODULE_SRC)
	$(FC) $(INCLUDE) $(FFLAGS) -c $<

# Compile main programs
$(MAIN_OBJ): $(MAIN_SRC) $(MODULE_OBJ)
	$(FC) $(INCLUDE) $(FFLAGS) -c $<


# Clean
clean:
	rm -f $(EXE) *.o *.mod
