FC=nvfortran
FFLAGS=-cuda 
INCLUDE=
LDFLAGS=-L /oscar/rt/9.2/software/0.20-generic/0.20.1/opt/spack/linux-rhel9-x86_64_v3/gcc-11.3.1/nvhpc-25.1-scvc2acaswih5fk5wowmbwgqstecutym/Linux_x86_64/25.1/math_libs/lib64 -lcufft  


# SOURCE Files
MODULE_SRC = cu_fft_module.cuf \
             cu_fft.cuf

MAIN_SRC   = main.cuf

# OBJECT Files
MODULE_OBJ = $(MODULE_SRC:.cuf=.o)
MAIN_OBJ   = $(MAIN_SRC:.cuf=.o)

# Executable
EXE = hit_gpu

# Default target
all: $(EXE)

# Link executable
$(EXE): $(MODULE_OBJ) $(MAIN_OBJ)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Compile modules first
%.o: %.cuf
	$(FC) $(FFLAGS) -c $<

# Clean
clean:
	rm -f $(EXE) *.o *.mod *.plt *.dat


