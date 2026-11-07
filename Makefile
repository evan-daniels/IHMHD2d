# IHMHD2d -- 2D Hall MHD with electron inertia (CentPack central scheme, CT EM solver)
# Needs MPI C++ (mpicxx, mpirun) and FFTW3 (static libfftw3.a) at $(FFTW_DIR).

FFTW_DIR ?= /opt/homebrew/opt/fftw
CXX      ?= mpicxx
CXXFLAGS ?= -std=c++11 -O3
NP       ?= 10

CORE    = $(wildcard src/core/*.cc)
PHYSICS = $(wildcard src/physics/*.cc)

.PHONY: all test test_helmholtz test_energy regress clean

all: IHMHD2d

IHMHD2d: $(CORE) $(PHYSICS) src/main.cc $(wildcard include/*.h)
	mpicxx $(CXXFLAGS) -Iinclude/ -I$(FFTW_DIR)/include $(CORE) $(PHYSICS) src/main.cc \
	  $(FFTW_DIR)/lib/libfftw3.a -o IHMHD2d

# Unit test: Helmholtz solve/apply with wall images (round-off residuals)
tests/test_helmholtz_walls: tests/test_helmholtz_walls.cc src/physics/helmholtz_walls.cc src/physics/boundary_conditions.cc
	mpicxx -std=c++11 -O2 -Iinclude/ -I$(FFTW_DIR)/include tests/test_helmholtz_walls.cc \
	  src/physics/helmholtz_walls.cc src/physics/boundary_conditions.cc \
	  src/core/doublearray*.cc $(FFTW_DIR)/lib/libfftw3.a -o tests/test_helmholtz_walls

# Unit test: energy identities (a) H symmetric with walls, (b) staggered product rule, (4) local field-energy law
tests/test_energy_identities: tests/test_energy_identities.cc src/physics/helmholtz_walls.cc src/physics/boundary_conditions.cc
	mpicxx -std=c++11 -O2 -Iinclude/ -I$(FFTW_DIR)/include tests/test_energy_identities.cc \
	  src/physics/helmholtz_walls.cc src/physics/boundary_conditions.cc \
	  src/core/doublearray*.cc $(FFTW_DIR)/lib/libfftw3.a -o tests/test_energy_identities

test_helmholtz: tests/test_helmholtz_walls
	cd tests && mpirun -np 1 ./test_helmholtz_walls && mpirun -np 4 ./test_helmholtz_walls

test_energy: tests/test_energy_identities
	cd tests && mpirun -np 1 ./test_energy_identities && mpirun -np 4 ./test_energy_identities

# d_e = 0 regression: GEM 1x to t = 1 on $(NP) ranks, bit-for-bit against tests/reference/ (CentPackL MHD2d)
regress: IHMHD2d
	tests/regress_de0.sh $(NP)

test: test_helmholtz test_energy regress

clean:
	rm -rf IHMHD2d IHMHD2d.dSYM tests/test_helmholtz_walls tests/test_energy_identities tests/*.dSYM tests/run_de0
