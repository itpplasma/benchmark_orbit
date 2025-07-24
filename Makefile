.PHONY: all build run plot clean help fortran run-fortran

# Default VMEC file - can be overridden with make run VMEC_FILE=path/to/wout.nc
VMEC_FILE ?= wout.nc

# Fortran compiler and flags for SIMPLE
FC = gfortran
SIMPLE_DIR = codes/SIMPLE
SIMPLE_BUILD = $(SIMPLE_DIR)/build
SIMPLE_INCLUDE = -I$(SIMPLE_BUILD)/include -I$(SIMPLE_BUILD)/libneo/include -I/opt/homebrew/include
SIMPLE_LIBS = $(SIMPLE_BUILD)/src/libsimple.a \
              $(SIMPLE_BUILD)/libneo/libneo.a \
              $(SIMPLE_BUILD)/libneo/src/field/libneo_field.a \
              $(SIMPLE_BUILD)/libneo/src/contrib/libCONTRIB.a \
              $(SIMPLE_BUILD)/libneo/src/contrib/librkf45.a \
              $(SIMPLE_BUILD)/libneo/src/hdf5_tools/libhdf5_tools.a
FFLAGS = $(SIMPLE_INCLUDE) -fopenmp -O2 -g
LDFLAGS = $(SIMPLE_LIBS) -L/opt/homebrew/opt/netcdf-fortran/lib -L/opt/homebrew/opt/netcdf/lib \
          -L/opt/homebrew/opt/hdf5/lib -L/opt/homebrew/opt/openblas/lib \
          -lnetcdff -lnetcdf -lhdf5_fortran -lhdf5_hl_fortran -lhdf5 -lhdf5_hl \
          -lopenblas -lgomp

all: build

build:
	@echo "Building codes..."
	$(MAKE) -C codes

fortran:
	@echo "Building Fortran orbit tracer..."
	@mkdir -p build
	$(FC) $(FFLAGS) src/trace_orbit_simple.f90 -o build/trace_orbit_simple $(LDFLAGS)

run: build
	@if [ ! -f "$(VMEC_FILE)" ]; then \
		echo "Error: VMEC file '$(VMEC_FILE)' not found!"; \
		echo "Usage: make run VMEC_FILE=path/to/wout.nc"; \
		echo ""; \
		echo "You can download a test file with:"; \
		echo "wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -O wout.nc"; \
		exit 1; \
	fi
	@echo "Running orbit trace with VMEC file: $(VMEC_FILE)"
	python scripts/trace_orbit_simple.py $(VMEC_FILE)

run-fortran: fortran
	@if [ ! -f "$(VMEC_FILE)" ]; then \
		echo "Error: VMEC file '$(VMEC_FILE)' not found!"; \
		echo "Usage: make run-fortran VMEC_FILE=path/to/wout.nc"; \
		exit 1; \
	fi
	@echo "Running Fortran orbit trace with VMEC file: $(VMEC_FILE)"
	cd build && ./trace_orbit_simple ../$(VMEC_FILE)

plot:
	@if [ ! -d "run" ] || [ -z "$$(ls -A run/*.nc 2>/dev/null)" ]; then \
		echo "No orbit data found in run/ directory!"; \
		echo "Run 'make run' first to generate orbit data."; \
		exit 1; \
	fi
	@echo "Creating orbit plots..."
	python scripts/plot_orbits.py

clean:
	@echo "Cleaning build artifacts..."
	$(MAKE) -C codes clean
	rm -rf run/* plot/*
	rm -f build/trace_orbit_simple

help:
	@echo "Benchmark Orbit Makefile"
	@echo "========================"
	@echo ""
	@echo "Targets:"
	@echo "  make build          - Build all codes (SIMPLE and firm3d)"
	@echo "  make run            - Run orbit tracing script (requires VMEC file)"
	@echo "  make run VMEC_FILE=path/to/wout.nc  - Run with specific VMEC file"
	@echo "  make plot           - Create orbit visualization plots from latest run"
	@echo "  make clean          - Clean all build artifacts, run outputs, and plots"
	@echo "  make help           - Show this help message"
	@echo ""
	@echo "Codes subdirectory targets (via 'make -C codes'):"
	@echo "  install-requirements - Install Python requirements"
	@echo "  clone-all           - Clone SIMPLE and firm3d repositories"
	@echo "  build-all           - Build all codes with pip install -e ."
	@echo "  clean-all           - Uninstall and remove all codes"
	@echo ""
	@echo "Example workflow:"
	@echo "  1. make build       - Build all codes"
	@echo "  2. make run         - Run orbit tracing (needs wout.nc)"
	@echo "  3. make plot        - Generate orbit visualizations"
