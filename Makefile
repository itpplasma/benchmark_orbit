.PHONY: all build run clean help

# Default VMEC file - can be overridden with make run VMEC_FILE=path/to/wout.nc
VMEC_FILE ?= wout.nc

all: build

build:
	@echo "Building codes..."
	$(MAKE) -C codes

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
	python scripts/trace_orbit.py $(VMEC_FILE)

clean:
	@echo "Cleaning build artifacts..."
	$(MAKE) -C codes clean
	rm -rf run/*

help:
	@echo "Benchmark Orbit Makefile"
	@echo "========================"
	@echo ""
	@echo "Targets:"
	@echo "  make build          - Build all codes (SIMPLE and firm3d)"
	@echo "  make run            - Run orbit tracing script (requires VMEC file)"
	@echo "  make run VMEC_FILE=path/to/wout.nc  - Run with specific VMEC file"
	@echo "  make clean          - Clean all build artifacts and run outputs"
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