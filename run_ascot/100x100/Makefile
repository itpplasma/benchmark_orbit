# Makefile for ASCOT5 orbit tracing with VMEC equilibrium
#
# Uses the same VMEC equilibrium as SIMPLE orbit tracer (wout.nc symlink)
#
# Targets:
#   make setup       - Setup ASCOT5 HDF5 file with VMEC field
#   make plot        - Plot equilibrium (PNG + interactive display)
#   make help        - Show this help message
#   make clean       - Remove generated files

PYTHON := python3
SETUP_SCRIPT := setup_ascot.py
PLOT_SCRIPT := plot_equilibrium.py
VMEC_FILE := wout.nc
H5_FILE := ascot.h5
PLOT_PNG := equilibrium.png

.PHONY: setup plot help clean check-vmec

# Check that VMEC file exists
check-vmec:
	@if [ ! -f "$(VMEC_FILE)" ]; then \
		echo "Error: VMEC file '$(VMEC_FILE)' not found"; \
		echo "Expected: $(shell pwd)/$(VMEC_FILE)"; \
		exit 1; \
	fi
	@echo "Found VMEC file: $(VMEC_FILE)"

# Setup ASCOT5 HDF5 file with VMEC field
setup: check-vmec $(SETUP_SCRIPT)
	@echo "=========================================="
	@echo "Setting up ASCOT5 HDF5 with VMEC Field"
	@echo "=========================================="
	$(PYTHON) $(SETUP_SCRIPT)
	@echo ""
	@echo "Done! HDF5 file: $(H5_FILE)"

# Plot equilibrium (depends on HDF5 file)
plot: $(H5_FILE) $(PLOT_SCRIPT)
	@echo "=========================================="
	@echo "Plotting ASCOT5 VMEC Equilibrium"
	@echo "=========================================="
	$(PYTHON) $(PLOT_SCRIPT)
	@echo ""
	@echo "Done! Check outputs:"
	@echo "  - PNG: $(PLOT_PNG)"

# Generate HDF5 file if it doesn't exist
$(H5_FILE): check-vmec $(SETUP_SCRIPT)
	$(PYTHON) $(SETUP_SCRIPT)

# Show help
help:
	@echo "ASCOT5 Orbit Tracing with VMEC Equilibrium"
	@echo ""
	@echo "Available targets:"
	@echo "  make setup  - Setup ASCOT5 HDF5 file with VMEC field"
	@echo "  make plot   - Plot equilibrium (PNG + interactive)"
	@echo "  make clean  - Remove generated files"
	@echo "  make help   - Show this help message"
	@echo ""
	@echo "Configuration:"
	@echo "  - VMEC file: $(VMEC_FILE)"
	@echo "  - HDF5 output: $(H5_FILE)"
	@echo "  - PNG output: $(PLOT_PNG)"
	@echo ""
	@echo "Workflow:"
	@echo "  1. make setup  - Create HDF5 with VMEC field"
	@echo "  2. make plot   - Visualize equilibrium (auto-runs setup if needed)"
	@echo ""
	@echo "Notes:"
	@echo "  - Uses same VMEC equilibrium as SIMPLE orbit tracer"
	@echo "  - Future: convert initial conditions and trace orbits in ASCOT"

# Clean generated files
clean:
	@echo "Cleaning generated files..."
	rm -f $(H5_FILE) $(PLOT_PNG)
	@echo "Clean complete"
