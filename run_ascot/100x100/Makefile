# Makefile for ASCOT5 orbit tracing with VMEC equilibrium
#
# Uses the same VMEC equilibrium as SIMPLE orbit tracer (wout.nc symlink)
# Traces orbits from same initial conditions as SIMPLE for comparison
#
# Targets:
#   make setup       - Setup ASCOT5 HDF5 file with VMEC field
#   make plot        - Plot equilibrium (PNG + interactive display)
#   make setup_mark  - Create markers from SIMPLE start.dat
#   make run         - Trace orbits (takes several minutes)
#   make extract     - Extract orbits to NetCDF
#   make plot_orbit  - Plot individual orbit (ORBIT=<id>)
#   make all         - Run full workflow
#   make help        - Show this help message
#   make clean       - Remove generated files

PYTHON := python3
SETUP_SCRIPT := setup_ascot.py
PLOT_SCRIPT := plot_equilibrium.py
MARKERS_SCRIPT := setup_markers.py
RUN_SCRIPT := run_orbits.py
EXTRACT_SCRIPT := extract_orbits.py
PLOT_ORBIT_SCRIPT := plot_orbit_ascot.py

VMEC_FILE := wout.nc
H5_FILE := ascot.h5
START_FILE := ../start.dat
MARKERS_H5 := $(H5_FILE)
RESULTS_H5 := orbits.h5
ORBITS_NC := orbits_ascot.nc
EQUILIBRIUM_PNG := equilibrium.png

.PHONY: setup plot setup_mark run extract plot_orbit all help clean check-vmec check-start

# Check that VMEC file exists
check-vmec:
	@if [ ! -f "$(VMEC_FILE)" ]; then \
		echo "Error: VMEC file '$(VMEC_FILE)' not found"; \
		exit 1; \
	fi
	@echo "Found VMEC file: $(VMEC_FILE)"

# Check that START file exists
check-start:
	@if [ ! -f "$(START_FILE)" ]; then \
		echo "Error: START file '$(START_FILE)' not found"; \
		echo "Run: cd .. && python setup_start.py"; \
		exit 1; \
	fi
	@echo "Found START file: $(START_FILE)"

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
	@echo "Done! Check output: $(EQUILIBRIUM_PNG)"

# Setup markers from SIMPLE start.dat
setup_mark: check-start $(H5_FILE) $(MARKERS_SCRIPT)
	@echo "=========================================="
	@echo "Setting up ASCOT5 Markers from SIMPLE"
	@echo "=========================================="
	$(PYTHON) $(MARKERS_SCRIPT)
	@echo ""
	@echo "Done! Markers added to: $(H5_FILE)"

# Run orbit tracing simulation
run: setup_mark $(RUN_SCRIPT)
	@echo "=========================================="
	@echo "Running ASCOT5 Orbit Tracing"
	@echo "=========================================="
	$(PYTHON) $(RUN_SCRIPT)
	@echo ""
	@echo "Done! Results saved to: $(RESULTS_H5)"

# Extract orbits to NetCDF
extract: $(RESULTS_H5) $(EXTRACT_SCRIPT)
	@echo "=========================================="
	@echo "Extracting ASCOT5 Orbits to NetCDF"
	@echo "=========================================="
	$(PYTHON) $(EXTRACT_SCRIPT)
	@echo ""
	@echo "Done! Orbits saved to: $(ORBITS_NC)"

# Plot individual orbit
plot_orbit: $(ORBITS_NC) $(PLOT_ORBIT_SCRIPT)
	@if [ -z "$(ORBIT)" ]; then \
		echo "Usage: make plot_orbit ORBIT=<particle_id>"; \
		echo ""; \
		echo "Example: make plot_orbit ORBIT=1"; \
		echo ""; \
		echo "To see available particles:"; \
		echo "  python3 -c \"import xarray as xr; ds = xr.open_dataset('$(ORBITS_NC)'); print('Particles:', sorted(ds.particle.values)[:20])\""; \
		exit 1; \
	fi
	@echo "Plotting orbit $(ORBIT)..."
	$(PYTHON) $(PLOT_ORBIT_SCRIPT) $(ORBIT)

# Generate HDF5 file if it doesn't exist
$(H5_FILE): check-vmec $(SETUP_SCRIPT)
	$(PYTHON) $(SETUP_SCRIPT)

# Full workflow
all: setup plot setup_mark run extract
	@echo ""
	@echo "=========================================="
	@echo "Full workflow complete!"
	@echo "=========================================="
	@echo "Try: make plot_orbit ORBIT=1"

# Show help
help:
	@echo "ASCOT5 Orbit Tracing with VMEC Equilibrium"
	@echo ""
	@echo "Available targets:"
	@echo "  make setup        - Setup ASCOT5 HDF5 with VMEC field"
	@echo "  make plot         - Plot equilibrium (interactive + PNG)"
	@echo "  make setup_mark   - Create markers from SIMPLE start.dat"
	@echo "  make run          - Trace orbits (several minutes)"
	@echo "  make extract      - Extract orbits to NetCDF"
	@echo "  make plot_orbit   - Plot orbit (use ORBIT=<id>)"
	@echo "  make all          - Run full workflow"
	@echo "  make clean        - Remove generated files"
	@echo "  make help         - Show this help message"
	@echo ""
	@echo "Quick start:"
	@echo "  1. make setup          # Initialize HDF5 with VMEC field"
	@echo "  2. make plot           # View equilibrium"
	@echo "  3. make setup_mark     # Add markers from SIMPLE"
	@echo "  4. make run            # Trace orbits (takes ~5-10 min)"
	@echo "  5. make extract        # Extract to NetCDF"
	@echo "  6. make plot_orbit ORBIT=1  # Visualize orbit 1"
	@echo ""
	@echo "Configuration:"
	@echo "  - VMEC file: $(VMEC_FILE)"
	@echo "  - START file: $(START_FILE)"
	@echo "  - HDF5 output: $(H5_FILE)"
	@echo "  - Results: $(RESULTS_H5)"
	@echo "  - Orbits NetCDF: $(ORBITS_NC)"
	@echo ""
	@echo "Notes:"
	@echo "  - Requires start.dat from SIMPLE (run ../setup_start.py)"
	@echo "  - Same initial conditions and runtime as SIMPLE (1 ms)"
	@echo "  - Uses guiding center (GC) approximation"

# Clean generated files
clean:
	@echo "Cleaning generated files..."
	rm -f $(H5_FILE) $(RESULTS_H5) $(ORBITS_NC) $(EQUILIBRIUM_PNG) orbit_ascot_*.png
	@echo "Clean complete"
