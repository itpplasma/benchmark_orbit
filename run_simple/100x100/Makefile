.PHONY: plot help clean

PYTHON := python3

# Default target
help:
	@echo "SIMPLE Orbit Analysis Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  make plot ORBIT=<id>   - Plot orbit <id> (requires orbits.nc)"
	@echo "  make plot              - Show help for plot target"
	@echo "  make clean             - Remove generated files"
	@echo ""
	@echo "Examples:"
	@echo "  make plot ORBIT=1"
	@echo "  make plot ORBIT=100"
	@echo ""
	@echo "Note: orbits.nc is now generated directly by SIMPLE"

# Plot a specific orbit
plot:
	@if [ -z "$(ORBIT)" ]; then \
		echo "Usage: make plot ORBIT=<particle_id>"; \
		echo ""; \
		echo "Example: make plot ORBIT=1"; \
		echo ""; \
		echo "To see available particles:"; \
		echo "  python3 -c \"import xarray as xr; ds = xr.open_dataset('orbits.nc'); print('Particles:', sorted(ds.particle.values))\""; \
		exit 1; \
	fi
	@if [ ! -f orbits.nc ]; then \
		echo "Error: orbits.nc not found. Run SIMPLE first."; \
		exit 1; \
	fi
	$(PYTHON) plot_orbit.py $(ORBIT)

# Clean up generated files
clean:
	rm -f orbits.nc orbit_*.png

.PHONY: help plot clean
