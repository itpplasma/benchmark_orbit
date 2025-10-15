.PHONY: plot help clean

PYTHON := python3

# Default target
help:
	@echo "FIRM3D Orbit Analysis Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  make plot FILE=<trajectory_file>  - Plot trajectory file"
	@echo "  make help                         - Show this help"
	@echo "  make clean                        - Remove generated plots"
	@echo ""
	@echo "Examples:"
	@echo "  make plot FILE=trajectory_data_vmec_tol_1e-06_resolution_96_tmax_0.001_trapped.txt"
	@echo "  make plot FILE=trajectory_data_tol_1e-06_resolution_48_tmax_0.001_trapped.txt"
	@echo ""
	@echo "Available trajectory files:"
	@ls -1 trajectory_*.txt 2>/dev/null | head -5 || echo "  (none found)"
	@echo "  ..."

# Plot a trajectory file
plot:
	@if [ -z "$(FILE)" ]; then \
		echo "Usage: make plot FILE=<trajectory_file>"; \
		echo ""; \
		echo "Example: make plot FILE=trajectory_data_vmec_tol_1e-06_resolution_96_tmax_0.001_trapped.txt"; \
		echo ""; \
		echo "Available files:"; \
		ls -1 trajectory_*.txt 2>/dev/null | head -10 || echo "  (none found)"; \
		exit 1; \
	fi
	@if [ ! -f "$(FILE)" ]; then \
		echo "Error: File not found: $(FILE)"; \
		exit 1; \
	fi
	$(PYTHON) plot_orbit_firm3d.py $(FILE)

# Clean up generated files
clean:
	rm -f trajectory_*.png
