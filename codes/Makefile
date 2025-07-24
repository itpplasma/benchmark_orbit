.PHONY: all clone-all pull-all build-all clean-all install-requirements
.PHONY: clone-simple pull-simple build-simple clean-simple
.PHONY: clone-firm3d pull-firm3d build-firm3d clean-firm3d

SIMPLE_REPO = https://github.com/itpplasma/SIMPLE.git
SIMPLE_DIR = SIMPLE

FIRM3D_REPO = https://github.com/ejpaul/firm3d.git
FIRM3D_DIR = firm3d

all: install-requirements clone-all build-all

install-requirements:
	@if [ -f "requirements.txt" ]; then \
		echo "Installing requirements from requirements.txt..."; \
		pip install -r requirements.txt; \
	else \
		echo "No requirements.txt found, skipping..."; \
	fi

clone-all: clone-simple clone-firm3d

pull-all: pull-simple pull-firm3d

build-all: build-simple build-firm3d

clean-all: clean-simple clean-firm3d

# SIMPLE targets
clone-simple:
	@if [ -d "$(SIMPLE_DIR)" ]; then \
		echo "SIMPLE directory already exists, pulling latest changes..."; \
		$(MAKE) pull-simple; \
	else \
		echo "Cloning SIMPLE repository..."; \
		git clone $(SIMPLE_REPO) $(SIMPLE_DIR); \
	fi

pull-simple:
	@if [ -d "$(SIMPLE_DIR)" ]; then \
		cd $(SIMPLE_DIR) && git pull; \
	else \
		echo "SIMPLE directory does not exist. Run 'make clone-simple' first."; \
		exit 1; \
	fi

build-simple:
	@if [ -d "$(SIMPLE_DIR)" ]; then \
		echo "Installing SIMPLE with pip..."; \
		cd $(SIMPLE_DIR) && pip install -e . --no-build-isolation; \
	else \
		echo "SIMPLE directory does not exist. Run 'make clone-simple' first."; \
		exit 1; \
	fi

clean-simple:
	@echo "Uninstalling SIMPLE..."; \
	pip uninstall -y simple 2>/dev/null || true; \
	if [ -d "$(SIMPLE_DIR)" ]; then \
		echo "Removing SIMPLE directory..."; \
		rm -rf $(SIMPLE_DIR); \
	fi

# firm3d targets
clone-firm3d:
	@if [ -d "$(FIRM3D_DIR)" ]; then \
		echo "firm3d directory already exists, pulling latest changes..."; \
		$(MAKE) pull-firm3d; \
	else \
		echo "Cloning firm3d repository..."; \
		git clone $(FIRM3D_REPO) $(FIRM3D_DIR); \
	fi

pull-firm3d:
	@if [ -d "$(FIRM3D_DIR)" ]; then \
		cd $(FIRM3D_DIR) && git pull; \
	else \
		echo "firm3d directory does not exist. Run 'make clone-firm3d' first."; \
		exit 1; \
	fi

build-firm3d:
	@if [ -d "$(FIRM3D_DIR)" ]; then \
		echo "Installing firm3d with pip..."; \
		cd $(FIRM3D_DIR) && pip install -e . --no-build-isolation; \
	else \
		echo "firm3d directory does not exist. Run 'make clone-firm3d' first."; \
		exit 1; \
	fi

clean-firm3d:
	@echo "Uninstalling firm3d..."; \
	pip uninstall -y firm3d 2>/dev/null || true; \
	if [ -d "$(FIRM3D_DIR)" ]; then \
		echo "Removing firm3d directory..."; \
		rm -rf $(FIRM3D_DIR); \
	fi
