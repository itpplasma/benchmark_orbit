.PHONY: all clone pull build clean

SIMPLE_REPO = https://github.com/itpplasma/SIMPLE.git
SIMPLE_DIR = SIMPLE

all: clone build

clone:
	@if [ -d "$(SIMPLE_DIR)" ]; then \
		echo "SIMPLE directory already exists, pulling latest changes..."; \
		$(MAKE) pull; \
	else \
		echo "Cloning SIMPLE repository..."; \
		git clone $(SIMPLE_REPO) $(SIMPLE_DIR); \
	fi

pull:
	@if [ -d "$(SIMPLE_DIR)" ]; then \
		cd $(SIMPLE_DIR) && git pull; \
	else \
		echo "SIMPLE directory does not exist. Run 'make clone' first."; \
		exit 1; \
	fi

build:
	@if [ -d "$(SIMPLE_DIR)" ]; then \
		echo "Building SIMPLE..."; \
		$(MAKE) -C $(SIMPLE_DIR); \
	else \
		echo "SIMPLE directory does not exist. Run 'make clone' first."; \
		exit 1; \
	fi

clean:
	@if [ -d "$(SIMPLE_DIR)" ]; then \
		$(MAKE) -C $(SIMPLE_DIR) clean; \
	fi