

PROJECT := ProteomeRiver

include config.mk

# ---------------------------------------------------------------
# Directories
# ---------------------------------------------------------------
BASE = $(shell /bin/pwd)
SCRIPS = $(BASE)/Source
TEST = $(BASE)/tests
SRCD = $(BASE)/R

ifeq ($(UNAME),Darwin)
BASH_PROFILE=$(HOME)/.profile
else
BASH_PROFILE=$(HOME)/.bashrc
endif

intro:
	@echo ''
	@echo '===> $(PROJECT)'
	@echo ''

setdirs:
	@mkdir -p $(PREFIX)/$(PROJECT)/Source
	@mkdir -p $(PREFIX)/$(PROJECT)/test
	@mkdir -p $(PREFIX)/$(PROJECT)/R



done:
	@echo ''
	@echo '============ All Done! ============'
	@echo ' '


# ---------------------------------------------------------------
# target : help
# ---------------------------------------------------------------

help: intro
	@echo '	'
	@echo '	---------------------------'
	@echo '	Makefile - HELP'
	@echo '	---------------------------'
	@echo '	'
	@echo '	make install			: install $(PROJECT)'
	@echo '	'
	@echo '	make setup			: set environment variable PATH '
	@echo '	'
	@echo '	make uninstall			: uninstall $(PROJECT)'
	@echo ' '


# ---------------------------------------------------------------
# target : install
# ---------------------------------------------------------------

install:intro setdirs
	@echo '	'
	@echo '-----------------------------------'
	@echo 'Installing $(PROJECT) on $(PREFIX)'
	@echo '-----------------------------------'
	@echo '	'
	@Rscript -e 'devtools::document()'
	@R CMD INSTALL ./
	@cp -r $(SCRIPS)/* $(PREFIX)/$(PROJECT)/Source/
	@chmod +x $(PREFIX)/$(PROJECT)/Source/R/*
	@chmod +x $(PREFIX)/$(PROJECT)/Source/bin/*
	@cp -r $(TEST)/* $(PREFIX)/$(PROJECT)/test/
	@cp -r $(SRCD)/* $(PREFIX)/$(PROJECT)/R/
	@echo '	'
	@echo 'append $(PREFIX)/$(PROJECT)/Source/R to your PATH environment variable'
	@echo "export PATH=$(PREFIX)/$(PROJECT)/Source/R:\$$PATH"
	@echo 'or run: make setup'
	@echo '	'
	@echo " === install done ! === "

# ---------------------------------------------------------------
# target : uninstall
# ---------------------------------------------------------------

uninstall:intro
	@echo '	'
	@echo '-----------------------------------'
	@echo 'Uninstalling $(PROJECT) from $(PREFIX)'
	@echo '-----------------------------------'
	@echo '	'
	@rm -rf $(PREFIX)/$(PROJECT)
	@echo '	'
	@echo " === uninstall done ! === "


# ---------------------------------------------------------------
# target : setup
# ---------------------------------------------------------------

setup:intro
	@echo "export PATH=$(PREFIX)/$(PROJECT)/Source/R:\$$PATH" >> $(BASH_PROFILE)
	@echo '	'
	@echo '-----------------------------------'
	@echo "run: source $(BASH_PROFILE)"
	@echo '-----------------------------------'
	@echo '	'
	@echo " === setup done ! === "
