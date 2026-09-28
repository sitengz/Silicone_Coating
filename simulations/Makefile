# Run with make -C simulations. Configuration output paths are relative to
# the repository root, so every generator invocation runs from there.
ROOT := ..
PYTHON ?= python3

V22_CASES := $(notdir $(patsubst %/,%,$(dir $(wildcard V22/*_*wt/model.conf))))
V35_CASES := $(notdir $(patsubst %/,%,$(dir $(wildcard V35/*_*wt/model.conf))))
V22_INFO := $(foreach case,$(V22_CASES),V22/$(case)/V22_$(case)/V22_$(case).info)
V35_INFO := $(foreach case,$(V35_CASES),V35/$(case)/V35_$(case)/V35_$(case).info)

.PHONY: all bulk bulk-v22 bulk-v35 films films-v22 films-v35 films-dry-run help

all: bulk

bulk: bulk-v22 bulk-v35

bulk-v22: $(V22_INFO)

bulk-v35: $(V35_INFO)

define V22_CASE_RULE
V22/$(1)/V22_$(1)/V22_$(1).info: V22/$(1)/model.conf $(ROOT)/build/v22_generator
	cd $(ROOT) && ./build/v22_generator --config simulations/V22/$(1)/model.conf
endef

define V35_CASE_RULE
V35/$(1)/V35_$(1)/V35_$(1).info: V35/$(1)/model.conf $(ROOT)/build/v35_generator
	cd $(ROOT) && ./build/v35_generator --config simulations/V35/$(1)/model.conf
endef

$(foreach case,$(V22_CASES),$(eval $(call V22_CASE_RULE,$(case))))
$(foreach case,$(V35_CASES),$(eval $(call V35_CASE_RULE,$(case))))

$(ROOT)/build/v22_generator $(ROOT)/build/v35_generator:
	$(MAKE) -C $(ROOT) $(@:$(ROOT)/%=%)

# Films require the final matched bulk data.<case>.npt_eq files.
# generate_films.py skips missing bulk results; it never runs LAMMPS or Slurm.
films: films-v22 films-v35

films-v22: $(ROOT)/build/v22_generator
	cd $(ROOT) && $(PYTHON) simulations/generate_films.py --formulation V22 --generator build/v22_generator

films-v35: $(ROOT)/build/v35_generator
	cd $(ROOT) && $(PYTHON) simulations/generate_films.py --formulation V35 --generator build/v35_generator

films-dry-run: $(ROOT)/build/v22_generator $(ROOT)/build/v35_generator
	cd $(ROOT) && $(PYTHON) simulations/generate_films.py --formulation V22 --generator build/v22_generator --dry-run
	cd $(ROOT) && $(PYTHON) simulations/generate_films.py --formulation V35 --generator build/v35_generator --dry-run

help:
	@printf '%s\n' 'make -C simulations                Generate all 28 bulk case packages' 'make -C simulations bulk-v22       Generate only the 14 V22 bulk packages' 'make -C simulations bulk-v35       Generate only the 14 V35 bulk packages' 'make -C simulations films-dry-run  Show films possible from completed bulk NPT data' 'make -C simulations films          Generate eligible V22 and V35 films (no submission)'
