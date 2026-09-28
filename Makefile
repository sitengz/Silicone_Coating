CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic
PYTHON ?= python3

BUILD_DIR := build
GENERATORS := $(BUILD_DIR)/v22_generator $(BUILD_DIR)/v35_generator $(BUILD_DIR)/oil_generator
ANALYZERS := $(addprefix $(BUILD_DIR)/,z_profile final_snapshot_analyzer msd_analyzer phase_separation_analyzer)

.PHONY: all generators analyzers test help

all: generators analyzers

generators: $(GENERATORS)

analyzers: $(ANALYZERS)

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/v22_generator: V22/v22_generator.cpp V35/v35_generator.cpp Formulation/config_input.hpp Formulation/silicone_oil_component.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BUILD_DIR)/v35_generator: V35/v35_generator.cpp Formulation/config_input.hpp Formulation/silicone_oil_component.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BUILD_DIR)/oil_generator: Oil/oil_generator.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BUILD_DIR)/%: Analysis/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

test:
	$(PYTHON) -m unittest discover -s Analysis -p 'test_*.py' -v

help:
	@printf '%s\n' 'make                 Build all generators and C++ analyzers in build/' 'make generators      Build V22, V35, and Oil generators' 'make analyzers       Build the four C++ analyzers' 'make test            Run Python analyzer tests' 'make -C simulations  Generate all 28 bulk case packages' 'make -C simulations films  Generate films after bulk NPT runs'
