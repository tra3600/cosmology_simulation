# Compilation : make        Tests : make test        Illustrations : make figures
CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -pedantic
BUILD    := build
LIB_SRC  := Cosmology.cpp SvgPlot.cpp Report.cpp
LIB_OBJ  := $(LIB_SRC:%.cpp=$(BUILD)/%.o)

all: $(BUILD)/cosmology_simulation

$(BUILD)/%.o: %.cpp $(wildcard *.h) | $(BUILD)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD)/cosmology_simulation: $(BUILD)/main.o $(LIB_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BUILD)/test_cosmology: tests/test_cosmology.cpp $(LIB_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BUILD):
	mkdir -p $(BUILD)

run: $(BUILD)/cosmology_simulation
	./$(BUILD)/cosmology_simulation

test: $(BUILD)/test_cosmology
	./$(BUILD)/test_cosmology

figures: $(BUILD)/cosmology_simulation
	./$(BUILD)/cosmology_simulation --svg docs/figures --no-table --no-plot > /dev/null

clean:
	rm -rf $(BUILD)

.PHONY: all run test figures clean
