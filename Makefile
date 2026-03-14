CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -O2
INCLUDES = -I include

SOURCES = src/main.cpp \
          src/grid/Grid.cpp \
          src/boundary/BoundaryConditions.cpp \
          src/solver/NavierStokesSolver.cpp

TARGET = NavierStockSolver

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SOURCES) -o $(TARGET)

clean:
	rm -f $(TARGET)

.PHONY: all clean
