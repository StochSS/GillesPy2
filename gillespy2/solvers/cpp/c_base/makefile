CC=g++
CFLAGS=-c -std=c++14 -Wall -O3
SIMFLAGS = -L. -std=c++14 -Wall -O3 -lsundials_cvode -lsundials_nvecserial
DEPS = model.h ODECSolver.h
OBJ = model.o ODECSolver.o
.PHONY: all

all: UserSimulation

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

UserSimulation.o:
	$(CC) -c -o UserSimulation.o UserSimulation.cpp $(CFLAGS)

UserSimulation: $(OBJ) UserSimulation.o
	$(CC) -o UserSimulation $(OBJ) UserSimulation.o $(SIMFLAGS)

ODESimulation:  $(OBJ) UserSimulation.o
	$(CC) -o ODESimulation model.o ODECSolver.o UserSimulation.o $(SIMFLAGS)

cleanSimulationTau:
	rm -f TauSimulation

cleanSimulation:
	rm -f UserSimulation

clean:
	rm -f *.o *~