CC=g++
CFLAGS=-c -std=c++14 -Wall -O3
SIMFLAGS = -std=c++14 -Wall -O3
#DEPS = model.h tau.h tau_leaper.h
#OBJ = model.o tau.o tau_leaper.o
DEPS = model.h ssa.h
OBJ = model.o ssa.o

all: UserSimulation

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

TauSimulation: model.o tau.o tau_leaper.o
	$(CC) TauSimulation.cpp $(SIMFLAGS) -o $@ $^

UserSimulation: $(OBJ)
	$(CC) UserSimulation.cpp $(SIMFLAGS) -o $@ $^

cleanSimulationTau:
	rm -f TauSimulation

cleanSimulation:
	rm -f UserSimulation

clean:
	rm -f *.o *~
