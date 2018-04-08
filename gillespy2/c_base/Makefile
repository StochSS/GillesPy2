CC=g++
CFLAGS=-std=c++14 -Wall -O3
SIMFLAGS = -std=c++14 -Wall -O3
DEPS = model.h ssa.h
OBJ = model.o ssa.o

all: MichaelisMenten

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

MichaelisMenten: $(OBJ)
	$(CC) MichaelisMenten.cpp $(SIMFLAGS) -o $@ $^

UserSimulation: $(OBJ)
	$(CC) UserSimulation.cpp $(SIMFLAGS) -o $@ $^

cleanSimulation:
	rm -f UserSimulation

clean:
	rm -f *.o *~ MichaelisMenten
