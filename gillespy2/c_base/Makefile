CC=g++
CFLAGS=-c -Wall -O3
SIMFLAGS = -Wall -O3
DEPS = model.h ssa.h
OBJ = model.o ssa.o

all: MichaelisMenten

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

MichaelisMenten: $(OBJ)
	$(CC) MichaelisMenten.cpp $(SIMFLAGS) -o $@ $^

UserSimulation: $(OBJ)
	$(CC) UserSimulation.cpp $(SIMFLAGS) -o $@ $^

clean:
	rm -f *.o *~ MichaelisMenten
