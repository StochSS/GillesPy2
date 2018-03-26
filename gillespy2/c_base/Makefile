CC=g++
CFLAGS=-c -Wall -O3
DEPS = model.h ssa.h
OBJ = model.o ssa.o

all: MichaelisMenten

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

MichaelisMenten: $(OBJ)
	$(CC) MichaelisMenten.cpp -o $@ $(CFLAGS) $^

UserSimulation: $(OBJ)
	$(CC) UserSimulation.cpp -o $@ $(CFLAGS) $^

clean:
	rm -f *.o *~ MichaelisMenten
