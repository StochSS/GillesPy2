CC=g++
CFLAGS=-c -std=c++14 -Wall -O3
LINKFLAGS = -L. -std=c++14 -Wall -O3
DEPS = model.h ssa.h
OBJ = model.o ssa.o
.PHONY: all

all: UserSimulation

model.o:
	$(CC) -c -o model.o $(CBASE_DIR)/model.cpp $(CFLAGS) -I$(CBASE_DIR)

ssa.o:
	$(CC) -c -o ssa.o $(GILLESPY_CPP_SSA_DIR)/ssa.cpp -I$(CBASE_DIR) $(CFLAGS)

UserSimulation.o:
	$(CC) -c -o UserSimulation.o UserSimulation.cpp -I$(GILLESPY_CPP_SSA_DIR) -I$(CBASE_DIR) $(CFLAGS)

UserSimulation: $(OBJ) UserSimulation.o
	$(CC) -o UserSimulation $(OBJ) UserSimulation.o -I$(CBASE_DIR) -I$(GILLESPY_CPP_SSA_DIR) $(LINKFLAGS)

cleanSimulation:
	rm -f UserSimulation

clean:
	rm -f *.o *~