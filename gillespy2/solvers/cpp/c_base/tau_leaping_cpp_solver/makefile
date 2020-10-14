CC=g++
CFLAGS = -c  -std=c++14 -Wall -O3
LINKFLAGS = -L. -std=c++14 -Wall -O3
DEPS = model.h tau_leaper.h
OBJ = tau_leaper.o model.o

model.o:
	$(CC) -c -o model.o $(CBASE_DIR)/model.cpp $(CFLAGS) -I$(CBASE_DIR)

tau_leaper.o:
	$(CC) -c -o tau_leaper.o $(GILLESPY_CPP_TAU_DIR)/tau_leaper.cpp $(CFLAGS) -I$(CBASE_DIR) -I$(GILLESPY_CPP_TAU_DIR)

TauSimulation.o:
	$(CC) -c -o TauSimulation.o TauSimulation.cpp $(CFLAGS) -I$(CBASE_DIR) -I$(GILLESPY_CPP_TAU_DIR)

TauSimulation: TauSimulation.o $(OBJ)
	$(CC) -o TauSimulation TauSimulation.o $(OBJ) $(LINKFLAGS)

cleanSimulationTau:
	rm -f TauSimulation

clean:
	rm -f *.o *~