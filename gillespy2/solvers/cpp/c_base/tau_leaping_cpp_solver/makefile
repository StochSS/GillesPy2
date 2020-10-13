CC=g++
CPPFLAGS = -c  -std=c++14 -Wall -O3
LINKFLAGS = -L. -std=c++14 -Wall -O3
DEPS = model.h tau_leaper.h
OBJ = tau_leaper.o model.o

model.o:
	$(CC) -c -o model.o $(CBASE_DIR)/model.cpp -I$(CBASE_DIR) $(CPPFLAGS)

tau_leaper.o:
	$(CC) -c -o tau_leaper.o $(GILLESPY_CPP_TAU_DIR)/tau_leaper.cpp -I$(CBASE_DIR) -I$(GILLESPY_CPP_TAU_DIR) $(CPPFLAGS)

TauSimulation.o:
	$(CC) -c -o TauSimulation.o TauSimulation.cpp -I$(CBASE_DIR) -I$(GILLESPY_CPP_TAU_DIR) $(CPPFLAGS)

TauSimulation: TauSimulation.o $(OBJ)
	$(CC) -o TauSimulation TauSimulation.o $(OBJ) $(LINKFLAGS)

cleanSimulationTau:
	rm -f TauSimulation

clean:
	rm -f *.o *~