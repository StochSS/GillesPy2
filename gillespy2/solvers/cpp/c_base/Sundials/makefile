CC=g++
CPPFLAGS = -c -I./include/ -I./src/ -I../ -std=c++14 -Wall -O3
SIMFLAGS = -L. -std=c++14 -Wall -O3

DEPS = model.h cvode_ls.h cvode_proj.h sundials_matrix.h sunmatrix_band.h sunmatrix_dense.h sundials_iterative.h sundials_nonlinearsolver.h \
sunnonlinsol_newton.h sundials_nvector.h nvector_serial.h sundials_linearsolver.h ODECSolver.h cvode.h cvode_spils.h \
sundials_types.h sundials_math.h sunlinsol_spgmr.h

OBJ = model.o cvode_nls.o cvode_io.o sundials_iterative.o cvode_proj.o sundials_matrix.o sunmatrix_band.o sunmatrix_dense.o cvode_ls.o \
sundials_linearsolver.o sundials_nonlinearsolver.o sundials_nvector_senswrapper.o sunnonlinsol_newton.o \
sundials_nvector.o nvector_serial.o cvode.o cvode_spils.o sundials_math.o sunlinsol_spgmr.o ODECSolver.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $ $

ODESimulation.o:
	$(CC) -c -o ODESimulation.o ODESimulation.cpp $(CPPFLAGS)

ODESimulation:  $(OBJ) UserSimulation.o
	$(CC) -o ODESimulation $(OBJ) ODECSolver.o UserSimulation.o $(SIMFLAGS)

cleanSimulationODE:
	rm -f ODESimulation

clean:
	rm -f *.o *~