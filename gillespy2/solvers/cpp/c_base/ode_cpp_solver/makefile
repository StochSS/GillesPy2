CC=g++
CPPFLAGS = -c  -std=c++14 -Wall -O3
SIMFLAGS = -L. -std=c++14 -Wall -O3
SRCDIR = ../Sundials/src
INCDIR = ../Sundials/include

SIMDEPS = model.h
SIMOBJ = model.o

SUNDEPS = cvode_ls.h cvode_proj.h sundials_matrix.h sunmatrix_band.h sunmatrix_dense.h sundials_iterative.h sundials_nonlinearsolver.h \
sunnonlinsol_newton.h sundials_nvector.h nvector_serial.h sundials_linearsolver.h ODECSolver.h cvode.h cvode_spils.h \
sundials_types.h sundials_math.h sunlinsol_spgmr.h

SUNOBJ = cvode_nls.o cvode_io.o sundials_iterative.o cvode_proj.o sundials_matrix.o sunmatrix_band.o sunmatrix_dense.o cvode_ls.o \
sundials_linearsolver.o sundials_nonlinearsolver.o sundials_nvector_senswrapper.o sunnonlinsol_newton.o \
sundials_nvector.o nvector_serial.o cvode.o cvode_spils.o sundials_math.o sunlinsol_spgmr.o ODECSolver.o


%.o: $(SRCDIR)/%.c
	$(CC) -c -o $@ $< -I$(INCDIR) $(CPPLAGS)

ODESimulation.o:
	$(CC) -c -o ODESimulation.o ODESimulation.cpp -I../ $(CPPFLAGS)

ODECSolver.o: $(SUNOBJ) $(SIMOBJ)
	$(CC) -c -o ODECSolver.o ODECSolver.cpp -I../ -I$(SRCDIR) -I$(INCDIR) $(OBJ) $(SIMOBJ)

ODESimulation: ODECSolver.o ODESimulation.o
	$(CC) -o ODESimulation ODECSolver.o ODESimulation.o $(SIMFLAGS)

cleanSimulationODE:
	rm -f ODESimulation

clean:
	rm -f *.o *~