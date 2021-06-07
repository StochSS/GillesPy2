CC=g++
CFLAGS = -c  -std=c++14 -Wall -O3
LINKFLAGS = -L. -std=c++14 -Wall -O3
SRCDIR = $(SUNDIALS_DIR)/src
INCDIR = $(SUNDIALS_DIR)/include

SIMDEPS = model.h
SIMOBJ = model.o

SUNDEPS = cvode_ls.h cvode_proj.h sundials_matrix.h sunmatrix_band.h sunmatrix_dense.h sundials_iterative.h sundials_nonlinearsolver.h \
sunnonlinsol_newton.h sundials_nvector.h nvector_serial.h sundials_linearsolver.h ODECSolver.h cvode.h cvode_spils.h \
sundials_types.h sundials_math.h sunlinsol_spgmr.h

SUNOBJ = cvode_nls.o cvode_io.o sundials_iterative.o cvode_proj.o sundials_matrix.o sunmatrix_band.o sunmatrix_dense.o cvode_ls.o \
sundials_linearsolver.o sundials_nonlinearsolver.o sundials_nvector_senswrapper.o sunnonlinsol_newton.o \
sundials_nvector.o nvector_serial.o cvode.o cvode_spils.o sundials_math.o sunlinsol_spgmr.o

model.o:
	$(CC) -c -o model.o $(CBASE_DIR)/model.cpp $(CFLAGS) -I$(CBASE_DIR)

%.o: $(SRCDIR)/%.c
	$(CC) -c -o $@ $< $(CFLAGS) -I$(INCDIR)

ODESimulation.o:
	$(CC) -c -o ODESimulation.o ODESimulation.cpp $(CFLAGS) -I$(CBASE_DIR) -I$(GILLESPY_C_ODE_DIR)

ODECSolver.o:
	$(CC) -c -o ODECSolver.o $(GILLESPY_C_ODE_DIR)/ODECSolver.cpp $(CFLAGS) -I$(CBASE_DIR) -I$(SRCDIR) -I$(INCDIR)

ODESimulation: ODECSolver.o ODESimulation.o $(SUNOBJ) model.o
	$(CC) -o ODESimulation ODECSolver.o ODESimulation.o $(SUNOBJ) model.o $(LINKFLAGS)

cleanSimulationODE:
	rm -f ODESimulation

clean:
	rm -f *.o *~