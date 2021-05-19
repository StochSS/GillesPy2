CC=g++
CFLAGS = -std=c++14 -Wall -O3
LINKFLAGS = -L. -std=c++14 -Wall -O3
SRCDIR = $(SUNDIALS_DIR)/src
INCDIR = $(SUNDIALS_DIR)/include

INCLUDES = -I$(CBASE_DIR) -I$(SUNDIALS_DIR) -I$(INCDIR) -I$(CBASE_DIR)/template
OBJDIR = obj
BINDIR = bin

SUNDEPS = cvode_ls.h cvode_proj.h sundials_matrix.h sunmatrix_band.h sunmatrix_dense.h sundials_iterative.h sundials_nonlinearsolver.h \
sunnonlinsol_newton.h sundials_nvector.h nvector_serial.h sundials_linearsolver.h ODECSolver.h cvode.h cvode_spils.h \
sundials_types.h sundials_math.h sunlinsol_spgmr.h

SUNOBJ = cvode_nls.o cvode_io.o sundials_iterative.o cvode_proj.o sundials_matrix.o sunmatrix_band.o sunmatrix_dense.o cvode_ls.o \
sundials_linearsolver.o sundials_nonlinearsolver.o sundials_nvector_senswrapper.o sunnonlinsol_newton.o \
sundials_nvector.o nvector_serial.o cvode.o cvode_spils.o sundials_math.o sunlinsol_spgmr.o
SUNOBJ_PATHS=$(SUNOBJ:%.o=$(OBJDIR)/%.o)

TEMPLATE_CPP = $(CBASE_DIR)/template/XBaseTemplate.cpp
ODE_SOLVER_PATH = $(CBASE_DIR)/ode_cpp_solver

model.o:
	$(CC) -c -o $(OBJDIR)/model.o $(CBASE_DIR)/model.cpp $(CFLAGS) -I$(CBASE_DIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) -c -o $@ $< $(CFLAGS) -I$(INCDIR)

ODECSolver.o:
	$(CC) -c -o $(OBJDIR)/ODECSolver.o $(ODE_SOLVER_PATH)/ODECSolver.cpp $(INCLUDES)

ODESimulation.o: ODECSolver.o $(CBASE_DIR)/ode_cpp_solver/ODETemplate.cpp
	$(CC) -c -o $(OBJDIR)/ODESimulation.o $(ODE_SOLVER_PATH)/ODETemplate.cpp $(INCLUDES)

ODESimulation: ODESimulation.o model.o ODECSolver.o $(SUNOBJ_PATHS)
	$(CC) -o $(BINDIR)/ODESimulation -pedantic $(SUNOBJ_PATHS) $(TEMPLATE_CPP) $(OBJDIR)/model.o $(OBJDIR)/ODESimulation.o $(OBJDIR)/ODECSolver.o $(INCLUDES)

clean:
	rm -f obj/*.o bin/*