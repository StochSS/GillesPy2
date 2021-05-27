### Compiler configuration/settings ###
CC=g++
CXXFLAGS = -std=c++14 -Wall -O3
#######################################

### Input directories ###
CBASE_DIR ?= .
TEMPLATE_DIR ?= $(CBASE_DIR)/template
TEMPLATE_CPP := $(TEMPLATE_DIR)/template.cpp
SUNDIALS_SRC := $(CBASE_DIR)/Sundials/src
SUNDIALS_INC := $(CBASE_DIR)/Sundials/include
ODE_SOLVER_PATH := $(CBASE_DIR)/ode_cpp_solver
#########################

### Output directories ###
OBJ_DIR ?= $(CBASE_DIR)
OUTPUT_DIR ?= $(CBASE_DIR)
OUTPUT_FILE ?= $(OUTPUT_DIR)/Simulation.out
SUNDIALS_OBJ ?= $(OBJ_DIR)
INCLUDES := -I$(CBASE_DIR) -I$(SUNDIALS_INC) -I$(TEMPLATE_DIR)
##########################

SUNDEPS = cvode_ls.h cvode_proj.h sundials_matrix.h sunmatrix_band.h sunmatrix_dense.h sundials_iterative.h sundials_nonlinearsolver.h \
sunnonlinsol_newton.h sundials_nvector.h nvector_serial.h sundials_linearsolver.h ODESolver.h cvode.h cvode_spils.h \
sundials_types.h sundials_math.h sunlinsol_spgmr.h

SUNOBJ = cvode_nls.o cvode_io.o sundials_iterative.o cvode_proj.o sundials_matrix.o sunmatrix_band.o sunmatrix_dense.o cvode_ls.o \
sundials_linearsolver.o sundials_nonlinearsolver.o sundials_nvector_senswrapper.o sunnonlinsol_newton.o \
sundials_nvector.o nvector_serial.o cvode.o cvode_spils.o sundials_math.o sunlinsol_spgmr.o
SUNOBJ_PATHS := $(SUNOBJ:%.o=$(SUNDIALS_OBJ)/%.o)

model.o:
	$(CC) $(CXXFLAGS) -c -o $(OBJ_DIR)/model.o $(CBASE_DIR)/model.cpp -I$(CBASE_DIR)

%.o: $(SUNDIALS_SRC)/%.c
	$(CC) -c -o $(SUNDIALS_OBJ)/$@ $< $(CXXFLAGS) -I$(SUNDIALS_INC)

ODESolver.o: $(ODE_SOLVER_PATH)/ODESolver.cpp 
	$(CC) $(CXXFLAGS) -c -o $(OBJ_DIR)/ODESolver.o $(ODE_SOLVER_PATH)/ODESolver.cpp  $(INCLUDES)

ODESimulation.o: ODESolver.o $(ODE_SOLVER_PATH)/ODESimulation.cpp
	$(CC) $(CXXFLAGS) -c -o $(OBJ_DIR)/ODESimulation.o $(ODE_SOLVER_PATH)/ODESimulation.cpp $(INCLUDES)

prebuild: ODESolver.o ODESimulation.o $(SUNOBJ)

ODESimulation: ODESimulation.o model.o ODESolver.o $(SUNOBJ_PATHS)
	$(CC) $(CXXFLAGS) -o $(OUTPUT_FILE) $(SUNOBJ_PATHS) $(TEMPLATE_CPP) $(OBJ_DIR)/model.o $(OBJ_DIR)/ODESimulation.o $(OBJ_DIR)/ODESolver.o $(INCLUDES)

clean:
	rm -rf $(OUTPUT_DIR)/*.out

clean_all:
	rm -rf $(OBJ_DIR)/*.o $(SUNDIALS_OBJ)/*.o $(OUTPUT_DIR)/*.out
