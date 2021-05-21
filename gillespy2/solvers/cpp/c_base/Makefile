### Compiler configuration/settings ###
CC=g++
CFLAGS = -std=c++14 -Wall -O3
#######################################

### Input directories ###
CBASE_DIR ?= .
TEMPLATE_DIR ?= $(CBASE_DIR)/template
TEMPLATE_CPP := $(TEMPLATE_DIR)/XBaseTemplate.cpp
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
sunnonlinsol_newton.h sundials_nvector.h nvector_serial.h sundials_linearsolver.h ODECSolver.h cvode.h cvode_spils.h \
sundials_types.h sundials_math.h sunlinsol_spgmr.h

SUNOBJ = cvode_nls.o cvode_io.o sundials_iterative.o cvode_proj.o sundials_matrix.o sunmatrix_band.o sunmatrix_dense.o cvode_ls.o \
sundials_linearsolver.o sundials_nonlinearsolver.o sundials_nvector_senswrapper.o sunnonlinsol_newton.o \
sundials_nvector.o nvector_serial.o cvode.o cvode_spils.o sundials_math.o sunlinsol_spgmr.o
SUNOBJ_PATHS := $(SUNOBJ:%.o=$(SUNDIALS_OBJ)/%.o)

model.o:
	$(CC) -c -o $(OBJ_DIR)/model.o $(CBASE_DIR)/model.cpp $(CFLAGS) -I$(CBASE_DIR)

%.o: $(SUNDIALS_SRC)/%.c
	$(CC) -c -o $(SUNDIALS_OBJ)/$@ $< $(CFLAGS) -I$(SUNDIALS_INC)

ODECSolver.o: $(ODE_SOLVER_PATH)/ODECSolver.cpp 
	$(CC) -c -o $(OBJ_DIR)/ODECSolver.o $(ODE_SOLVER_PATH)/ODECSolver.cpp $(INCLUDES)

ODESimulation.o: ODECSolver.o $(ODE_SOLVER_PATH)/ODETemplate.cpp
	$(CC) -c -o $(OBJ_DIR)/ODESimulation.o $(ODE_SOLVER_PATH)/ODETemplate.cpp $(INCLUDES)

ODESimulation: ODESimulation.o model.o ODECSolver.o $(SUNOBJ)
	$(CC) -o $(OUTPUT_FILE) $(SUNOBJ_PATHS) $(TEMPLATE_CPP) $(OBJ_DIR)/model.o $(OBJ_DIR)/ODESimulation.o $(OBJ_DIR)/ODECSolver.o $(INCLUDES)

clean:
	rm -rf $(OBJ_DIR)/*.o $(SUNDIALS_OBJ)/*.o $(OUTPUT_DIR)/*.out
