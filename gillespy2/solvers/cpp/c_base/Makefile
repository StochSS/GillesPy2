### Compiler configuration/settings ###
CC=g++
CXXFLAGS = -std=c++14 -Wall -O3
GPY_SOLVER ?= ssa
#######################################

### Input directories ###
CBASE_DIR ?= .
TEMPLATE_DIR ?= $(CBASE_DIR)/template
TEMPLATE_CPP := $(TEMPLATE_DIR)/template.cpp
SUNDIALS_SRC := $(CBASE_DIR)/Sundials/src
SUNDIALS_INC := $(CBASE_DIR)/Sundials/include
TAU_DIR := $(CBASE_DIR)/Tau
ODE_SOLVER_PATH := $(CBASE_DIR)/ode_cpp_solver
SSA_SOLVER_PATH := $(CBASE_DIR)/ssa_cpp_solver
TAU_LEAPING_SOLVER_PATH := $(CBASE_DIR)/tau_leaping_cpp_solver
TAU_HYBRID_SOLVER_PATH := $(CBASE_DIR)/tau_hybrid_cpp_solver
#########################

### Output directories ###
OBJ_DIR ?= $(CBASE_DIR)
OUTPUT_DIR ?= $(CBASE_DIR)
OUTPUT_FILE ?= $(OUTPUT_DIR)/Simulation.out
EXE_CMD ?= ./
EXE_ARGS ?= -trajectories 1 -end 10 -timesteps 11 -increment 1
SUNDIALS_OBJ ?= $(OBJ_DIR)
INCLUDES := -I$(CBASE_DIR) -I$(SUNDIALS_INC) -I$(TEMPLATE_DIR) -I$(TAU_DIR)
##########################

SUNOBJ = cvode_nls.o cvode_io.o sundials_iterative.o cvode_proj.o sundials_matrix.o sunmatrix_band.o sunmatrix_dense.o cvode_ls.o \
sundials_linearsolver.o sundials_nonlinearsolver.o sundials_nvector_senswrapper.o sunnonlinsol_newton.o \
sundials_nvector.o nvector_serial.o cvode.o cvode_spils.o sundials_math.o sunlinsol_spgmr.o
SUNOBJ_PATHS := $(SUNOBJ:%.o=$(SUNDIALS_OBJ)/%.o)

###################################
### SOLVER DEPENDENCIES COMPILE ###
$(OBJ_DIR)/model.o: $(CBASE_DIR)/model.cpp $(CBASE_DIR)/model.h
	$(CC) $(CXXFLAGS) -c -o $(OBJ_DIR)/model.o $(CBASE_DIR)/model.cpp -I$(CBASE_DIR)
model: $(OBJ_DIR)/model.o ;

$(SUNOBJ_PATHS): $(SUNDIALS_OBJ)/%.o: $(SUNDIALS_SRC)/%.c
	$(CC) -c -o $@ $< $(CXXFLAGS) -I$(SUNDIALS_INC)
sundials: $(SUNOBJ_PATHS) ;

Tau.o: $(TAU_DIR)/tau.cpp $(TAU_DIR)/tau.h
	$(CC) $(CXXFLAGS) -c -o Tau.o $(TAU_DIR)/tau.cpp $(INCLUDES)

#################################
### SOLVER ALGORITHM COMPILE ####
$(OBJ_DIR)/ODESolver.o: $(ODE_SOLVER_PATH)/ODESolver.cpp 
	$(CC) $(CXXFLAGS) -c -o $(OBJ_DIR)/ODESolver.o $(ODE_SOLVER_PATH)/ODESolver.cpp  $(INCLUDES)

$(OBJ_DIR)/SSASolver.o: $(SSA_SOLVER_PATH)/SSASolver.cpp
	$(CC) $(CXXFLAGS) -c -o $(OBJ_DIR)/SSASolver.o $(SSA_SOLVER_PATH)/SSASolver.cpp $(INCLUDES)

$(OBJ_DIR)/TauLeapingSolver.o: $(TAU_LEAPING_SOLVER_PATH)/TauLeapingSolver.cpp
	$(CC) $(CXXFLAGS) -c -o $(OBJ_DIR)/TauLeapingSolver.o $(TAU_LEAPING_SOLVER_PATH)/TauLeapingSolver.cpp $(INCLUDES)

$(OBJ_DIR)/TauHybridSolver.o: $(TAU_HYBRID_SOLVER_PATH)/TauHybridSolver.cpp
	$(CC) $(CXXFLAGS) -c -o $(OBJ_DIR)/TauHybridSolver.o $(TAU_HYBRID_SOLVER_PATH)/TauHybridSolver.cpp $(INCLUDES)

%Solver.o: $(OBJ_DIR)/%Solver.o ;

#################################
### SOLVER SIMULATION COMPILE ###

$(OBJ_DIR)/ODESimulation.o: $(ODE_SOLVER_PATH)/ODESimulation.cpp
	$(CC) $(CXXFLAGS) -c -o $(OBJ_DIR)/ODESimulation.o $(ODE_SOLVER_PATH)/ODESimulation.cpp $(INCLUDES)

$(OBJ_DIR)/SSASimulation.o: $(SSA_SOLVER_PATH)/SSASimulation.cpp
	$(CC) $(CXXFLAGS) -c -o $(OBJ_DIR)/SSASimulation.o $(SSA_SOLVER_PATH)/SSASimulation.cpp $(INCLUDES)

$(OBJ_DIR)/TauLeapingSimulation.o: $(TAU_LEAPING_SOLVER_PATH)/TauLeapingSimulation.cpp
	$(CC) $(CXXFLAGS) -c -o $(OBJ_DIR)/TauLeapingSimulation.o $(TAU_LEAPING_SOLVER_PATH)/TauLeapingSimulation.cpp $(INCLUDES)

$(OBJ_DIR)/TauHybridSimulation.o: $(TAU_HYBRID_SOLVER_PATH)/TauHybridSimulation.cpp
	$(CC) $(CXXFLAGS) -c -o $(OBJ_DIR)/TauHybridSimulation.o $(TAU_HYBRID_SOLVER_PATH)/TauHybridSimulation.cpp $(INCLUDES)

%Simulation.o: $(OBJ_DIR)/%Simulation.o ;

#########################
### PRE-COMPILE RULES ###
prebuild_solvers: ODESolver.o SSASolver.o TauLeapingSolver.o TauHybridSolver.o ;

prebuild_simulations: ODESimulation.o SSASimulation.o TauLeapingSimulation.o TauHybridSolver.o ;

prebuild: prebuild_solvers prebuild_simulations sundials ;

#########################

##########################
### FINAL COMPILATIONS ###
COMPILATION_ARGS = $(CXXFLAGS) -o $(OUTPUT_FILE) $(TEMPLATE_CPP) $(OBJ_DIR)/model.o $(INCLUDES)

ode: ODESimulation.o ODESolver.o sundials model $(TEMPLATE_CPP)
	$(CC) $(COMPILATION_ARGS) $(SUNOBJ_PATHS) $(OBJ_DIR)/ODESimulation.o $(OBJ_DIR)/ODESolver.o

ssa: SSASimulation.o SSASolver.o model $(TEMPLATE_CPP)
	$(CC) $(COMPILATION_ARGS) $(OBJ_DIR)/SSASimulation.o $(OBJ_DIR)/SSASolver.o

tau_leap: TauLeapingSimulation.o TauLeapingSolver.o model $(TEMPLATE_CPP)
	$(CC) $(COMPILATION_ARGS) $(OBJ_DIR)/TauLeapingSimulation.o $(OBJ_DIR)/TauLeapingSolver.o

hybrid: TauHybridSimulation.o TauHybridSolver.o model sundials $(TAU_DIR)/tau.cpp $(TEMPLATE_CPP)
	$(CC) $(COMPILATION_ARGS) $(SUNOBJ_PATHS) $(TAU_DIR)/tau.cpp $(OBJ_DIR)/TauHybridSimulation.o $(OBJ_DIR)/TauHybridSolver.o

build: $(GPY_SOLVER);
debug: CXXFLAGS = -std=c++14 -Wall -g
debug: clean_all build;

run: build
	./$(OUTPUT_FILE) $(EXE_ARGS)
gdb: debug
	gdb --args $(OUTPUT_FILE) $(EXE_ARGS)

clean:
	rm -rf $(OUTPUT_DIR)/*.out

clean_all:
	rm -rf $(OBJ_DIR)/*.o $(SUNDIALS_OBJ)/*.o $(OUTPUT_DIR)/*.out
