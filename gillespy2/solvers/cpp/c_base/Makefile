### Compiler configuration/settings ###
CXX := g++
CFLAGS := -Wall -O3
CXXFLAGS := -std=c++14 -Wall -O3
SOLVER ?= ssa
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
EXE_ARGS ?= --trajectories 1 --end 10 --timesteps 11 --increment 1
SUNDIALS_OBJ ?= $(OBJ_DIR)
INCLUDES := -I$(CBASE_DIR) -I$(SUNDIALS_INC) -I$(TEMPLATE_DIR) -I$(TAU_DIR)
##########################

SUNOBJ = cvode_nls.o cvode_io.o sundials_iterative.o cvode_proj.o sundials_matrix.o sunmatrix_band.o sunmatrix_dense.o cvode_ls.o \
sundials_linearsolver.o sundials_nonlinearsolver.o sundials_nvector_senswrapper.o sunnonlinsol_newton.o \
sundials_nvector.o nvector_serial.o cvode.o cvode_spils.o sundials_math.o sunlinsol_spgmr.o
SUNOBJ_PATHS := $(SUNOBJ:%.o=$(SUNDIALS_OBJ)/%.o)

TAU_DEPENDENCIES = tau.o
TAU_DEPS_PATHS := $(TAU_DEPENDENCIES:%.o=$(OBJ_DIR)/%.o)

TAU_HYBRID_DEPENDENCIES = HybridModel.o rhs.o integrator.o hybrid_template.o hybridutils.o
TAU_HYBRID_DEPS_PATHS := $(TAU_HYBRID_DEPENDENCIES:%.o=$(OBJ_DIR)/%.o)

###################################
### SOLVER DEPENDENCIES COMPILE ###
$(OBJ_DIR)/model.o: $(CBASE_DIR)/model.cpp $(CBASE_DIR)/model.h 
	$(CXX) $(CXXFLAGS) -c -o $(OBJ_DIR)/model.o $(CBASE_DIR)/model.cpp -I$(CBASE_DIR)
model: $(OBJ_DIR)/model.o ;

$(OBJ_DIR)/arg_parser.o: $(CBASE_DIR)/arg_parser.cpp $(CBASE_DIR)/arg_parser.h
	$(CXX) $(CXXFLAGS) -c -o $(OBJ_DIR)/arg_parser.o $(CBASE_DIR)/arg_parser.cpp $(INCLUDES)
arg_parser: $(OBJ_DIR)/arg_parser.o ;

$(SUNOBJ_PATHS): $(SUNDIALS_OBJ)/%.o: $(SUNDIALS_SRC)/%.c
	$(CXX) -c -o $@ $< $(CXXFLAGS) -I$(SUNDIALS_INC)
sundials: $(SUNOBJ_PATHS);

$(TAU_DEPS_PATHS): $(OBJ_DIR)/%.o: $(TAU_DIR)/%.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(INCLUDES)
tau_deps: $(TAU_DEPS_PATHS);

$(TAU_HYBRID_DEPS_PATHS): $(OBJ_DIR)/%.o: $(TAU_HYBRID_SOLVER_PATH)/%.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(INCLUDES)
hybrid_deps: $(TAU_HYBRID_DEPS_PATHS);

#################################
### SOLVER ALGORITHM COMPILE ####
$(OBJ_DIR)/ODESolver.o: $(ODE_SOLVER_PATH)/ODESolver.cpp 
	$(CXX) $(CXXFLAGS) -c -o $(OBJ_DIR)/ODESolver.o $(ODE_SOLVER_PATH)/ODESolver.cpp  $(INCLUDES)

$(OBJ_DIR)/SSASolver.o: $(SSA_SOLVER_PATH)/SSASolver.cpp 
	$(CXX) $(CXXFLAGS) -c -o $(OBJ_DIR)/SSASolver.o $(SSA_SOLVER_PATH)/SSASolver.cpp $(INCLUDES)

$(OBJ_DIR)/TauLeapingSolver.o: $(TAU_LEAPING_SOLVER_PATH)/TauLeapingSolver.cpp
	$(CXX) $(CXXFLAGS) -c -o $(OBJ_DIR)/TauLeapingSolver.o $(TAU_LEAPING_SOLVER_PATH)/TauLeapingSolver.cpp $(INCLUDES)

$(OBJ_DIR)/TauHybridSolver.o: $(TAU_HYBRID_SOLVER_PATH)/TauHybridSolver.cpp
	$(CXX) $(CXXFLAGS) -c -o $(OBJ_DIR)/TauHybridSolver.o $(TAU_HYBRID_SOLVER_PATH)/TauHybridSolver.cpp $(INCLUDES)

%Solver.o: $(OBJ_DIR)/%Solver.o ;

#################################
### SOLVER SIMULATION COMPILE ###

$(OBJ_DIR)/ODESimulation.o: $(ODE_SOLVER_PATH)/ODESimulation.cpp
	$(CXX) $(CXXFLAGS) -c -o $(OBJ_DIR)/ODESimulation.o $(ODE_SOLVER_PATH)/ODESimulation.cpp $(INCLUDES)

$(OBJ_DIR)/SSASimulation.o: $(SSA_SOLVER_PATH)/SSASimulation.cpp
	$(CXX) $(CXXFLAGS) -c -o $(OBJ_DIR)/SSASimulation.o $(SSA_SOLVER_PATH)/SSASimulation.cpp $(INCLUDES)

$(OBJ_DIR)/TauLeapingSimulation.o: $(TAU_LEAPING_SOLVER_PATH)/TauLeapingSimulation.cpp
	$(CXX) $(CXXFLAGS) -c -o $(OBJ_DIR)/TauLeapingSimulation.o $(TAU_LEAPING_SOLVER_PATH)/TauLeapingSimulation.cpp $(INCLUDES)

$(OBJ_DIR)/TauHybridSimulation.o: $(TAU_HYBRID_SOLVER_PATH)/TauHybridSimulation.cpp
	$(CXX) $(CXXFLAGS) -c -o $(OBJ_DIR)/TauHybridSimulation.o $(TAU_HYBRID_SOLVER_PATH)/TauHybridSimulation.cpp $(INCLUDES)

%Simulation.o: $(OBJ_DIR)/%Simulation.o ;

#########################
### PRE-COMPILE RULES ###
prebuild_solvers: ODESolver.o SSASolver.o TauLeapingSolver.o TauHybridSolver.o ;

prebuild_simulations: ODESimulation.o SSASimulation.o TauLeapingSimulation.o TauHybridSolver.o ;

prebuild: prebuild_solvers prebuild_simulations sundials arg_parser ;

#########################

##########################
### FINAL COMPILATIONS ###
COMPILATION_ARGS := $(CXXFLAGS) -o $(OUTPUT_FILE) $(TEMPLATE_CPP) $(OBJ_DIR)/model.o $(OBJ_DIR)/arg_parser.o $(INCLUDES)

ode: ODESimulation.o ODESolver.o sundials model arg_parser $(TEMPLATE_CPP)
	$(CXX) $(COMPILATION_ARGS) $(SUNOBJ_PATHS) $(OBJ_DIR)/ODESimulation.o $(OBJ_DIR)/ODESolver.o

ssa: SSASimulation.o SSASolver.o model arg_parser $(TEMPLATE_CPP)
	$(CXX) $(COMPILATION_ARGS) $(OBJ_DIR)/SSASimulation.o $(OBJ_DIR)/SSASolver.o

tau_leap: TauLeapingSimulation.o TauLeapingSolver.o tau_deps model arg_parser $(TEMPLATE_CPP)
	$(CXX) $(COMPILATION_ARGS) $(OBJ_DIR)/TauLeapingSimulation.o $(OBJ_DIR)/tau.o $(OBJ_DIR)/TauLeapingSolver.o

hybrid: TauHybridSimulation.o TauHybridSolver.o model arg_parser sundials tau_deps hybrid_deps $(TEMPLATE_CPP)
	$(CXX) $(COMPILATION_ARGS) $(SUNOBJ_PATHS) $(TAU_HYBRID_DEPS_PATHS) $(TAU_DEPS_PATHS) $(OBJ_DIR)/TauHybridSimulation.o  $(OBJ_DIR)/TauHybridSolver.o

build: $(SOLVER);

run: build
	./$(OUTPUT_FILE) $(EXE_ARGS)

debug: CXXFLAGS = -std=c++14 -Wall -g
debug: build;
	gdb --args $(OUTPUT_FILE) $(EXE_ARGS)

all: prebuild

clean:
	rm -rf $(OUTPUT_DIR)/*.out

clean_all:
	rm -rf $(OBJ_DIR)/*.o $(SUNDIALS_OBJ)/*.o $(OUTPUT_DIR)/*.out
