### Compiler configuration/settings ###
CXX := g++
CFLAGS ?= -Wall -O3
CXXFLAGS ?= -std=c++14 -Wall -O3
GPY_SOLVER ?= ssa
GPY_EXE_ARGS ?= --trajectories 1 --timesteps 21 --end 20
#######################################
.PHONY: run clean build profile ode ssa tau_leap hybrid

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
INCLUDES := -I$(CBASE_DIR) -I$(SUNDIALS_INC) -I$(TEMPLATE_DIR) -I$(TAU_DIR)
SUNDIALS_OBJ ?= $(OBJ_DIR)
OUTPUT_DIR ?= $(CBASE_DIR)
ifeq ($(OS),Windows_NT)
	OUTPUT_FILE ?= $(OUTPUT_DIR)/Simulation.exe
else
	OUTPUT_FILE ?= $(OUTPUT_DIR)/Simulation.out
endif
##########################

SUNOBJ = cvode_nls.o cvode_io.o sundials_iterative.o cvode_proj.o sundials_matrix.o sunmatrix_band.o sunmatrix_dense.o cvode_ls.o \
sundials_linearsolver.o sundials_nonlinearsolver.o sundials_nvector_senswrapper.o sunnonlinsol_newton.o \
sundials_nvector.o nvector_serial.o cvode.o cvode_spils.o sundials_math.o sunlinsol_spgmr.o
SUNOBJ_PATHS := $(SUNOBJ:%.o=$(SUNDIALS_OBJ)/%.o)

###################################
### SOLVER DEPENDENCIES COMPILE ###

### DEPENDENCIES FOR ALL SOLVERS ###
GPY_SRC = model.cpp arg_parser.cpp
GPY_OBJ := $(GPY_SRC:%.cpp=$(OBJ_DIR)/%.o)

### DEPENDENCIES FOR TAU-LEAPING BASED SOLVERS ###
GPY_SRC_TAUBASE = tau.cpp
GPY_OBJ_TAUBASE := $(GPY_SRC_TAUBASE:%.cpp=$(OBJ_DIR)/%.o)

###    DEPENDENCIES FOR ODE SOLVER    ###
GPY_SRC_ODE = ODESimulation.cpp ODESolver.cpp
GPY_OBJ_ODE := $(GPY_SRC_ODE:%.cpp=$(OBJ_DIR)/%.o)

### DEPENDENCIES FOR HYBRID SOLVER ###
GPY_SRC_HYBRID = TauHybridSimulation.cpp TauHybridSolver.cpp HybridModel.cpp \
hybrid_template.cpp integrator.cpp
GPY_OBJ_HYBRID := $(GPY_SRC_HYBRID:%.cpp=$(OBJ_DIR)/%.o)

### DEPENDENCIES FOR SSA SOLVER ###
GPY_SRC_SSA = SSASimulation.cpp SSASolver.cpp
GPY_OBJ_SSA := $(GPY_SRC_SSA:%.cpp=$(OBJ_DIR)/%.o)

### DEPENDENCIES FOR TAU LEAPING SOLVER ###
GPY_SRC_TAU = TauLeapingSimulation.cpp TauLeapingSolver.cpp
GPY_OBJ_TAU := $(GPY_SRC_TAU:%.cpp=$(OBJ_DIR)/%.o)

$(GPY_OBJ): $(OBJ_DIR)/%.o: $(CBASE_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^ $(INCLUDES)

$(GPY_OBJ_HYBRID): $(OBJ_DIR)/%.o: $(TAU_HYBRID_SOLVER_PATH)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^ $(INCLUDES)

$(GPY_OBJ_ODE): $(OBJ_DIR)/%.o: $(ODE_SOLVER_PATH)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^ $(INCLUDES)

$(GPY_OBJ_SSA): $(OBJ_DIR)/%.o: $(SSA_SOLVER_PATH)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^ $(INCLUDES)

$(GPY_OBJ_TAU): $(OBJ_DIR)/%.o: $(TAU_LEAPING_SOLVER_PATH)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^ $(INCLUDES)

$(GPY_OBJ_TAUBASE): $(OBJ_DIR)/%.o: $(TAU_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^ $(INCLUDES)

$(SUNOBJ_PATHS): $(SUNDIALS_OBJ)/%.o: $(SUNDIALS_SRC)/%.c
	$(CXX) -c -o $@ $< $(CFLAGS) -I$(SUNDIALS_INC)
sundials: $(SUNOBJ_PATHS) ;

##########################
### FINAL COMPILATIONS ###
GPY_ALL_DEPS := $(GPY_OBJ) $(TEMPLATE_CPP)

ode: $(GPY_ALL_DEPS) $(GPY_OBJ_ODE) $(SUNOBJ_PATHS)
	$(CXX) $(CXXFLAGS) -o $(OUTPUT_FILE) $^ $(INCLUDES)

ssa: $(GPY_ALL_DEPS) $(GPY_OBJ_SSA)
	$(CXX) $(CXXFLAGS) -o $(OUTPUT_FILE) $^ $(INCLUDES)

tau_leap: $(GPY_ALL_DEPS) $(GPY_OBJ_TAU) $(GPY_OBJ_TAUBASE)
	$(CXX) $(CXXFLAGS) -o $(OUTPUT_FILE) $^ $(INCLUDES)

hybrid: $(GPY_ALL_DEPS) $(GPY_OBJ_HYBRID) $(GPY_OBJ_TAUBASE) $(SUNOBJ_PATHS)
	$(CXX) $(CXXFLAGS) -o $(OUTPUT_FILE) $^ $(INCLUDES)

build: $(GPY_ALL_DEPS) $(GPY_OBJ_TAUBASE) $(GPY_OBJ_SSA) $(GPY_OBJ_ODE) $(GPY_OBJ_TAU) $(GPY_OBJ_HYBRID)) ;

run: $(GPY_SOLVER)
	$(OUTPUT_FILE) $(GPY_EXE_ARGS)

debug: CXXFLAGS = -std=c++14 -Wall -g
debug: $(GPY_SOLVER)
	gdb --args $(OUTPUT_FILE) $(GPY_EXE_ARGS)

profile: CXXFLAGS = -std=c++14 -Wall -pg
profile: $(GPY_SOLVER)
	$(OUTPUT_FILE) $(GPY_EXE_ARGS)

clean:
	rm -rf $(OUTPUT_DIR)/*.out $(SUNDIALS_OBJ)/*.o $(OBJ_DIR)/*.o
