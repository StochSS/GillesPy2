CC=gcc
CFLAGS= -O -Wall
LFLAGS= -pthread

all: c_multi_solver

c_multi_solver: c_multi_solver.c
	$(CC) $(LFLAGS) -o c_multi_solver c_multi_solver.c

test: 
	./c_multi_solver ../MichaelisMenten 4
sim:
	./c_multi_solver ../UserSimulation 4
	
clean:
	rm -f c_multi_solver
