CC=gcc
CFLAGS= -O -Wall
LFLAGS=

all: c_multi_solver

c_multi_solver: c_multi_solver.c semaphore_lib.h semaphore_lib.c
	$(CC) $(CFLAGS) -c -o semaphore_lib.o semaphore_lib.c
	$(CC) $(CFLAGS) -c -o c_multi_solver.o c_multi_solver.c
	$(CC) $(LFLAGS) -o c_multi_solver c_multi_solver.o semaphore_lib.o

test: 
	./c_multi_solver ../MichaelisMenten 4
	
clean:
	rm -f c_multi_solver c_multi_solver.o semaphore_lib.o
