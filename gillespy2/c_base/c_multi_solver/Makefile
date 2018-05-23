CC=gcc
CFLAGS= -O -Wall
LFLAGS= -pthread

all: c_multi_solver

c_multi_solver: c_multi_solver.c linked_list.c linked_list.h
	$(CC) $(LFLAGS) -c -o linked_list.o linked_list.c
	$(CC) $(LFLAGS) -c -o c_multi_solver.o c_multi_solver.c
	$(CC) $(LFLAGS) -o c_multi_solver c_multi_solver.o linked_list.o

test: 
	./c_multi_solver ../MichaelisMenten 4 28 101
sim:
	./c_multi_solver ../UserSimulation 4 4 100
	
clean:
	rm -f c_multi_solver c_multi_solver.exe c_multi_solver.o linked_list.o
