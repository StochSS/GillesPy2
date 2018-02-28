CC=g++
CFLAGS=-c -Wall
DEPS = model.h ssa.h
OBJ = model.o ssa.o

all: MichaelisMenten


%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)


MichaelisMenten: $(OBJ)
	$(CC) MichaelisMenten.cpp -o $@ $^


clean:
	rm -f *.o *~ MichaelisMenten
