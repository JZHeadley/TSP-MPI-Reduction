INCLUDES=-I. 
CC=mpic++
CFLAGS=-g -Wall $(INCLUDES) -std=c++11
LINKARGS=-lm 
LIBS=-lm  
DEPS=assignment2.h

all: tsp

tsp: tsp.cpp $(DEPS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
	rm -f tsp

%.o: %.cpp $(DEPS)
	$(CC) -c $(CFLAGS) $< -o $@

run: all
	# mpirun -np 6 tsp 	
	./tsp 15 4 500 500