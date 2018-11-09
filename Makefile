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
	mpirun -np 13 ./tsp 15 16 500 500 	
	# ./tsp 15 4 500 500