CC=g++
CXXFLAGS= -Wall -g -std=c++11

lsh: lsh.o curve.o multi_hash_table.o
	$(CC) $(CXXFLAGS) lsh.o curve.o multi_hash_table.o -o lsh

lsh.o: lsh.cpp lsh.h 
	$(CC) $(CXXFLAGS) -c lsh.cpp

curve.o: curve.cpp curve.h
	$(CC) $(CXXFLAGS) -c curve.cpp

multi_hash_table.o: multi_hash_table.cpp multi_hash_table.h 
	$(CC) $(CXXFLAGS) -c multi_hash_table.cpp

.PHONY: clean

clean:
	rm -f lsh *.o
