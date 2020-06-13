CC=g++
CXXFLAGS= -Wall -g -std=c++11

cluster: main.o lsh.o curve.o cluster.o multi_hash_table.o bst.o
	$(CC) $(CXXFLAGS) main.o lsh.o curve.o cluster.o multi_hash_table.o bst.o -o cluster

main.o: main.cpp lsh.h bst.h
	$(CC) $(CXXFLAGS) -c main.cpp
	
lsh.o: lsh.cpp lsh.h 
	$(CC) $(CXXFLAGS) -c lsh.cpp

curve.o: curve.cpp curve.h cluster.h
	$(CC) $(CXXFLAGS) -c curve.cpp
	
cluster.o: cluster.cpp cluster.h
	$(CC) $(CXXFLAGS) -c cluster.cpp

multi_hash_table.o: multi_hash_table.cpp multi_hash_table.h 
	$(CC) $(CXXFLAGS) -c multi_hash_table.cpp
	
bst.o: bst.cpp bst.h lsh.h
	$(CC) $(CXXFLAGS) -c bst.cpp

.PHONY: clean

clean:
	rm -f cluster *.o
