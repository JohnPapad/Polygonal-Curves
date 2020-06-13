CC=g++
CXXFLAGS= -Wall -g -std=c++11

curves_apps: main.o distances.o curve.o cluster.o clustering.o protein.o roads.o
	$(CC) $(CXXFLAGS) main.o distances.o curve.o cluster.o clustering.o protein.o roads.o -o curves_apps

main.o: main.cpp header.h cluster.h
	$(CC) $(CXXFLAGS) -c main.cpp
	
distances.o: distances.cpp header.h
	$(CC) $(CXXFLAGS) -c -I ./ distances.cpp

curve.o: curve.cpp curve.h cluster.h
	$(CC) $(CXXFLAGS) -c curve.cpp
	
cluster.o: cluster.cpp cluster.h
	$(CC) $(CXXFLAGS) -c cluster.cpp
	
clustering.o: clustering.cpp header.h
	$(CC) $(CXXFLAGS) -c clustering.cpp
    
protein.o: protein.cpp header.h
	$(CC) $(CXXFLAGS) -c protein.cpp

roads.o: roads.cpp header.h
	$(CC) $(CXXFLAGS) -c roads.cpp


.PHONY: clean

clean:
	rm -f curves_apps main.o distances.o curve.o cluster.o clustering.o protein.o roads.o
