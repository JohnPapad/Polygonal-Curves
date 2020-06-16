# Polygonal Curves Project

A collection of libraries implementing Locality Sensitive Hashing (LSH), Clustering, and Applications of it;  
 developed for the class of ***"Software Development for Algorithmic Problems"*** in the Informatics Department.  
The whole project was implemented with **C++** using **STL** (following the C++11 standard) and **Object Oriented Programming** practices.


### Table of Contents

[Part 1: Locality Sensitive Hashing (LSH)](#part1)

[Part 2: Clustering](#part2)

[Part 3: Applications](#part3)

<a name="part1"/>

## Part 1: Locality Sensitive Hashing (LSH)

An implementation of LSH  with classic hash function for polygonal curves.

### Description

LSH is implemented for the purpose of finding the approximate nearest neighboring curve *q* from a set of polygonal curves.  

A polygonal curve is a sequence of vertices in the euclidean space &#8476;<sup>2</sup>. Each curve can have a different number of vertices.  

The distance metric used is the Discrete Frechet Distance.  

The LSH will map every polygonal curve to a grid curve which can be represented as a vector. This grid curve will then be hashed.

*L* hash tables will be created in total. Every input polygonal curve will be hashed *L* times (differently each time) and saved in a different hash table (in a hash table bucket).

A query curve will be hashed once and will search for all curves that have the same grid curve, in all hash tables. From those found, the nearest neighboring one will be found, by calculating the corresponding Discrete Frechet Distances.

In case that a search in the hash tables' stored grid curves is not successful, all curves from the corresponding buckets will be examined.

Apart from the approximate nearest neighbor, found by the LSH, the true nearest neighbor will be also found, through exhaustive search.


### Inputs

- A text tab-separated dataset with the following format:

    ```
    curve_id1   m1  (x11, y11)  (x12, y12) ... (x1m1, y1m1)
    .           .       .                  ...
    .           .       .                  ...
    .           .       .                  ...
    curve_idN   mN  (xN1, yN1)  (xN2, yN2) ... (xNmN, y1mN)
    ```

    where ```(xij, yij)``` are the coordinates (in double representation) of point ```j``` of curve ```i```, ```j``` &#8804; ```mi``` and ```mi``` the number of points of curve ```i```. The dataset's filename is provided as a command line argument with the *-d* flag.

- A text file that includes all the query curves. A query curve is a polygonal curve, for which a nearest neighboring curve search will be performed. This file also follows the same format as the previous input. The queries' filename is provided as a command line argument with the *-q* flag.
   

### Outputs

- A text file that for every query curve includes the approximate nearest neighboring curve ID, found by the LSH, the true nearest neighboring curve ID, found by exhaustive search and the corresponding distances, as well. The file has the following format:
    ```
    Query: curveJ
    DistanceFunction: DFT
    HashFunction: Classic
    FoundGridCurve: {True, False}
    LSH Nearest neighbor: curveY
    True Nearest neighbor: curveX
    distanceLSH: <double>
    distanceTrue: <double>
    ```
- In case that the command line *-stats* argument is given the program will be executed 100 times. Each time the query curve's approximate nearest neighbor will be found, by recreating the *L* hash tables, and the following additional statistics will be calculated and depicted in this format:
    ```
    Query: curveJ
    DistanceFunction: DFT
    HashFunction: Classic
    |minDistanceLSH – distanceTrue|: <double>
    |maxDistanceLSH – distanceTrue|: <double>
    |avgDistanceLSH – distanceTrue|: <double>
    tLSHmin: <double>
    tLSHmax: <double>
    tLSHavg: <double>
    tTrue: <double>
    ```

The output's filename is provided as a command line argument with the *-o* flag.


### Command line arguments
- ***-d*** : input dataset's filename *(if omitted will be asked from the standard input)*
- ***-q*** : query's filename *(if omitted will be asked from the standard input)*
- ***-o*** : output's filename *(if omitted will be asked from the standard input)*
- ***-L*** : number of hash tables *(default value: 3)*
- ***-k*** : number of locality-sensitive functions *hi* that define a grid curve *(default value: 2)*
- ***-stats*** : *as mentioned above*

### How to run

```
$ make (to compile)
$ ./lsh -d <input file> -q <query file> -k <int> -L <int> -ο <output file>
-stats [optional] 
$ make clean (to delete all .o files)
```

<a name="part2"/>

## Part 2: Clustering

An implementation of clustering using various combinations of methods, such as K-medoids, PAM, LSH etc.

## Description

Clustering algorithms for polygonal curves are implemented using 8 combinations of the following methods:

**Initialization**
1. K-means++
2. Random selection of k-points (*simplest*)

**Assignment**
1. Lloyd's assignment (*simplest approach*)
2. Assignment by Range search (*LSH*)

**Update**
1. Calculate the Mean Discrete Frechet curve
2. Partitioning Around Medoids (*PAM*) - Improved Update


### Inputs

- A text tab-separated dataset with the following format:

    ```
    curve_id1   m1  (x11, y11)  (x12, y12) ... (x1m1, y1m1)
    .           .       .                  ...
    .           .       .                  ...
    .           .       .                  ...
    curve_idN   mN  (xN1, yN1)  (xN2, yN2) ... (xNmN, y1mN)
    ```

    where ```(xij, yij)``` are the coordinates (in double representation) of point ```j``` of curve ```i```, ```j``` &#8804; ```mi``` and ```mi``` the number of points of curve ```i```. The dataset's filename is provided as a command line argument with the *-i* flag.

- A configuration file for the various parameters with the following format:
    ```
    number_of_clusters:<int>    // k 
    number_of_grid_curves:<int> // if omitted: default value = 2
    number_of_hash_tables:<int> // if omitted: default value L = 3
    ```

    The configuration's filename is provided as a command line argument with the *-c* flag.

### Outputs

A text file that includes the produced clusters from every combination. The *Silhouette* is an eternal evaluation metric.  
The output file will have the following format:

```
    Algorithm: ΙxAxUx
    Metric: Frechet
    CLUSTER-1 {size: <int>, centroid: <curve_id> or an array with the centroid's vertices in the case of the Mean Discrete Frechet curve update}
        .               .
        .               .
        .               .
    CLUSTER-k {size: <int>, centroid: <curve_id>}
    clustering_time: <double> // in seconds
    Silhouette: [s1, ..., si, ..., sk, stotal]
    /* si=average s(p) of points in cluster i, stotal=average s(p) of points in dataset */
    
    /* Optionally with command line parameter -complete */
    CLUSTER-1 {curve_idA, curve_idB, ..., curve_idC}
    .                   .
    .                   .
    .                   .
    CLUSTER-k {curve_idR, curve_idT, ..., curve_idZ}
```
The output's filename is provided as a command line argument with the *-o* flag.


### Command line arguments

- ***-i*** : input dataset's filename *(if omitted will be asked from the standard input)*
- ***-c*** : configuration's filename *(if omitted will be asked from the standard input)*
- ***-o*** : output's filename *(if omitted will be asked from the standard input)*
- ***-complete*** : *see output's format above*


### How to run

```
$ make (to compile)
$ ./cluster -i <input file> -c <configuration file> -o <output file> 
$ make clean (to delete all .o files)
```

<a name="part3"/>

## Part 3: Applications

Part 3 is 2 applications of the above in molecular conformations and roads of Athens, Greece with data taken from OpenStreetMap.

### Application 1: Molecular Conformations Clustering

### Description

A set of molecular conformations is given as input. Each conformation consists of a specific sequence of *N* points in the three-dimensional Euclidean space.

The molecular conformations will be clustered in *k* clusters, using *K-means++* for initialization, *Lloyd’s* method for assignment, *Partitioning Around Medoids (PAM) - Improved* for update and the *c-RMSD* function as the distance function.

There will also be another clustering using the *Discrete Frechét* distance function. In this case, before calculating the distance, the conformations are shifted and rotated based on minimizing the *c-RMSD* distance.

The *c-RMSD* distance function is implemented, using the external linear algebra library *Eigen*. The necessary files for its operation have been included.

Furthermore, the optimal *k* (number of cluster) is calculated through the *silhouette* evaluation measure.

### Inputs

The input is an *<input.dat>* text file. The first two lines contain the number of conformations and variable *N*. The rest of the file is separated by tabs in 3 columns that
contain the *N* triads of coordinates: ```x y z``` of the 1st conformation, then *N* triads of the 2nd conformation and so on, i.e. ```numConform * N``` lines, in total, after the first 2:

```
numConform: <Integer>
N: <Integer>
-32.5   91.2    11.7
12.8    -18.3   79.1
...
```

### Outputs

- A *<crmsd.dat>* text file, which depicts the clustering with the optimal *silhouette* value, separated by tabs with the following format:

    ```
    k: <Integer>
    s: <real in [-1,1]>
    1   9   11  12 ...
    2   3   17 [the 2nd cluster contains 3 elements]
    ```
    where the 1st line contains the *k*, the 2nd the *silhouette value*, and the next *k* lines contain the indicators of the corresponding cluster elements in ascending order.

- A *<frechet.dat>* text file with the same format as above.

### Application 2: Road Segments Clustering

### Description

A number of road segments from Athens, Greece (extracted from OpenStreetMap) are given as input.

The road segments will be clustered in *k* clusters, using *K-means++* for initialization, *Lloyd’s* method for assignment, *Partitioning Around Medoids (PAM) - Improved* for update and the *Discrete Frechét* distance metric. Before calculating the distance, the road segments are shifted and rotated based on minimizing the *c-RMSD* distance.

Furthermore, the optimal *k* (number of cluster) is calculated through the *silhouette* evaluation measure.

### Inputs

The input is a csv file with the following format:
```
nο, wayID, num_of_segs, node1_lat, node1_lon, . . ., last_node_lat, last_node_lon
0, 15918864, 2, 40.6157885, 22.9575427, . . ., 40.6163788, 22.9574743
... ... ... ... ... ... ... ... ... ... ...
```

### Outputs

The output is a *<kmeans_ways_frechet.dat>* text file, which depicts the clustering with the optimal *silhouette* value,
separated by tabs with the following format:
```
k: <Integer>
s: <real in [-1,1]>
clustering_time: <number of msec>
1   9   11  12 ...
2   3   17  [the 2nd cluster contains 3 elements]
```
where the 1st line contains the *k*, the 2nd the *silhouette* value, and the next *k* lines contain the indicators *(no)* of the road segments, belonging to the respective clusters in
ascending order.

### Command line arguments

- ***-i1*** : input filename for the 1st application *(if omitted will be asked from the standard input)*
- ***-i2*** : input filename for the 2nd application *(if omitted will be asked from the standard input)*
- ***-d*** : distance function option. Can be either *DFD* or *c-RMSD* *(if omitted c-RMSD will be chosen as default)*

### How to run
The 2 applications will be serially executed.  
(First the *molecular conformations* and then the *road segments*)
```
$ make (to compile)
$ ./curves_apps -i1 <input file1> -i2 <input file2> -d {DFT, CRMSD} 
$ make clean (to delete all .o files)
```
