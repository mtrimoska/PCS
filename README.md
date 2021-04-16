# PCS
This is a C-code implementation of the Parallel Collision Search algorithm by van Oorschot and Wiener. It is adapted for one-collision and multi-collision search on elliptic curves. 

This implementation is a result of joint work with [Sorina Ionica](https://home.mis.u-picardie.fr/~ionica/) and [Gilles Dequen](https://home.mis.u-picardie.fr/~dequen/doku.php). See our paper: [Time-Memory Analysis for Parallel Collision Search Algorithms](https://eprint.iacr.org/2017/581.pdf).

Distinguished points are stored in a Packed Radix-Tree-List (PRTL) structure, introduced in our research paper. A traditional hash table is implemented as well, as the goal is to compare the different approaches. The hash function we used for this comparison is the ElfHash function, which is used in the UNIX ELF format for object files. This code is modeled in a way that makes it very easy to add and compare other storage alternatives. 

### Dependencies
GMP - The GNU Multiple Precision Arithmetic Library

OpenMP for parallel programming

### Compiling in Linux and Mac OS
```bash
mkdir build && cd build
cmake ..
make
```
An executable ```pcs_exec``` will be created in the project's home directory.

### Command-line arguments
By default, the program solves the ECDLP for a random point P and a random secret key x. The code has several configuration options:
```
-f : choose an f-bit elliptic curve (default is 35 and currently possible values are 5k with k=7,...,23)
-t : number of threads to use (default is the number of cores avaliable)
-n : number of runs with different random secret keys (default is 10)
-s : storage structure (PRTL - default or hash_unix)
-l : level of the absract radix tree (default is 7 - see paper on how to choose this optimally)
-d : number of trailling zero bits in a distinguished point (default is floor(f/4))
-c : number of collisions that need to be found (default is one - for solving the ECDLP)
```

### Setting the value of the __DATA_SIZE_IN_BITS__ constant for optimal memory use
The PRTL structure stores all relevant data for one entry in one byte-vector. Since byte-vectors are statically allocated, we use a constant __DATA_SIZE_IN_BITS__ to define the size of byte-vectors. For optimal memory use, this constant should be set to the minimum required for a specific attack. The constant is set in the ```pcs_vect_bin.h``` file and should be equal to the maximum number of bits you need to store your data in the structure, which can be calculated as per the parameters used for your attack. For example, for the PCS we store the x-coordinate of the distinguished point and a coefficient 'a'. Don't forget to subtract the trailling zero bits and the used prefix (which is equal to l). If we solve on an f-bit curve and we use level l and d trailling zero bits, the number of bits we need is : f - d - l (for the x-coordinate) + f (for the a coefficient). There is also a __DATA_SIZE_IN_BYTES__ constant that is calculated as \ceil{__DATA_SIZE_IN_BITS__/8}. Depending on the compilator, you may need to set this constant manually as well. 

### Experimental results
The data from the experimental results is written in the ```results``` directory. A line in one of the ```*.all``` files corresponds to one run and states the configuration options used followed by the specific measurement. For example, a line in the ```time.all``` file has the form

``` f s t d l time_measured ```.

The script ```refresh_avg.sh``` computes average values for each existing configuration and stores them in corresponding ```*.avg``` files.

### Organization of the source code
The main execution file of the source code is ```pcs_exec.c```. It also contains code for management of experimental results. The following is a brief description of the other files:

```pcs_elliptic_curve_operations.c``` - Functions for initializing the Point and Curve structures and performing elliptic curve operations.

```pcs_pollard_rho.c``` - Computing the random walk function and the classical Pollard's rho algorithm.

```pcs_storage.c``` - Modelization of the storage functionality. Transfers calls of all data-management related functions to the respective data structures.

```pcs_struct_PRTL.c``` - Implementation of the PRTL structure.

```pcs_struct_hash.c``` - Implementation of a hash table.

```pcs_struct_hash_UNIX.c``` - Computing the ElfHash hash table function.

```pcs_vect_bin.c``` - A byte-vector implementation used for the 'packed' property of the PRTL structure.

```pcs.c``` - Functions relative to the Parallel Collision Search algorithm. 

### Adding other data structures for storing points
To add an implementation of a new data structure you need to create a new C file and its corresponding header. For consistency, you can name the files ```pcs_struct_XX.c``` and ```pcs_struct_XX.h```, replacing XX with the name of your structure. Then, include ```pcs_struct_XX.h``` in ```pcs_storage.c```. Your structure needs to implement four required functions:

```struct_init_XX``` - initializes the distinguished-point-storing structure.

```struct_add_XX``` - looks for a point in the structure. If the point is not found it is added with the corresponding a coefficient.

```struct_free_XX``` - frees the distinguished-point-storing structure.

```struct_memory_XX``` - gets the memory occupation of the distinguished-point-storing structure. This is required if you need to make experimental comparisons on memory use between the different structures. 

To know which parameters are available for each of these functions, see the Doxygen comments in ```pcs_storage.c```.

If you need a hash table, you can use our classical implementation of a hash table with another function. In this case, you just need to create a ```pcs_struct_hash_XX.c``` file implementing the ```get_hash_XX```function that computes the hash value.

Secondly, you need to add a ```case``` for your structure in all four functions in ```pcs_storage.c```.

Finally, in the file ```pcs_exec.c```, on line 102 of the current version, you need to add a key-word for your structure in the list of available structures, making sure that the position of your structure in the list corresponds to the ```case```number that you chose in the previous step. For instance, if you add a structre called ```binary-tree```, and you used ```case: 2``` in the ```switch``` in ```pcs_storage.c```, line 102 should look as follows:

```char *struct_i_str[] = {"PRTL", "hash_unix", "binary-tree"};```. 

To use your structure, you need to execute the program using the parameter ```-s binary-tree```. 


### Adding curves and points
The elliptic curves that are currently available for experiments are defined over \mathbb{F}_p, with p prime. There are f-bit curves for f=35,40,45,...,115. There are 10 points available for each curve, each of order equal to the cardinality of the group of points on the curve. To add curves and points without modifying the source code, the following rules need to be respected. Each line in the 'curves' file corresponds to one curve starting from a 35-bit field and growing in 5-bit increments. A line is composed of 5 arguments separated by spaces. For an f-bit curve E: y^2 = x^3 + Ax + B, defined over \mathbb{F}_p and of cardinality n, the arguments are as follows
``` f A B p n ```.
Each line is 83 characters long (add spaces to complete). Similarly, points are stored in the 'points' file. Each line in this file is 79 characters long and hold one point P(x,y) with two arguments
``` x y ```.

### Parameter choices for reproducing the results of the paper
The following are some examples of executions with appropriate command-line arguments, which correspond to the experiments performed for our paper. Command-line arguments are written as \[optional\] when the default value is the same as the specified value. Our experiments were performed on a 28-core processor and running times, as well as the default value of the -t parameter, may vary on different machines. 

* Figure 1.
  * The __DATA_SIZE_IN_BYTES__ constant does not need to be specified, as these experiments use only the hash finction structure.
  * Execution commands  
     ./pcs_exec -f 65 -s hash_unix -n 100 -d 6 \[-t 28 -c 1\]  
     ./pcs_exec -f 65 -s hash_unix -n 100 -d 8 \[-t 28 -c 1\]  
     ./pcs_exec -f 65 -s hash_unix -n 100 -d 10 \[-t 28 -c 1\]  
     ./pcs_exec -f 65 -s hash_unix -n 100 -d 12 \[-t 28 -c 1\]  
     ./pcs_exec -f 65 -s hash_unix -n 100 -d 14 \[-t 28 -c 1\]  
     ./pcs_exec -f 65 -s hash_unix -n 100 -d 16 \[-t 28 -c 1\]  
     ./pcs_exec -f 65 -s hash_unix -n 100 -d 18 \[-t 28 -c 1\]  
     ./pcs_exec -f 65 -s hash_unix -n 100 -d 20 \[-t 28 -c 1\]  
     ./pcs_exec -f 65 -s hash_unix -n 100 -d 22 \[-t 28 -c 1\]  
     ./pcs_exec -f 65 -s hash_unix -n 100 -d 24 \[-t 28 -c 1\]  
     ./pcs_exec -f 65 -s hash_unix -n 100 -d 26 \[-t 28 -c 1\]  
  * The results/time.avg file would contain:  
      65 hash_unix 28 6 7 :2963275486: (100 tests)  
      65 hash_unix 28 8 7 :1741456996: (100 tests)  
      65 hash_unix 28 10 7 :1361757547: (100 tests)  
      65 hash_unix 28 12 7 :1321082279: (100 tests)  
      65 hash_unix 28 14 7 :1226721996: (100 tests)  
      65 hash_unix 28 16 7 :1272840757: (100 tests)  
      65 hash_unix 28 18 7 :1358805440: (100 tests)  
      65 hash_unix 28 20 7 :1315671062: (100 tests)  
      65 hash_unix 28 22 7 :1385919515: (100 tests)  
      65 hash_unix 28 24 7 :1808912024: (100 tests)  
      65 hash_unix 28 26 7 :3009994583: (100 tests)  
   * The results/memory.avg file would contain:  
      65 hash_unix 28 6 7 :5236174321: (100 tests)  
      65 hash_unix 28 8 7 :1344970051: (100 tests)  
      65 hash_unix 28 10 7 :330122598: (100 tests)  
      65 hash_unix 28 12 7 :84414141: (100 tests)  
      65 hash_unix 28 14 7 :20087148: (100 tests)  
      65 hash_unix 28 16 7 :5129036: (100 tests)  
      65 hash_unix 28 18 7 :1350889: (100 tests)  
      65 hash_unix 28 20 7 :325659: (100 tests)  
      65 hash_unix 28 22 7 :84388: (100 tests)  
      65 hash_unix 28 24 7 :24932: (100 tests)  
      65 hash_unix 28 26 7 :8543: (100 tests)  
   * The results/points.avg file would contain:  
      65 hash_unix 28 6 7 :81741548: (100 tests)  
      65 hash_unix 28 8 7 :21462189: (100 tests)  
      65 hash_unix 28 10 7 :5266027: (100 tests)  
      65 hash_unix 28 12 7 :1375345: (100 tests)  
      65 hash_unix 28 14 7 :324774: (100 tests)  
      65 hash_unix 28 16 7 :84714: (100 tests)  
      65 hash_unix 28 18 7 :22674: (100 tests)  
      65 hash_unix 28 20 7 :5511: (100 tests)  
      65 hash_unix 28 22 7 :1446: (100 tests)  
      65 hash_unix 28 24 7 :449: (100 tests)  
      65 hash_unix 28 26 7 :162: (100 tests)  
   * The results/rate.avg file would contain:  
      65 hash_unix 28 6 7 :85.67 (50.3): (100 tests)  
      65 hash_unix 28 8 7 :87.26 (56.33): (100 tests)  
      65 hash_unix 28 10 7 :87.27 (57.32): (100 tests)  
      65 hash_unix 28 12 7 :88.95 (60.12): (100 tests)  
      65 hash_unix 28 14 7 :87.5 (56.98): (100 tests)  
      65 hash_unix 28 16 7 :88.1 (59.5): (100 tests)  
      65 hash_unix 28 18 7 :90.36 (62.55): (100 tests)  
      65 hash_unix 28 20 7 :89.28 (61.7): (100 tests)  
      65 hash_unix 28 22 7 :90.55 (63.28): (100 tests)  
      65 hash_unix 28 24 7 :91.6 (51.77): (100 tests)  
      65 hash_unix 28 26 7 :94.38 (57.10): (100 tests)  
* Table 1.
  * The __DATA_SIZE_IN_BYTES__ is different for each execution and specified after the corresponding command.
  * Execution commands  
    ./pcs_exec -f 55 -s PRTL -n 100 -d 13 -c 100 -l 15 \[-t 28\] (__DATA_SIZE_IN_BYTES__ should be set to 11)  
    ./pcs_exec -f 55 -s PRTL -n 100 -d 13 -c 500 -l 16 \[-t 28\] (__DATA_SIZE_IN_BYTES__ should be set to 11)   
    ./pcs_exec -f 55 -s PRTL -n 100 -d 13 -c 1000 -l 16 \[-t 28\] (__DATA_SIZE_IN_BYTES__ should be set to 11)  
    ./pcs_exec -f 55 -s PRTL -n 100 -d 13 -c 2000 -l 17 \[-t 28\] (__DATA_SIZE_IN_BYTES__ should be set to 10)    
    ./pcs_exec -f 55 -s PRTL -n 100 -d 13 -c 5000 -l 17 \[-t 28\] (__DATA_SIZE_IN_BYTES__ should be set to 10)    
    ./pcs_exec -f 55 -s PRTL -n 100 -d 13 -c 7000 -l 18 \[-t 28\] (__DATA_SIZE_IN_BYTES__ should be set to 10)  
* Experiments in Tables 2, 3 and 4 consist in adding random points to the PRTL and hash table structures, instead of performing an ECDLP attack, so the code was modified perform these experiments. 
* Table 5.
  * Execution commands  
    ./pcs_exec -f 55 -s PRTL -n 100 -l 12 \[-d 13 -t 28 -c 1\] (__DATA_SIZE_IN_BYTES__ should be set to 11)  
    ./pcs_exec -f 55 -s hash_unix -n 100 \[-d 13 -t 28 -c 1\]  
    ./pcs_exec -f 60 -s PRTL -n 100 -l 12 \[-d 15 -t 28 -c 1\] (__DATA_SIZE_IN_BYTES__ should be set to 12)  
    ./pcs_exec -f 60 -s hash_unix -n 100 \[-d 15 -t 28 -c 1\]  
    ./pcs_exec -f 65 -s PRTL -n 100 -l 14 \[-d 16 -t 28 -c 1\] (__DATA_SIZE_IN_BYTES__ should be set to 13)  
    ./pcs_exec -f 65 -s hash_unix -n 100 \[-d 16 -t 28 -c 1\]  
* Table 6.
  * The __DATA_SIZE_IN_BYTES__ constant should be set to 12
  * Execution commands  
    ./pcs_exec -f 60 -s PRTL -n 100 -l 12 -t 1 \[-d 15 -c 1\]  
    ./pcs_exec -f 60 -s PRTL -n 100 -l 12 -t 2 \[-d 15 -c 1\]  
    ./pcs_exec -f 60 -s PRTL -n 100 -l 12 -t 7 \[-d 15 -c 1\]  
    ./pcs_exec -f 60 -s PRTL -n 100 -l 12 -t 14 \[-d 15 -c 1\]  
    ./pcs_exec -f 60 -s PRTL -n 100 -l 12 \[-t 28 -d 15 -c 1\]  
    
    
