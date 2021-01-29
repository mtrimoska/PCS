# PCS
This is a C-code implementation of the Parallel Collision Search algorithm by van Oorschot and Wiener. It is adapted for one-collision and multi-collision search on elliptic curves.  

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
