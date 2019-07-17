# PCS
C-code implementation of the Parallel Collision Search algorithm by van Oorschot and Wiener. It is adapted for one-collision and multi-collision search on elliptic curves.  

This implementation is a result of joint work with [Sorina Ionica](https://home.mis.u-picardie.fr/~ionica/) and [Gilles Dequen](https://home.mis.u-picardie.fr/~dequen/doku.php). See our paper: [Time-Memory Trade-offs for Parallel Collision Search Algorithms](https://eprint.iacr.org/2017/581.pdf).

Distinguished points are stored in a Packed Radix-Tree-List (PRTL) structure, introduced in our research paper. A traditional hash table is implemented as well, as the goal is to compare the different approaches. The hash function we used for this comparison is the ElfHash function, which is used in the UNIX ELF format for object files. This code is modeled in a way that makes it very easy to add and compare other storage alternatives (How-to coming soon). 

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
-f : choose an f-bit elliptic curve (default is 35)
-t : number of threads to use (default is the number of cores avaliable)
-n : number of runs with different random secret keys (default is 10)
-s : storage structure (PRTL - default or hash_unix)
-l : level of the absract radix tree (default is 7 - see paper on how to choose this optimally)
-d : number of trailling zero bits in a distinguished point (default is floor(f/4))
-c : number of collisions that need to be found (default is one - for solving the ECDLP)
```

### Setting the value of the __DATA_SIZE_IN_BITS__ constant for optimal memory use
This constant is set in the ```pcs_vect_bin.h``` file and should be equal to the maximum number of bits you need to store your data in the structure. For example, for the PCS we store the x-coordinate of the distinguished point and a coefficient 'a'. Don't forget to subtract the trailling zero bits and the used prefix (which is equal to l). If we solve on a f-bit curve and we use level l and d trailling zero bits, the number of bits we need is : f - d - l (for the x-coordinate) + f (for the a coefficient). All this makes more sense in the [paper](https://eprint.iacr.org/2017/581.pdf).

### Experimental results
The data from the experimental results is written in the ```results``` directory. A line in one of the ```*.all``` files corresponds to one run and states the configuration options used followed by the specific measurement. For example, a line in the ```time.all``` file has the form

``` f s t d l time_measured ```.

The script ```refresh_avg.sh``` computes average values for each existing configuration and stores them in corresponding ```*.avg``` files.
