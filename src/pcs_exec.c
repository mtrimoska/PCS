/** @file pcs_exec.c
 *  @brief The main execution file. Contains code for experiments management. 
 *
 *	Created by Monika Trimoska on 03/12/2015.
 *	Copyright © 2015 Monika Trimoska. All rights reserved.
 */
#include <stdio.h>
#include <gmp.h>
#include <omp.h>
#include <inttypes.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "pcs_elliptic_curve_operations.h"
#include "pcs_pollard_rho.h"
#include "pcs_storage.h"
#include "pcs.h"
#include "pcs_vect_bin.h"

#define RESULTS_PATH "./results/"
#define __NB_STRUCTURES__ 2

/** Generates random number of EXACTLY nb_bits bits stored as an mpz_t type.
 * 	
 * 	@param[out]	s			Will hold the resulting number.
 * 	@param[in]	nb_bits		The number of bits (size of the number). 
 */
void generate_random_key(mpz_t s, int nb_bits)
{
	time_t t;
	mpz_t min;
	mpz_t max;
	mpz_t interval;
	unsigned long int seed;
	gmp_randstate_t r_state;
	srand((unsigned) time(&t));
	seed = rand();
	gmp_randinit_default(r_state);
	gmp_randseed_ui(r_state, seed);
	
	mpz_inits(min, max, interval, NULL);
	mpz_ui_pow_ui(min, 2, (nb_bits - 1));
	mpz_ui_pow_ui(max, 2, nb_bits);
	mpz_sub(interval, max, min);
	
	mpz_urandomm(s, r_state, max);
	if(mpz_cmp(s, min) < 0)
	{
		mpz_add(s, s, interval);
	}
	mpz_clears(min, max, interval, NULL);
	gmp_randclear(r_state);
}

/** Generates 20 sets for the adding walks.
 * 
 * 	@param[out]		A	The A coefficient set.
 * 	@param[out]		B	The B coefficient set.
 * 	@param[in]		max	The maximum value of a coefficient.
 */
void generate_adding_sets(mpz_t A[20], mpz_t B[20], mpz_t max)
{
	gmp_randstate_t r_state;
	gmp_randinit_default(r_state);
	gmp_randseed_ui(r_state, time(NULL));
	uint8_t i;
	for(i = 0; i < 20; i++)
	{
		mpz_urandomm(A[i], r_state, max);
		mpz_urandomm(B[i], r_state, max);
	}
	gmp_randclear(r_state);
}

/** Print out executable usage.
 */
void print_usage() {
    printf("Usage: \n-f : choose an f-bit elliptic curve (default is 35 and currently possible values are 5k with k=7,...,23)\n-t : number of threads to use (default is the number of cores avaliable)\n-n : number of runs with different random secret keys (default is 10)\n-s : storage structure (PRTL - default or hash_unix)\n-l : level of the absract radix tree (default is 7 - see paper on how to choose this optimally)\n-d : number of trailling zero bits in a distinguished point (default is floor(f/4))\n-c : number of collisions that need to be found (default is one - for solving the ECDLP)\n");
}

/**	Add a structure to the list of structures to be used.
 */
void add_to_struct_options(uint8_t structs[], char **struct_i_str, char *opt, uint8_t *struct_chosen)
{
	uint8_t i;
	for(i = 0; i < __NB_STRUCTURES__; i++)
	{
		if(strncmp(opt, struct_i_str[i], strlen(struct_i_str[i]) + 1) == 0)
		{
			structs[i] = 1;
			*struct_chosen = 1;
		}
	}
}

int main(int argc,char * argv[])
{	
	elliptic_curve_t E;
	char str_A[4], str_B[4], str_p[40], str_large_prime[40], str_X[40],str_Y[40];
	char *struct_i_str[] = {"PRTL", "hash_unix"};
	point_t P;
	point_t Q;
	mpz_t large_prime;
	mpz_t A[__NB_ENSEMBLES__];
	mpz_t B[__NB_ENSEMBLES__];
	FILE *file_res;
	FILE *file_curves;
	FILE *file_points;
	FILE *file_conf;
	char conf_str[1000] = {0};
	char conf_str_cpy[1000] = {0};
	char *conf_value;
	uint8_t update = 0;
	uint8_t struct_chosen = 0;
	struct timeval tv1;
	struct timeval tv2;
	unsigned long long int time, time1, time2;
	unsigned long long int memory;
	unsigned long int nb_points;
	float rate_of_use, rate_slots;
	mpz_t key;
	mpz_t x;
	char option;
	uint8_t nb_bits, trailling_bits, nb_curve, line_file_curves, line_file_points, nb_points_file, nb_point, j, struct_i, level;
	int test_i, nb_tests, nb_threads;
	int nb_collisions = 1;
	int trailling_bits_is_set = 0;
	int correct_data_size_in_bytes;
	uint8_t structs[__NB_STRUCTURES__] = {0};
    rate_slots = 0.0;
    level = 7;
	nb_bits = 35;
    trailling_bits = 0;
	nb_threads =  omp_get_max_threads();
	nb_tests = 10;
	line_file_curves = 84;
	line_file_points = 80;
	nb_points_file = 10;

	while ((option = getopt(argc, argv,"f:t:n:s:l:d:c:h")) != -1) {
        switch (option) {
			case 'f' : nb_bits = atoi(optarg);
				break;
			case 't' : nb_threads = atoi(optarg);
				break;
			case 'n' : nb_tests = atoi(optarg);
				break;
			case 's' : add_to_struct_options(structs, struct_i_str, optarg, &struct_chosen);
				break;
            case 'l' : level = atoi(optarg);
				break;
			case 'c' : nb_collisions = atoi(optarg);
				break;
			case 'd' : trailling_bits = atoi(optarg); trailling_bits_is_set = 1;
				break;
			case 'h' : {print_usage();exit(0);}
				break;
		}
	}
	
	/*** BEGIN: check input parameters boundary conditions */
	if(nb_bits < 35 || nb_bits > 115 || nb_bits%5 != 0)
	{
		fprintf(stdout, "\n********\n\033[0;31mWarning:\033[0m There is no example curve of %2" SCNu8 " bits available in the current implementation.\n", nb_bits);
		fprintf(stdout, "Available choices for the -f parameter are: 35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115.\n");
		nb_bits = (nb_bits / 5) * 5;
		if(nb_bits < 35) nb_bits = 35;
		if(nb_bits > 115) nb_bits = 115;
		fprintf(stdout, "\033[0;31mThis execution will use a %2" SCNu8 "-bit curve.\033[0m For a different choice, please restart the program.\n********\n\n", nb_bits);
		
	}
	
	if(nb_threads < 1 || nb_threads > 2000)
	{
		fprintf(stderr, "Can not use %d threads. Choose a value in the [1;2000] interval.\n", nb_threads);
		exit(1);
	}
	
	if(nb_threads > omp_get_max_threads())
	{
		fprintf(stdout, "\n********\n\033[0;31mWarning:\033[0m Using %d threads on this machine may result in slower running times. Recommended number of threads: %d. Continuing execution with %d...\n********\n\n", nb_threads, omp_get_max_threads(), nb_threads);
	}
    
	if(nb_tests < 1)
	{
		fprintf(stderr, "Invalid number of tests: %d.\n", nb_tests);
		exit(1);
	}
	
    if(trailling_bits_is_set == 0)
	   trailling_bits = nb_bits / 4;
	if(trailling_bits > nb_bits)
	{
		fprintf(stderr, "The number of trailling zero bits can not be greater than the number of bits of the x-coordinate.\n");
		exit(1);
	}
	
	if(!struct_chosen)
	{
		fprintf(stdout, "\n********\n\033[0;31mWarning:\033[0m There are no chosen storage structures. Adding PRTL by default.\n********\n\n");
		add_to_struct_options(structs, struct_i_str, "PRTL", &struct_chosen);
	}
	
	if(structs[0] == 1) //checks the boundary conditions of level l only if the PRTL structure is used
	{
		if(level < 0)
		{
			fprintf(stderr, "Invalid level: %d.\n", level);
			exit(1);
		}
		if(level > nb_bits - trailling_bits)
		{
			fprintf(stderr, "The level (prefix) can not be greater than the length of a stored word: %d.\n", nb_bits - trailling_bits);
			exit(1);
		}
	}
	
	if(nb_collisions < 1)
	{
		fprintf(stderr, "Invalid number of collisions: %d.\n", nb_collisions);
		exit(1);
	}
	
	/*** END: check input parameters boundary conditions */
	
	/***BEGIN: check if the __DATA_SIZE_IN_BYTES__ constant is properly set for the chosen parameters */
	if(structs[0] == 1) //these checks apply only if the PRTL structure is used
	{
		correct_data_size_in_bytes = 2*nb_bits - trailling_bits - level;
		correct_data_size_in_bytes = (((correct_data_size_in_bytes / (sizeof(_vect_bin_t) << 3)) + ((correct_data_size_in_bytes % (sizeof(_vect_bin_t) << 3)) ? 1 : 0)));
		if(__DATA_SIZE_IN_BYTES__ < correct_data_size_in_bytes)
		{
			fprintf(stderr, "Not enough memory per point allocated. The __DATA_SIZE_IN_BYTES__ constant in the pcs_vect_bin.h file has to be set to %d.\n", correct_data_size_in_bytes);
			exit(1);
		}
		if(__DATA_SIZE_IN_BYTES__ > correct_data_size_in_bytes)
		{
			fprintf(stdout, "\n********\n\033[0;31mWarning:\033[0m The choice of compile-time parameters for the PRTL structure is not optimal for these parameters and will result in higher memory requirements. \nIf these experiments assess the memory reqirements of the PRTL structure, the __DATA_SIZE_IN_BYTES__ constant in the pcs_vect_bin.h file should be set to %d.\n********\n\n", correct_data_size_in_bytes);
		}
	}
	
	/***END: check if the __DATA_SIZE_IN_BYTES__ constant is properly set for the chosen parameters */
	
	/*** set the number of threads for preallocation***/
	set_nb_threads(nb_threads);
	
	curve_init(&E);
	point_init(&P);
	point_init(&Q);
	mpz_inits(x,large_prime, key, NULL);
	for(j=0;j<__NB_ENSEMBLES__;j++)
	{
		mpz_inits(A[j],B[j],NULL);
	}

	nb_curve = nb_bits / 5 - 3;

	/*** read curve ***/
	file_curves = fopen("curves","r");
	if (file_curves == NULL) 
	{
		fprintf(stderr, "Can not open file curves.\n");
		exit(1);
	}
	fseek(file_curves, nb_curve * line_file_curves, SEEK_SET);
   	if(fscanf(file_curves, "%2" SCNu8 "%s %s %s %s", &nb_bits, str_A, str_B, str_p, str_large_prime) < 5)
   	{
		fprintf(stderr, "Can not read file curves.\n");
		exit(1);
	}
	fclose(file_curves);
	mpz_set_str(E.A, str_A, 10);
	mpz_set_str(E.B, str_B, 10);
	mpz_set_str(E.p, str_p, 10);
	mpz_set_str(large_prime, str_large_prime, 10);
		
	generate_adding_sets(A, B, large_prime);
    
	test_i = 0;
	while(test_i < nb_tests)
	{
		/*** read point P ***/
		nb_point = test_i % 10 + 1;
		file_points = fopen("points","r");
		if (file_points == NULL) 
		{
			fprintf(stderr, "Can not open file points.\n");
			exit(1);
		}
		fseek(file_points, nb_curve * (nb_points_file + 1) * line_file_points + (nb_point * line_file_points), SEEK_SET);
		if(fscanf(file_points, "%s %s",str_X, str_Y) < 2)
		{
			fprintf(stderr, "Can not read file points.\n");
			exit(1);
		}
		fclose(file_points);
		mpz_set_str(P.x, str_X, 10);
		mpz_set_str(P.y, str_Y, 10);
		mpz_set_ui(P.z, 1);
		
		//choose a key of size: nb_bits
		generate_random_key(key, nb_bits - 1);
		//compute Q
		double_and_add(&Q, P, key, E);
		
		
		/******* BEGIN: Setting up environment for experiment running and statistics *******/
		
		/*** Update possible values of argument f (field/nb_bits) in experiments ***/
		file_conf = fopen(RESULTS_PATH"conf_avg/f.conf","r");
		if (file_conf == NULL) 
		{
			fprintf(stderr, "Can not open configuration file (see constant RESULTS_PATH in main.c)\n");
			exit(1);
		}
		memset(conf_str, 0, 1000);
		if(fgets(conf_str, 1000, file_conf) != NULL)
		{
			strncpy(conf_str_cpy, conf_str, 1000);
			conf_value = strtok(conf_str_cpy, " ");
		}
		fclose(file_conf);
		update = 1;
		while(conf_value != NULL)
		{
			if(atoi(conf_value) == nb_bits)
			{
				update = 0;
			}
			conf_value = strtok(NULL, " ");
		}
		if(update)
		{
			file_conf=fopen(RESULTS_PATH"conf_avg/f.conf","w");
			if (file_conf == NULL) 
			{
				fprintf(stderr, "Can not open configuration file (see constant RESULTS_PATH in main.c)\n");
				exit(1);
			}
			fprintf(file_conf, "%s %2" SCNu8, conf_str, nb_bits);
			fclose(file_conf);
		}
		
		/*** Update possible values of argument s (storage structure) in experiments ***/
		for(struct_i = 0; struct_i < __NB_STRUCTURES__; struct_i++)
		{
			if(structs[struct_i] == 1)
			{
				file_conf = fopen(RESULTS_PATH"conf_avg/s.conf","r");
				if (file_conf == NULL) 
				{
					fprintf(stderr, "Can not open configuration file (see constant RESULTS_PATH in main.c)\n");
					exit(1);
				}
				memset(conf_str, 0, 1000);
				if(fgets(conf_str, 1000, file_conf) != NULL)
				{
					strncpy(conf_str_cpy, conf_str, 1000);
					conf_value = strtok(conf_str_cpy, " ");
				}
				fclose(file_conf);
				update = 1;
				while(conf_value != NULL)
				{
					if(strncmp(conf_value, struct_i_str[struct_i], strlen(conf_value)) == 0 && strlen(conf_value) == strlen(struct_i_str[struct_i]))
					{
						update = 0;
					}
					conf_value = strtok(NULL, " ");
				}
				if(update)
				{
					file_conf=fopen(RESULTS_PATH"conf_avg/s.conf","w");
					if (file_conf == NULL) 
					{
						fprintf(stderr, "Can not open configuration file (see constant RESULTS_PATH in main.c)\n");
						exit(1);
					}
					fprintf(file_conf, "%s %s", conf_str, struct_i_str[struct_i]);
					fclose(file_conf);
				}
			}
		}
		
		/*** Update possible values of argument t (thread number) in experiments ***/
		file_conf = fopen(RESULTS_PATH"conf_avg/t.conf","r");
		if (file_conf == NULL) 
		{
			fprintf(stderr, "Can not open configuration file (see constant RESULTS_PATH in main.c)\n");
			exit(1);
		}
		memset(conf_str, 0, 1000);
		if(fgets(conf_str, 1000, file_conf) != NULL)
		{
			strncpy(conf_str_cpy, conf_str, 1000);
			conf_value = strtok(conf_str_cpy, " ");
		}
		fclose(file_conf);
		update = 1;
		while(conf_value != NULL)
		{
			if(atoi(conf_value) == nb_threads)
			{
				update = 0;
			}
			conf_value = strtok(NULL, " ");
		}
		if(update)
		{
			file_conf=fopen(RESULTS_PATH"conf_avg/t.conf","w");
			if (file_conf == NULL) 
			{
				fprintf(stderr, "Can not open configuration file (see constant RESULTS_PATH in main.c)\n");
				exit(1);
			}
			fprintf(file_conf, "%s %d", conf_str, nb_threads);
			fclose(file_conf);
		}
		
		/*** Update possible values of argument d (number of trailling bits to zero) in experiments ***/
		file_conf = fopen(RESULTS_PATH"conf_avg/theta.conf","r");
		if (file_conf == NULL) 
		{
			fprintf(stderr, "Can not open configuration file (see constant RESULTS_PATH in main.c)\n");
			exit(1);
		}
		memset(conf_str, 0, 1000);
		if(fgets(conf_str, 1000, file_conf) != NULL)
		{
			strncpy(conf_str_cpy, conf_str, 1000);
			conf_value = strtok(conf_str_cpy, " ");
		}
		fclose(file_conf);
		update = 1;
		while(conf_value != NULL)
		{
			if(atoi(conf_value) == trailling_bits)
			{
				update = 0;
			}
			conf_value = strtok(NULL, " ");
		}
		if(update)
		{
			file_conf=fopen(RESULTS_PATH"conf_avg/theta.conf","w");
			if (file_conf == NULL) 
			{
				fprintf(stderr, "Can not open configuration file (see constant RESULTS_PATH in main.c)\n");
				exit(1);
			}
			fprintf(file_conf, "%s %2" SCNu8, conf_str, trailling_bits);
			fclose(file_conf);
		}
		
		/*** Update possible values of argument l (level of the abstract radix tree) in experiments ***/
		file_conf = fopen(RESULTS_PATH"conf_avg/l.conf","r");
		if (file_conf == NULL)
		{
			fprintf(stderr, "Can not open configuration file (see constant RESULTS_PATH in main.c)\n");
			exit(1);
		}
		memset(conf_str, 0, 1000);
		if(fgets(conf_str, 1000, file_conf) != NULL)
		{
			strncpy(conf_str_cpy, conf_str, 1000);
			conf_value = strtok(conf_str_cpy, " ");
		}
		fclose(file_conf);
		update = 1;
		while(conf_value != NULL)
		{
			if(atoi(conf_value) == level)
			{
				update = 0;
			}
			conf_value = strtok(NULL, " ");
		}
		if(update)
		{
			file_conf=fopen(RESULTS_PATH"conf_avg/l.conf","w");
			if (file_conf == NULL)
			{
				fprintf(stderr, "Can not open configuration file (see constant RESULTS_PATH in main.c)\n");
				exit(1);
			}
			fprintf(file_conf, "%s %2" SCNu8, conf_str, level);
			fclose(file_conf);
		}
		
		/******* END: Setting up environment for experiment running and statistics *******/
		
		
		/* Test different structures */
		printf("*** Test %d ***\n", test_i + 1);
		for(struct_i = 0; struct_i < __NB_STRUCTURES__; struct_i++)
		{
			if(structs[struct_i] == 1)
			{
				printf("\t**Structure %s\n", struct_i_str[struct_i]);
				pcs_init(P, Q, E, large_prime, A, B, nb_bits, trailling_bits, struct_i, nb_threads, level);
				gettimeofday(&tv1,NULL);
				pcs_run(x, nb_threads, nb_collisions);
				gettimeofday(&tv2, NULL);
				time1=(tv1.tv_sec) * 1000000 + tv1.tv_usec;
				time2 = (tv2.tv_sec) * 1000000 + tv2.tv_usec;
				time = time2 - time1;
                memory = struct_memory(&nb_points, &rate_of_use, &rate_slots, nb_threads);
				pcs_clear();
                
				if(mpz_cmp(x, key)!=0)
				{
					fprintf(stderr, "Error in PCS computation.\n");
					//exit(2);
					break;
				}
				
				/*** Write execution time ***/
				file_res=fopen(RESULTS_PATH"time.all","a");
				if (file_res == NULL) 
				{
					fprintf(stderr, "Can not open file time.all (see constant RESULTS_PATH in main.c)\n");
					exit(1);
				}
				fprintf(file_res,"%d %s %d %d %"SCNu8" %llu\n", nb_bits, struct_i_str[struct_i], nb_threads, trailling_bits, level, time);
				fclose(file_res);
				
				/*** Write memory usage ***/
				file_res=fopen(RESULTS_PATH"memory.all","a");
				if (file_res == NULL) 
				{
					fprintf(stderr, "Can not open file memory.all (see constant RESULTS_PATH in main.c)\n");
					exit(1);
				}
				fprintf(file_res,"%d %s %d %d %"SCNu8" %llu\n", nb_bits, struct_i_str[struct_i], nb_threads, trailling_bits, level, memory);
				fclose(file_res);
				
				/*** Write number of stored points ***/
				file_res=fopen(RESULTS_PATH"points.all","a");
				if (file_res == NULL) 
				{
					fprintf(stderr, "Can not open file points.all (see constant RESULTS_PATH in main.c)\n");
					exit(1);
				}
				fprintf(file_res,"%d %s %d %d %"SCNu8" %lu\n", nb_bits, struct_i_str[struct_i], nb_threads, trailling_bits, level, nb_points);
				fclose(file_res);
				
				/*** Write rate of memory use ***/
				file_res=fopen(RESULTS_PATH"rate.all","a");
				if (file_res == NULL) 
				{
					fprintf(stderr, "Can not open file rate.all (see constant RESULTS_PATH in main.c)\n");
					exit(1);
				}
				fprintf(file_res,"%d %s %d %d %"SCNu8" %.2f (%.2f)\n", nb_bits, struct_i_str[struct_i], nb_threads, trailling_bits, level, rate_of_use, rate_slots);
				fclose(file_res);
			}
		}
		test_i++;
	}
	curve_clear(&E);
	point_clear(&P);
	point_clear(&Q);
	mpz_clears(x, large_prime, key, NULL);
	for(j=0;j<__NB_ENSEMBLES__;j++)
	{
		mpz_clears(A[j],B[j],NULL);
	}
	preallocation_clear();
}
