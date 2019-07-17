/** @file pcs_exec.c
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
    printf("Usage: [\n-f nb_bits \n-t nb_threads \n-n nb_tests \n-s structure (PRTL or hash_unix) \n-l tree_level \n-d nb_trailling_bits_to_0 \n-c nb_collisions\n] \n");
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
	char conf_str[1001] = {0};
	char conf_str_cpy[1001] = {0};
	char *conf_value;
    char joker[20] = "_";
	uint8_t update = 0;
	uint8_t struct_chosen = 0;
	struct timeval tv1;
	struct timeval tv2;
	unsigned long long int time, time1, time2;
	unsigned long long int sum_time[__NB_STRUCTURES__] = {0};
	unsigned long long int sum_memory[__NB_STRUCTURES__] = {0};
	unsigned long long int sum_nb_points[__NB_STRUCTURES__] = {0};
	unsigned long long int memory;
	unsigned long int nb_points;
	float rate_of_use, rate_slots;
	mpz_t key;
	mpz_t x;
	char option;
	uint8_t nb_bits, trailling_bits, nb_curve, line_file_curves, line_file_points, nb_points_file, nb_point, j, struct_i, level;
	int test_i, nb_tests, nb_threads;
	int nb_collisions = 1;
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
			case 'f' : nb_bits = atoi(optarg) / 5 * 5;
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
			case 'd' : trailling_bits = atoi(optarg);
				break;
			case 'h' : {print_usage();exit(0);}
				break;
		}
	}
    
    if(trailling_bits == 0)
	   trailling_bits = nb_bits / 4;
	
	if(!struct_chosen)
		add_to_struct_options(structs, struct_i_str, "PRTL", &struct_chosen);
	
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
}
