/** @file pcs_struct_hash.c
 *  @brief Implementation of a hash table.
 *
 *	Created by Monika Trimoska on 03/12/2015.
 *	Copyright © 2015 Monika Trimoska. All rights reserved.
 */

#include <omp.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pcs_struct_hash_UNIX.h"
#include "pcs_struct_hash.h"

#define __PI_NUMERATOR__ 355  	// correct to three digits
#define __PI_DENOMINATOR__ 113	// correct to three digits

static uint8_t hash_type;
static hashUNIX_t **table;
static unsigned long int table_size;
static omp_lock_t *table_locks;
static omp_lock_t memory_alloc_lock;
static unsigned long long int memory_alloc;
/***Memory limiting feature is turned off
static unsigned long long int memory_limit;
 ***/


/** Calculate the hash table recommended size.
 *
 *	@brief The recommended size corresponds to the expected number
 *	of stored points. This function is used only if the size is not specified
 *	(using the level parameter).
 *
 */
void set_table_size(mpz_t n, uint8_t trailling_bits)
{
	unsigned long long int distinguished;
	mpz_t table_size_inter;
	mpz_init(table_size_inter);
	mpz_mul_ui(table_size_inter, n, __PI_NUMERATOR__);
	mpz_tdiv_q_ui(table_size_inter, table_size_inter, 2 * __PI_DENOMINATOR__);
	mpz_sqrt(table_size_inter, table_size_inter);
	distinguished = (unsigned long int)pow(2, trailling_bits);
	mpz_tdiv_q_ui(table_size_inter, table_size_inter, distinguished);
	table_size = mpz_get_ui(table_size_inter);
	mpz_clear(table_size_inter);
}

/** Initialize the hash table and allocate memory.
 *
 */
void struct_init_hash(uint8_t hash_type_init, mpz_t n, uint8_t trailling_bits, uint8_t level)
{
	unsigned long int i;
	hash_type = hash_type_init;
    if(level != 7)
    {
        table_size = pow(2, level);
    }
    else
    {
        set_table_size(n, trailling_bits);
    }

	printf("\t\ttable_size: %lu\n",table_size);
	table = malloc(sizeof(*table) * table_size); //i.e. sizeof(hashUNIX_t *)
	table_locks = malloc(sizeof(omp_lock_t) * table_size);
    omp_init_lock(&memory_alloc_lock);
	/***Memory limiting feature is turned off
    memory_limit = 100000000;
	 ***/
    memory_alloc = 0LL;
    memory_alloc += sizeof(*table) * table_size;
	memory_alloc += sizeof(omp_lock_t) * table_size;
	for(i = 0; i < table_size; i++)
	{
		table[i] = NULL;
		omp_init_lock(&table_locks[i]);
	}
}

unsigned long int get_hash(char *xDist)
{
	switch(hash_type)
	{
		default: 
			return (get_hash_UNIX(xDist) % table_size);
	}
}

int struct_add_hash(mpz_t a_out, mpz_t a_in, char xDist[])
{
	unsigned long int h = 0;
	hashUNIX_t *new;
	hashUNIX_t *last;
	hashUNIX_t *next;
	uint8_t retval = 0;
	
	h = get_hash(xDist);
	omp_set_lock(&table_locks[h]);
	next = table[h];
    while(next != NULL && next-> key != NULL && strncmp(xDist, next->key, strlen(xDist)) > 0)
	{
		last = next;
		next = next->next;
	}
	
	if(next != NULL && next->key != NULL && strcmp(xDist, next->key) == 0 ) //collision
	{
		mpz_set_str(a_out, next->a_, 62);
		retval = 1;
	}
	else
	{
		/***Memory limiting feature is turned off
		if(memory_alloc < memory_limit)
        {
		 ***/
            new = malloc(sizeof(hashUNIX_t));
            new->key = strdup(xDist);
            new->a_ = NULL;
            new->a_ = mpz_get_str(new->a_, 62, a_in);
            new->next = NULL;

            if(next == table[h]) //add at the beginning
            {
                new->next = next;
                table[h] = new;
            }
            else
            {
                if(next != NULL) //add in the middle
                {	
                    new->next = next;
                }
                last->next = new;
            }
            omp_set_lock(&memory_alloc_lock);
            memory_alloc += strlen(new->key) + 1;
            memory_alloc += strlen(new->a_) + 1;
            memory_alloc += sizeof(hashUNIX_t);
            omp_unset_lock(&memory_alloc_lock);
        //}
	}
	omp_unset_lock(&table_locks[h]);
	return retval;
}
void struct_free_hash(void)
{
	unsigned long int i;
	hashUNIX_t *last;
	hashUNIX_t *next;
    omp_destroy_lock(&memory_alloc_lock);
	for(i = 0; i < table_size; i++)
	{
		next = table[i];
		omp_destroy_lock(&table_locks[i]);
		while(next != NULL)
		{
			free(next->key);
			free(next->a_);
			last = next;
			next = next->next;
			free(last);
		}
	}
	free(table);
	free(table_locks);
}

unsigned long long int struct_memory_hash_rec(hashUNIX_t *it, unsigned long int *nb_points, int *link)
{
	unsigned long long int sum = 0;
	
	if (it == NULL)
		return 0;
		
    (*link)++;
	sum += strlen(it->key) + 1;
	sum += strlen(it->a_) + 1;
	sum += sizeof(hashUNIX_t);
	(*nb_points)++;
	return (sum + struct_memory_hash_rec(it->next, nb_points, link));
}

unsigned long long int struct_memory_hash(unsigned long int *nb_points, float *rate_of_use, float *rate_slots) 
{
    int link_long[1000] = {0};
	unsigned long long int sum = 0;
	unsigned long long int lost = 0;
	unsigned long int i;
	unsigned long int empty_slots = 0;
    hashUNIX_t *next;
	*nb_points = 0;
    int link;
    sum += sizeof(*table) * table_size;
	sum += sizeof(omp_lock_t) * table_size;
	for(i = 0; i < table_size; i++)
	{
        link = 0;
		if(table[i] == NULL)
		{
			lost += sizeof(*table);
			lost += sizeof(omp_lock_t);
			empty_slots++;
		}
		else
		{
            link = 1;
			(*nb_points)++;
			sum += sizeof(hashUNIX_t);
			sum += strlen(table[i]->key) + 1;
			sum += strlen(table[i]->a_) + 1;
			if(table[i]->next != NULL)
				sum += struct_memory_hash_rec(table[i]->next, nb_points, &link);
		}
        //link_long[link]++;
	}
     //for(i = 0; i < 60; i++)
        //printf("long %d: %d\n",i,link_long[i]);
	*rate_of_use = (1.0 - ((float)lost) / ((float)sum)) * 100.0;
    *rate_slots = (1.0 - ((float)empty_slots) / ((float)table_size)) * 100.0;
    printf("\t\tPoints: %lu\n", *nb_points);
	printf("\t\tEmpty slots: %lu\n", empty_slots);
	return sum;
}
