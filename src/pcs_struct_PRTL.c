/** @file pcs_struct_PRTL.c
 *  @brief Implementation of the PRTL structure
 *
 *	Created by Monika Trimoska on 03/12/2016.
 *	Copyright © 2016 Monika Trimoska. All rights reserved.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <inttypes.h>
#include "pcs_vect_bin.h"
#include "pcs_struct_PRTL.h"

static uint8_t nb_bits;
static uint8_t level;
static _vect_bin_chain_t *chain_array;
static int chain_array_size;
static omp_lock_t *locks;
static int xDist_start;
static int xDist_end;
static int a_start;
static int a_end;
static int suffix_len;
/***Memory limiting feature is turned off
static unsigned long long int memory_limit;
 ***/
static unsigned long long int memory_alloc;
static mpz_t mask;

/** Initialize the Packed Radix-Tree-List.
 *
 *	@brief Initialize the PRTL, allocate memory and
 *	create mask which will be used to map a stored point
 *	to an index of the chain array.
 *
 */
void struct_init_PRTL(uint8_t _nb_bits, uint8_t trailling_bits, int nb_threads, uint8_t _level)
{
	int i;
	nb_bits = _nb_bits;
    uint8_t c = nb_bits - trailling_bits;
	level = _level;
    suffix_len = c - level;
    unsigned int __vect_bin_size = nb_bits + suffix_len; //already init as a constant
    xDist_start = 0;
    xDist_end = suffix_len - 1;
    a_start = suffix_len;
    a_end = xDist_end + nb_bits;
    chain_array_size = pow(2, level);
     _vect_bin_t_initiate;
    
    /* allocate chain table */
	
    /***Memory limiting feature is turned off
	memory_limit = 100000000;
	 ***/
    chain_array = (_vect_bin_chain_t *) malloc(sizeof(_vect_bin_chain_t) * chain_array_size);
    _vect_bin_t_count_memory(chain_array_size);
    locks = malloc(sizeof(omp_lock_t) * chain_array_size);
    memory_alloc = sizeof(omp_lock_t) * chain_array_size;
    for(i = 0; i < chain_array_size; i++)
	{
		vect_bin_t_reset(chain_array[i].v);
        chain_array[i].nxt = NULL;
        omp_init_lock(&locks[i]);
	}
    
    /* create mask */
    mpz_inits(mask, NULL);
    mpz_set_ui(mask, 1);
    mpz_mul_2exp(mask, mask, level); //left shift
    mpz_sub_ui(mask, mask, 1);
}

/** Search and insert function for the PRTL structure.
 *
 *  @brief Look for a point in the structure. If the point is not found
 *  it is added with the corresponding a coefficient.
 *
 *  @param[out]	a_out	The a coefficient of the found point.
 *  @param[in]	a_in	The a coefficient of the newly added point.
 *  @param[in]	xDist	The x-coordinate, without the trailling zeros.
 *  @return 	1 if the point was found, 0 otherwise.
 */
int struct_add_PRTL(mpz_t a_out, mpz_t a_in, mpz_t xDist)
{
    uint8_t retval = 0;
	_vect_bin_chain_t *new;
	_vect_bin_chain_t *last;
	_vect_bin_chain_t *next;
    mpz_t key_mpz;
    int key;
    
    mpz_inits(key_mpz, NULL);
	
    mpz_and(key_mpz, xDist, mask);
    key = mpz_get_ui(key_mpz);
    omp_set_lock(&locks[key]);
    next = &chain_array[key];
    if(vect_bin_is_empty(next->v))
    {
        vect_bin_set_mpz(next->v, xDist_start, suffix_len, xDist, level);
        vect_bin_set_mpz(next->v, a_start, nb_bits, a_in, 0);
        next->nxt=NULL;
    }
    else
    {
        while(next != NULL && vect_bin_cmp_mpz(next->v, xDist_start, suffix_len, xDist, level) < 0)
        {
            last = next;
            next = next->nxt;
        }
        if(next != NULL && vect_bin_cmp_mpz(next->v, xDist_start, suffix_len, xDist, level) == 0 ) //collision
        {
            vect_bin_get_mpz(next->v, a_start, nb_bits, a_out);
            retval = 1;
        }
        else
        {
			/***Memory limiting feature is turned off
            if(memory_alloc + _vect_bin_alloc_size < memory_limit)
            {
			 ***/
                _vect_bin_chain_t_new(new);
                vect_bin_t_reset(new->v);
                if(next == &chain_array[key]) //add at the beginning
                {
                    vect_bin_cpy(new->v, next->v);
                    new->nxt = next->nxt;

                    vect_bin_t_reset(next->v);
                    vect_bin_set_mpz(next->v, xDist_start, suffix_len, xDist, level);
                    vect_bin_set_mpz(next->v, a_start, nb_bits, a_in, 0);

                    next->nxt = new;   
                }
                else
                {
                    vect_bin_set_mpz(new->v, xDist_start, suffix_len, xDist, level);
                    vect_bin_set_mpz(new->v, a_start, nb_bits, a_in, 0);
                    if(next != NULL) //add in the middle
                    {	
                        new->nxt = next;
                    }
                    last->nxt = new;  
                }
            //}
        }
    }
	omp_unset_lock(&locks[key]);
    mpz_clears(key_mpz, NULL);
	return retval;
}

/** Free the allocated memory for the Packed Radix-Tree-List.
 *
 */
void struct_free_PRTL(void)
{
    int i;
	_vect_bin_chain_t *last;
	_vect_bin_chain_t *next;
    omp_destroy_lock(&_vect_alloc_size_lock);
	for(i = 0; i < chain_array_size; i++)
	{
        omp_destroy_lock(&locks[i]);
		next = chain_array[i].nxt;
		while(next != NULL)
		{
			last = next;
			next = next->nxt;
			_vect_bin_chain_t_free(last);
		}
	}
	free(chain_array);
    _vect_bin_t_count_memory(-chain_array_size);
	free(locks);
    mpz_clears(mask, NULL);
}

/** Recursive function used in struct_memory_PRTL_rec.
 *
 */
void struct_memory_PRTL_rec(_vect_bin_chain_t *it, unsigned long int *nb_points)
{
	if (it != NULL)
    {
        (*nb_points)++;
        struct_memory_PRTL_rec(it->nxt, nb_points);
    }
}

/** Get the memory occupation of the PRTL structure.
 *
 *	@brief Calculates the total memory occupation,
 *	the total number of stored points,
 *	the number of empty slots and the rate of memory use
 *	(ratio between the allocated empty and allocated used memory).
 *
 *  @return	The memory occupation in bytes.
 */
unsigned long long int struct_memory_PRTL(unsigned long int *nb_points, float *rate_of_use, float *rate_slots)
{
	int i = 0;
	unsigned long long int sum = 0;
	unsigned long long int lost = 0;
    int empty_slots = 0;
    *nb_points = 0;
	sum += sizeof(omp_lock_t) * chain_array_size;
	for(i = 0; i < chain_array_size; i++)
	{
		if(vect_bin_is_empty(chain_array[i].v))
		{
			lost += sizeof(omp_lock_t) + sizeof(_vect_bin_chain_t);
			empty_slots++;
		}
		else
		{
			(*nb_points)++;
            if(chain_array[i].nxt != NULL)
				struct_memory_PRTL_rec(chain_array[i].nxt, nb_points);
		}
	}
    
    sum += _vect_bin_alloc_size;
    *rate_of_use = (1.0 - ((float)lost) / ((float)sum)) * 100.0;
    *rate_slots = (1.0 - ((float)empty_slots) / ((float)chain_array_size)) * 100.0;
	printf("\t\tPoints: %lu\n", *nb_points);
    printf("\t\tEmpty slots: %d\n", empty_slots);
    return sum;
}
