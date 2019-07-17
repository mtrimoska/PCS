/** @file pcs_struct_hash.h
 *
 *	Created by Monika Trimoska on 03/12/2015.
 *	Copyright © 2015 Monika Trimoska. All rights reserved.
 */

#include <omp.h>
#include <gmp.h>
#include <inttypes.h>

typedef struct hashUNIX
{
	char *key;
	char *a_;
	struct hashUNIX *next;
	
}hashUNIX_t;

void struct_init_hash(uint8_t hash_type_init, mpz_t n, uint8_t trailling_bits, uint8_t level);
int struct_add_hash(mpz_t a_out, mpz_t a_in, char xDist[]);
void struct_free_hash(void);
unsigned long long int struct_memory_hash(unsigned long int *nb_points, float *rate_of_use, float *rate_slots);
