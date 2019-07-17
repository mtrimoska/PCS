/** @file pcs_struct_PRTL.h
 *
 *	Created by Monika Trimoska on 03/12/2016.
 *	Copyright © 2016 Monika Trimoska. All rights reserved.
 */

#include <gmp.h>
#include <omp.h>

void struct_init_PRTL(uint8_t nb_bits, uint8_t trailling_bits, int nb_threads, uint8_t _level);
int struct_add_PRTL(mpz_t a_out, mpz_t a_in, mpz_t xDist);
void struct_free_PRTL(void);
unsigned long long int struct_memory_PRTL(unsigned long int *nb_points, float *rate_of_use, float *rate_slots);

