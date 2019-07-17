/** @file pcs_storage.h
 *
 *	Created by Monika Trimoska on 03/12/2015.
 *	Copyright © 2015 Monika Trimoska. All rights reserved.
 */

#include <inttypes.h>

void struct_init(uint8_t type, mpz_t n, uint8_t trailling_bits, uint8_t nb_bits, int nb_threads, uint8_t level);
int struct_add(mpz_t a_out, mpz_t a_in, mpz_t xDist, char xDist_str[]);
void struct_free();
unsigned long long int struct_memory(unsigned long int *nb_points, float *rate_of_use, float *rate_slots, int nb_threads);
