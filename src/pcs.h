/** @file pcs.h
 *
 *	Created by Monika Trimoska on 03/12/2015.
 *	Copyright © 2015 Monika Trimoska. All rights reserved.
 */

#include <gmp.h>
#include <omp.h>
#include <inttypes.h>

#define __NB_ENSEMBLES__ 20

void combLin(point_t * R, mpz_t a, mpz_t b);
void pcs_init(point_t P_init, point_t Q_init, elliptic_curve_t E_init, mpz_t n_init, mpz_t *A_init, mpz_t *B_init, uint8_t nb_bits_init, uint8_t trailling_bits_init, int type_struct, int nb_threads, uint8_t level);
long long int pcs_run(mpz_t x_res, int nb_threads, int nb_collisions);
void pcs_clear();
