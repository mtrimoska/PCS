/** @file pcs.c
 *
 *	Created by Monika Trimoska on 03/12/2015.
 *	Copyright © 2015 Monika Trimoska. All rights reserved.
 */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <gmp.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include "pcs_elliptic_curve_operations.h"
#include "pcs_pollard_rho.h"
#include "pcs_storage.h"
#include "pcs.h"

elliptic_curve_t E;
point_t P;
point_t Q;
mpz_t n;
mpz_t *A; 
mpz_t *B;
point_t M[__NB_ENSEMBLES__];
uint8_t trailling_bits;
uint8_t nb_bits;

/** Determines whether a point is a distinguished one.
 *
 *  @param[in]	R				A point on an elliptic curve.
 *  @param[in]	trailling_bits	Number of trailling zero bits in a ditinguished point.
 *  @param[out]	q				The x-coordinate, without the trailling zeros.
 *  @return 	1 if the point is distinguished, 0 otherwise.
 */
int is_distinguished(point_t R, int trailling_bits, mpz_t *q)
{
	int res;
	mpz_t r;
	mpz_inits(r, NULL);
	mpz_tdiv_qr_ui(*q, r, R.x, (unsigned long int)pow(2, trailling_bits));
	res=(mpz_sgn(r) == 0);
	mpz_clears(r, NULL);
	return (res);
}

/** Checks if the linear combination aP+bQ is equal to R or its inverse.
 *
 *  @param[in]	R	A point on an elliptic curve.
 *  @param[in]	a	a coefficient.
 *  @param[in]	b	b coefficient.
 *  @return 	1 if aP+bQ = R, 0 if aP+bQ = -R.
 */
int same_point(point_t R, mpz_t a, mpz_t b)
{
	int res;
	point_t S1, S2, S;
	mpz_inits(S1.x, S1.y, S1.z, S2.x, S2.y, S2.z, S.x, S.y, S.z, NULL);
	double_and_add(&S1, P, a, E);
	double_and_add(&S2, Q, b, E);
	add(&S, S1, S2, E);
	res=(mpz_cmp(R.y, S.y) == 0);
	mpz_clears(S1.x, S1.y, S1.z, S2.x, S2.y, S2.z, S.x, S.y, S.z, NULL);
	return res;
}

/** Computes the linear combination aP+bQ on E.
 *
 *  @param[out]	R	Resulting point.
 *  @param[in]	a	a coefficient.
 *  @param[in]	b	b coefficient.
 */
void lin_comb(point_t * R, mpz_t a, mpz_t b)
{
	point_t S1, S2;
	mpz_inits(S1.x, S1.y, S1.z, S2.x, S2.y, S2.z, NULL);
	double_and_add(&S1, P, a, E);
	double_and_add(&S2, Q, b, E);
	add(R, S1, S2, E);
	mpz_clears(S1.x, S1.y, S1.z, S2.x, S2.y, S2.z, NULL);
}

/** Checks if there is a collision.
 *
 */
int is_collision(mpz_t x, mpz_t a1, mpz_t a2, int trailling_bits)
{
	uint8_t r;
	mpz_t xDist_;
	int retval = 0;
	mpz_t b1, b2;
	point_t R;
	point_init(&R);
	mpz_inits(b1, b2, xDist_, NULL);
	
	mpz_set_ui(b2, 0);
	mpz_set_ui(b1, 0);
	double_and_add(&R, P, a1, E);
	//recompute first a,b pair
	while(!is_distinguished(R, trailling_bits, &xDist_))
	{
		r = hash(R.y);
		compute_a(a1, A[r], n);
		compute_b(b1, B[r], n);
		f(R, M[r], &R, E);
	}
	
	//recompute second a,b pair
	double_and_add(&R, P, a2, E);
	while(!is_distinguished(R, trailling_bits, &xDist_))
	{
		r = hash(R.y);
		compute_a(a2, A[r], n);
		compute_b(b2, B[r], n);
		f(R, M[r], &R, E);
	}
	if(mpz_cmp(b1, b2) != 0) //we found two different pairs, so collision
	{
		if(!same_point(R, a1, b1)) //it's the inverse point
		{	
			mpz_neg(a2, a2); 
			mpz_mmod(a2, a2, n);
			mpz_neg(b2, b2);
			mpz_mmod(b2, b2, n);
		}
		compute_x(x, a1, a2, b1, b2, n);
		retval = 1;
	}
	point_clear(&R);
	mpz_clears(b1, b2, xDist_, NULL);
	return retval;
}

/** Initialize all variables needed to do a PCS algorithm.
 *
 */
void pcs_init(point_t P_init, point_t Q_init, elliptic_curve_t E_init, mpz_t n_init, mpz_t *A_init, mpz_t *B_init, uint8_t nb_bits_init, uint8_t trailling_bits_init, int type_struct, int nb_threads, uint8_t level)
{
	uint8_t i;
	
	point_init(&P);
	point_init(&Q);
	curve_init(&E);
	mpz_init(n);
	
	mpz_set(P.x, P_init.x);
	mpz_set(P.y, P_init.y);
	mpz_set(P.z, P_init.z);
	
	mpz_set(Q.x, Q_init.x);
	mpz_set(Q.y, Q_init.y);
	mpz_set(Q.z, Q_init.z);
	
	mpz_set(E.A, E_init.A);
	mpz_set(E.B, E_init.B);
	mpz_set(E.p, E_init.p);
	
	mpz_set(n, n_init);
	
	A = A_init;
	B = B_init;
	
	for(i=0; i<__NB_ENSEMBLES__; i++)
	{
		mpz_inits(M[i].x,M[i].y,M[i].z,NULL);
		lin_comb(&M[i],A[i],B[i]);
	}
	
	trailling_bits = trailling_bits_init;
	nb_bits = nb_bits_init;
	
	struct_init(type_struct, n, trailling_bits, nb_bits, nb_threads, level);
}

/** Run the PCS algorithm.
 *
 */
long long int pcs_run(mpz_t x_res, int nb_threads, int nb_collisions)
{
	point_t R;
	mpz_t a, a2;
	mpz_t x, xDist;
	uint8_t r;
    int trail_length;
    int trail_length_max = pow(2, trailling_bits) * 20;
    int collision_count = 0;
	char xDist_str[50];
	#pragma omp parallel private(R, a, a2, x, r, xDist, xDist_str, trail_length) shared(collision_count, x_res, trail_length_max) num_threads(nb_threads)
	{
		point_init(&R);
		mpz_inits(x, a2, a, xDist, NULL);
		
		//Initialize a starting point
		gmp_randstate_t r_state;
		gmp_randinit_default(r_state);
		gmp_randseed_ui(r_state, time(NULL) * (omp_get_thread_num() + 1));
		mpz_urandomb(a, r_state, nb_bits);
		double_and_add(&R, P, a, E);
        trail_length = 0;
		
		while(collision_count < nb_collisions)
		{
            if(is_distinguished(R, trailling_bits, &xDist))
			{
                if(struct_add(a2, a, xDist, xDist_str))
				{
					if(is_collision(x, a, a2, trailling_bits))
					{
						#pragma omp critical
						{
							collision_count++;
                            mpz_set(x_res, x);
						}
					}
				}
				mpz_urandomb(a, r_state, nb_bits);
				double_and_add(&R, P, a, E);
                trail_length = 0;
			}
			else
			{
				r=hash(R.y);
				f(R, M[r], &R, E);
                trail_length++;
                if(trail_length > trail_length_max)
                {
                    mpz_urandomb(a, r_state, nb_bits);
                    double_and_add(&R, P, a, E);
                    trail_length = 0;
                }
			}
		}
		point_clear(&R);
		mpz_clears(a, a2, x, xDist, NULL);
		gmp_randclear(r_state);
	}
	return 0;
}

/** Free all variables used in the previous PCS run.
 *
 */
void pcs_clear()
{
	uint8_t i;
	point_clear(&P);
	point_clear(&Q);
	curve_clear(&E);
	mpz_clear(n);
	for(i = 0; i < __NB_ENSEMBLES__; i++)
	{
		mpz_clears(M[i].x, M[i].y, M[i].z, NULL);
	}
	struct_free();
}
