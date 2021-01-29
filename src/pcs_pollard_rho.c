/** @file pcs_pollard_rho.c
 *  @brief Computing the random walk function and the classical Pollard's rho algorithm
 *
 *	Created by Monika Trimoska on 03/12/2015.
 *	Copyright © 2015 Monika Trimoska. All rights reserved.
 */

#include<stdio.h>
#include<math.h>
#include<gmp.h>
#include "pcs_elliptic_curve_operations.h"
#include "pcs_pollard_rho.h"
#include "pcs.h"

/** Chooses adding walk set.
 *
 * 	@param[in]		data	We choose the set according to this value.
 *	@return			The chosen set.
 */
int hash(mpz_t data)
{
	return mpz_mmod_ui(NULL, data, 20);
}

/** Compute coefficient a.
 *
 * 	@param[in,out]	a	Current coefficient a.
 * 	@param[in]		a2	Added coefficient a.
 * 	@param[in]		n	Group order.
 */
void compute_a(mpz_t a, mpz_t a2, mpz_t n)
{
	mpz_add(a, a, a2);
	mpz_mmod(a, a, n);
}

/** Compute coefficient b.
 *
 * 	@param[in,out]	b	Current coefficient b.
 * 	@param[in]		b2	Added coefficient b.
 * 	@param[in]		n	Group order.
 */
void compute_b(mpz_t b, mpz_t b2, mpz_t n)
{
	mpz_add(b, b, b2);
	mpz_mmod(b, b, n);
}

/** Compute the next step on the adding walk.
 *
 * 	@param[in]	oldR	Current point on the adding walk.
 * 	@param[in]	M		Adding step.
 * 	@param[out]	newR	Next point on the adding walk.
 * 	@param[in]	e		The elliptic curve.
 */
void f(point_t oldR, point_t M, point_t * newR, elliptic_curve_t e)
{
	add(newR, oldR, M, e);
}

/** Compute the discrete log using a and b coefficients.
 *
 * 	@param[out]	x	The discrete log.
 * 	@param[in]	a1	First a coefficient.
 * 	@param[in]	a2	Second a coefficient.
 * 	@param[in]	b1	First b coefficient.
 * 	@param[in]	b2	Second b coefficient.
 * 	@param[in]	n	Group order.
 */
void compute_x(mpz_t x, mpz_t a1, mpz_t a2, mpz_t b1, mpz_t b2, mpz_t n)
{
	mpz_t xUP, xDOWN;
	mpz_inits(xUP, xDOWN, NULL);	
	mpz_sub(xUP, a2, a1);
	mpz_mmod(xUP, xUP, n);
	mpz_sub(xDOWN, b1, b2);
	mpz_mmod(xDOWN, xDOWN, n);
	mpz_invert(xDOWN, xDOWN, n);
	mpz_mul(x, xUP, xDOWN);
	mpz_mmod(x, x, n);
	mpz_clears(xUP, xDOWN, NULL);
}

/** The classic Pollard's rho method.
 *
 * 	@param[out]	x	The discrete log.
 * 	@param[in]	P	Point P.
 * 	@param[in]	Q	Point Q = xP.
 * 	@param[in]	E	Elliptic curve.
 * 	@param[in]	n	Group order.
 * 	@param[in]	A	a coefficients of the adding walk sets.
 * 	@param[in]	B	b coefficients of the adding walk sets.
 */
void pollard_rho(mpz_t x, point_t P, point_t Q, elliptic_curve_t E, mpz_t n, mpz_t A[20], mpz_t B[20])
{
	mpz_t xUP, xDOWN;
	point_t R[2];
	mpz_t a[2];
	mpz_t b[2];
	mpz_init_set(R[0].x, P.x);
	mpz_init_set(R[0].y, P.y);
	mpz_init_set(R[0].z, P.z);
	mpz_init_set(R[1].x, P.x);
	mpz_init_set(R[1].y, P.y);
	mpz_init_set(R[1].z, P.z);
	mpz_init_set_ui(a[0], 1);
	mpz_init_set_ui(a[1], 1);
	mpz_init_set_ui(b[0], 0);
	mpz_init_set_ui(b[1], 0);
	int r;
	
	point_t M[20];
	for(r = 0; r < 20; r++)
	{
		mpz_inits(M[r].x,M[r].y,M[r].z,NULL);
	}
	
	r = hash(R[1].y);
	compute_a(a[1], A[r], n);
	compute_b(b[1], B[r], n);
	f(R[1], M[r], &R[1], E);
	while(!equal(R[0], R[1]))
	{
		r = hash(R[0].y);
		compute_a(a[0], A[r], n);
		compute_b(b[0], B[r], n);
		f(R[0], M[r], &R[0], E);
		
		r = hash(R[1].y);
		compute_a(a[1], A[r], n);
		compute_b(b[1], B[r], n);
		f(R[1], M[r], &R[1], E);  
		
		r = hash(R[1].y);
		compute_a(a[1], A[r], n);
		compute_b(b[1], B[r], n);
		f(R[1], M[r], &R[1], E);
		
	}
	compute_x(x, a[0], a[1], b[0], b[1], n);
	mpz_clears(a[0], a[1], b[0], b[1], NULL);
}
