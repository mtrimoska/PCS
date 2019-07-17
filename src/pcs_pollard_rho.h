/** @file pcs_pollard_rho.h
 *
 *	Created by Monika Trimoska on 03/12/2015.
 *	Copyright © 2015 Monika Trimoska. All rights reserved.
 */

void compute_a(mpz_t a, mpz_t a2, mpz_t n);
void compute_b(mpz_t b, mpz_t b2, mpz_t n);
void f(point_t oldR, point_t M, point_t * newR, elliptic_curve_t e);
void pollard_rho(mpz_t x, point_t P, point_t Q, elliptic_curve_t E, mpz_t n, mpz_t A[20], mpz_t B[20]);
void compute_x(mpz_t x, mpz_t a1, mpz_t a2, mpz_t b1, mpz_t b2, mpz_t n);
int hash(mpz_t donnee);
