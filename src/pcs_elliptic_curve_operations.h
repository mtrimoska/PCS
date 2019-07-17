/** @file pcs_elliptic_curve_operations.h
 * 
 *	Created by Monika Trimoska on 03/12/2015.
 *	Copyright © 2015 Monika Trimoska. All rights reserved.
 */
#include<gmp.h>

/** Elliptic curve point structure
 */
typedef struct
{
	mpz_t x;
	mpz_t y;
	mpz_t z;
}point_t;

/** Elliptic curve structure
 *  @brief Elliptic curve defined by an equation of the form y^2 = x^3 + Ax + B
 */
typedef struct
{
	mpz_t A;
	mpz_t B;
	mpz_t p;
}elliptic_curve_t;

void point_init(point_t *P);
void curve_init(elliptic_curve_t *E);
void point_clear(point_t *P);
void curve_clear(elliptic_curve_t *E);
void mod(mpz_t a, mpz_t p);
void point_to_string(point_t P);
void curve_to_string(elliptic_curve_t E);
int is_elliptic_curve(elliptic_curve_t E);
int P_is_on_E(point_t P, elliptic_curve_t E);
int equal(point_t P1, point_t P2);
int add(point_t *P3, point_t P1, point_t P2, elliptic_curve_t E);
int _double(point_t * R, point_t P, elliptic_curve_t E);
int double_and_add(point_t *R, point_t P, mpz_t s, elliptic_curve_t E);
