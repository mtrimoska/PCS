/** @file pcs_elliptic_curve_operations.c
 *  @brief Functions for initializing the Point and Curve structures and performing elliptic curve operations.
 * 
 *	Created by Monika Trimoska on 03/12/2015.
 *	Copyright © 2015 Monika Trimoska. All rights reserved.
 */
 
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gmp.h>
#include "pcs_elliptic_curve_operations.h"

/*** BEGIN: Preallocation for GMP objects*/
char preallocation_init_done = 0;
mpz_t temp_obj[__NB_TEMP_MPZ_OBJ__];
point_t temp_point[__NB_TEMP_POINTS__];

/** Allocates the temp objects that are used in all functions.
 *
 */
void preallocation_init()
{
	int i;
	for(i = 0; i < __NB_TEMP_MPZ_OBJ__; i++)
	{
		mpz_init(temp_obj[i]);
	}
	for(i = 0; i < __NB_TEMP_MPZ_OBJ__; i++)
	{
		point_init(&temp_point[i]);
	}
	preallocation_init_done = 1;
}

/** Clears the temp objects that are used in all functions.
 *
 */
void preallocation_clear()
{
	int i;
	for(i = 0; i < __NB_TEMP_MPZ_OBJ__; i++)
	{
		mpz_clear(temp_obj[i]);
	}
	for(i = 0; i < __NB_TEMP_MPZ_OBJ__; i++)
	{
		point_clear(&temp_point[i]);
	}
	preallocation_init_done = 0;
}

/*** END: Preallocation for GMP objects*/

/** Initializes a point.
 * 
 * 	@param[in,out]	P	The point to be initialized.
 */
void point_init(point_t *P)
{
	mpz_inits(P->x, P->y, P->z, NULL);
}

/** Initializes a curve.
 * 
 * 	@param[in,out]	E	The curve to be initialized.
 */
void curve_init(elliptic_curve_t *E)
{
	mpz_inits(E->A, E->B, E->p, NULL);
}

/** Clears a point structure.
 * 
 * 	@param[in,out]	P	The point to be cleared.
 */
void point_clear(point_t *P)
{
	mpz_clears(P->x, P->y, P->z, NULL);
}

/** Clears a curve structure.
 * 
 * 	@param[in,out]	P	The curve to be cleared.
 */
void curve_clear(elliptic_curve_t *E)
{
	mpz_clears(E->A, E->B, E->p, NULL);
}

/**	The modulo operation.
 * 	@brief Performs the modulo operation on mpz numbers. The result is always non-negative.
 * 
 * 	@param[in,out]	a	The dividend. The value of a is changed in the function.
 * 	@param[in]		p	The divisor.
 */
void mod(mpz_t a, mpz_t p)
{
	mpz_mod(a, a, p);
	if(mpz_sgn(a) < 0)
	{
		mpz_add(a, a, p);
	}
}

/** Prints out a point to stdout in the form [x,y,z].
 * 
 * 	@param[in]	P	The point to be printed.
 */
void point_to_string(point_t P)
{
	printf("[%s,%s,%s]\n", mpz_get_str(NULL, 10, P.x), mpz_get_str(NULL, 10, P.y), mpz_get_str(NULL, 10, P.z));
}


/** Prints out an elliptic curve to stdout in the form y^2 = x^3 + Ax + B (mod p).
 * 
 * 	@param[in]	E	The curve to be printed.
 */
void curve_to_string(elliptic_curve_t E)
{
	printf("y^2 = x^3 + %s*x + %s (mod %s)\n", mpz_get_str(NULL, 10, E.A), mpz_get_str(NULL, 10, E.B), mpz_get_str(NULL, 10, E.p));
}

/** Checks if the structure repsents a nonsingular elliptic curve of the form y^2 = x^3 + Ax + B.
 * 
 * 	@param[in]	E	The curve to be checked.
 * 	@return 	Returns 1 for yes and 0 for no.
 */
int is_elliptic_curve(elliptic_curve_t E)
{
	int result;
	mpz_t *discriminant, *k1, *k2;
	if(!preallocation_init_done)
	{
		preallocation_init();
	}
	discriminant = &(temp_obj[0]);
	k1 = &(temp_obj[1]);
	k2 = &(temp_obj[2]);
	//discriminant=4*pow(A,3)+27*pow(B,2);
	mpz_mul(*k1, E.A, E.A);
	mpz_mul(*k1, *k1, E.A);
	mpz_mul_ui(*k1, *k1, 4);
	mpz_mul(*k2, E.B, E.B);
	mpz_mul_ui(*k2, *k2, 27);
	mpz_add(*discriminant, *k1, *k2);
	
	result = (mpz_sgn(*discriminant) != 0);
	return result;
}

/** Checks if P is a point on the curve E.
 * 	
 * 	@param[in]	P	The point.
 * 	@param[in]	E	The curve.
 * 	@return		Returns 1 for yes and 0 for no.
 */
int P_is_on_E(point_t P, elliptic_curve_t E)
{
	int result;
	if(mpz_get_ui(P.z) != 1) //P is the identity element or an invalid point
	{
		if(mpz_sgn(P.x) !=0 || mpz_get_ui(P.y) != 1)//P is an invalid point
		{
			return 0;
		}
		else //P is the identity element
		{
			return 1;
		}
	}
	
	mpz_t *left, *right, *sub_right;
	if(!preallocation_init_done)
	{
		preallocation_init();
	}
	left = &(temp_obj[3]);
	right = &(temp_obj[4]);
	sub_right = &(temp_obj[5]);
	mpz_mul(*left, P.y, P.y);
	mod(*left, E.p);
	mpz_mul(*sub_right, E.A, P.x);
	mpz_mul(*right, P.x, P.x);
	mpz_mul(*right, *right, P.x);
	mpz_add(*right, *right, *sub_right);
	mpz_add(*right, *right, E.B);
	mpz_mmod(*right, *right, E.p);
	
	result = (mpz_cmp(*left, *right) == 0);
	return result;
}

/** Compares two points on an elliptic curve.
 * 	@brief	The function assumes that P1 and P2 are valid points on a curve.
 * 
 * 	@param[in]	P1	The first point.
 * 	@param[in]	P2	The second point.
 * 	@return		Returns 1 is the two points are equal, 0 otherwise.
 */
int equal(point_t P1, point_t P2)
{
	return ((mpz_cmp(P1.x, P2.x) == 0) && (mpz_cmp(P1.y, P2.y) == 0) && (mpz_cmp(P1.z, P2.z) == 0));
}

/** Adds two points on a curve.
 * 
 * 	@param[out]	P3	P3 will be set to P1+P2. P3 should be initialized beforehand.
 * 	@param[in]	P1	The first point.
 * 	@param[in]	P2	The second point.
 * 	@param[in]	E	The curve.
 * 	@return		Returns 0 if the addition was performed, and 1 if a problem occured, exemple one of the two points is not on the curve.
 */
int add(point_t *P3, point_t P1, point_t P2, elliptic_curve_t E)
{
	if((P_is_on_E(P1, E) == 0) || (P_is_on_E(P2, E) == 0))
		return 1;
	
	//Addition with the identity element
	if(mpz_get_ui(P1.z) != 1)
	{
		mpz_set(P3->x, P2.x);
		mpz_set(P3->y, P2.y);
		mpz_set_ui(P3->z, 1);
		return 0;
	}
	if(mpz_get_ui(P2.z) != 1)
	{
		mpz_set(P3->x, P1.x);
		mpz_set(P3->y, P1.y);
		mpz_set_ui(P3->z, 1);
		return 0;
	}
	
	//Other cases
	mpz_t *l, *up, *down, *v, *up_bis, *x3, *y3;
	if(!preallocation_init_done)
	{
		preallocation_init();
	}
	l = &(temp_obj[6]);
	up = &(temp_obj[7]);
	down = &(temp_obj[8]);
	v = &(temp_obj[9]);
	up_bis = &(temp_obj[10]);
	x3 = &(temp_obj[11]);
	y3 = &(temp_obj[12]);
	
	if(equal(P1, P2))
	{
		if (mpz_sgn(P1.y) == 0)
		{
			mpz_set_ui(P3->x, 0);
			mpz_set_ui(P3->y, 1);
			mpz_set_ui(P3->z, 0);
			return 0;
		}
		else
		{	
			mpz_mul(*up, P1.x, P1.x),
			mpz_mul_ui(*up, *up, 3);
			mpz_add(*up, *up, E.A);
			mpz_mmod(*up, *up, E.p);
			
			mpz_mul_ui(*down, P1.y, 2);
			mpz_invert(*down, *down, E.p);
			
			mpz_mul(*l, *up, *down);
			mpz_mmod(*l, *l, E.p);
			
			mpz_mul(*up, P1.x, P1.x);
			mpz_mul(*up, *up, P1.x);
			mpz_neg(*up, *up);
			mpz_mul(*up_bis, P1.x, E.A);
			mpz_add(*up, *up, *up_bis);
			mpz_mul_ui(*up_bis, E.B, 2);
			mpz_add(*up, *up, *up_bis);
			
			mpz_mul(*v, *up, *down);
			mpz_mmod(*v, *v, E.p);
            
		}
	}
	else
	{
		if(mpz_cmp(P1.x, P2.x) == 0)
		{
			mpz_set_ui(P3->x, 0);
			mpz_set_ui(P3->y, 1);
			mpz_set_ui(P3->z, 0);
			return 0;
		}
		else
		{
			mpz_sub(*up, P2.y, P1.y);
			mpz_sub(*down, P2.x, P1.x);
			mpz_invert(*down, *down, E.p);
			mpz_mul(*l, *up, *down);
			mpz_mmod(*l, *l, E.p);
			
			mpz_mul(*up_bis, P2.y, P1.x);
			mpz_mul(*up, P2.x, P1.y);
			mpz_sub(*up, *up, *up_bis);
			mpz_mul(*v, *up, *down);
			mpz_mmod(*v, *v, E.p);
		}
	}
	
	//Compute x3
	mpz_mul(*x3, *l, *l);
	mpz_sub(*x3, *x3, P1.x);
	mpz_sub(*x3, *x3, P2.x);
	mpz_mmod(*x3, *x3, E.p);
	//Compute y3
	mpz_neg(*y3, *l);
	mpz_mul(*y3, *y3, *x3);
	mpz_sub(*y3, *y3, *v);
	mpz_mmod(*y3, *y3, E.p);
	//Initialize P3
	mpz_set(P3->x, *x3);
	mpz_set(P3->y, *y3);
	mpz_set_ui(P3->z, 1);
	return 0;
}

/** Local function used in square_and_multiply. 
 * 	Returns the square of the point P on the curve E.
 * 
 * 	@param[out]	R	The result point.
 * 	@param[in]	P	The point.
 * 	@param[in]	E	The curve.
 * 	@return		Returns 0 if the operation was successful, 1 otherwise.
 */
int _double(point_t * R, point_t P, elliptic_curve_t E)
{
	return add(R, P, P, E);
}

/**	Multiplies point with a scalar using the double-and-add method.
 * 
 * 	@param[out]	R	The point result.
 * 	@param[in]	P	The point to be multiplied.
 * 	@param[in]	s	The scalar.
 * 	@param[in]	E	The curve.
 * 	@return		Returns 0 if the operation was successful, 1 otherwise.
 */
int double_and_add(point_t *R, point_t P, mpz_t s, elliptic_curve_t E)
{
	if(P_is_on_E(P, E) == 0) 
		return 1;
		
	mpz_t *remainder, *s_cpy;
	point_t *temp;
	if(!preallocation_init_done)
	{
		preallocation_init();
	}
	remainder = &(temp_obj[13]);
	s_cpy = &(temp_obj[14]);
	temp = &(temp_point[0]);
	mpz_set(*s_cpy, s);
	mpz_set(temp->x, P.x);
	mpz_set(temp->y, P.y);
	mpz_set(temp->z, P.z);
	
	//Set result to the identity element at first
	mpz_set_ui(R->x, 0);
	mpz_set_ui(R->y, 1);
	mpz_set_ui(R->z, 0); 
	
	while(mpz_sgn(*s_cpy) != 0)
	{
		mpz_divmod_ui(*s_cpy, *remainder, *s_cpy, 2);
		if(mpz_sgn(*remainder) != 0)
		{
			add(R, *temp, *R, E);
		}
		_double(temp, *temp, E);
	}
	return 0;
}
