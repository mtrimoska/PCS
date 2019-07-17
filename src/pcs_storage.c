/** @file pcs_storage.c
 *
 *	Created by Monika Trimoska on 03/12/2015.
 *	Copyright © 2015 Monika Trimoska. All rights reserved.
 */

#include<gmp.h>
#include "pcs_storage.h"
#include "pcs_struct_hash.h"
#include "pcs_struct_PRTL.h"

uint8_t struct_type;

/** Initialize the distinguished-point-storing structure.
 * 
 */
void struct_init(uint8_t type, mpz_t n, uint8_t trailling_bits, uint8_t nb_bits, int nb_threads, uint8_t level)
{
    struct_type = type;
	switch(struct_type)
	{
		case 0: struct_init_PRTL(nb_bits, trailling_bits, nb_threads, level);
			break;
        default:
			struct_init_hash(struct_type, n, trailling_bits, level);
	}
}

/** Search and insert.
 *  
 *  @brief Look for a point in the structure. If the point is not found 
 *  it is added with the corresponding a coefficient.  
 *  
 *  @param[out]	a_out	The a coefficient of the found point.
 *  @param[in]	a_in	The a coefficient of the newly added point.
 *  @param[in]	xDist	The x coordinate, without the trailling zeros.
 *  @return 	1 if the point was found, 0 otherwise.
 */
int struct_add(mpz_t a_out, mpz_t a_in, mpz_t xDist, char xDist_str[])
{
	switch(struct_type)
	{
		case 0: return struct_add_PRTL(a_out, a_in, xDist);
			break;
        default:
			{mpz_get_str(xDist_str, 16, xDist); return struct_add_hash(a_out, a_in, xDist_str);}
	}
}

/** Free the distinguished-point-storing structure.
 * 
 */
void struct_free()
{
	switch(struct_type)
	{
		case 0: struct_free_PRTL();
			break;
        default: 
			struct_free_hash();
	}
}

/** Get the memory occupation of the distinguished-point-storing structure.
 *  
 *  @return 	The memory occupation in bytes.
 */
unsigned long long int struct_memory(unsigned long int *nb_points, float *rate_of_use, float *rate_slots, int nb_threads)
{
	switch(struct_type)
	{
		case 0: return struct_memory_PRTL(nb_points, rate_of_use, rate_slots);
			break;
        default:
			return struct_memory_hash(nb_points, rate_of_use, rate_slots);
	}
}
