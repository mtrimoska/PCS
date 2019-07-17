/** @file pcs_struct_hash_UNIX.c
 *
 *	Created by Monika Trimoska on 03/12/2015.
 *	Copyright © 2015 Monika Trimoska. All rights reserved.
 */

#include <gmp.h>
#include <string.h>
#include <inttypes.h>

/** Implementation of the ElfHash function.
 *
 *  @brief This hash function is used in the UNIX ELF format for object files.
 *
 *  @param[in]	xDist	x-coordinate of a distinguished point.
 *  @return 	Hash value.
 */
unsigned long int get_hash_UNIX(char *xDist)
{
	unsigned long int hash_val = 0;
	unsigned long int hi_bits = 0;
	uint8_t pos;
	for(pos = 0; pos < strlen(xDist); pos++)
	{
		hash_val = (hash_val << 4) + xDist[pos];
		hi_bits = hash_val & 0xF0000000;
		if(hi_bits !=0)
		{
			hash_val ^= hi_bits >> 24;
		}
		hash_val &= ~hi_bits;
	}
	return hash_val;
}
