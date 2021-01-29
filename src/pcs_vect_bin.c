/** @file pcs_vect_bin.c
 *  @brief A byte-vector implementation used for the 'packed' property of the PRTL structure.
 *
 *	Created by Gilles Dequen on 10/11/2016.
 *	Copyright © 2016 Gilles Dequen. All rights reserved.
 */

#include<assert.h>
#include<stdlib.h>
#include <stdio.h>
#include<string.h>
#include "pcs_vect_bin.h"

/// Adressing bytes (and bits) is like this in _vect_bin_t type
///      0         1        2         3         4
/// [N-------][--------][--------][--------][-------0]
/// [--------][--------][--------][--------][--------]


/// get bit at rank return _true or _false respectively to 1 and 0
inline _bool_t vect_bin_get_bit(_vect_bin_t *t, int rank) {
//  printf(" Rank(%d) : t[%lu] @ Rank(%lu)\n", rank, _vect_bin_array_size - 1 - (rank / (sizeof(_vect_bin_t) << 3)), (rank % (sizeof(_vect_bin_t) << 3)));
  return((t[_vect_bin_array_size - 1 - (rank / (sizeof(_vect_bin_t) << 3))] & ((_vect_bin_t) 1 << (rank % (sizeof(_vect_bin_t) << 3)))) ? 1 : 0);
}

/// set 'rank' bit to 1
inline void vect_bin_set_1(_vect_bin_t *t, int rank) {
  t[_vect_bin_array_size - 1 - (rank / (sizeof(_vect_bin_t) << 3))] |= ((_vect_bin_t) 1 << (rank % (sizeof(_vect_bin_t) << 3)));
}

/// set 'rank' bit to 0
inline void vect_bin_set_0(_vect_bin_t *t, int rank) {
  t[_vect_bin_array_size - 1 - (rank / (sizeof(_vect_bin_t) << 3))] &= (~((_vect_bin_t) 1 << (rank % (sizeof(_vect_bin_t) << 3))));
}

/// set an mpz value to a _vect_bin_t type starting at a given bit number
_vect_bin_t *vect_bin_set_mpz(_vect_bin_t *t, int from_bit_vect, int nb_bits, mpz_t value, int from_bit_mpz) {
  for(int i = 0; i < nb_bits; ++i)
    if(mpz_tstbit(value, from_bit_mpz + i)) vect_bin_set_1(t, i + from_bit_vect);
  return(t);
}

/// get an mpz value from a _vect_bin_t type starting at a given bit number
void vect_bin_get_mpz(_vect_bin_t *t, int from_bit, int nb_bits, mpz_t value) {
  mpz_set_ui(value, 0);
  for(int i = 0; i < nb_bits; ++i)
    if(vect_bin_get_bit(t, from_bit + i)) mpz_setbit(value, i);
}

///compare an mpz with a _vect_bin_t type starting at a given bit number 
int vect_bin_cmp_mpz(_vect_bin_t *t, int from_bit_vect, int nb_bits, mpz_t value, int from_bit_mpz) {
    for(int i = nb_bits - 1; i >= 0; --i)
        if(vect_bin_get_bit(t, from_bit_vect + i) > mpz_tstbit(value, from_bit_mpz + i)) return 1;
        else if(vect_bin_get_bit(t, from_bit_vect + i) < mpz_tstbit(value, from_bit_mpz + i)) return -1;
    return 0;
}

/// set an int value to a _vect_bin_t type starting at a given bit number
_vect_bin_t *vect_bin_set_int(_vect_bin_t *t, int from_bit, int value) {
  const int int_bitsize = (sizeof(int) << 3);
  for(int i = 0; i < int_bitsize; ++i)
    if(value & (1 << i)) vect_bin_set_1(t, i + from_bit);
  return(t);
}

/// get an int value from a _vect_bin_t type starting at a given bit number
int vect_bin_get_int(_vect_bin_t *t, int from_bit) {
  int out = 0;
  const int int_bitsize = (sizeof(int) << 3);
  for(int i = 0; i < int_bitsize; ++i)
    if(vect_bin_get_bit(t, from_bit + i)) out |= (1 << i);
  return(out);
}

void print_vect_bin(_vect_bin_t *v) {
  for(int i = 0; i < _vect_bin_array_size; ++i) printf(" %d", v[i]);
  printf("\n");
}

//check if _vect_bin_t is zero
int vect_bin_is_empty(_vect_bin_t *v)
{
    for(int i = 0; i < _vect_bin_array_size; ++i) if(v[i] != 0) return 0;
    return 1;
}

//deep clone _vect_bin_t
void vect_bin_cpy(_vect_bin_t *out, _vect_bin_t *in)
{
    memcpy(out, in, _vect_bin_array_size);
}

/// return a string that represents the vect_bin_t number in decimal format
char *vect_bin_to_binary_string(_vect_bin_t *t, char *s) {
  char *out = ((s == NULL) ? (char *)malloc(sizeof(char) * (_vect_bin_size + 1)) : s);
  for(int i = _vect_bin_size - 1; i >= 0; --i)
    out[_vect_bin_size - 1 - i] = (vect_bin_get_bit(t, i) ? '1' : '0');
  out[_vect_bin_size] = '\0';
  return(out);
}

/// _vect_bin_t assigned to 0
_vect_bin_t *vect_bin_t_reset(_vect_bin_t *_v) {
  if(_v == NULL) return(NULL);
  for(int i = 0; i < _vect_bin_array_size; ++i) _v[i] = 0;
  return _v;
}

/// return a substring from t from from_bit to to_bit. o could be NULL. It's then allocated
void vect_bin_get_vect_bin_from(_vect_bin_t *t, int from_bit, int to_bit, _vect_bin_t *o) {
  int i, j;
  for(i = from_bit, j = 0; i < to_bit; ++i, ++j)
    if(vect_bin_get_bit(t, i)) vect_bin_set_1(o, j);
    else vect_bin_set_0(o, j);
}
