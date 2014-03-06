//! \file
//! \brief Module provides Base64 encoding functionality
//!
//! Valid XML files are pure text files. One possibility to achieve this in
//! the case of binary data is to encode them. This is the identical approach
//! used in email traffic for non-ascii attachments. We employ Base64 encoding
//! for this. This type of encoding takes a tuple of 3 bytes of the binary
//! data and converts it into 4 symbols from a subset of the ASCII character
//! set. Thus, the amount of data increases by approximately a factor of 4/3.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vtkw_types.h"
#include "vtkw_encode.h"

/*
0       A       16      Q       32      g       48      w
1       B       17      R       33      h       49      x
2       C       18      S       34      i       50      y
3       D       19      T       35      j       51      z
4       E       20      U       36      k       52      0
5       F       21      V       37      l       53      1
6       G       22      W       38      m       54      2
7       H       23      X       39      n       55      3
8       I       24      Y       40      o       56      4
9       J       25      Z       41      p       57      5
10      K       26      a       42      q       58      6
11      L       27      b       43      r       59      7
12      M       28      c       44      s       60      8
13      N       29      d       45      t       61      9
14      O       30      e       46      u       62      +
15      P       31      f       47      v       63      /
*/

/* We need three masks for sextet generation */

//! Bitmask for generation of bit sextet from byte triple (OOOO OOLL)
#define MASK1 0x03

//! Bitmask for generation of bit sextet from byte triple (OOOO LLLL)
#define MASK2 0x0F

//! Bitmask for generation of bit sextet from byte triple (OOLL LLLL)
#define MASK3 0x3F


//! \brief This is the translation table for converting numerical sextet to
//! ASCII character symbol
static const unsigned char map[] = { 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
                                     'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
                                     'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
                                     'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
                                     'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
                                     'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
                                     'w', 'x', 'y', 'z', '0', '1', '2', '3',
                                     '4', '5', '6', '7', '8', '9', '+', '/' };


//! Base64 encode data from input to output buffer

//! Take data from input buffer and write them to output buffer performing
//! base64 encoding on the fly.
//! \param in        input buffer; contains data to be encoded
//! \param numBytes  number of bytes in input buffer to be encoded
//! \param out       output buffer; contains encoded data
void vtkwEncodeBase64Buffer2Buffer( unsigned char *in, Int32 numBytes,
                                    unsigned char *out ) {

  Int32 i = 0, k = 0;
  Int32 aux1, aux2, aux3, aux4;

  /* safety check */
  if ( sizeof( unsigned char ) != 1 ) {
    printf( "\n\n Severe problem: unsigned char has %d bytes!!!\n\n",
            (Int32)sizeof(unsigned char) );
  }

  /* treat complete triples of bytes first */
  for ( i = 0; i < numBytes - numBytes%3; i += 3 ) {

    /* extract first sextet */
    aux1 = in[i] >> 2;

    /* extract second sextet */
    aux2 = ( in[i] & MASK1 ) << 4 | in[i+1] >> 4;

    /* extract third sextet */
    aux3 = ( in[i+1] & MASK2 ) << 2 | in[i+2] >> 6;

    /* extract fourth sextet */
    aux4 = in[i+2] & MASK3;

    /* encode the four sextets */
    out[k  ] = map[aux1];
    out[k+1] = map[aux2];
    out[k+2] = map[aux3];
    out[k+3] = map[aux4];
    k += 4;
  }

  /* determine excess, i.e. number of extra bytes not treated */
  switch ( numBytes%3 ) {

  case 1:

    /* first sextet comes from execess byte */
    aux1 = in[i] >> 2;
    
    /* second sextet is incomplete, just take two bits from execess byte */
    aux2 = ( in[i] & MASK1 ) << 4;

    /* encode two sextets and add two '='s to mark the missing bytes */
    out[k  ] = map[aux1];
    out[k+1] = map[aux2];
    out[k+2] = '=';
    out[k+3] = '=';
    break;

  case 2:

    /* first sextet comes from excess byte 1 */
    aux1 = in[i] >> 2;
    
    /* second sextet comes from excess bytes 1 and 2 */
    aux2 = ( in[i] & MASK1 ) << 4 | in[i+1] >> 4;

    /* third sextet is incomplete, just take bits from excess byte 2 */
    aux3 = ( in[i+1] & MASK2 ) << 2;

    /* encode three sextets and add one '=' to mark missing byte */
    out[k  ] = map[aux1];
    out[k+1] = map[aux2];
    out[k+2] = map[aux3];
    out[k+3] = '=';
    break;
  }

}
