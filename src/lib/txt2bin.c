
#include <assert.h>
#include "txt2bin.h"

/**fill with "1" or "0" until byte-aligned.
 */
static int fill_bit( int data, int residual_bit, int nonzero )
{
	unsigned int rest_len = 8 - residual_bit;
	data = data << rest_len;
	if ( nonzero ) {
		int mask = 1;
		int i = rest_len;

		while ( i-- ) {
			data |= mask;
			mask <<= 1;
		}
	}
	return data;
}

/**Convert binary digits in text expression to raw data.
 *'0', '1', and blank letter are allowed.
 *\return -1 on error. Otherwise, size of raw data in bit units.
 */
int txt2bin( const char *src, unsigned char *dst )
{
	unsigned int bitlen;
	unsigned int residual;
	int c;
	int d;

	d = 0;
	bitlen = 0;
	residual = 0;
	while ( ( c = *src++ ) != '\0' ) {
		if ( c == '0' ) {
			d <<= 1;
			bitlen++;
		} else if ( c == '1' ) {
			d = ( d << 1 ) | 1;
			bitlen++;
		} else if ( ( c == ' ' ) || ( c == '\t' ) ) {
			continue;
		} else if ( c == 'p' ) {
			/* fill '1' until next boundary */
			if ( ( bitlen & 7 ) == 0 ) {
				continue;
			}
			d = fill_bit( d, residual, 1 );
			bitlen = ( bitlen + 7 ) & ~7;
		} else if ( c == 'z' ) {
			/* fill '0' until next boundary */
			if ( ( bitlen & 7 ) == 0 ) {
				continue;
			}
			d = fill_bit( d, residual, 0 );
			bitlen = ( bitlen + 7 ) & ~7;
		} else {
			/* error */
			assert( 1 );
			return -1;
		}
		residual = bitlen & 7;
		if ( residual == 0 ) {
			*dst++ = d;
		}
	}
	if ( residual ) {
		unsigned int rest_len = 8 - residual;
		*dst = d << rest_len;
/*		bitlen += rest_len; */
	}
	return bitlen;
}

