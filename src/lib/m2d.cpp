/** Yet Another Video decoder
 *  Copyright 2011 Takayuki Minegishi
 *
 *  Permission is hereby granted, free of charge, to any person
 *  obtaining a copy of this software and associated documentation
 *  files (the "Software"), to deal in the Software without
 *  restriction, including without limitation the rights to use, copy,
 *  modify, merge, publish, distribute, sublicense, and/or sell copies
 *  of the Software, and to permit persons to whom the Software is
 *  furnished to do so, subject to the following conditions:
 *  
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the Software.
 *  
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *  HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *  WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 *  DEALINGS IN THE SOFTWARE.
 */

#include <assert.h>
#include "m2d.h"
#include "bitio.h"

/**Decode Variable Length Code which has one argument.
 */
int m2d_dec_vld_unary(dec_bits *stream, const vlc_t *vld_tab, int bitlen)
{
	const vlc_t *code;
	int idx;
	int len;

	code = &vld_tab[show_bits(stream, bitlen)];
	len = code->length;
	idx = 0;
	while (len <= 0) {
		int rest_len;
		if (len == 0) {
			/* invalid code */
			dec_bits_tell_error(stream);
		}
		skip_bits(stream, bitlen);
		rest_len = -len;
		rest_len = rest_len < bitlen ? rest_len : bitlen;
		idx += code->pattern;
		code = &vld_tab[show_bits(stream, rest_len) + idx];
		len = code->length;
	}
	skip_bits(stream, len);
	return code->pattern;
}

/**Search "00 00 01 xx" pattern.
 */
int m2d_next_start_code(const byte_t *org_src, int byte_len)
{
	const signed char *src;

	src = (const signed char *)org_src + 2;
	byte_len = byte_len - 2;

	while (0 < byte_len) {
		int d = *src++;
		if (1 < (unsigned)d) {
			/* most probable case */
			src += 2;
			byte_len -= 3;
		} else if (d == 0) {
			byte_len -= 1;
		} else {
			/* d == 1 */
			src -= 3;
			if ((src[0] == 0) && (src[1] == 0)) {
				/* found "00 00 01" */
				src += 3;
				break;
			} else {
				src += 5;
				byte_len -= 3;
			}
		}
	}
	return (byte_len <= 0) ? -1 : (int)((const byte_t *)src - org_src);
}

/** Search start code for block(s) of input data.
 */
int m2d_find_mpeg_data(dec_bits *stream)
{
	byte_align(stream);
	show_bits(stream, 8);
	const byte_t *indata = dec_bits_current(stream);
	const byte_t *tail = dec_bits_tail(stream);
	int indata_bytes = (int)(tail - indata);
	int read_bytes;

	while ((read_bytes = m2d_next_start_code(indata, indata_bytes)) < 0) {
		int d1 = indata[indata_bytes - 2];
		int d0 = indata[indata_bytes - 1];
		indata = dec_bits_load_next(stream, &indata_bytes);
		if (indata == 0) {
			return -1;
		} else if (d1 == 0 && d0 == 0 && indata[0] == 1) {
			read_bytes = 1;
			break;
		} else if (d0 == 0 && indata[0] == 0 && indata[1] == 1) {
			read_bytes = 2;
			break;
		}
	}
	skip_bytes(stream, read_bytes);
	return 0;
}

