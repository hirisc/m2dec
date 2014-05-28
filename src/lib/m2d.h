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

#ifndef __M2D_H__
#define __M2D_H__

#include "m2types.h"
#include "bitio.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	uint8_t *luma;
	uint8_t *chroma;
	void *id;
	int32_t cnt;
	int16_t width, height;
	int16_t crop[4];
} m2d_frame_t;

typedef struct {
	int16_t pattern;
	int8_t length;
} vlc_t;

typedef struct {
	int16_t src_width, src_height;
	int16_t disp_width, disp_height;
	int16_t frame_num;
	int16_t crop[4];
	int additional_size;
} m2d_info_t;

typedef struct {
	size_t context_size;
	int (*init)(void *, int, int (*)(void *, void *), void *);
	dec_bits *(*stream_pos)(void *);
	int (*get_info)(void *, m2d_info_t *);
	int (*set_frames)(void *, int, m2d_frame_t *, uint8_t *, int);
	int (*decode_picture)(void *);
	int (*peek_decoded_frame)(void *, m2d_frame_t *, int);
	int (*get_decoded_frame)(void *, m2d_frame_t *, int);
} m2d_func_table_t;

int m2d_dec_vld_unary(dec_bits *stream, const vlc_t *vld_tab, int bitlen);
void m2d_load_bytes_skip03(dec_bits *ths, int read_bytes);
int m2d_find_mpeg_data(dec_bits *stream);
int m2d_next_start_code(const byte_t *org_src, int byte_len);

static inline uint32_t get_bits32(dec_bits *ths, int bit_len)
{
	if (bit_len <= 24) {
		return get_bits(ths, bit_len);
	} else {
		int rest = bit_len - 24;
		return (get_bits(ths, 24) << rest) | get_bits(ths, rest);
	}
}

static inline uint32_t ue_golomb(dec_bits *str)
{
	int bits, rest;
	int i;

	if (get_onebit_inline(str)) {
		return 0;
	}
	rest = 0;
	i = 16;
	do {
		bits = get_bits(str, 2);
		switch (bits) {
		case 0:
			rest += 2;
			break;
		case 1:
			return static_cast<uint32_t>(get_bits32(str, rest + 2) + ((bits << 2) << rest) - 1);
			/* NOTREACHED */
			break;
		case 2:
			/* FALLTHROUGH */
		case 3:
			return static_cast<uint32_t>((rest ? get_bits32(str, rest) : 0) + (bits << rest) - 1);
			/* NOTREACHED */
			break;
		}
	} while (--i);
	return 0;
}

static inline int32_t se_golomb(dec_bits *stream)
{
	int32_t ue = ue_golomb(stream);
	int32_t t = (ue + 1) >> 1;
	return (ue & 1) ? t : -t;
}

#ifdef __cplusplus
}
#endif

#endif /* __M2D_H__ */
