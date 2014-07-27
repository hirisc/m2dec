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
	uint32_t range;
	uint32_t offset;
	int8_t* context;
} m2d_cabac_t;

typedef struct {
	int8_t m, n;
} m2d_cabac_init_mn_t;

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

static inline void init_cabac_engine(m2d_cabac_t *cb, dec_bits *st)
{
	cb->range = 0x1fe;
	cb->offset = get_bits(st, 9);
}

static inline void init_cabac_context(m2d_cabac_t *cabac, int slice_qp, const m2d_cabac_init_mn_t *lut, int lut_len)
{
	int8_t *ctx = cabac->context;
	do {
		int pre_state = ((lut->m * slice_qp) >> 4) + lut->n;
		if (pre_state < 64) {
			/* valMPS == 0 */
			pre_state = (pre_state <= 0) ? 1 : pre_state;
			*ctx = (63 - pre_state) * 2; /* -1..-63 ==> 124, 122,..0 */
		} else {
			/* valMPS == 1 */
			pre_state = (126 < pre_state) ? 126 : pre_state;
			*ctx = (pre_state - 64) * 2 + 1; /* 0..62 ==> 1, 3,..125 */
		}
		ctx++;
		lut++;
	} while (--lut_len);
}

static inline void cabac_renorm(m2d_cabac_t *cb, dec_bits *st, int range, int offset)
{
	static const uint8_t num_leading_zeros_plus1[256] = {
		9, 8, 7, 7, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
		2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
		2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
		2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	};
	uint32_t bits = num_leading_zeros_plus1[range];
	cb->range = range << bits;
	cb->offset = (offset << bits) | get_bits(st, bits);
}

static inline int cabac_decode_decision_raw(m2d_cabac_t *cb, dec_bits *st, int8_t *ctx)
{
	static const uint8_t rangeTabLPS[64][4] = {
		{128, 176, 208, 240}, {128, 167, 197, 227},
		{128, 158, 187, 216}, {123, 150, 178, 205},
		{116, 142, 169, 195}, {111, 135, 160, 185},
		{105, 128, 152, 175}, {100, 122, 144, 166},
		{95, 116, 137, 158}, {90, 110, 130, 150},
		{85, 104, 123, 142}, {81, 99, 117, 135},
		{77, 94, 111, 128}, {73, 89, 105, 122},
		{69, 85, 100, 116}, {66, 80, 95, 110},
		{62, 76, 90, 104}, {59, 72, 86, 99},
		{56, 69, 81, 94}, {53, 65, 77, 89},
		{51, 62, 73, 85}, {48, 59, 69, 80},
		{46, 56, 66, 76}, {43, 53, 63, 72},
		{41, 50, 59, 69}, {39, 48, 56, 65},
		{37, 45, 54, 62}, {35, 43, 51, 59},
		{33, 41, 48, 56}, {32, 39, 46, 53},
		{30, 37, 43, 50}, {29, 35, 41, 48},
		{27, 33, 39, 45}, {26, 31, 37, 43},
		{24, 30, 35, 41}, {23, 28, 33, 39},
		{22, 27, 32, 37}, {21, 26, 30, 35},
		{20, 24, 29, 33}, {19, 23, 27, 31},
		{18, 22, 26, 30}, {17, 21, 25, 28},
		{16, 20, 23, 27}, {15, 19, 22, 25},
		{14, 18, 21, 24}, {14, 17, 20, 23},
		{13, 16, 19, 22}, {12, 15, 18, 21},
		{12, 14, 17, 20}, {11, 14, 16, 19},
		{11, 13, 15, 18}, {10, 12, 15, 17},
		{10, 12, 14, 16}, {9 ,11, 13, 15},
		{9, 11, 12, 14}, {8, 10, 12, 14},
		{8, 9, 11, 13}, {7, 9, 11, 12},
		{7, 9, 10, 12}, {7, 8, 10, 11},
		{6, 8, 9, 11}, {6, 7, 9, 10},
		{6, 7, 8, 9}, {2, 2, 2, 2}
	};
	static const int8_t state_trans[64] = {
		1, 0, 2, 4, 4, 8, 8, 10,
		12, 14, 16, 18, 18, 22, 22, 24,
		26, 26, 30, 30, 32, 32, 36, 36,
		38, 38, 42, 42, 44, 44, 46, 48,
		48, 50, 52, 52, 54, 54, 56, 58,
		58, 60, 60, 60, 62, 64, 64, 66,
		66, 66, 68, 68, 70, 70, 70, 72,
		72, 72, 74, 74, 74, 76, 76, 126
	};
	int pStateIdx = *ctx;
	uint32_t valMPS = pStateIdx & 1;
	uint32_t range, offset, lps;
	pStateIdx >>= 1;
	range = cb->range;
	offset = cb->offset;
	lps = rangeTabLPS[pStateIdx][(range >> 6) & 3];
	range = range - lps;
	if (offset < range) {
		*ctx = ((pStateIdx + (pStateIdx < 62)) * 2) | valMPS;
		if (256 <= range) {
			cb->range = range;
			return valMPS;
		}
	} else {
		offset = offset - range;
		range = lps;
		*ctx = state_trans[pStateIdx] ^ valMPS;
		valMPS ^= 1;
	}
	cabac_renorm(cb, st, range, offset);
	return valMPS;
}

static inline int cabac_decode_multibypass(m2d_cabac_t *cb, dec_bits *st, uint32_t num)
{
	int range = cb->range;
	int offset = cb->offset;
	uint32_t bin = 0;
	assert (num <= 16);
	offset = (offset << num) | get_bits(st, num);
	do {
		bin = bin * 2;
		if (range <= (offset >> (num - 1))) {
			offset -= range << (num - 1);
			bin |= 1;
		}
	} while (--num);
	cb->offset = offset;
	return bin;
}

static inline int cabac_decode_bypass(m2d_cabac_t *cb, dec_bits *st)
{
	int range, offset;
	range = cb->range;
	offset = (cb->offset << 1) | get_onebit_inline(st);
	if (offset < range) {
		cb->offset = offset;
		return 0;
	} else {
		cb->offset = offset - range;
		return 1;
	}
}

static int header_dummyfunc(void *arg, void *seq_id) {return 0;}

#ifdef __cplusplus
}

template <typename T>
struct AddSaturate {
	T operator()(T x, T y) const {
		T msk;
		msk = ((x & y) + (((x ^ y) >> 1) & 0x7f7f7f7f7f7f7f7fULL)) & ~0x7f7f7f7f7f7f7f7fULL;
		msk = (msk << 1) - (msk >> 7);
		return ((x + y) - msk) | msk;
	}
};

template <typename T>
struct SubSaturate {
	T operator()(T x, T y) const {
		T msk;
		msk = ((~x & y) + (((~x ^ y) >> 1) & 0x7f7f7f7f7f7f7f7fULL)) & ~0x7f7f7f7f7f7f7f7fULL;
		msk = (msk << 1) - (msk >> 7);
		return (x | msk) - (y | msk);
	}
};

template<int N, int colour, typename T, typename F>
static inline void acNxNtransform_dconly_base(uint8_t *dst, T dc, int stride, F saturate)
{
	int y = N;
	dc *= static_cast<T>((colour == 0) ? 0x0101010101010101ULL : ((colour == 1) ? 0x0001000100010001ULL : 0x0100010001000100ULL));
	const int max = ((colour == 0 ? N : N * 2) / sizeof(T)) ? ((colour == 0 ? N : N * 2) / sizeof(T)) : 1;
	do {
		for (int x = 0; x < max; ++x) {
			((T*)dst)[x] = saturate(((T*)dst)[x], dc);
		}
		dst += stride;
	} while (--y);
}

template <int SHIFT>
struct IdctAdjust {
	int operator()(int dc) const {
		return (dc + (1 << (SHIFT - 1))) >> SHIFT;
	}
};

template<int N, int colour, typename T, typename F>
static void acNxNtransform_dconly_generic(uint8_t *dst, int dc, int stride, F Adjust)
{
	dc = Adjust(dc);
	if (dc < 0) {
		acNxNtransform_dconly_base<N, colour>(dst, static_cast<T>(-dc), stride, SubSaturate<T>());
	} else {
		acNxNtransform_dconly_base<N, colour>(dst, static_cast<T>(dc), stride, AddSaturate<T>());
	}
}

template<int N, int SHIFT, int colour, typename T>
static void acNxNtransform_dconly(uint8_t *dst, int dc, int stride) {
	acNxNtransform_dconly_generic<N, colour, T>(dst, dc, stride, IdctAdjust<SHIFT>());
}

#endif

#endif /* __M2D_H__ */
