/** Yet Another H.265 decoder
 *  Copyright 2015 Takayuki Minegishi
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

#if (defined(__GNUC__) && defined(__SSE2__)) || defined(_M_IX86) || defined(_M_AMD64)
#define X86ASM
#include <emmintrin.h>

#include <assert.h>
#include <string.h>
#include <setjmp.h>
#include <limits.h>
#include <cstdlib>
#include <algorithm>
#include "h265.h"
#include "m2d_macro.h"

#if defined(_M_IX86) || defined(_M_AMD64)
#define ALIGNVC(x) __declspec(align(x))
#else
#define ALIGNVC(x) __attribute__((aligned(x)))

#endif

ALIGNVC(16) static const int rounding[4] = {
	64, 64, 64, 64
};

ALIGNVC(16) static const int ones32[4] = {
	1, 1, 1, 1
};

ALIGNVC(16) static const uint8_t ofs[16] = {
	128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128
};

ALIGNVC(16) static const int zeros[] = {
	0, 0, 0, 0
};

static inline void transform4x4vert(const int16_t* src0, int gap, const int16_t coef[][8], __m128i& out0, __m128i& out1, __m128i& out2, __m128i& out3) {
	__m128i a0 = _mm_loadu_si128((const __m128i*)src0);
	__m128i a2 = _mm_loadu_si128((const __m128i*)(src0 + gap * 2));
	__m128i b0 = _mm_unpacklo_epi16(a0, a2);
	__m128i b2 = _mm_unpackhi_epi16(a0, a2);
	out0 = _mm_add_epi32(_mm_madd_epi16(b0, *(const __m128i*)coef[0]), _mm_madd_epi16(b2, *(const __m128i*)coef[1]));
	out1 = _mm_add_epi32(_mm_madd_epi16(b0, *(const __m128i*)coef[2]), _mm_madd_epi16(b2, *(const __m128i*)coef[3]));
	out2 = _mm_add_epi32(_mm_madd_epi16(b0, *(const __m128i*)coef[4]), _mm_madd_epi16(b2, *(const __m128i*)coef[5]));
	out3 = _mm_add_epi32(_mm_madd_epi16(b0, *(const __m128i*)coef[6]), _mm_madd_epi16(b2, *(const __m128i*)coef[7]));
}

static inline void transform4x4horiz_half(const int16_t coef[][8], const __m128i& src, __m128i& dst0, __m128i& dst1) {
	__m128i a0 = _mm_unpacklo_epi32(src, src);
	__m128i a1 = _mm_shuffle_epi32(src, 0x05);
	__m128i a2 = _mm_unpackhi_epi32(src, src);
	__m128i a3 = _mm_shuffle_epi32(src, 0xaf);
	__m128i b0 = _mm_madd_epi16(a0, *(const __m128i*)*((coef) + 0));
	__m128i b1 = _mm_madd_epi16(a1, *(const __m128i*)*((coef) + 1));
	__m128i b2 = _mm_madd_epi16(a2, *(const __m128i*)*((coef) + 0));
	__m128i b3 = _mm_madd_epi16(a3, *(const __m128i*)*((coef) + 1));
	dst0 = _mm_add_epi32(b0, b1);
	dst1 = _mm_add_epi32(b2, b3);
}

static inline void transform4x4horiz(const int16_t coef[][8], const __m128i& src0, const __m128i& src1, const __m128i& src2, const __m128i& src3, __m128i& dst0, __m128i& dst1, __m128i& dst2, __m128i& dst3) {
	__m128i p0 = _mm_srai_epi32(_mm_add_epi32(src0, *(const __m128i*)rounding), 7);
	__m128i p1 = _mm_srai_epi32(_mm_add_epi32(src1, *(const __m128i*)rounding), 7);
	__m128i p2 = _mm_srai_epi32(_mm_add_epi32(src2, *(const __m128i*)rounding), 7);
	__m128i p3 = _mm_srai_epi32(_mm_add_epi32(src3, *(const __m128i*)rounding), 7);
	__m128i l0 = _mm_packs_epi32(p0, p1);
	__m128i l2 = _mm_packs_epi32(p2, p3);
	transform4x4horiz_half(coef, l0, dst0, dst1);
	transform4x4horiz_half(coef, l2, dst2, dst3);
}

static inline void transform4x4core(const int16_t* in, const int16_t coeffs[][8], int bitdepth, __m128i& dst0, __m128i& dst1, __m128i& dst2, __m128i& dst3) {
	__m128i d0, d1, d2, d3;
	__m128i o0, o1, o2, o3;
	__m128i rnd;
	int shift;
	transform4x4vert(in, 4, coeffs, d0, d1, d2, d3);
	transform4x4horiz(coeffs + 8, d0, d1, d2, d3, o0, o1, o2, o3);
	shift = 20 - bitdepth;
	rnd = _mm_slli_epi32(_mm_load_si128((const __m128i*)ones32), shift - 1);
	dst0 = _mm_srai_epi32(_mm_add_epi32(o0, rnd), shift);
	dst1 = _mm_srai_epi32(_mm_add_epi32(o1, rnd), shift);
	dst2 = _mm_srai_epi32(_mm_add_epi32(o2, rnd), shift);
	dst3 = _mm_srai_epi32(_mm_add_epi32(o3, rnd), shift);
}

static inline void store4x4(uint8_t* dst, int stride, const __m128i& d) {
	*(uint32_t*)dst = _mm_cvtsi128_si32(d);
	*(uint32_t*)(dst + stride) = _mm_cvtsi128_si32(_mm_srli_si128(d, 4));
	*(uint32_t*)(dst + stride * 2) = _mm_cvtsi128_si32(_mm_srli_si128(d, 8));
	*(uint32_t*)(dst + stride * 3) = _mm_cvtsi128_si32(_mm_srli_si128(d, 12));
}

static inline void transform4x4(uint8_t* dst, const int16_t* coeff, int stride, const int16_t matrices[][8]) {
	__m128i p0, p1, p2, p3;

	transform4x4core(coeff, matrices, 8, p0, p1, p2, p3);
	__m128i s0 = _mm_set_epi32(*(uint32_t*)(dst + stride * 3), *(uint32_t*)(dst + stride * 2), *(uint32_t*)(dst + stride), *(uint32_t*)dst);
	s0 = _mm_sub_epi8(s0, *(const __m128i*)ofs);
	__m128i d0 = _mm_packs_epi16(_mm_packs_epi32(p0, p1), _mm_packs_epi32(p2, p3));
	__m128i e0 = _mm_adds_epi8(s0, d0);
	e0 = _mm_add_epi8(e0, *(const __m128i*)ofs);
	store4x4(dst, stride, e0);
}

void transformdst_ac4x4(uint8_t* dst, int16_t* coeff, int stride) {
	ALIGNVC(16) static const int16_t coeffs[][8] = {
		{29, 84, 29, 84, 29, 84, 29, 84},
		{74, 55, 74, 55, 74, 55, 74, 55},
		{55, -29, 55, -29, 55, -29, 55, -29},
		{74, -84, 74, -84, 74, -84, 74, -84},
		{74, -74, 74, -74, 74, -74, 74, -74},
		{0, 74, 0, 74, 0, 74, 0, 74},
		{84, 55, 84, 55, 84, 55, 84, 55},
		{-74, -29, -74, -29, -74, -29, -74, -29},

		{29, 74, 55, 74, -74, 74, 55, -29},
		{84, 55, -29, -84, 74, 0, 84, -74},
	};
	transform4x4(dst, coeff, stride, coeffs);
}

ALIGNVC(16) static const signed short idct8coeffs[][8] = {
	{89, 50, 89, 50, 89, 50, 89, 50},
	{75, 18, 75, 18, 75, 18, 75, 18},
	{75, -89, 75, -89, 75, -89, 75, -89},
	{-18, -50, -18, -50, -18, -50, -18, -50},
	{50, 18, 50, 18, 50, 18, 50, 18},
	{-89, 75, -89, 75, -89, 75, -89, 75},
	{18, 75, 18, 75, 18, 75, 18, 75},
	{-50, -89, -50, -89, -50, -89, -50, -89},
};

static inline void idct4x4vert4(const __m128i& t0, const __m128i& t1, __m128i& out0, __m128i& out1, __m128i& out2, __m128i& out3) {
	ALIGNVC(16) static const int16_t coef[][8] = {
		{64, 64, 64, 64, 64, 64, 64, 64},
		{83, 36, 83, 36, 83, 36, 83, 36},
		{64, -64, 64, -64, 64, -64, 64, -64},
		{36, -83, 36, -83, 36, -83, 36, -83},
	};
	__m128i a0 = _mm_madd_epi16(t0, *(const __m128i*)coef[0]);
	__m128i b0 = _mm_madd_epi16(t1, *(const __m128i*)coef[1]);
	__m128i a1 = _mm_madd_epi16(t0, *(const __m128i*)coef[2]);
	__m128i b1 = _mm_madd_epi16(t1, *(const __m128i*)coef[3]);
	out0 = _mm_add_epi32(a0, b0);
	out1 = _mm_add_epi32(a1, b1);
	out2 = _mm_sub_epi32(a1, b1);
	out3 = _mm_sub_epi32(a0, b0);
}

static inline void idct4x4horiz_half(const __m128i& i0, const __m128i& i1, __m128i& w0, __m128i& w1, int shift, const __m128i& rnd) {
	ALIGNVC(16) static const signed short coef[][8] = {
		{64, 83, -64, -83, 64, 83, -64, -83},
		{64, -36, 64, -36, 64, -36, 64, -36},
		{64, 36, 64, 36, 64, 36, 64, 36},
		{-64, 83, 64, -83, -64, 83, 64, -83},
	};
	__m128i a0 = _mm_madd_epi16(i0, *(const __m128i*)coef[0]);
	__m128i a2 = _mm_madd_epi16(i0, *(const __m128i*)coef[1]);
	__m128i a1 = _mm_madd_epi16(i1, *(const __m128i*)coef[2]);
	__m128i a3 = _mm_madd_epi16(i1, *(const __m128i*)coef[3]);
	__m128i bl = _mm_add_epi32(a0, a1);
	__m128i bh = _mm_add_epi32(a2, a3);
	__m128i cl = _mm_add_epi32(bl, rnd);
	__m128i ch = _mm_add_epi32(bh, rnd);
	__m128i dl = _mm_srai_epi32(cl, shift);
	__m128i dh = _mm_srai_epi32(ch, shift);
	w0 = _mm_unpacklo_epi64(dl, dh);
	w1 = _mm_unpackhi_epi64(dl, dh);
}

static inline void idct4x4horiz(int bitdepth, const __m128i& i0, const __m128i& i2, __m128i& w0, __m128i& w1, __m128i& w2, __m128i& w3) {
	__m128i i1 = _mm_shuffle_epi32(i0, 0xb1); // 2, 3, 0, 1, 6, 7, 4, 5
	__m128i i3 = _mm_shuffle_epi32(i2, 0xb1);
	int shift = 20 - bitdepth;
	__m128i rnd = _mm_slli_epi32(_mm_load_si128((const __m128i*)ones32), shift - 1);
	idct4x4horiz_half(i0, i1, w0, w1, shift, rnd);
	idct4x4horiz_half(i2, i3, w2, w3, shift, rnd);
}

void transform_ac4x4(uint8_t* dst, int16_t* coeff, int stride) {
	__m128i s0 = _mm_loadu_si128((const __m128i*)coeff);
	__m128i s2 = _mm_loadu_si128((const __m128i*)(coeff + 8));
	__m128i t0 = _mm_unpacklo_epi16(s0, s2);
	__m128i t2 = _mm_unpackhi_epi16(s0, s2);
	__m128i u0, u1, u2, u3;
	idct4x4vert4(t0, t2, u0, u1, u2, u3);
	u0 = _mm_srai_epi32(_mm_add_epi32(u0, *(const __m128i*)rounding), 7);
	u1 = _mm_srai_epi32(_mm_add_epi32(u1, *(const __m128i*)rounding), 7);
	u2 = _mm_srai_epi32(_mm_add_epi32(u2, *(const __m128i*)rounding), 7);
	u3 = _mm_srai_epi32(_mm_add_epi32(u3, *(const __m128i*)rounding), 7);
	__m128i v0 = _mm_packs_epi32(u0, u1);
	__m128i v2 = _mm_packs_epi32(u2, u3);
	__m128i w0, w1, w2, w3;
	idct4x4horiz(8, v0, v2, w0, w1, w2, w3);
	__m128i x0 = _mm_sub_epi8(_mm_set_epi32(*(uint32_t*)(dst + stride * 3), *(uint32_t*)(dst + stride * 2), *(uint32_t*)(dst + stride), *(uint32_t*)dst), *(const __m128i*)ofs);
	__m128i y0 = _mm_packs_epi16(_mm_packs_epi32(w0, w1), _mm_packs_epi32(w2, w3));
	__m128i e0 = _mm_add_epi8(_mm_adds_epi8(y0, x0), *(const __m128i*)ofs);
	store4x4(dst, stride, e0);
}

static inline void idct8x8vertodd4(const __m128i& t0, const __m128i& t1, __m128i& out0, __m128i& out1, __m128i& out2, __m128i& out3) {
	out0 = _mm_add_epi32(_mm_madd_epi16(t0, *(const __m128i*)idct8coeffs[0]), _mm_madd_epi16(t1, *(const __m128i*)idct8coeffs[1]));
	out1 = _mm_add_epi32(_mm_madd_epi16(t0, *(const __m128i*)idct8coeffs[2]), _mm_madd_epi16(t1, *(const __m128i*)idct8coeffs[3]));
	out2 = _mm_add_epi32(_mm_madd_epi16(t0, *(const __m128i*)idct8coeffs[4]), _mm_madd_epi16(t1, *(const __m128i*)idct8coeffs[5]));
	out3 = _mm_add_epi32(_mm_madd_epi16(t0, *(const __m128i*)idct8coeffs[6]), _mm_madd_epi16(t1, *(const __m128i*)idct8coeffs[7]));
}

static inline void idct8x8horiz1core(const __m128i& m0, __m128i& dst0, __m128i& dst1) {
	/*
	64 * s0 + 83 * s2 + 64 * s4 + 36 * s6
	64 * s0 + 36 * s2 - 64 * s4 - 83 * s6
	64 * s0 - 36 * s2 - 64 * s4 + 83 * s6
	64 * s0 - 83 * s2 + 64 * s4 - 36 * s6
         [89, 75, 50, 18],
         [75, -18, -89, -50],
         [50, -89, 18, 75],
         [18, -50, 75, -89],
	 */
	ALIGNVC(16) static const int16_t coeff[][8] = {
		{64, 83, -64, -83, 64, -36, 64, -36},
		{64, 36, 64, 36, -64, 83, 64, -83},

		{89, 75, -89, -50, 50, -89, 75, -89},
		{50, 18, 75, -18, 18, 75, 18, -50}
	};
	__m128i m1 = _mm_shuffle_epi32(m0, 0xb1); // 4, 6, 0, 2, 5, 7, 1, 3
	__m128i even = _mm_add_epi32(_mm_madd_epi16(_mm_unpacklo_epi64(m0, m0), *(const __m128i*)coeff[0]), _mm_madd_epi16(_mm_unpacklo_epi64(m1, m1), *(const __m128i*)coeff[1]));
	__m128i odd = _mm_add_epi32(_mm_madd_epi16(_mm_unpackhi_epi64(m0, m0), *(const __m128i*)coeff[2]), _mm_madd_epi16(_mm_unpackhi_epi64(m1, m1), *(const __m128i*)coeff[3]));
	dst1 = _mm_shuffle_epi32(_mm_sub_epi32(even, odd), 0x1b);
	dst0 = _mm_add_epi32(even, odd);
}

static inline __m128i idct8x8preswap(const __m128i& d0l, const __m128i& d0h) {
	__m128i d0 = _mm_unpacklo_epi32(d0l, d0h); // 0, 4, 1, 5
	__m128i d1 = _mm_unpackhi_epi32(d0l, d0h); // 2, 6, 3, 7
	__m128i e0 = _mm_srai_epi32(_mm_add_epi32(_mm_unpacklo_epi32(d0, d1), *(const __m128i*)rounding), 7); // 0, 2, 4, 6
	__m128i e1 = _mm_srai_epi32(_mm_add_epi32(_mm_unpackhi_epi32(d0, d1), *(const __m128i*)rounding), 7); // 1, 3, 5, 7
	return _mm_packs_epi32(e0, e1); // 0, 2, 4, 6, 1, 3, 5, 7
}

static inline void idct8x8horiz1(uint8_t* dst, int shift, __m128i rnd, const __m128i& d0l, const __m128i& d0h) {
	__m128i a0, a1;
	idct8x8horiz1core(idct8x8preswap(d0l, d0h), a0, a1);
//	__m128i s0 = _mm_sub_epi8(_mm_loadu_si128((const __m128i*)dst), *(const __m128i*)ofs);
	__m128i s0 = _mm_sub_epi8(_mm_loadl_epi64((const __m128i*)dst), *(const __m128i*)ofs);
	__m128i t0 = _mm_packs_epi16(_mm_packs_epi32(_mm_srai_epi32(_mm_add_epi32(a0, rnd), shift), _mm_srai_epi32(_mm_add_epi32(a1, rnd), shift)), *(const __m128i*)zeros);
	__m128i u0 = _mm_add_epi8(_mm_adds_epi8(s0, t0), *(const __m128i*)ofs);
//	_mm_storeu_si128((__m128i*)dst, u0);
	_mm_storel_epi64((__m128i*)dst, u0);
}

static inline void transform_ac8x8bit(uint8_t* dst, const int16_t* src, int stride, int bitdepth) {
	__m128i s0 = _mm_loadu_si128((const __m128i*)src);
	__m128i s1 = _mm_loadu_si128((const __m128i*)(src + 8 * 2));
	__m128i s2 = _mm_loadu_si128((const __m128i*)(src + 8 * 4));
	__m128i s3 = _mm_loadu_si128((const __m128i*)(src + 8 * 6));
	__m128i e0l, e1l, e2l, e3l;
	__m128i o0l, o1l, o2l, o3l;
	__m128i e0h, e1h, e2h, e3h;
	__m128i o0h, o1h, o2h, o3h;
	idct4x4vert4(_mm_unpacklo_epi16(s0, s2), _mm_unpacklo_epi16(s1, s3), e0l, e1l, e2l, e3l);
	__m128i t0 = _mm_loadu_si128((const __m128i*)(src + 8));
	__m128i t1 = _mm_loadu_si128((const __m128i*)(src + 8 * 3));
	__m128i t2 = _mm_loadu_si128((const __m128i*)(src + 8 * 5));
	__m128i t3 = _mm_loadu_si128((const __m128i*)(src + 8 * 7));
	idct8x8vertodd4(_mm_unpacklo_epi16(t0, t2), _mm_unpacklo_epi16(t1, t3), o0l, o1l, o2l, o3l);
	__m128i u0l = _mm_add_epi32(e0l, o0l);
	__m128i u7l = _mm_sub_epi32(e0l, o0l);
	__m128i u1l = _mm_add_epi32(e1l, o1l);
	__m128i u6l = _mm_sub_epi32(e1l, o1l);
	__m128i u2l = _mm_add_epi32(e2l, o2l);
	__m128i u5l = _mm_sub_epi32(e2l, o2l);
	__m128i u3l = _mm_add_epi32(e3l, o3l);
	__m128i u4l = _mm_sub_epi32(e3l, o3l);
	idct4x4vert4(_mm_unpackhi_epi16(s0, s2), _mm_unpackhi_epi16(s1, s3), e0h, e1h, e2h, e3h);
	idct8x8vertodd4(_mm_unpackhi_epi16(t0, t2), _mm_unpackhi_epi16(t1, t3), o0h, o1h, o2h, o3h);
	int shift = 20 - bitdepth;
	__m128i rnd = _mm_slli_epi32(_mm_load_si128((const __m128i*)ones32), shift - 1);
	__m128i u0h = _mm_add_epi32(e0h, o0h);
	idct8x8horiz1(dst, shift, rnd, u0l, u0h);
	__m128i u7h = _mm_sub_epi32(e0h, o0h);
	idct8x8horiz1(dst + stride * 7, shift, rnd, u7l, u7h);
	__m128i u1h = _mm_add_epi32(e1h, o1h);
	idct8x8horiz1(dst + stride, shift, rnd, u1l, u1h);
	__m128i u6h = _mm_sub_epi32(e1h, o1h);
	idct8x8horiz1(dst + stride * 6, shift, rnd, u6l, u6h);
	__m128i u2h = _mm_add_epi32(e2h, o2h);
	idct8x8horiz1(dst + stride * 2, shift, rnd, u2l, u2h);
	__m128i u5h = _mm_sub_epi32(e2h, o2h);
	idct8x8horiz1(dst + stride * 5, shift, rnd, u5l, u5h);
	__m128i u3h = _mm_add_epi32(e3h, o3h);
	idct8x8horiz1(dst + stride * 3, shift, rnd, u3l, u3h);
	__m128i u4h = _mm_sub_epi32(e3h, o3h);
	idct8x8horiz1(dst + stride * 4, shift, rnd, u4l, u4h);
}

void transform_ac8x8(uint8_t* dst, int16_t* src, int stride) {
	transform_ac8x8bit(dst, src, stride, 8);
}

ALIGNVC(16) static const int16_t idct16coeffs[][8] = {
	{90, 87, 90, 87, 90, 87, 90, 87},
	{80, 70, 80, 70, 80, 70, 80, 70},
	{57, 43, 57, 43, 57, 43, 57, 43},
	{25, 9, 25, 9, 25, 9, 25, 9},
	{87, 57, 87, 57, 87, 57, 87, 57},
	{9, -43, 9, -43, 9, -43, 9, -43},
	{-80, -90, -80, -90, -80, -90, -80, -90},
	{-70, -25, -70, -25, -70, -25, -70, -25},
	{80, 9, 80, 9, 80, 9, 80, 9},
	{-70, -87, -70, -87, -70, -87, -70, -87},
	{-25, 57, -25, 57, -25, 57, -25, 57},
	{90, 43, 90, 43, 90, 43, 90, 43},
	{70, -43, 70, -43, 70, -43, 70, -43},
	{-87, 9, -87, 9, -87, 9, -87, 9},
	{90, 25, 90, 25, 90, 25, 90, 25},
	{-80, -57, -80, -57, -80, -57, -80, -57},
	{57, -80, 57, -80, 57, -80, 57, -80},
	{-25, 90, -25, 90, -25, 90, -25, 90},
	{-9, -87, -9, -87, -9, -87, -9, -87},
	{43, 70, 43, 70, 43, 70, 43, 70},
	{43, -90, 43, -90, 43, -90, 43, -90},
	{57, 25, 57, 25, 57, 25, 57, 25},
	{-87, 70, -87, 70, -87, 70, -87, 70},
	{9, -80, 9, -80, 9, -80, 9, -80},
	{25, -70, 25, -70, 25, -70, 25, -70},
	{90, -80, 90, -80, 90, -80, 90, -80},
	{43, 9, 43, 9, 43, 9, 43, 9},
	{-57, 87, -57, 87, -57, 87, -57, 87},
	{9, -25, 9, -25, 9, -25, 9, -25},
	{43, -57, 43, -57, 43, -57, 43, -57},
	{70, -80, 70, -80, 70, -80, 70, -80},
	{87, -90, 87, -90, 87, -90, 87, -90}
};

static inline void idct16x16horizodd(const __m128i& r0, __m128i& dst0, __m128i& dst1) {
	ALIGNVC(16) static const int16_t coef[8][8] = {
/*          [90, 87, 80, 70, 57, 43, 25, 9],
          [87, 57, 9, -43, -80, -90, -70, -25],
          [80, 9, -70, -87, -25, 57, 90, 43],
          [70, -43, -87, 9, 90, 25, -80, -57],
          [57, -80, -25, 90, -9, -87, 43, 70],
          [43, -90, 57, 25, -87, 70, 9, -80],
          [25, -70, 90, -80, 43, 9, -57, 87],
          [9, -25, 43, -57, 70, -80, 87, -90]
*/
		{90, 87, 9, -43, -25, 57, -80, -57},
		{57, -80, 57, 25, 43, 9, 87, -90},
		{25, 9, 87, 57, -70, -87, 90, 25},
		{43, 70, 43, -90, 90, -80, 70, -80},
		{57, 43, -70, -25, 80, 9, -87, 9},
		{-9, -87, 9, -80, 25, -70, 43, -57},
		{80, 70, -80, -90, 90, 43, 70, -43},
		{-25, 90, -87, 70, -57, 87, 9, -25}
	};
	__m128i r1 = _mm_shuffle_epi32(r0, 0x93); // 13, 15, 1, 3, 5, 7, 9, 11
	__m128i r2 = _mm_shuffle_epi32(r0, 0x4e); // 9, 11, 13, 15, 1, 3, 5, 7
	__m128i r3 = _mm_shuffle_epi32(r0, 0x39); // 5, 7, 9, 11, 13, 15, 1, 3
	__m128i d0 = _mm_madd_epi16(r0, *(const __m128i*)coef[0]);
	__m128i e0 = _mm_madd_epi16(r0, *(const __m128i*)coef[1]);
	__m128i d1 = _mm_madd_epi16(r1, *(const __m128i*)coef[2]);
	__m128i e1 = _mm_madd_epi16(r1, *(const __m128i*)coef[3]);
	__m128i d2 = _mm_madd_epi16(r2, *(const __m128i*)coef[4]);
	__m128i e2 = _mm_madd_epi16(r2, *(const __m128i*)coef[5]);
	__m128i d3 = _mm_madd_epi16(r3, *(const __m128i*)coef[6]);
	__m128i e3 = _mm_madd_epi16(r3, *(const __m128i*)coef[7]);
	dst0 = _mm_add_epi32(_mm_add_epi32(d0, d1), _mm_add_epi32(d2, d3));
	dst1 = _mm_add_epi32(_mm_add_epi32(e0, e1), _mm_add_epi32(e2, e3));
}

static inline void idct16x16horizcore(int shift, const __m128i& v0, const __m128i& v1, __m128i& a0, __m128i& a1, __m128i& a2, __m128i& a3) {
	__m128i even0, even1;
	idct8x8horiz1core(v0, even0, even1);
	__m128i rnd = _mm_slli_epi32(_mm_load_si128((const __m128i*)ones32), shift - 1);
	even0 = _mm_add_epi32(even0, rnd);
	even1 = _mm_add_epi32(even1, rnd);
	__m128i odd0, odd1;
	idct16x16horizodd(v1, odd0, odd1);
	a0 = _mm_add_epi32(even0, odd0);
	a1 = _mm_add_epi32(even1, odd1);
	a2 = _mm_shuffle_epi32(_mm_sub_epi32(even1, odd1), 0x1b);
	a3 = _mm_shuffle_epi32(_mm_sub_epi32(even0, odd0), 0x1b);
}

static inline __m128i round_pack16(int shift, const __m128i& d0, const __m128i& d1, const __m128i& d2, const __m128i& d3) {
	return _mm_packs_epi16(_mm_packs_epi32(_mm_srai_epi32(d0, shift), _mm_srai_epi32(d1, shift)), _mm_packs_epi32(_mm_srai_epi32(d2, shift), _mm_srai_epi32(d3, shift)));
}

static inline void idctwrite16(uint8_t* dst, int shift, const __m128i& d0, const __m128i& d1, const __m128i& d2, const __m128i& d3) {
	__m128i s0 = _mm_sub_epi8(_mm_loadu_si128((const __m128i*)dst), *(const __m128i*)ofs);
	__m128i f0 = round_pack16(shift, d0, d1, d2, d3);
	__m128i g0 = _mm_add_epi8(_mm_adds_epi8(s0, f0), *(const __m128i*)ofs);
	_mm_storeu_si128((__m128i*)dst, g0);
}

static inline void idctwrite16chroma(uint8_t* dst, int shift, const __m128i& d0, const __m128i& d1, const __m128i& d2, const __m128i& d3) {
	__m128i s0 = _mm_sub_epi8(_mm_loadu_si128((const __m128i*)dst), *(const __m128i*)ofs);
	__m128i s1 = _mm_sub_epi8(_mm_loadu_si128((const __m128i*)(dst + 16)), *(const __m128i*)ofs);
	__m128i f0 = round_pack16(shift, d0, d1, d2, d3);
	__m128i f0l = _mm_unpacklo_epi8(f0, *(const __m128i*)zeros);
	__m128i f0h = _mm_unpackhi_epi8(f0, *(const __m128i*)zeros);
	__m128i g0 = _mm_add_epi8(_mm_adds_epi8(s0, f0l), *(const __m128i*)ofs);
	__m128i g1 = _mm_add_epi8(_mm_adds_epi8(s1, f0h), *(const __m128i*)ofs);
	_mm_storeu_si128((__m128i*)dst, g0);
	_mm_storeu_si128((__m128i*)(dst + 16), g1);
}

struct idct16x16write {
	void operator()(uint8_t* dst, int shift, const __m128i& d0, const __m128i& d1, const __m128i& d2, const __m128i& d3) const {
		idctwrite16(dst, shift, d0, d1, d2, d3);
	}
};

struct idct16x16writechroma {
	void operator()(uint8_t* dst, int shift, const __m128i& d0, const __m128i& d1, const __m128i& d2, const __m128i& d3) const {
		idctwrite16chroma(dst, shift, d0, d1, d2, d3);
	}
};

template <typename F>
static inline void idct16x16horiz(const int16_t* in, uint8_t* dst, int bitdepth, F Write) {
	__m128i s0 = _mm_loadu_si128((const __m128i*)(in + 0));
	__m128i s1 = _mm_loadu_si128((const __m128i*)(in + 8));
	__m128i t0 = _mm_unpacklo_epi16(s0, s1); // 0, 8, 1, 9, 2, 10, 3, 11
	__m128i t1 = _mm_unpackhi_epi16(s0, s1); // 4, 12, 5, 13, 6, 14, 7, 15
	__m128i u0 = _mm_unpacklo_epi16(t0, t1); // 0, 4, 8, 12, 1, 5, 9, 13
	__m128i u1 = _mm_unpackhi_epi16(t0, t1); // 2, 6, 10, 14, 3, 7, 11, 15
	__m128i v0 = _mm_unpacklo_epi64(u0, u1); // 0, 4, 8, 12, 2, 6, 10, 14
	__m128i v1 = _mm_unpackhi_epi16(u0, u1); // 1, 3, 5, 7, 9, 11, 13, 15
	int shift = 20 - bitdepth;
	__m128i a0, a1, a2, a3;
	idct16x16horizcore(shift, v0, v1, a0, a1, a2, a3);
	Write(dst, shift, a0, a1, a2, a3);
}

static inline __m128i idct16x16vertodd4col(const int16_t coef[][8], const __m128i& a0, const __m128i& a1, const __m128i& a2, const __m128i& a3) {
	__m128i b0 = _mm_add_epi32(_mm_madd_epi16(a0, *(const __m128i*)coef[0]), _mm_madd_epi16(a1, *(const __m128i*)coef[1]));
	__m128i b1 = _mm_add_epi32(_mm_madd_epi16(a2, *(const __m128i*)coef[2]), _mm_madd_epi16(a3, *(const __m128i*)coef[3]));
	return _mm_add_epi32(b0, b1);
}

static inline void idct16x16vertoddhalf(const int16_t coef[][8], __m128i dst[], const __m128i& l0, const __m128i& l1, const __m128i& l2, const __m128i& l3) {
	for (int i = 0; i < 8; ++i) {
		_mm_store_si128(dst + i, idct16x16vertodd4col(coef + i * 4, l0, l1, l2, l3));
	}
}

static inline void idct16x16vertodd(const int16_t* src, int gap, const int16_t coef[][8], __m128i odd[]) {
	__m128i i0 = _mm_loadu_si128((const __m128i*)(src + gap * 1));
	__m128i i1 = _mm_loadu_si128((const __m128i*)(src + gap * 3));
	__m128i i2 = _mm_loadu_si128((const __m128i*)(src + gap * 5));
	__m128i i3 = _mm_loadu_si128((const __m128i*)(src + gap * 7));
	__m128i i4 = _mm_loadu_si128((const __m128i*)(src + gap * 9));
	__m128i i5 = _mm_loadu_si128((const __m128i*)(src + gap * 11));
	__m128i i6 = _mm_loadu_si128((const __m128i*)(src + gap * 13));
	__m128i i7 = _mm_loadu_si128((const __m128i*)(src + gap * 15));
	idct16x16vertoddhalf(coef, odd, _mm_unpackhi_epi16(i0, i1), _mm_unpackhi_epi16(i2, i3), _mm_unpackhi_epi16(i4, i5), _mm_unpackhi_epi16(i6, i7));
	idct16x16vertoddhalf(coef, odd + 8, _mm_unpacklo_epi16(i0, i1), _mm_unpacklo_epi16(i2, i3), _mm_unpacklo_epi16(i4, i5), _mm_unpacklo_epi16(i6, i7));
}

static inline void idct16x16vert_half_pack(int pos, const __m128i odd[], const __m128i& e0, const __m128i& o0, __m128i& dst0, __m128i& dst1) {
	__m128i e0a = _mm_add_epi32(e0, *(const __m128i*)rounding);
	__m128i u0 = _mm_add_epi32(e0a, o0);
	__m128i u7 = _mm_sub_epi32(e0a, o0);
	dst0 = _mm_packs_epi32(_mm_srai_epi32(_mm_add_epi32(u0, *(odd + pos)), 7), _mm_srai_epi32(_mm_sub_epi32(u0, *(odd + pos)), 7));
	dst1 = _mm_packs_epi32(_mm_srai_epi32(_mm_add_epi32(u7, *(odd + 7 - pos)), 7), _mm_srai_epi32(_mm_sub_epi32(u7, *(odd + 7 - pos)), 7));
}

static inline void idct16x16vert_storehalf(int pos, __m128i tmp[], const __m128i& e0, const __m128i& o0) {
	__m128i u0, u7;
	idct16x16vert_half_pack(pos, tmp, e0, o0, u0, u7);
	_mm_store_si128(tmp + pos, u0);
	_mm_store_si128(tmp + 7 - pos, u7);
}

struct idct16x16vert_former {
	void operator()(__m128i half[], const __m128i& e0, const __m128i& e1, const __m128i& e2, const __m128i& e3, const __m128i& o0, const __m128i& o1, const __m128i& o2, const __m128i& o3) const {
		idct16x16vert_storehalf(0, half, e0, o0);
		idct16x16vert_storehalf(1, half, e1, o1);
		idct16x16vert_storehalf(2, half, e2, o2);
		idct16x16vert_storehalf(3, half, e3, o3);
	}
};

void idct16x16vert_latter4(int pos, int16_t* out, const __m128i half[], const __m128i& e0, const __m128i& o0) {
	__m128i u0, u7;
	idct16x16vert_half_pack(pos, half + 8, e0, o0, u0, u7);
	_mm_storeu_si128((__m128i*)(out + 16 * pos), _mm_unpacklo_epi64(u0, *(half + pos)));
	_mm_storeu_si128((__m128i*)(out + 16 * (15 - pos)), _mm_unpackhi_epi64(u0, *(half + pos)));
	_mm_storeu_si128((__m128i*)(out + 16 * (7 - pos)), _mm_unpacklo_epi64(u7, *(half + 7 - pos)));
	_mm_storeu_si128((__m128i*)(out + 16 * (8 + pos)), _mm_unpackhi_epi64(u7, *(half + 7 - pos)));
}

struct idct16x16vert_latter {
	void operator()(int16_t* out, const __m128i half[], const __m128i& e0, const __m128i& e1, const __m128i& e2, const __m128i& e3, const __m128i& o0, const __m128i& o1, const __m128i& o2, const __m128i& o3) const {
		idct16x16vert_latter4(0, out, half, e0, o0);
		idct16x16vert_latter4(1, out, half, e1, o1);
		idct16x16vert_latter4(2, out, half, e2, o2);
		idct16x16vert_latter4(3, out, half, e3, o3);
	}
};

template <typename T, typename F0, typename F1>
static inline void idct16x16verteven(const int16_t* src, int gap, T* out, __m128i odd[], F0 Former, F1 Latter) {
	__m128i s0 = _mm_loadu_si128((const __m128i*)src);
	__m128i s1 = _mm_loadu_si128((const __m128i*)(src + gap * 2));
	__m128i s2 = _mm_loadu_si128((const __m128i*)(src + gap * 4));
	__m128i s3 = _mm_loadu_si128((const __m128i*)(src + gap * 6));
	__m128i e0l, e1l, e2l, e3l;
	__m128i o0l, o1l, o2l, o3l;
	__m128i e0, e1, e2, e3;
	__m128i o0, o1, o2, o3;
	idct4x4vert4(_mm_unpackhi_epi16(s0, s2), _mm_unpackhi_epi16(s1, s3), e0l, e1l, e2l, e3l);
	__m128i t0 = _mm_loadu_si128((const __m128i*)(src + gap));
	__m128i t1 = _mm_loadu_si128((const __m128i*)(src + gap * 3));
	__m128i t2 = _mm_loadu_si128((const __m128i*)(src + gap * 5));
	__m128i t3 = _mm_loadu_si128((const __m128i*)(src + gap * 7));
	idct8x8vertodd4(_mm_unpackhi_epi16(t0, t2), _mm_unpackhi_epi16(t1, t3), o0l, o1l, o2l, o3l);
	idct16x16vert_storehalf(0, odd, e0l, o0l);
	idct16x16vert_storehalf(1, odd, e1l, o1l);
	idct16x16vert_storehalf(2, odd, e2l, o2l);
	idct16x16vert_storehalf(3, odd, e3l, o3l);
	idct4x4vert4(_mm_unpacklo_epi16(s0, s2), _mm_unpacklo_epi16(s1, s3), e0, e1, e2, e3);
	idct8x8vertodd4(_mm_unpacklo_epi16(t0, t2), _mm_unpacklo_epi16(t1, t3), o0, o1, o2, o3);
	Latter(out, odd, e0, e1, e2, e3, o0, o1, o2, o3);
}

static inline void idct16x16vert(const int16_t* src, int gap, int16_t* out) {
	ALIGNVC(16) __m128i odd[16];
	for (int x = 0; x < 2; ++x) {
		idct16x16vertodd(src + x * 8, gap, idct16coeffs, odd);
		idct16x16verteven(src + x * 8, gap * 2, out + x * 8, odd, idct16x16vert_former(), idct16x16vert_latter());
	}
}

template <typename F>
static void transform_ac16x16bit(uint8_t* out, const int16_t* in, int16_t* tmp, int stride, int bitdepth, F Write) {
	idct16x16vert(in, 16, tmp);
	for (int y = 0; y < 16; y++) {
		idct16x16horiz(tmp + y * 16, out + y * stride, bitdepth, Write);
	}
}

void transform_ac16x16(uint8_t* dst, int16_t* src, int stride) {
	transform_ac16x16bit(dst, src, src + 16 * 16, stride, 8, idct16x16write());
}

void transform_ac16x16chroma(uint8_t* dst, int16_t* src, int stride) {
	transform_ac16x16bit(dst, src, src + 16 * 16, stride, 8, idct16x16writechroma());
}

static inline void idct32x32vertodd4col(__m128i odd[], const __m128i& i0, const __m128i& i1, const __m128i& i2, const __m128i& i3, const __m128i& i4, const __m128i& i5, const __m128i& i6, const __m128i& i7) {
	static const int16_t coef[][8] = {
		{90, 90, 90, 90, 90, 90, 90, 90},
		{88, 85, 88, 85, 88, 85, 88, 85},
		{82, 78, 82, 78, 82, 78, 82, 78},
		{73, 67, 73, 67, 73, 67, 73, 67},
		{61, 54, 61, 54, 61, 54, 61, 54},
		{46, 38, 46, 38, 46, 38, 46, 38},
		{31, 22, 31, 22, 31, 22, 31, 22},
		{13, 4, 13, 4, 13, 4, 13, 4},
		{90, 82, 90, 82, 90, 82, 90, 82},
		{67, 46, 67, 46, 67, 46, 67, 46},
		{22, -4, 22, -4, 22, -4, 22, -4},
		{-31, -54, -31, -54, -31, -54, -31, -54},
		{-73, -85, -73, -85, -73, -85, -73, -85},
		{-90, -88, -90, -88, -90, -88, -90, -88},
		{-78, -61, -78, -61, -78, -61, -78, -61},
		{-38, -13, -38, -13, -38, -13, -38, -13},
		{88, 67, 88, 67, 88, 67, 88, 67},
		{31, -13, 31, -13, 31, -13, 31, -13},
		{-54, -82, -54, -82, -54, -82, -54, -82},
		{-90, -78, -90, -78, -90, -78, -90, -78},
		{-46, -4, -46, -4, -46, -4, -46, -4},
		{38, 73, 38, 73, 38, 73, 38, 73},
		{90, 85, 90, 85, 90, 85, 90, 85},
		{61, 22, 61, 22, 61, 22, 61, 22},
		{85, 46, 85, 46, 85, 46, 85, 46},
		{-13, -67, -13, -67, -13, -67, -13, -67},
		{-90, -73, -90, -73, -90, -73, -90, -73},
		{-22, 38, -22, 38, -22, 38, -22, 38},
		{82, 88, 82, 88, 82, 88, 82, 88},
		{54, -4, 54, -4, 54, -4, 54, -4},
		{-61, -90, -61, -90, -61, -90, -61, -90},
		{-78, -31, -78, -31, -78, -31, -78, -31},
		{82, 22, 82, 22, 82, 22, 82, 22},
		{-54, -90, -54, -90, -54, -90, -54, -90},
		{-61, 13, -61, 13, -61, 13, -61, 13},
		{78, 85, 78, 85, 78, 85, 78, 85},
		{31, -46, 31, -46, 31, -46, 31, -46},
		{-90, -67, -90, -67, -90, -67, -90, -67},
		{4, 73, 4, 73, 4, 73, 4, 73},
		{88, 38, 88, 38, 88, 38, 88, 38},
		{78, -4, 78, -4, 78, -4, 78, -4},
		{-82, -73, -82, -73, -82, -73, -82, -73},
		{13, 85, 13, 85, 13, 85, 13, 85},
		{67, -22, 67, -22, 67, -22, 67, -22},
		{-88, -61, -88, -61, -88, -61, -88, -61},
		{31, 90, 31, 90, 31, 90, 31, 90},
		{54, -38, 54, -38, 54, -38, 54, -38},
		{-90, -46, -90, -46, -90, -46, -90, -46},
		{73, -31, 73, -31, 73, -31, 73, -31},
		{-90, -22, -90, -22, -90, -22, -90, -22},
		{78, 67, 78, 67, 78, 67, 78, 67},
		{-38, -90, -38, -90, -38, -90, -38, -90},
		{-13, 82, -13, 82, -13, 82, -13, 82},
		{61, -46, 61, -46, 61, -46, 61, -46},
		{-88, -4, -88, -4, -88, -4, -88, -4},
		{85, 54, 85, 54, 85, 54, 85, 54},
		{67, -54, 67, -54, 67, -54, 67, -54},
		{-78, 38, -78, 38, -78, 38, -78, 38},
		{85, -22, 85, -22, 85, -22, 85, -22},
		{-90, 4, -90, 4, -90, 4, -90, 4},
		{90, 13, 90, 13, 90, 13, 90, 13},
		{-88, -31, -88, -31, -88, -31, -88, -31},
		{82, 46, 82, 46, 82, 46, 82, 46},
		{-73, -61, -73, -61, -73, -61, -73, -61},
		{61, -73, 61, -73, 61, -73, 61, -73},
		{-46, 82, -46, 82, -46, 82, -46, 82},
		{31, -88, 31, -88, 31, -88, 31, -88},
		{-13, 90, -13, 90, -13, 90, -13, 90},
		{-4, -90, -4, -90, -4, -90, -4, -90},
		{22, 85, 22, 85, 22, 85, 22, 85},
		{-38, -78, -38, -78, -38, -78, -38, -78},
		{54, 67, 54, 67, 54, 67, 54, 67},
		{54, -85, 54, -85, 54, -85, 54, -85},
		{-4, 88, -4, 88, -4, 88, -4, 88},
		{-46, -61, -46, -61, -46, -61, -46, -61},
		{82, 13, 82, 13, 82, 13, 82, 13},
		{-90, 38, -90, 38, -90, 38, -90, 38},
		{67, -78, 67, -78, 67, -78, 67, -78},
		{-22, 90, -22, 90, -22, 90, -22, 90},
		{-31, -73, -31, -73, -31, -73, -31, -73},
		{46, -90, 46, -90, 46, -90, 46, -90},
		{38, 54, 38, 54, 38, 54, 38, 54},
		{-90, 31, -90, 31, -90, 31, -90, 31},
		{61, -88, 61, -88, 61, -88, 61, -88},
		{22, 67, 22, 67, 22, 67, 22, 67},
		{-85, 13, -85, 13, -85, 13, -85, 13},
		{73, -82, 73, -82, 73, -82, 73, -82},
		{4, 78, 4, 78, 4, 78, 4, 78},
		{38, -88, 38, -88, 38, -88, 38, -88},
		{73, -4, 73, -4, 73, -4, 73, -4},
		{-67, 90, -67, 90, -67, 90, -67, 90},
		{-46, -31, -46, -31, -46, -31, -46, -31},
		{85, -78, 85, -78, 85, -78, 85, -78},
		{13, 61, 13, 61, 13, 61, 13, 61},
		{-90, 54, -90, 54, -90, 54, -90, 54},
		{22, -82, 22, -82, 22, -82, 22, -82},
		{31, -78, 31, -78, 31, -78, 31, -78},
		{90, -61, 90, -61, 90, -61, 90, -61},
		{4, 54, 4, 54, 4, 54, 4, 54},
		{-88, 82, -88, 82, -88, 82, -88, 82},
		{-38, -22, -38, -22, -38, -22, -38, -22},
		{73, -90, 73, -90, 73, -90, 73, -90},
		{67, -13, 67, -13, 67, -13, 67, -13},
		{-46, 85, -46, 85, -46, 85, -46, 85},
		{22, -61, 22, -61, 22, -61, 22, -61},
		{85, -90, 85, -90, 85, -90, 85, -90},
		{73, -38, 73, -38, 73, -38, 73, -38},
		{-4, 46, -4, 46, -4, 46, -4, 46},
		{-78, 90, -78, 90, -78, 90, -78, 90},
		{-82, 54, -82, 54, -82, 54, -82, 54},
		{-13, -31, -13, -31, -13, -31, -13, -31},
		{67, -88, 67, -88, 67, -88, 67, -88},
		{13, -38, 13, -38, 13, -38, 13, -38},
		{61, -78, 61, -78, 61, -78, 61, -78},
		{88, -90, 88, -90, 88, -90, 88, -90},
		{85, -73, 85, -73, 85, -73, 85, -73},
		{54, -31, 54, -31, 54, -31, 54, -31},
		{4, 22, 4, 22, 4, 22, 4, 22},
		{-46, 67, -46, 67, -46, 67, -46, 67},
		{-82, 90, -82, 90, -82, 90, -82, 90},
		{4, -13, 4, -13, 4, -13, 4, -13},
		{22, -31, 22, -31, 22, -31, 22, -31},
		{38, -46, 38, -46, 38, -46, 38, -46},
		{54, -61, 54, -61, 54, -61, 54, -61},
		{67, -73, 67, -73, 67, -73, 67, -73},
		{78, -82, 78, -82, 78, -82, 78, -82},
		{85, -88, 85, -88, 85, -88, 85, -88},
		{90, -90, 90, -90, 90, -90, 90, -90},
	};
	for (int y = 0; y < 16; ++y) {
		__m128i d0 = idct16x16vertodd4col(coef + y * 8, i0, i1, i2, i3);
		__m128i d1 = idct16x16vertodd4col(coef + y * 8 + 4, i4, i5, i6, i7);
		_mm_store_si128(odd + y, _mm_add_epi32(d0, d1));
	}
}

template <typename F>
static inline void unpack4(const int16_t* src, const __m128i& s1, const __m128i& s5, const __m128i& s9, const __m128i& s13, __m128i& d0, __m128i& d1, __m128i& d2, __m128i& d3, F unpack) {
	d0 = unpack(s1, *(const __m128i*)(src + 32 * 3));
	d1 = unpack(s5, *(const __m128i*)(src + 32 * 7));
	d2 = unpack(s9, *(const __m128i*)(src + 32 * 11));
	d3 = unpack(s13, *(const __m128i*)(src + 32 * 15));
}

struct unpack16low {
	__m128i operator()(const __m128i& i0, const __m128i& i1) const {
		return _mm_unpacklo_epi16(i0, i1);
	}
};

struct unpack16high {
	__m128i operator()(const __m128i& i0, const __m128i& i1) const {
		return _mm_unpackhi_epi16(i0, i1);
	}
};

static inline void idct32x32vertodd(const int16_t* src, __m128i odd[]) {
	__m128i i0, i1, i2, i3, i4, i5, i6, i7;
	__m128i s1 = _mm_load_si128((const __m128i*)(src + 32 * 1));
	__m128i s5 = _mm_load_si128((const __m128i*)(src + 32 * 5));
	__m128i s9 = _mm_load_si128((const __m128i*)(src + 32 * 9));
	__m128i s13 = _mm_load_si128((const __m128i*)(src + 32 * 13));
	unpack4(src, s1, s5, s9, s13, i0, i1, i2, i3, unpack16high());
	__m128i s17 = _mm_load_si128((const __m128i*)(src + 32 * 17));
	__m128i s21 = _mm_load_si128((const __m128i*)(src + 32 * 21));
	__m128i s25 = _mm_load_si128((const __m128i*)(src + 32 * 25));
	__m128i s29 = _mm_load_si128((const __m128i*)(src + 32 * 29));
	unpack4(src + 32 * 16, s17, s21, s25, s29, i4, i5, i6, i7, unpack16high());
	idct32x32vertodd4col(odd, i0, i1, i2, i3, i4, i5, i6, i7);
	unpack4(src, s1, s5, s9, s13, i0, i1, i2, i3, unpack16low());
	unpack4(src + 32 * 16, s17, s21, s25, s29, i4, i5, i6, i7, unpack16low());
	idct32x32vertodd4col(odd + 16, i0, i1, i2, i3, i4, i5, i6, i7);
}

struct store_former {
	void operator()(int16_t* out0, int16_t* out31, __m128i p0, __m128i odd16[]) const {
		_mm_store_si128(odd16, p0);
	}
};

struct store_latter {
	void operator()(int16_t* out0, int16_t* out31, __m128i p0, __m128i odd16[]) const {
		_mm_store_si128((__m128i*)out0, _mm_unpacklo_epi64(p0, *odd16));
		_mm_store_si128((__m128i*)out31, _mm_unpackhi_epi64(p0, *odd16));
	}
};

template <typename F>
static inline void idct32x32vert_half(int16_t* out, int pos, const __m128i odd8[], const __m128i odd16[], __m128i odd16w[], const __m128i& e0, const __m128i& o0, F Write) {
	__m128i ea = _mm_add_epi32(e0, *(const __m128i*)rounding);
	__m128i u0 = _mm_add_epi32(ea, o0);
	__m128i u7 = _mm_sub_epi32(ea, o0);

	__m128i i0 = _mm_add_epi32(u0, *(odd8 + pos));
	__m128i i15 = _mm_sub_epi32(u0, *(odd8 + pos));

	__m128i t0 = _mm_add_epi32(i0, *(odd16 + pos));
	__m128i t31 = _mm_sub_epi32(i0, *(odd16 + pos));
	__m128i p0 = _mm_packs_epi32(_mm_srai_epi32(t0, 7), _mm_srai_epi32(t31, 7));
	Write(out + pos * 32, out + (31 - pos) * 32, p0, odd16w + pos);

	__m128i t15 = _mm_add_epi32(i15, *(odd16 + 15 - pos));
	__m128i t16 = _mm_sub_epi32(i15, *(odd16 + 15 - pos));
	__m128i p1 = _mm_packs_epi32(_mm_srai_epi32(t15, 7), _mm_srai_epi32(t16, 7));
	Write(out + (15 - pos) * 32, out + (16 + pos) * 32, p1, odd16w + 15 - pos);

	__m128i i7 = _mm_add_epi32(u7, *(odd8 + 7 - pos));
	__m128i i8 = _mm_sub_epi32(u7, *(odd8 + 7 - pos));

	__m128i t7 = _mm_add_epi32(i7, *(odd16 + 7 - pos));
	__m128i t24 = _mm_sub_epi32(i7, *(odd16 + 7 - pos));
	__m128i p2 = _mm_packs_epi32(_mm_srai_epi32(t7, 7), _mm_srai_epi32(t24, 7));
	Write(out + (7 - pos) * 32, out + (24 + pos) * 32, p2, odd16w + 7 - pos);

	__m128i t8 = _mm_add_epi32(i8, *(odd16 + 8 + pos));
	__m128i t23 = _mm_sub_epi32(i8, *(odd16 + 8 + pos));
	__m128i p3 = _mm_packs_epi32(_mm_srai_epi32(t8, 7), _mm_srai_epi32(t23, 7));
	Write(out + (8 + pos) * 32, out + (23 - pos) * 32, p3, odd16w + 8 + pos);
}

static inline void idct32x32verteven(const int16_t* src, int16_t* out, __m128i odd16[], const __m128i odd8[]) {
	int gap = 32 * 4;
	__m128i s0 = _mm_loadu_si128((const __m128i*)src);
	__m128i s1 = _mm_loadu_si128((const __m128i*)(src + gap * 2));
	__m128i s2 = _mm_loadu_si128((const __m128i*)(src + gap * 4));
	__m128i s3 = _mm_loadu_si128((const __m128i*)(src + gap * 6));
	__m128i e0l, e1l, e2l, e3l;
	__m128i o0l, o1l, o2l, o3l;
	idct4x4vert4(_mm_unpackhi_epi16(s0, s2), _mm_unpackhi_epi16(s1, s3), e0l, e1l, e2l, e3l);
	__m128i t0 = _mm_loadu_si128((const __m128i*)(src + gap));
	__m128i t1 = _mm_loadu_si128((const __m128i*)(src + gap * 3));
	__m128i t2 = _mm_loadu_si128((const __m128i*)(src + gap * 5));
	__m128i t3 = _mm_loadu_si128((const __m128i*)(src + gap * 7));
	idct8x8vertodd4(_mm_unpackhi_epi16(t0, t2), _mm_unpackhi_epi16(t1, t3), o0l, o1l, o2l, o3l);
	idct32x32vert_half(out, 0, odd8, odd16, odd16, e0l, o0l, store_former());
	idct32x32vert_half(out, 1, odd8, odd16, odd16, e1l, o1l, store_former());
	idct32x32vert_half(out, 2, odd8, odd16, odd16, e2l, o2l, store_former());
	idct32x32vert_half(out, 3, odd8, odd16, odd16, e3l, o3l, store_former());
	__m128i e0, e1, e2, e3;
	__m128i o0, o1, o2, o3;
	idct4x4vert4(_mm_unpacklo_epi16(s0, s2), _mm_unpacklo_epi16(s1, s3), e0, e1, e2, e3);
	idct8x8vertodd4(_mm_unpacklo_epi16(t0, t2), _mm_unpacklo_epi16(t1, t3), o0, o1, o2, o3);
	idct32x32vert_half(out, 0, odd8 + 8, odd16 + 16, odd16, e0, o0, store_latter());
	idct32x32vert_half(out, 1, odd8 + 8, odd16 + 16, odd16, e1, o1, store_latter());
	idct32x32vert_half(out, 2, odd8 + 8, odd16 + 16, odd16, e2, o2, store_latter());
	idct32x32vert_half(out, 3, odd8 + 8, odd16 + 16, odd16, e3, o3, store_latter());
}

static inline void idct32x32vert(const int16_t* src, int16_t* out) {
	ALIGNVC(16) __m128i odd[32];
	ALIGNVC(16) __m128i odd2[16];
	for (int x = 0; x < 4; ++x) {
		idct32x32vertodd(src + x * 8, odd);
		idct16x16vertodd(src + x * 8, 64, idct16coeffs, odd2);
		idct32x32verteven(src + x * 8, out + x * 8, odd, odd2);
	}
}

static inline void maccum4(__m128i& a0, __m128i& a1, __m128i& a2, __m128i& a3, const __m128i& d0, const int16_t coef[][8]) {
	a0 = _mm_add_epi32(a0, _mm_madd_epi16(d0, *(const __m128i*)coef[0]));
	a1 = _mm_add_epi32(a1, _mm_madd_epi16(d0, *(const __m128i*)coef[1]));
	a2 = _mm_add_epi32(a2, _mm_madd_epi16(d0, *(const __m128i*)coef[2]));
	a3 = _mm_add_epi32(a3, _mm_madd_epi16(d0, *(const __m128i*)coef[3]));
}

static inline void idct32x32horizodd(const __m128i& dl0, const __m128i& dh0, __m128i& b0, __m128i& b1, __m128i& b2, __m128i& b3) {
	ALIGNVC(16) static const int16_t coefs[][8] = {
		// assume dl0, dh0 were ordered as:
		// dl0 = [1, 5, 9, 13, 17, 21, 25, 29]
		// dh0 = [3, 7, 11, 15, 19, 23, 27, 31]
		{90, 88, 22, -31, -46, 38, -61, -78},
		{82, -54, 13, 67, -13, 61, 82, -73},
		{61, -46, -46, 82, 22, -85, -90, 22},
		{31, 90, 73, -4, 54, 4, 85, 90},
		{90, 85, -4, -54, -4, 73, -90, -31},
		{22, -90, 85, -22, 82, -46, 46, -61},
		{-73, 82, -61, 13, 67, 13, 54, -82},
		{-78, -61, -38, 46, -31, 22, -88, -90},
		{82, 73, -73, -90, 90, 61, 85, -13},
		{-61, 78, -88, 31, -88, 85, 67, -78},
		{31, -13, -90, 67, 73, 4, 38, 73},
		{4, -88, -78, -82, -46, -82, 4, 22},
		{78, 67, -85, -88, 85, 22, 46, -67},
		{13, 85, -61, 90, -4, 54, -54, 38},
		{-88, 90, 38, -78, -82, 78, -88, -4},
		{54, 82, 90, 54, 67, 90, -13, -31},
		{61, 46, -78, -38, 88, 31, -90, -22},
		{31, -90, 54, -90, 73, -90, 85, -90},
		{-4, 22, -22, -31, 46, 38, -67, -46},
		{-38, 73, -13, 67, 13, 61, 38, 54},
		{54, 38, -61, -13, 67, -13, -73, 38},
		{-46, -67, -38, -46, -31, -22, -22, 4},
		{-90, 85, 90, -73, -90, 54, 90, -31},
		{-22, -90, -31, -88, -38, -78, -46, -61},
		{31, 13, 90, 67, -54, -90, 82, 54},
		{4, 88, 78, -82, 78, -38, 90, -88},
		{-38, 54, 54, -4, -90, 61, 85, 13},
		{67, -46, 22, 85, 88, 85, 67, 78},
		{22, 4, 82, 46, -82, -78, 88, -4},
		{73, 38, -4, -73, 67, -90, 13, -31},
		{-78, 67, -85, 88, 31, -88, -78, 61},
		{-13, 85, -61, -90, -90, -73, -73, -82}
	};
	__m128i a0 = _mm_madd_epi16(dl0, *(const __m128i*)coefs[0]);
	__m128i a1 = _mm_madd_epi16(dl0, *(const __m128i*)coefs[1]);
	__m128i a2 = _mm_madd_epi16(dl0, *(const __m128i*)coefs[2]);
	__m128i a3 = _mm_madd_epi16(dl0, *(const __m128i*)coefs[3]);
	__m128i dl1 = _mm_shuffle_epi32(dl0, 0x39);
	maccum4(a0, a1, a2, a3, dh0, coefs + 4);
	__m128i dh1 = _mm_shuffle_epi32(dh0, 0x39);
	maccum4(a0, a1, a2, a3, dl1, coefs + 8);
	__m128i dl2 = _mm_shuffle_epi32(dl0, 0x4e);
	maccum4(a0, a1, a2, a3, dh1, coefs + 12);
	__m128i dh2 = _mm_shuffle_epi32(dh0, 0x4e);
	maccum4(a0, a1, a2, a3, dl2, coefs + 16);
	__m128i dl3 = _mm_shuffle_epi32(dl0, 0x93);
	maccum4(a0, a1, a2, a3, dh2, coefs + 20);
	__m128i dh3 = _mm_shuffle_epi32(dh0, 0x93);
	maccum4(a0, a1, a2, a3, dl3, coefs + 24);
	b0 = _mm_add_epi32(a0, _mm_madd_epi16(dh3, *(const __m128i*)coefs[28]));
	b1 = _mm_add_epi32(a1, _mm_madd_epi16(dh3, *(const __m128i*)coefs[29]));
	b2 = _mm_add_epi32(a2, _mm_madd_epi16(dh3, *(const __m128i*)coefs[30]));
	b3 = _mm_add_epi32(a3, _mm_madd_epi16(dh3, *(const __m128i*)coefs[31]));
}

struct idct32x32write {
	void operator()(uint8_t* dst, int shift, const __m128i& d0, const __m128i& d1, const __m128i& d2, const __m128i& d3, const __m128i& d4, const __m128i& d5, const __m128i& d6, const __m128i& d7) const {
		idctwrite16(dst, shift, _mm_add_epi32(d0, d4), _mm_add_epi32(d1, d5), _mm_add_epi32(d2, d6), _mm_add_epi32(d3, d7));
		idctwrite16(dst + 16, shift, _mm_shuffle_epi32(_mm_sub_epi32(d3, d7), 0x1b), _mm_shuffle_epi32(_mm_sub_epi32(d2, d6), 0x1b), _mm_shuffle_epi32(_mm_sub_epi32(d1, d5), 0x1b), _mm_shuffle_epi32(_mm_sub_epi32(d0, d4), 0x1b));
	}
};

struct idct32x32writechroma {
	void operator()(uint8_t* dst, int shift, const __m128i& d0, const __m128i& d1, const __m128i& d2, const __m128i& d3, const __m128i& d4, const __m128i& d5, const __m128i& d6, const __m128i& d7) const {
		idctwrite16chroma(dst, shift, _mm_add_epi32(d0, d4), _mm_add_epi32(d1, d5), _mm_add_epi32(d2, d6), _mm_add_epi32(d3, d7));
		idctwrite16chroma(dst + 32, shift, _mm_shuffle_epi32(_mm_sub_epi32(d3, d7), 0x1b), _mm_shuffle_epi32(_mm_sub_epi32(d2, d6), 0x1b), _mm_shuffle_epi32(_mm_sub_epi32(d1, d5), 0x1b), _mm_shuffle_epi32(_mm_sub_epi32(d0, d4), 0x1b));
	}
};

template <typename F>
static inline void idct32x32horiz(const int16_t* in, uint8_t* dst, int bitdepth, F Write) {
	__m128i s0 = _mm_load_si128((const __m128i*)(in + 0));
	__m128i s1 = _mm_load_si128((const __m128i*)(in + 8));
	__m128i s2 = _mm_load_si128((const __m128i*)(in + 16));
	__m128i s3 = _mm_load_si128((const __m128i*)(in + 24));

	__m128i t0 = _mm_unpacklo_epi16(s0, s1); // 0, 8, 1, 9, 2, 10, 3, 11
	__m128i t1 = _mm_unpackhi_epi16(s0, s1); // 4, 12, 5, 13, 6, 14, 7, 15
	__m128i t2 = _mm_unpacklo_epi16(s2, s3); // 16, 24, 17, 25, 18, 26, 19, 27
	__m128i t3 = _mm_unpackhi_epi16(s2, s3); // 20, 28, 21, 29, 22, 30, 23, 31
	__m128i u0 = _mm_unpacklo_epi32(t0, t2); // 0, 8, 16, 24, 1, 9, 17, 25
	__m128i u1 = _mm_unpacklo_epi32(t1, t3); // 4, 12, 20, 28, 5, 13, 21, 29
	__m128i u2 = _mm_unpackhi_epi16(t0, t1); // 2, 6, 10, 14, 3, 7, 11, 15
	__m128i u3 = _mm_unpackhi_epi16(t2, t3); // 18, 22, 26, 30, 19, 23, 27, 31
	__m128i v0 = _mm_unpacklo_epi64(u0, u1); // 0, 8, 16, 24, 4, 12, 20, 28
	__m128i v1 = _mm_unpacklo_epi64(u2, u3); // 2, 6, 10, 14, 18, 22, 26, 30
	int shift = 20 - bitdepth;
	__m128i a0, a1, a2, a3;
	idct16x16horizcore(shift, v0, v1, a0, a1, a2, a3);

	// odd part
	__m128i w0 = _mm_unpackhi_epi16(u0, u1); // 1, 5, 9, 13, 17, 21, 25, 29
	__m128i w1 = _mm_unpackhi_epi64(u2, u3); // 3, 7, 11, 15, 19, 23, 27, 31
	__m128i b0, b1, b2, b3;
	idct32x32horizodd(w0, w1, b0, b1, b2, b3);
#if 1
	Write(dst, shift, a0, a1, a2, a3, b0, b1, b2, b3);
#else
	idct32x32write16(dst, shift, _mm_add_epi32(a0, b0), _mm_add_epi32(a1, b1), _mm_add_epi32(a2, b2), _mm_add_epi32(a3, b3));
	idct32x32write16(dst + 16, shift, _mm_shuffle_epi32(_mm_sub_epi32(a3, b3), 0x1b), _mm_shuffle_epi32(_mm_sub_epi32(a2, b2), 0x1b), _mm_shuffle_epi32(_mm_sub_epi32(a1, b1), 0x1b), _mm_shuffle_epi32(_mm_sub_epi32(a0, b0), 0x1b));
#endif
}

template <typename F>
static void transform_ac32x32bit(uint8_t* out, const int16_t* in, int16_t* tmp, int stride, int bitdepth, F Write) {
	idct32x32vert(in, tmp);
	for (int y = 0; y < 32; y++) {
		idct32x32horiz(tmp + y * 32, out + y * stride, bitdepth, Write);
	}
}

void transform_ac32x32(uint8_t* dst, int16_t* src, int stride) {
	transform_ac32x32bit(dst, src, src + 32 * 32, stride, 8, idct32x32write());
}

void transform_ac32x32chroma(uint8_t* dst, int16_t* src, int stride) {
	transform_ac32x32bit(dst, src, src + 32 * 32, stride, 8, idct32x32writechroma());
}

#endif
