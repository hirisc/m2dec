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

ALIGNVC(16) static const int16_t idctcoef4x4[][8] = {
	{64, 64, 64, 64, 64, 64, 64, 64},
	{83, 36, 83, 36, 83, 36, 83, 36},
	{64, -64, 64, -64, 64, -64, 64, -64},
	{36, -83, 36, -83, 36, -83, 36, -83},
	{64, -64, 64, -64, 64, -64, 64, -64},
	{-36, 83, -36, 83, -36, 83, -36, 83},
	{64, 64, 64, 64, 64, 64, 64, 64},
	{-83, -36, -83, -36, -83, -36, -83, -36},

	{64, 83, 64, 36, -64, 83, 64, -36},
	{64, 36, -64, -83, 64, -36, 64, -83},
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

ALIGNVC(16) static const uint8_t ofs[16] = {
	128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128
};

static inline void transform4x4(uint8_t* dst, const int16_t* coeff, int stride, const int16_t matrices[][8]) {
	__m128i p0, p1, p2, p3;

	transform4x4core(coeff, matrices, 8, p0, p1, p2, p3);
	__m128i s0 = _mm_set_epi32(*(uint32_t*)(dst + stride * 3), *(uint32_t*)(dst + stride * 2), *(uint32_t*)(dst + stride), *(uint32_t*)dst);
	s0 = _mm_sub_epi8(s0, *(const __m128i*)ofs);
	__m128i d0 = _mm_packs_epi16(_mm_packs_epi32(p0, p1), _mm_packs_epi32(p2, p3));
	__m128i e0 = _mm_adds_epi8(s0, d0);
	e0 = _mm_add_epi8(e0, *(const __m128i*)ofs);
	*(uint16_t*)dst = _mm_extract_epi16(e0, 0);
	*(uint16_t*)(dst + 2) = _mm_extract_epi16(e0, 1);
	*(uint16_t*)(dst + stride) = _mm_extract_epi16(e0, 2);
	*(uint16_t*)(dst + stride + 2) = _mm_extract_epi16(e0, 3);
	*(uint16_t*)(dst + stride * 2) = _mm_extract_epi16(e0, 4);
	*(uint16_t*)(dst + stride * 2 + 2) = _mm_extract_epi16(e0, 5);
	*(uint16_t*)(dst + stride * 3) = _mm_extract_epi16(e0, 6);
	*(uint16_t*)(dst + stride * 3 + 2) = _mm_extract_epi16(e0, 7);
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

void transform_ac4x4(uint8_t* dst, int16_t* coeff, int stride) {
	transform4x4(dst, coeff, stride, idctcoef4x4);
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
	__m128i a0 = _mm_madd_epi16(t0, *(const __m128i*)idctcoef4x4[0]);
	__m128i b0 = _mm_madd_epi16(t1, *(const __m128i*)idctcoef4x4[1]);
	__m128i a1 = _mm_madd_epi16(t0, *(const __m128i*)idctcoef4x4[2]);
	__m128i b1 = _mm_madd_epi16(t1, *(const __m128i*)idctcoef4x4[3]);
	out0 = _mm_add_epi32(a0, b0);
	out1 = _mm_add_epi32(a1, b1);
	out2 = _mm_sub_epi32(a1, b1);
	out3 = _mm_sub_epi32(a0, b0);
}

static inline void idct8x8vertodd4(const __m128i& t0, const __m128i& t1, __m128i& out0, __m128i& out1, __m128i& out2, __m128i& out3) {
	out0 = _mm_add_epi32(_mm_madd_epi16(t0, *(const __m128i*)idct8coeffs[0]), _mm_madd_epi16(t1, *(const __m128i*)idct8coeffs[1]));
	out1 = _mm_add_epi32(_mm_madd_epi16(t0, *(const __m128i*)idct8coeffs[2]), _mm_madd_epi16(t1, *(const __m128i*)idct8coeffs[3]));
	out2 = _mm_add_epi32(_mm_madd_epi16(t0, *(const __m128i*)idct8coeffs[4]), _mm_madd_epi16(t1, *(const __m128i*)idct8coeffs[5]));
	out3 = _mm_add_epi32(_mm_madd_epi16(t0, *(const __m128i*)idct8coeffs[6]), _mm_madd_epi16(t1, *(const __m128i*)idct8coeffs[7]));
}

static inline __m128i idct8x8preswap(const __m128i& d0l, const __m128i& d0h) {
	__m128i d0 = _mm_unpacklo_epi32(d0l, d0h); // 0, 4, 1, 5
	__m128i d1 = _mm_unpackhi_epi32(d0l, d0h); // 2, 6, 3, 7
	__m128i e0 = _mm_srai_epi32(_mm_add_epi32(_mm_unpacklo_epi32(d0, d1), *(const __m128i*)rounding), 7); // 0, 2, 4, 6
	__m128i e1 = _mm_srai_epi32(_mm_add_epi32(_mm_unpackhi_epi32(d0, d1), *(const __m128i*)rounding), 7); // 1, 3, 5, 7
	return _mm_packs_epi32(e0, e1); // 0, 2, 4, 6, 1, 3, 5, 7
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

static inline void idct8x8horiz1(uint8_t* dst, int shift, __m128i rnd, const __m128i& d0l, const __m128i& d0h) {
	ALIGNVC(16) static const int zeros[] = {
		0, 0, 0, 0
	};
	__m128i a0, a1;
	idct8x8horiz1core(idct8x8preswap(d0l, d0h), a0, a1);
	__m128i s0 = _mm_sub_epi8(_mm_loadu_si128((const __m128i*)dst), *(const __m128i*)ofs);
	__m128i t0 = _mm_packs_epi16(_mm_packs_epi32(_mm_srai_epi32(_mm_add_epi32(a0, rnd), shift), _mm_srai_epi32(_mm_add_epi32(a1, rnd), shift)), *(const __m128i*)zeros);
	__m128i u0 = _mm_add_epi8(_mm_adds_epi8(s0, t0), *(const __m128i*)ofs);
	_mm_storeu_si128((__m128i*)dst, u0);
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

static inline void idct16x16vertodd4(int idx, const __m128i& a0, const __m128i& a1, const __m128i& a2, const __m128i& a3, const __m128i& e0, const __m128i& e1, __m128i& dst0, __m128i& dst1, __m128i& dst2, __m128i& dst3) {
	__m128i b0 = _mm_add_epi32(_mm_madd_epi16(a0, *(const __m128i*)idct16coeffs[idx * 8]), _mm_madd_epi16(a1, *(const __m128i*)idct16coeffs[idx * 8 + 1]));
	__m128i b1 = _mm_add_epi32(_mm_madd_epi16(a2, *(const __m128i*)idct16coeffs[idx * 8 + 2]), _mm_madd_epi16(a3, *(const __m128i*)idct16coeffs[idx * 8 + 3]));
	__m128i d0 = _mm_add_epi32(b0, b1);
	dst0 = _mm_add_epi32(e0, d0);
	dst2 = _mm_sub_epi32(e0, d0);
	__m128i b2 = _mm_add_epi32(_mm_madd_epi16(a0, *(const __m128i*)idct16coeffs[idx * 8 + 4]), _mm_madd_epi16(a1, *(const __m128i*)idct16coeffs[idx * 8 + 5]));
	__m128i b3 = _mm_add_epi32(_mm_madd_epi16(a2, *(const __m128i*)idct16coeffs[idx * 8 + 6]), _mm_madd_epi16(a3, *(const __m128i*)idct16coeffs[idx * 8 + 7]));
	__m128i d1 = _mm_add_epi32(b2, b3);
	dst1 = _mm_add_epi32(e1, d1);
	dst3 = _mm_sub_epi32(e1, d1);
}

static inline void load16xNpack(const uint8_t* src, int gap, int pos, __m128i& dst0, __m128i& dst1, __m128i& dst2, __m128i& dst3) {
	__m128i s0 = _mm_packs_epi32(_mm_loadu_si128((const __m128i*)(src + gap * (pos + 1))), _mm_loadu_si128((const __m128i*)(src + gap * (pos + 5))));
	__m128i s1 = _mm_packs_epi32(_mm_loadu_si128((const __m128i*)(src + gap * (pos + 3))), _mm_loadu_si128((const __m128i*)(src + gap * (pos + 7))));
	__m128i s2 = _mm_packs_epi32(_mm_loadu_si128((const __m128i*)(src + gap * (pos + 9))), _mm_loadu_si128((const __m128i*)(src + gap * (pos + 13))));
	__m128i s3 = _mm_packs_epi32(_mm_loadu_si128((const __m128i*)(src + gap * (pos + 11))), _mm_loadu_si128((const __m128i*)(src + gap * (pos + 15))));
	dst0 = _mm_unpacklo_epi16(s0, s1);
	dst1 = _mm_unpackhi_epi16(s0, s1);
	dst2 = _mm_unpacklo_epi16(s2, s3);
	dst3 = _mm_unpackhi_epi16(s2, s3);
}

static inline void idct16x16vertodd4pack(int idx, int16_t* dst, int pos0, int pos1, const __m128i& a0, const __m128i& a1, const __m128i& a2, const __m128i& a3, const __m128i& e0, const __m128i& e1) {
	__m128i l0, l1, h0, h1;
	idct16x16vertodd4(idx, a0, a1, a2, a3, e0, e1, l0, l1, h0, h1);
	__m128i r0 = _mm_packs_epi32(_mm_srai_epi32(_mm_add_epi32(l0, *(const __m128i*)rounding), 7), _mm_srai_epi32(_mm_add_epi32(l1, *(const __m128i*)rounding), 7));
	__m128i r1 = _mm_packs_epi32(_mm_srai_epi32(_mm_add_epi32(h1, *(const __m128i*)rounding), 7), _mm_srai_epi32(_mm_add_epi32(h0, *(const __m128i*)rounding), 7));
	_mm_storeu_si128((__m128i*)(dst + 16 * pos0), r0);
	_mm_storeu_si128((__m128i*)(dst + 16 * pos1), r1);
}

static inline void idct16x16vertodd(const uint8_t* src, int gap, int16_t* out, const __m128i& even0, const __m128i& even1, const __m128i& even2, const __m128i& even3, const __m128i& even4, const __m128i& even5, const __m128i& even6, const __m128i& even7) {
	__m128i s0, s1, s2, s3;
	load16xNpack(src, gap, 0, s0, s1, s2, s3);
	idct16x16vertodd4pack(0, out, 0, 7, s0, s1, s2, s3, even0, even1);
	idct16x16vertodd4pack(1, out, 1, 6, s0, s1, s2, s3, even2, even3);
	idct16x16vertodd4pack(2, out, 2, 5, s0, s1, s2, s3, even4, even5);
	idct16x16vertodd4pack(3, out, 3, 4, s0, s1, s2, s3, even6, even7);
}

#define madd8(cf, r0, r1, r2, r3) _mm_add_epi32(_mm_add_epi32(_mm_madd_epi16(r0, *(const __m128i*)*((cf) + 0)), _mm_madd_epi16(r1, *(const __m128i*)*((cf) + 1))), _mm_add_epi32(_mm_madd_epi16(r2, *(const __m128i*)*((cf) + 2)), _mm_madd_epi16(r3, *(const __m128i*)*((cf) + 3))))

#define IDCT16x16HORIZODD1LINE(reg0, dst0, dst1) do {			\
	__m128i reg1 = _mm_shuffle_epi32((reg0), 0xb1);			\
	__m128i reg2 = _mm_shuffle_epi32((reg0), 0x4e);			\
	__m128i reg3 = _mm_shuffle_epi32((reg0), 0x1b);			\
	dst0 = madd8(coeff16oddh + 0, (reg0), reg1, reg2, reg3);	\
	dst1 = madd8(coeff16oddh + 4, (reg0), reg1, reg2, reg3);	\
} while (0)

#define IDCT16x16_HorizontalOdd(odd0, odd1, dst0, dst1, dst2, dst3) do { \
	__m128i m0 = _mm_unpacklo_epi64(odd0, odd1);		\
	__m128i n0 = _mm_unpackhi_epi64(odd0, odd1);		\
	IDCT16x16HORIZODD1LINE(m0, dst0, dst1);				\
	IDCT16x16HORIZODD1LINE(n0, dst2, dst3);				\
} while (0)

#define IDCT16x16_Add1line(dst, even0, even1, odd0, odd1, rnd, shft)  do { \
	__m128i d0 = _mm_add_epi32(even0, odd0);				\
	__m128i d1 = _mm_add_epi32(_mm_shuffle_epi32(even1, 0x1b), odd1);	\
	__m128i d3 = _mm_shuffle_epi32(_mm_sub_epi32(even0, odd0), 0x1b);	\
	__m128i d2 = _mm_sub_epi32(even1, _mm_shuffle_epi32(odd1, 0x1b));	\
	_mm_storeu_si128((__m128i*)dst, _mm_srai_epi32(_mm_add_epi32(d0, rnd), shft)); \
	_mm_storeu_si128((__m128i*)(dst + 4), _mm_srai_epi32(_mm_add_epi32(d1, rnd), shft)); \
	_mm_storeu_si128((__m128i*)(dst + 8), _mm_srai_epi32(_mm_add_epi32(d2, rnd), shft)); \
	_mm_storeu_si128((__m128i*)(dst + 12), _mm_srai_epi32(_mm_add_epi32(d3, rnd), shft)); \
} while (0)

ALIGNVC(16) static const int16_t coeff16oddh[][8] = {
	{90, 57, 9, -70, 9, 57, 9, -57},
	{80, 25, 87, -80, -87, 43, -43, 25},
	{87, 43, -43, -25, 80, -25, -87, -80},
	{70, 9, 57, -90, -70, 90, 70, 90},
	{57, -9, 57, 9, -70, 9, -57, -90},
	{-25, 43, 43, -87, -80, 87, -25, -80},
	{-80, -87, 25, -80, 25, 43, 43, 87},
	{90, 70, -90, 70, 90, -57, 9, 70}
};

static inline void idct16x16horizodd(const __m128i& r0, __m128i& dst0, __m128i& dst1) {
	ALIGNVC(16) static const int16_t coef[][8] = {
/*          [90, 87, 80, 70, 57, 43, 25, 9],
          [87, 57, 9, -43, -80, -90, -70, -25],
          [80, 9, -70, -87, -25, 57, 90, 43],
          [70, -43, -87, 9, 90, 25, -80, -57],
          [57, -80, -25, 90, -9, -87, 43, 70],
          [43, -90, 57, 25, -87, 70, 9, -80],
          [25, -70, 90, -80, 43, 9, -57, 87],
          [9, -25, 43, -57, 70, -80, 87, -90]
*/
	};
	__m128i r1 = _mm_shuffle_epi32(r0, 0xb1); // 13, 15, 1, 3, 5, 7, 9, 11
	__m128i t2 = _mm_shuffle_epi32(r0, 0x4e); // 9, 11, 13, 15, 1, 3, 5, 7
	__m128i t3 = _mm_shuffle_epi32(r0, 0x4e); // 5, 7, 9, 11, 13, 15, 1, 3
	__m128i d0 = _mm_madd_epi16(r0, *(const __m128i*)coef[0]);
	__m128i e0 = _mm_add_epi32(d0, d1);
}

static inline void IDCT16x16_Horizontal(const int16_t* in, int* out, int bitdepth) {
	int shift;
	__m128i round;
	__m128i s0 = _mm_loadu_si128((const __m128i*)(in + 0));
	__m128i s1 = _mm_loadu_si128((const __m128i*)(in + 8));
	__m128i t0 = _mm_unpacklo_epi16(s0, s1); // 0, 8, 1, 9, 2, 10, 3, 11
	__m128i t1 = _mm_unpackhi_epi16(s0, s1); // 4, 12, 5, 13, 6, 14, 7, 15
	__m128i u0 = _mm_unpacklo_epi16(t0, t1); // 0, 4, 8, 12, 1, 5, 9, 13
	__m128i u1 = _mm_unpackhi_epi16(t0, t1); // 2, 6, 10, 14, 3, 7, 11, 15
	__m128i v0 = _mm_unpacklo_epi64(u0, u1); // 0, 4, 8, 12, 2, 6, 10, 14
	__m128i even0, even1;
	idct8x8horiz1core(v0, even0, even1);
	__m128i v1 = _mm_unpackhi_epi16(u0, u1); // 1, 3, 5, 7, 9, 11, 13, 15
	idct16x16horizodd(v1);
#if 0
	__m128i t1 = _mm_unpackhi_epi16(s0, s2); // 16, 24, 17, 25, 18, 26, 19, 27
	__m128i t2 = _mm_unpacklo_epi16(s1, s3); // 4, 12, 5, 13, 6, 14, 7, 15
	__m128i t3 = _mm_unpackhi_epi16(s1, s3); // 20, 28, 21, 29, 22, 30, 23, 31
	__m128i u0 = _mm_unpacklo_epi32(t0, t2); // 0, 8, 4, 12, 1, 9, 5, 13
	__m128i u1 = _mm_unpackhi_epi32(t0, t2); // 2, 10, 6, 14, 3, 11, 7, 15
	__m128i u2 = _mm_unpacklo_epi32(t1, t3); // 16, 24, 20, 28, 17, 25, 21, 29
	__m128i u3 = _mm_unpackhi_epi32(t1, t3); // 18, 26, 22, 30, 19, 27, 23, 31
	__m128i even0 = _mm_unpacklo_epi64(u0, u2);
	__m128i odd0 = _mm_unpackhi_epi64(u0, u2); // 1, 9, 5, 13, 17, 25, 21, 29
	__m128i even1 = _mm_unpacklo_epi64(u1, u3);
	__m128i odd1 = _mm_unpackhi_epi64(u1, u3); // 3, 11, 7, 15, 19, 27, 23, 31
	__m128i e0, e1, e2, e3;
	__m128i o0, o1, o2, o3;
	/*
	  0, 4, 2, 6, 16, 20, 18, 22
	  1, 5, 3, 7, 17, 21, 19, 23
	  8, 12, 10, 14, 24, 28, 26, 30
	  9, 13, 11, 15, 25, 29, 27, 31
	 */
	IDCT8x8_Horiz2core(even0, even1, e0, e1, e2, e3);
	IDCT16x16_HorizontalOdd(odd0, odd1, o0, o1, o2, o3);
	shift = 20 - bitdepth;
	round = _mm_sll_epi32(_mm_load_si128((const __m128i*)ones32), shift - 1);
	IDCT16x16_Add1line(out, e0, e1, o0, o1, round, shift);
	IDCT16x16_Add1line(out + 16, e2, e3, o2, o3, round, shift);
#endif
}

static inline void idct16x16even8l(const int16_t* src, int gap, __m128i& dst0, __m128i& dst1, __m128i& dst2, __m128i& dst3, __m128i& dst4, __m128i& dst5, __m128i& dst6, __m128i& dst7) {
	__m128i s0 = _mm_loadu_si128((const __m128i*)src);
	__m128i s1 = _mm_loadu_si128((const __m128i*)(src + gap * 2));
	__m128i s2 = _mm_loadu_si128((const __m128i*)(src + gap * 4));
	__m128i s3 = _mm_loadu_si128((const __m128i*)(src + gap * 6));
	__m128i e0l, e1l, e2l, e3l;
	__m128i o0l, o1l, o2l, o3l;
	__m128i e0h, e1h, e2h, e3h;
	__m128i o0h, o1h, o2h, o3h;
	idct4x4vert4(_mm_unpacklo_epi16(s0, s2), _mm_unpacklo_epi16(s1, s3), e0l, e1l, e2l, e3l);
	__m128i t0 = _mm_loadu_si128((const __m128i*)(src + gap));
	__m128i t1 = _mm_loadu_si128((const __m128i*)(src + gap * 3));
	__m128i t2 = _mm_loadu_si128((const __m128i*)(src + gap * 5));
	__m128i t3 = _mm_loadu_si128((const __m128i*)(src + gap * 7));
	idct8x8vertodd4(_mm_unpacklo_epi16(t0, t2), _mm_unpacklo_epi16(t1, t3), o0l, o1l, o2l, o3l);
	dst0 = _mm_add_epi32(e0l, o0l);
	dst7 = _mm_sub_epi32(e0l, o0l);
	dst1 = _mm_add_epi32(e1l, o1l);
	dst6 = _mm_sub_epi32(e1l, o1l);
	dst2 = _mm_add_epi32(e2l, o2l);
	dst5 = _mm_sub_epi32(e2l, o2l);
	dst3 = _mm_add_epi32(e3l, o3l);
	dst4 = _mm_sub_epi32(e3l, o3l);
}

static void DCT16x16_2D_InverseTransform(const int16_t* in, int16_t* tmp, uint8_t* out, int bitdepth) {
	int x, y;
	for (x = 0; x < 16; x += 8) {
		__m128i t0, t1, t2, t3, t4, t5, t6, t7;
		idct16x16even8l(in + x, 16, t0, t1, t2, t3, t4, t5, t6, t7);
		idct16x16vertodd(in + x, 16, tmp + x, t1, t2, t3, t4, t5, t6, t7);
}
	for (y = 0; y < 16; y += 2) {
		IDCT16x16_Horizontal(tmp + 8 * y, out + 16 * y, bitdepth);
	}
}
#endif
