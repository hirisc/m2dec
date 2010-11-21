/** Motion Compensation for MPEG-2.
 *  Copyright 2008 Takayuki Minegishi
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

#include "bitio.h"
#include "motioncomp.h"

#define HALFPEL(x, y) ((((y) & 1) * 2) | ((x) & 1))

#define X01 0x01010101
#define X3F 0x3F3F3F3F
#define X03 0x03030303
#define X02 0x02020202

#if defined(__RENESAS_VERSION__)
#pragma inline(AVERAGE2)
#endif

static inline uint32_t AVERAGE2(uint32_t s1, uint32_t s2)
{
	uint32_t x = s1 ^ s2;
	return (s1 & s2) + ((x & ~X01) >> 1) + (x & X01);
}
#define AVERAGE2_NORND(s1, s2) (((s1) & (s2)) + ((((s1) ^ (s2)) & ~X01) >> 1))



#ifndef FAST_DECODE


/** Functor that just stores passed value.
 * This is for mono-directional motion compensation.
 */
template <typename T>
class Store {
public:
	void operator ()(T *dst, int val) {
		*dst = val;
	}
};

/** Functor that stores average value of passed value and destination.
 * This is for bi-directional motion compensation.
 */
template <typename T>
class AveStore {
public:
	void operator ()(T *dst, int val) {
		*dst = (*dst + val + 1) >> 1;
	}
};

template <>
void AveStore<uint32_t>::operator ()(uint32_t *dst, int val) {
		*dst = AVERAGE2(*dst, val);
}


#if defined(__RENESAS_VERSION__) && (defined(_SH4ALDSP) || defined(_SH4A))
extern "C" {
void m2d_copy16xn(const uint8_t *src, uint8_t *dst, int stride, int height);
void m2d_bilinear16_horiz_22_rnd(const uint8_t *src, uint8_t *dst, int stride, int height);
void m2d_bilinear16_vert_22_rnd(const uint8_t *src, uint8_t *dst, int stride, int height);
void m2d_bilinear_chroma_horiz_rnd(const uint8_t *src, uint8_t *dst, int stride, int height);

void m2d_copy16xn_add(const uint8_t *src, uint8_t *dst, int stride, int height);
void m2d_bilinear16_horiz_22_rnd_add(const uint8_t *src, uint8_t *dst, int stride, int height);
void m2d_bilinear16_vert_22_rnd_add(const uint8_t *src, uint8_t *dst, int stride, int height);
void m2d_bilinear_chroma_horiz_rnd_add(const uint8_t *src, uint8_t *dst, int stride, int height);
}
#ifdef UNIT_TEST
#define m2d_copy16xn org_copy16xn
#define m2d_bilinear16_horiz_22_rnd org_bilinear16_horiz_22_rnd
#define m2d_bilinear16_vert_22_rnd org_bilinear16_vert_22_rnd
#define m2d_bilinear_chroma_horiz_rnd org_bilinear_chroma_horiz_rnd
#define m2d_copy16xn_add org_copy16xn_add
#define m2d_bilinear16_horiz_22_rnd_add org_bilinear16_horiz_22_rnd_add
#define m2d_bilinear16_vert_22_rnd_add org_bilinear16_vert_22_rnd_add
#define m2d_bilinear_chroma_horiz_rnd_add org_bilinear_chroma_horiz_rnd_add
#endif /* UNIT_TEST */
#endif /* __RENESAS_VERSION__ */

#if !defined(__RENESAS_VERSION__) || !(defined(_SH4ALDSP) || defined(_SH4A)) || defined(UNIT_TEST)

/**
 *\brief Dead-copy interpolation in case src is in four bytes alignment.
 *\param src Pointer to source region in reference frame.
 *\param dst Pointer to destination region in current frame.
 *\param stride Width of frame.
 *\param height Height of the block.
 */

template <typename F>
static inline void m2d_copy16xn_align4_base(const uint8_t *src, uint8_t *dst, int stride, int height, F _store)
{
	do {
		_store((uint32_t *)dst, *(const uint32_t *)src);
		_store((uint32_t *)(dst + 4), *(const uint32_t *)(src + 4));
		_store((uint32_t *)(dst + 8), *(const uint32_t *)(src + 8));
		_store((uint32_t *)(dst + 12), *(const uint32_t *)(src + 12));
		src += stride;
		dst += stride;
	} while (--height);
}

template <typename F>
static inline void m2d_copy16xn_align1_base(const uint8_t *src, uint8_t *dst, int stride, int height, F _store)
{
	do {
		const uint32_t *src32 = (const uint32_t *)src;
		uint32_t s0, s1, s2, s3;
		s0 = read4_unalign(src32++);
		s1 = read4_unalign(src32++);
		s2 = read4_unalign(src32++);
		s3 = read4_unalign(src32);
		_store((uint32_t *)dst, s0);
		_store((uint32_t *)(dst + 4), s1);
		_store((uint32_t *)(dst + 8), s2);
		_store((uint32_t *)(dst + 12), s3);
		src += stride;
		dst += stride;
	} while (--height);
}

void m2d_copy16xn(const uint8_t *src, uint8_t *dst, int stride, int height)
{
#if defined(__GNUC__) && defined(__i386__)
	asm volatile ("\n\t"
		"push	%%ebx\n\t"
		"push	%%esi\n\t"
		"movl	%2, %%ecx\n\t"
		"movl	%0, %%eax\n\t"
		"movl	%1, %%ebx\n\t"
		"movl	%%ecx, %%esi\n\t"
		"movl	%3, %%edx\n\t"
		"addl	%%ecx, %%esi\n\t"
	"L01_copy:\n\t"
		"movdqu	(%%eax), %%xmm0\n\t"
		"movdqu	(%%eax, %%ecx), %%xmm1\n\t"
		"movntdq	%%xmm0, (%%ebx)\n\t"
		"movntdq	%%xmm1, (%%ebx, %%ecx)\n\t"
		"addl	%%esi, %%eax\n\t"
		"addl	%%esi, %%ebx\n\t"
		"addl	$-2, %%edx\n\t"
		"jnz	L01_copy\n\t"
		"pop	%%esi\n\t"
		"pop	%%ebx\n\t"
		:
		: "m"(src), "m"(dst), "m"(stride), "m"(height));
#elif defined(_M_IX86)
	__asm {
		push	ebx
		push	esi
		mov	ecx, stride
		mov	eax, src
		mov	ebx, dst
		mov	esi, ecx
		mov	edx, height
		add	esi, ecx
L01_copy:
		movdqu	xmm0, [eax]
		movdqu	xmm1, [eax][ecx]
		movntdq	[ebx], xmm0
		movntdq	[ebx][ecx], xmm1
		add	eax, esi
		add	ebx, esi
		add	edx, -2
		jnz	L01_copy
		pop	esi
		pop	ebx
	};
#else
	if (((intptr_t)src & 3) == 0) {
		m2d_copy16xn_align4_base(src, dst, stride, height, Store<uint32_t>());
	} else {
		m2d_copy16xn_align1_base(src, dst, stride, height, Store<uint32_t>());
	}
#endif
}

static inline void m2d_copy16xn_add(const uint8_t *src, uint8_t *dst, int stride, int height)
{
	if (((intptr_t)src & 3) == 0) {
		m2d_copy16xn_align4_base(src, dst, stride, height, AveStore<uint32_t>());
	} else {
		m2d_copy16xn_align1_base(src, dst, stride, height, AveStore<uint32_t>());
	}
}

/**
 *\brief Bilinear 1:1 interpolation for whole block with rounding. Vertical only.
 */
template <typename F>
static inline void m2d_bilinear16_vert_22_rnd_base(const uint8_t *src, uint8_t *dst, int stride, int height, F _store)
{
	uint32_t s0, s1, s2, s3;

	s0 = read4_unalign((const uint32_t *)src);
	s1 = read4_unalign((const uint32_t *)(src + 4));
	s2 = read4_unalign((const uint32_t *)(src + 8));
	s3 = read4_unalign((const uint32_t *)(src + 12));
	src = src + stride;
	do {
		uint32_t s_bot;

		s_bot = read4_unalign((const uint32_t *)src);
		_store((uint32_t *)dst, AVERAGE2(s0, s_bot));
		s0 = s_bot;
		s_bot = read4_unalign((const uint32_t *)(src + 4));
		_store((uint32_t *)(dst + 4), AVERAGE2(s1, s_bot));
		s1 = s_bot;
		s_bot = read4_unalign((const uint32_t *)(src + 8));
		_store((uint32_t *)(dst + 8), AVERAGE2(s2, s_bot));
		s2 = s_bot;
		s_bot = read4_unalign((const uint32_t *)(src + 12));
		_store((uint32_t *)(dst + 12), AVERAGE2(s3, s_bot));
		s3 = s_bot;
		src += stride;
		dst += stride;
	} while (--height);
}

template <typename F>
static inline void m2d_bilinear16_horiz_22_rnd_base(const uint8_t *src, uint8_t *dst, int stride, int height, F _store)
{
	const uint8_t *sr2 = src + 1 + (height & 1);
	height &= ~1; // LSB of height is abused to share this function with luma/chroma.
	do {
		uint32_t t0, t1;

		t0 = read4_unalign((const uint32_t *)src);
		t1 = read4_unalign((const uint32_t *)sr2);
		_store((uint32_t *)dst, AVERAGE2(t0, t1));
		t0 = read4_unalign((const uint32_t *)(src + 4));
		t1 = read4_unalign((const uint32_t *)(sr2 + 4));
		_store((uint32_t *)(dst + 4), AVERAGE2(t0, t1));
		t0 = read4_unalign((const uint32_t *)(src + 8));
		t1 = read4_unalign((const uint32_t *)(sr2 + 8));
		_store((uint32_t *)(dst + 8), AVERAGE2(t0, t1));
		t0 = read4_unalign((const uint32_t *)(src + 12));
		t1 = read4_unalign((const uint32_t *)(sr2 + 12));
		_store((uint32_t *)(dst + 12), AVERAGE2(t0, t1));

		src += stride;
		sr2 += stride;
		dst += stride;
	} while (--height);
}

static void m2d_bilinear16_vert_22_rnd(const uint8_t *src, uint8_t *dst, int stride, int height)
{
	m2d_bilinear16_vert_22_rnd_base(src, dst, stride, height, Store<uint32_t>());
}

static void m2d_bilinear16_horiz_22_rnd(const uint8_t *src, uint8_t *dst, int stride, int height)
{
	m2d_bilinear16_horiz_22_rnd_base(src, dst, stride, height, Store<uint32_t>());
}

static void m2d_bilinear16_vert_22_rnd_add(const uint8_t *src, uint8_t *dst, int stride, int height)
{
	m2d_bilinear16_vert_22_rnd_base(src, dst, stride, height, AveStore<uint32_t>());
}

static void m2d_bilinear16_horiz_22_rnd_add(const uint8_t *src, uint8_t *dst, int stride, int height)
{
	m2d_bilinear16_horiz_22_rnd_base(src, dst, stride, height, AveStore<uint32_t>());
}

static void m2d_bilinear_chroma_horiz_rnd(const uint8_t *src, uint8_t *dst, int stride, int height)
{
	m2d_bilinear16_horiz_22_rnd_base(src, dst, stride, height + 1, Store<uint32_t>());
}

static void m2d_bilinear_chroma_horiz_rnd_add(const uint8_t *src, uint8_t *dst, int stride, int height)
{
	m2d_bilinear16_horiz_22_rnd_base(src, dst, stride, height + 1, AveStore<uint32_t>());
}

#endif /* __RENESAS_VERSION__ */

const int MB_LEN = 16;

/**
 *\brief WMV8 bilinear 1:1 for both horizontal/vertical direction.
 *
 * Rounding is to be done by addition of 2.
 *\param *src Pointer to MB of Reference frame. 
 *\param *dst Pointer to MB of Prediction frame.
 *\param stride Width of frame.
 *\param height Number of lines.
 */
template <typename F>
static inline void m2d_bilinear_22_22_rnd_base(const uint8_t *src, uint8_t *dst, int stride, int height, F _store)
{
	int y;
	unsigned int upper_line[MB_LEN];
	unsigned int *upper;
	unsigned int s0, s1;

	s0 = *src++;
	upper = upper_line;
	y = MB_LEN / 2;
	do {
		s1 = *src++ + 1;
		upper[0] = s0 + s1;
		s0 = *src++;
		upper[1] = s0 + s1;
		upper += 2;
	} while (--y);

	stride -= MB_LEN;
	src = src + stride - 1;
	do {
		int width;
		upper = upper_line;
		s0 = *src++;
		width = MB_LEN / 2;
		do {
			unsigned int knl_bot;
			unsigned int s1;
			s1 = *src++ + 1;
			knl_bot = s0 + s1;
			s0 = *src++;
			_store(dst, (knl_bot + upper[0]) >> 2);
			upper[0] = knl_bot;
			knl_bot = s0 + s1;
			_store(dst + 1, (knl_bot + upper[1]) >> 2);
			upper[1] = knl_bot;
			upper += 2;
			dst += 2;
		} while (--width);
		src = src + stride - 1;
		dst += stride;
	} while (--height);
}

static void m2d_bilinear_22_22_rnd(const uint8_t *src, uint8_t *dst, int stride, int height)
{
	m2d_bilinear_22_22_rnd_base(src, dst, stride, height, Store<uint8_t>());
}

static void m2d_bilinear_22_22_rnd_add(const uint8_t *src, uint8_t *dst, int stride, int height)
{
	m2d_bilinear_22_22_rnd_base(src, dst, stride, height, AveStore<uint8_t>());
}


template <typename F>
static inline void m2d_bilinear_chroma_22_22_rnd_base(const uint8_t *src, uint8_t *dst, int src_stride, int height, F _store)
{
	uint16_t upper_line[8 * 2];
	uint16_t *upper;
	int y;
	uint32_t cr0, cr1, cb0, cb1;

	upper = upper_line;
	y = 8 / 2; // Not height but width actuallly
	cb0 = *src++ + 1;
	cr0 = *src++ + 1;
	do {
		cb1 = *src++;
		cr1 = *src++;
		upper[0] = cb0 + cb1;
		cb0 = *src++ + 1;
		upper[1] = cr0 + cr1;
		cr0 = *src++ + 1;
		upper[2] = cb1 + cb0;
		upper[3] = cr1 + cr0;
		upper += 4;
	} while (--y);

	src_stride -= 8 * 2;
	src = src + src_stride - 2;
	do {
		int width;
		upper = upper_line;
		cb0 = *src++;
		cr0 = *src++;
		width = 8 / 2;
		do {
			unsigned int knl_bot;
			cb1 = *src++ + 1;
			cr1 = *src++ + 1;
			knl_bot = cb0 + cb1;
			_store(dst, (knl_bot + upper[0]) >> 2);
			upper[0] = knl_bot;
			knl_bot = cr0 + cr1;
			_store(dst + 1, (knl_bot + upper[1]) >> 2);
			upper[1] = knl_bot;

			cb0 = *src++;
			cr0 = *src++;
			knl_bot = cb0 + cb1;
			_store(dst + 2, (knl_bot + upper[2]) >> 2);
			upper[2] = knl_bot;
			knl_bot = cr0 + cr1;
			_store(dst + 3, (knl_bot + upper[3]) >> 2);
			upper[3] = knl_bot;

			upper += 4;
			dst += 4;
		} while (--width);
		src = src + src_stride - 2;
		dst += src_stride;
	} while (--height);
}

static void m2d_bilinear_chroma_22_22_rnd(const uint8_t *src, uint8_t *dst, int stride, int height)
{
	m2d_bilinear_chroma_22_22_rnd_base(src, dst, stride, height, Store<uint8_t>());
}


static void m2d_bilinear_chroma_22_22_rnd_add(const uint8_t *src, uint8_t *dst, int stride, int height)
{
	m2d_bilinear_chroma_22_22_rnd_base(src, dst, stride, height, AveStore<uint8_t>());
}

#if defined(__RENESAS_VERSION__) && (defined(_SH4ALDSP) || defined(_SH4A))
#undef m2d_copy16xn
#undef m2d_bilinear16_horiz_22_rnd
#undef m2d_bilinear16_vert_22_rnd
#undef m2d_bilinear_chroma_horiz_rnd
#undef m2d_copy16xn_add
#undef m2d_bilinear16_horiz_22_rnd_add
#undef m2d_bilinear16_vert_22_rnd_add
#undef m2d_bilinear_chroma_horiz_rnd_add
#endif

void (* const m2d_motion_comp_luma[4])(const uint8_t *src, uint8_t *dst, int stride, int height) = {
		m2d_copy16xn,
		m2d_bilinear16_horiz_22_rnd,
		m2d_bilinear16_vert_22_rnd,
		m2d_bilinear_22_22_rnd
};

void (* const m2d_motion_comp_chroma[4])(const uint8_t *src, uint8_t *dst, int stride, int height) = {
		m2d_copy16xn,
		m2d_bilinear_chroma_horiz_rnd,
		m2d_bilinear16_vert_22_rnd,
		m2d_bilinear_chroma_22_22_rnd
};

void (* const m2d_motion_comp_luma_add[4])(const uint8_t *src, uint8_t *dst, int stride, int height) = {
		m2d_copy16xn_add,
		m2d_bilinear16_horiz_22_rnd_add,
		m2d_bilinear16_vert_22_rnd_add,
		m2d_bilinear_22_22_rnd_add
};

void (* const m2d_motion_comp_chroma_add[4])(const uint8_t *src, uint8_t *dst, int stride, int height) = {
		m2d_copy16xn_add,
		m2d_bilinear_chroma_horiz_rnd_add,
		m2d_bilinear16_vert_22_rnd_add,
		m2d_bilinear_chroma_22_22_rnd_add
};

/**
 *\brief Bilinear motion compensation for 16x(16|8) block in case round-bit is active.
 *
 *\param src Pointer to MB of Reference frame.
 *\param dst Pointer to MB of Prediction frame.
 *\param stride Width of frame.
 *\param mvxy pointer to array of motion vector.
 *\param height Number of lines.
 */
void m2d_motion_compensation_luma(const uint8_t *src, uint8_t *dst, int stride, int *mvxy, int height)
{
	src = src + stride * (mvxy[1] >> 1) + (mvxy[0] >> 1);
	m2d_motion_comp_luma[HALFPEL(mvxy[0], mvxy[1])](src, dst, stride, height);
}


/**
 *\brief WMV8 Bilinear motion compensation for 8x8 block in case round-bit is active.
 *
 *\param dst Pointer to MB of Prediction frame.
 *\param src Pointer to MB of Reference frame.
 *\param stride Width of frame.
 *\param mvxy pointer to array of motion vector.
 *\param height Number of lines.
 */
void m2d_motion_compensation_chroma(const uint8_t *src, uint8_t *dst, int stride, int *mvxy, int height)
{
	int mvx = mvxy[0] / 2;
	int mvy = mvxy[1] / 2;
	src = src + stride * (mvy >> 1) + (mvx & ~1);
	m2d_motion_comp_chroma[HALFPEL(mvx, mvy)](src, dst, stride, height);
}


/**
 *\brief Bilinear motion compensation for 16x(16|8) block in case round-bit is active.
 * This is for backword prediction of bi-directional MC.
 *
 *\param src Pointer to MB of Reference frame.
 *\param dst Pointer to MB of Prediction frame.
 *\param stride Width of frame.
 *\param mvxy pointer to array of motion vector.
 *\param height Number of lines.
 */
void m2d_motion_compensation_luma_add(const uint8_t *src, uint8_t *dst, int stride, int *mvxy, int height)
{
	src = src + stride * (mvxy[1] >> 1) + (mvxy[0] >> 1);
	m2d_motion_comp_luma_add[HALFPEL(mvxy[0], mvxy[1])](src, dst, stride, height);
}


/**
 *\brief Bilinear motion compensation for chroma 8x(8|4) block in case round-bit is active.
 * This is for backword prediction of bi-directional MC.
 *
 *\param dst Pointer to MB of Prediction frame.
 *\param src Pointer to MB of Reference frame.
 *\param stride Width of frame.
 *\param mvxy pointer to array of motion vector.
 *\param height Number of lines.
 */
void m2d_motion_compensation_chroma_add(const uint8_t *src, uint8_t *dst, int stride, int *mvxy, int height)
{
	int mvx = mvxy[0] / 2;
	int mvy = mvxy[1] / 2;
	src = src + stride * (mvy >> 1) + (mvx & ~1);
	m2d_motion_comp_chroma_add[HALFPEL(mvx, mvy)](src, dst, stride, height);
}

#else /* FAST_DECODE */

void m2d_motion_compensation_luma(const uint8_t *src, uint8_t *dst, int stride, int *mvxy, int height)
{
	src = src + stride * ((mvxy[1] + 8) >> 4) + ((mvxy[0] + 8) >> 4);
	dst[0] = src[0];
	dst[1] = src[1];
	src += stride;
	dst += stride;
	dst[0] = src[0];
	dst[1] = src[1];
}

void m2d_motion_compensation_chroma(const uint8_t *src, uint8_t *dst, int stride, int *mvxy, int height)
{
	int mvx = mvxy[0] / 2;
	int mvy = mvxy[1] / 2;
	src = src + stride * ((mvy + 8) >> 4) + (((mvx + 4) >> 3) & ~1);
	*(sint16_t *)dst = *(sint16_t *)src;
}

void m2d_copy16xn(const uint8_t *src, uint8_t *dst, int stride, int height)
{
	*(sint16_t *)dst = *(sint16_t *)src;
	if (--height != 0) {
		src += stride;
		dst += stride;
		*(sint16_t *)dst = *(sint16_t *)src;
	}
}

void m2d_motion_compensation_luma_add(const uint8_t *src, uint8_t *dst, int stride, int *mvxy, int height)
{
	src = src + stride * ((mvxy[1] + 8) >> 4) + ((mvxy[0] + 8) >> 4);
	dst[0] = (dst[0] + src[0]) >> 1;
	dst[1] = (dst[1] + src[1]) >> 1;
	src += stride;
	dst += stride;
	dst[0] = (dst[0] + src[0]) >> 1;
	dst[1] = (dst[1] + src[1]) >> 1;
}

void m2d_motion_compensation_chroma_add(const uint8_t *src, uint8_t *dst, int stride, int *mvxy, int height)
{
	int mvx = mvxy[0] / 2;
	int mvy = mvxy[1] / 2;
	src = src + stride * ((mvy + 8) >> 4) + (((mvx + 4) >> 3) & ~1);
	dst[0] = (dst[0] + src[0]) >> 1;
	dst[1] = (dst[1] + src[1]) >> 1;
}

#endif /* FAST_DECODE */

void (* const m2d_motion_compensation[2][2])(const uint8_t *src, uint8_t *dst, int stride, int *mvxy, int height) = {
	{m2d_motion_compensation_luma, m2d_motion_compensation_chroma},
	{m2d_motion_compensation_luma_add, m2d_motion_compensation_chroma_add},
};

#if defined(__RENESAS_VERSION__) && (defined(_SH4ALDSP) || defined(_SH4A)) && defined(UNIT_TEST)

static void (* const *m2d_motion_comp[4])(const uint8_t *src, uint8_t *dst, int stride, int height) = {
	m2d_motion_comp_luma, m2d_motion_comp_chroma,
	m2d_motion_comp_luma_add, m2d_motion_comp_chroma_add
};

static void (* const org_motion_comp[4][3])(const uint8_t *src, uint8_t *dst, int stride, int height) = {
		{org_copy16xn, org_bilinear16_horiz_22_rnd, org_bilinear16_vert_22_rnd},
		{org_copy16xn, org_bilinear_chroma_horiz_rnd, org_bilinear16_vert_22_rnd},
		{org_copy16xn_add, org_bilinear16_horiz_22_rnd_add, org_bilinear16_vert_22_rnd_add},
		{org_copy16xn_add, org_bilinear_chroma_horiz_rnd_add, org_bilinear16_vert_22_rnd_add}
};

#define MB_SQ (MB_LEN * MB_LEN * 2)
static uint8_t blockbuf[3][MB_SQ];
#include <string.h>

static void generate_pattern(uint8_t *src)
{
	memset(src, 0, sizeof(blockbuf));
	int k = 0;
	for (int i = 1; i < MB_LEN; ++i) {
		for (int j = 0; j < MB_LEN * 2; ++j) {
			*src++ = k;
			k = k + i;
		}
	}
}

static int compare_block(const uint8_t *src0, const uint8_t *src1)
{
	int err = 0;
	int i = MB_SQ;
	do {
		if (*src0++ != *src1++) {
			err++;
			break;
		}
	} while (--i); 
	return err;
}

int test_motioncomp()
{
	generate_pattern(blockbuf[0]);
	int err_num = 0;
	for (int block_type = 0; block_type < 4; ++block_type) {
		for (int mc = 0; mc < 3; ++mc) {
			for (int align = 0; align < 4; ++align) {
				m2d_motion_comp[block_type][mc](blockbuf[0] + align, blockbuf[1], MB_LEN * 2, MB_LEN);
				org_motion_comp[block_type][mc](blockbuf[0] + align, blockbuf[2], MB_LEN * 2, MB_LEN);
				err_num += compare_block(blockbuf[1], blockbuf[2]);
			}
		}
	}
	return err_num;
}
#endif
