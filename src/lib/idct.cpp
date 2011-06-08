/** iDCT for MPEG-2 or WMV8.
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

#ifndef FAST_DECODE

/**iDCT 8x8 for MPEG-2.
 * This implementation is same as fast idct of reference decoder.
 */

#include "config.h"
#include "m2types.h"
#include "idct.h"

#define W1 2841 /* 2048*sqrt(2)*cos(1*pi/16) */
#define W2 2676 /* 2048*sqrt(2)*cos(2*pi/16) */
#define W3 2408 /* 2048*sqrt(2)*cos(3*pi/16) */
#define W5 1609 /* 2048*sqrt(2)*cos(5*pi/16) */
#define W6 1108 /* 2048*sqrt(2)*cos(6*pi/16) */
#define W7 565  /* 2048*sqrt(2)*cos(7*pi/16) */

#ifdef __RENESAS_VERSION__
#pragma inline(m2d_idct_dconly_horiz8)
#endif

void m2d_idct_add_dconly8xN(int32_t dc, int height, uint8_t *dst_base, int32_t dst_stride);

#if defined(__RENESAS_VERSION__) && (defined(_SH4ALDSP) || defined(_SH4A))
extern "C" void m2d_idct_horizontal(int16_t *src_coef, uint32_t coef_exist);
#ifdef UNIT_TEST
#define m2d_idct_horizontal org_idct_horizontal
#endif
#endif

#if !defined(__RENESAS_VERSION__) || !(defined(_SH4ALDSP) || defined(_SH4A)) || defined(UNIT_TEST)

extern "C" const unsigned long long
m2d_idct_coef[] = {
	0x08000b190a740968LL,
	0x0800064904540235LL,
	0x0800023504540649LL,
	0x080009680a740b19LL,
	0x0080000000000000LL,
};

/**
 *\brief iDCT horizontally for 8x8.
 */
static void m2d_idct_horizontal(int16_t *src_coef, uint32_t coef_exist)
{
//#if defined(__GNUC__) && defined(__i386__)
#if 0 //defined(_M_IX86)
	__asm {
		push	ebp
		push	ebx
		push	edi
		push	esi
		mov	ebx,coef_exist
		mov	ecx,src_coef
		mov	edx,8
L01_horiz:
		shr	ebx, 1
		jc	H01_fullcalc
		movd	mm1, [ecx]
		punpcklwd	mm1, mm1
		punpckldq	mm1, mm1
		psllw	mm1, 3
		movq	[ecx], mm1
		movq	[ecx + 8], mm1
		add	ecx, 16
		add	edx, -1
		jnz	L01_horiz
		jmp	H02_fin
H01_fullcalc:
		movq	mm0, [ecx]
		movq	mm1, [ecx + 8]
		movq	mm2, mm0
		movq	mm3, mm1
		pmullw	mm0, m2d_idct_coef	; W3s3 W2s2 W1s1 W0s0
		pmullw	mm1, m2d_idct_coef + 8	; W7s7 W6s6 W5s5 W0s4
		pmullw	mm2, m2d_idct_coef + 16	; W5s3 W6s2 W7s1 W0s0
		pmullw	mm3, m2d_idct_coef + 24	; W1s7 W2s6 W3s5 W0s4
		paddw	mm0, m2d_idct_coef + 32
		paddw	mm2, m2d_idct_coef + 32
		pshufw	mm0, mm0, 01101100B	; W1s1 W2s2 W3s3 W0s0
		pshufw	mm2, mm2, 01101100B	; W7s1 W6s2 W5s3 W0s0
		paddsw	mm0, mm1	; (W1s1+W7s7)(W2s2+W6s6)(W3s3+W5s5)(W0s0+W0s4)
		psubsw	mm2, mm3	; (W7s1-W1s7)(W6s2-W2s6)(W5s3-W3s5)(W0s0-W0s4)
		movq	mm1, mm0
		movq	mm3, mm2

		punpcklwd	mm0, mm2	; (W5s3-W3s5)(W3s3+W5s5)(W0s0-W0s4)(W0s0+W0s4)
		punpckhwd	mm1, mm2	; (W7s1-W1s7)(W1s1+W7s7)(W6s2-W2s6)(W2s2+W6s6)

;	x0 = 2048s0 - 2048s4 + 128
;	x1 = 2048s0 + 2048s4 + 128
;	x4 = x4 - x6 = ((W7s7+W1s1)-(W3s3+W5s5));
;	x6 = x4 + x6 = (W7s7 + W1s1) + (W3s3 + W5s5);
;	x5 = x5 - x7 = ((W7s1-W1s7)+(W5s3-W3s5));
;	x7 = x5 + x7 = (W7s1 - W1s7) - (W5s3 - W3s5);
;	x5 = (int32_t) ((((W7s7+W1s1)-(W3s3+W5s5)) + ((W7s1-W1s7)+(W5s3-W3s5))) * 181 + 128) >> 8;
;	x4 = (int32_t) ((((W7s7+W1s1)-(W3s3+W5s5)) - ((W7s1-W1s7)+(W5s3-W3s5))) * 181 + 128) >> 8;
;	src_coef[0] = (int32_t) (((2048s0+2048s4+128)+(W6s6+W2s2)) + x6) >> 8;
;	src_coef[1] = (int32_t) (((2048s0-2048s4+128)+(W6s2-W2s6)) + x5) >> 8;
;	src_coef[2] = (int32_t) (((2048s0-2048s4+128)-(W6s2-W2s6)) + x4) >> 8;
;	src_coef[3] = (int32_t) (((2048s0+2048s4+128)-(W6s6+W2s2)) + x7) >> 8;
;	src_coef[4] = (int32_t) (((2048s0+2048s4+128)-(W6s6+W2s2)) - x7) >> 8;
;	src_coef[5] = (int32_t) (((2048s0-2048s4+128)-(W6s2-W2s6)) - x4) >> 8;
;	src_coef[6] = (int32_t) (((2048s0-2048s4+128)+(W6s2-W2s6)) - x5) >> 8;
;	src_coef[7] = (int32_t) (((2048s0+2048s4+128)+(W6s6+W2s2)) - x6) >> 8;
		movq	[ecx], mm1
		movq	[ecx + 8], mm2
		add	ecx,16
		add	edx,-1
		jnz	L01_horiz
H02_fin:
		pop	esi
		pop	edi
		pop	ebx
		pop	ebp
	}
	__asm EMMS;
#else
	int column = 8;
	do {
		int32_t x0;
		if ((coef_exist & 1) == 0) {
			x0 = src_coef[0];
			if (x0) {
				x0 *= 8;
				src_coef[0] = x0;
				src_coef[1] = x0;
				src_coef[2] = x0;
				src_coef[3] = x0;
				src_coef[4] = x0;
				src_coef[5] = x0;
				src_coef[6] = x0;
				src_coef[7] = x0;
			}
		} else {
			int32_t x1, x2, x3, x4, x5, x6, x7, t;
			x0 = *src_coef++ * 2048 + 128;
			x4 = *src_coef++;
			x5 = src_coef[5];
			x3 = *src_coef++;
			x7 = *src_coef++;
			x1 = *src_coef++ * 2048;
			t = x0;
			x0 = x0 - x1;
			x1 = x1 + t;
			if (x4 || x5) {
				t = W7 * (x4 + x5);
				x4 = t + (W1 - W7) * x4;
				x5 = t - (W1 + W7) * x5;
			}
			x6 = *src_coef++;
			if (x7 || x6) {
				t = W3 * (x6 + x7);
				x6 = t - (W3 - W5) * x6;
				x7 = t - (W3 + W5) * x7;
			}
			t = x4;
			x4 -= x6;
			x6 += t;
			t = x5;
			x5 -= x7;
			x7 += t;
			t = x5;
			x5 = (int32_t) ((x4 + x5) * 181 + 128) >> 8;
			x4 = (int32_t) ((x4 - t) * 181 + 128) >> 8;
			x2 = *src_coef;
			src_coef -= 6;
			if (x3 || x2) {
				t = W6 * (x3 + x2);
				x2 = t - (W2 + W6) * x2;
				x3 = t + (W2 - W6) * x3;
			}
			t = x0;
			x0 = x0 - x2;
			x2 = t + x2;
			t = x1;
			x1 = x1 - x3;
			x3 = t + x3;
//	x0 = 2048s0 - 2048s4 + 128
//	x1 = 2048s0 + 2048s4 + 128
//	x4 = W1s1 + W7s7;
//	x5 = W7s1 - W1s7;
//	x6 = W5s5 + W3s3
//	x7 = W3s5 - W5s3
//	x2 = W6s2 - W2s6;
//	x3 = W2s2 + W6s6;
//	x4 = x4 - x6 = (W7s7 + W1s1) - (W3s3 + W5s5);
//	x6 = x4 + x6 = (W7s7 + W1s1) + (W3s3 + W5s5);
//	x5 = x5 - x7 = (W7s1 - W1s7) - (W3s5 - W5s3);
//	x7 = x5 + x7 = (W7s1 - W1s7) + (W3s5 - W5s3);
//	x5 = (int32_t) ((x4 + x5) * 181 + 128) >> 8;
//	x4 = (int32_t) ((x4 - x5) * 181 + 128) >> 8;
//	x0 = x0 - x2
//	x2 = x0 + x2
//	x1 = x1 + x3
//	x3 = x1 - x3
			/* fourth stage */
			src_coef[0] = (int32_t) ((x3 + x6) >> 8);
			src_coef[1] = (int32_t) ((x2 + x5) >> 8);
			src_coef[2] = (int32_t) ((x0 + x4) >> 8);
			src_coef[3] = (int32_t) ((x1 + x7) >> 8);
			src_coef[4] = (int32_t) ((x1 - x7) >> 8);
			src_coef[5] = (int32_t) ((x0 - x4) >> 8);
			src_coef[6] = (int32_t) ((x2 - x5) >> 8);
			src_coef[7] = (int32_t) ((x3 - x6) >> 8);
		}
		src_coef += 8;
		coef_exist >>= 1;
	} while (--column);
#endif
}

#if defined(__RENESAS_VERSION__) && (defined(_SH4ALDSP) || defined(_SH4A)) && defined(UNIT_TEST)
#undef m2d_idct_horizontal

static int16_t coef[2][8 * 8];

static void gen_pattern(int16_t *coef)
{
	for (int y = 0; y < 8; ++y) {
		int k = y * 256 - 1024;
		for (int x = 0; x < 8; ++x) {
			coef[x] = k;
			k = (k * 2) / 3;
		}
		coef += 8;
	}
}

static int compare_coef(const int16_t *src0, const int16_t *src1)
{
	int err = 0;
	int i = 64;
	do {
		if (*src0++ != *src1++) {
			err++;
		}
	} while (--i); 
	return err;
}

#include <string.h>

int test_idct_horizontal()
{
	int err = 0;
	for (int coef_exist = 0; coef_exist < 64; ++coef_exist) {
		gen_pattern(coef[0]);
		memcpy(coef[1], coef[0], sizeof(coef[0]));
		org_idct_horizontal(coef[0], coef_exist);
		m2d_idct_horizontal(coef[1], coef_exist);
		err += compare_coef(coef[0], coef[1]);
	}
	return err;
}

#endif

#endif /* __RENESAS_VERSION__ */

template <int _N, typename _F>
static inline void m2d_idct_vertical(uint8_t *dst_base, int32_t dst_stride, int16_t *src_coef, _F _Store)
{
	int column = 8;
	do {
		int32_t x0, x1, x2, x3, x4, x5, x6, x7, x8;
		uint8_t *dst8;
		int t;

		/* first stage */
		if ((x7 = src_coef[8 * 3]) | (x6 = src_coef[8 * 5])) {
			x8 = W3 * (x6 + x7) + 4;
			x6 = (x8 - (W3 - W5) * x6) >> 3;
			x7 = (x8 - (W3 + W5) * x7) >> 3;
		}
		if ((x4 = src_coef[8 * 1]) | (x5 = src_coef[8 * 7])) {
			x8 = W7 * (x4 + x5) + 4;
			x4 = (x8 + (W1 - W7) * x4) >> 3;
			x5 = (x8 - (W1 + W7) * x5) >> 3;
		}

		/* second stage */
		if ((x3 = src_coef[8 * 2]) | (x2 = src_coef[8 * 6])) {
			x1 = W6 * (x3 + x2) + 4;
			x2 = (x1 - (W2 + W6) * x2) >> 3;
			x3 = (x1 + (W2 - W6) * x3) >> 3;
		}
		x1 = x4 + x6;
		x4 = x4 - x6;
		x6 = x5 + x7;
		x5 = x5 - x7;

		/* third stage */
		x0 = (int32_t)src_coef[8 * 0] * 256 + 8192;
		x7 = (int32_t)src_coef[8 * 4] * 256;
		src_coef++;
		x8 = x0 + x7;
		x0 = x0 - x7;
		x7 = x8 + x3;
		x8 -= x3;
		x3 = x0 + x2;
		x0 -= x2;
		x2 = ((x4 + x5) * 181 + 128) >> 8;
		x4 = ((x4 - x5) * 181 + 128) >> 8;

		/* fourth stage */
		dst8 = dst_base;
		dst_base += _N;
		t = (x7 + x1) >> 14;
		_Store(dst8, t);
		t = (x3 + x2) >> 14;
		dst8 += dst_stride;
		_Store(dst8, t);
		t = (x0 + x4) >> 14;
		dst8 += dst_stride;
		_Store(dst8, t);
		t = (x8 + x6) >> 14;
		dst8 += dst_stride;
		_Store(dst8, t);
		t = (x8 - x6) >> 14;
		dst8 += dst_stride;
		_Store(dst8, t);
		t = (x0 - x4) >> 14;
		dst8 += dst_stride;
		_Store(dst8, t);
		t = (x3 - x2) >> 14;
		dst8 += dst_stride;
		_Store(dst8, t);
		t = (x7 - x1) >> 14;
		dst8 += dst_stride;
		_Store(dst8, t);
	} while (--column);
}


/** Functor that just stores passed value.
 * This is for mono-directional motion compensation.
 */
template <typename T>
class ClipStore {
public:
	void operator ()(T *dst, int val) {
		*dst = CLIP255C(val);
	}
};

/** Functor that stores average value of passed value and destination.
 * This is for bi-directional motion compensation.
 */
template <typename T>
class AddStore {
public:
	void operator ()(T *dst, int val) {
		int t = *dst + val;
		*dst = CLIP255C(t);
	}
};

/**
 *\brief Integer iDCT for intra 8x8 block. Result will be stored directly
 * (without adding) into frame.
 */
void m2d_idct_intra_luma(uint8_t *dst_base, int32_t dst_stride, int16_t *src_coef, uint32_t coef_exist)
{
	m2d_idct_horizontal(src_coef, coef_exist);
	m2d_idct_vertical<1>(dst_base, dst_stride, src_coef, ClipStore<uint8_t>());
}

/**
 *\brief Integer iDCT for intra 8x8 block. Result will be stored directly
 * (without adding) into frame.
 */
void m2d_idct_intra_chroma(uint8_t *dst_base, int32_t dst_stride, int16_t *src_coef, uint32_t coef_exist)
{
	m2d_idct_horizontal(src_coef, coef_exist);
	m2d_idct_vertical<2>(dst_base, dst_stride, src_coef, ClipStore<uint8_t>());
}

/**
 *\brief Integer iDCT for intra 8x8 block. Result will be stored directly
 * (without adding) into frame.
 */
void m2d_idct_inter_luma(uint8_t *dst_base, int32_t dst_stride, int16_t *src_coef, uint32_t coef_exist)
{
	m2d_idct_horizontal(src_coef, coef_exist);
	m2d_idct_vertical<1>(dst_base, dst_stride, src_coef, AddStore<uint8_t>());
}

/**
 *\brief Integer iDCT for intra 8x8 block. Result will be stored directly
 * (without adding) into frame.
 */
void m2d_idct_inter_chroma(uint8_t *dst_base, int32_t dst_stride, int16_t *src_coef, uint32_t coef_exist)
{
	m2d_idct_horizontal(src_coef, coef_exist);
	m2d_idct_vertical<2>(dst_base, dst_stride, src_coef, AddStore<uint8_t>());
}


#ifdef __RENESAS_VERSION__
#pragma inline_asm(m2d_idct_add_dconly8xN)
void m2d_idct_add_dconly8xN(int32_t dc, int height, uint8_t *dst_base, int32_t dst_stride)
{
	MOV.L	R8,@-R15
	MOV.L	R9,@-R15
	MOV.L	R10,@-R15
	MOV.L	R11,@-R15
	.align	4
	MOV.L	R12,@-R15
	MOVA	@(?Hdata - ?Hcommon,PC),R0
?Hcommon:
	CMP/PZ	R4
	MOV.L	@R0+,R1
	MOVT	R8
	MOV	#-1,R3
	BT/S	?Hadd0
	EXTU.B	R3,R3	; 255
	NEG	R4,R4
?Hadd0:
	CMP/HI	R4,R3
	MOVT	R2	; 1 if (dc <= 255)
	ADD	R3,R2
	AND	R3,R2	; 0 or 255
	OR	R4,R2	; saturated dc
	MOV.L	@R0,R0
	TST	R8,R8
	MUL.L	R1,R2
	MOV	#-7,R12
	NOT	R0,R1	; H'80808080
	BF/S	?Ladd
	STS	MACL,R11
?Lsub:
	MOV.L	@R6,R4
			MOV.L	@(4,R6),R10
	NOT	R4,R2
			NOT	R10,R8
	MOV	R2,R3
			MOV	R8,R9
	XOR	R11,R2
	AND	R11,R3
	SHLR	R2
			XOR	R11,R8
			AND	R11,R9
	AND	R0,R2
			SHLR	R8
	ADD	R3,R2
			AND	R0,R8
	AND	R1,R2
			ADD	R9,R8
	MOV	R2,R3
			AND	R1,R8
	SHLD	R12,R2
			MOV	R8,R9
	ADD	R3,R3
			SHLD	R12,R8
	SUB	R2,R3	; mask
			ADD	R9,R9
	OR	R3,R4
			SUB	R8,R9
	OR	R11,R3
			OR	R9,R10
	SUB	R3,R4
			OR	R11,R9
	MOV.L	R4,@R6
			SUB	R9,R10
	DT	R5
	MOV.L	R10,@(4,R6)
	BF/S	?Lsub
	ADD	R7,R6
	BRA	?Hend
	NOP
	.align	4
?Hdata:
	.data.l	H'01010101,H'7f7f7f7f
?Ladd:
	MOV.L	@R6,R2
			MOV.L	@(4,R6),R8
	MOV	R2,R3
			MOV	R8,R9
	MOV	R2,R4
			MOV	R8,R10
	XOR	R11,R2
	AND	R11,R3
	SHLR	R2
			XOR	R11,R8
			AND	R11,R9
			SHLR	R8
	AND	R0,R2
			AND	R0,R8
	ADD	R3,R2
			ADD	R9,R8
	AND	R1,R2
			AND	R1,R8
	MOV	R2,R3
			MOV	R8,R9
	SHLD	R12,R2
			SHLD	R12,R8
	ADD	R3,R3
			ADD	R9,R9
	SUB	R2,R3	; mask
			SUB	R8,R9
	ADD	R11,R4
			ADD	R11,R10
	SUB	R3,R4
			SUB	R9,R10
	OR	R3,R4
			OR	R9,R10
	MOV.L	R4,@R6
	DT	R5
	MOV.L	R10,@(4,R6)
	BF/S	?Ladd
	ADD	R7,R6
?Hend:
	MOV.L	@R15+,R12
	MOV.L	@R15+,R11
	MOV.L	@R15+,R10
	MOV.L	@R15+,R9
	MOV.L	@R15+,R8
}
#else /* __RENESAS_VERSION__ */

/**Saturated addition of unsigned packed four bytes.
 *Same as "paddb" in Intel MMX instruction.
 */
static inline uint32_t ADD_CLIP(uint32_t a, uint32_t b)
{
	uint32_t mask = ((a & b) + (((a ^ b) >> 1) & 0x7f7f7f7f)) & ~0x7f7f7f7f;
	mask = (mask << 1) - (mask >> 7);
	return ((a + b) - mask) | mask;
}

/**Saturated subtraction of unsigned packed four bytes.
 *Same as "psubb" in Intel MMX instruction.
 */
static inline uint32_t SUB_CLIP(uint32_t a, uint32_t b)
{
	uint32_t mask = ((~a & b) + (((~a ^ b) >> 1) & 0x7f7f7f7f)) & ~0x7f7f7f7f;
	mask = (mask << 1) - (mask >> 7);
	return (a | mask) - (b | mask);
}

void m2d_idct_add_dconly8xN(int32_t dc, int height, uint8_t *dst_base, int32_t dst_stride)
{
	if (dc <= 0) {
		dc = -dc;
		dc = CLIP255I(dc);
		dc = EXTEND_BYTE(dc);
		do {
			uint32_t s0, s1;
			s0 = *(uint32_t *)dst_base;
			s1 = *(uint32_t *)(dst_base + 4);
			*(uint32_t *)dst_base = SUB_CLIP(s0, dc);
			*(uint32_t *)(dst_base + 4) = SUB_CLIP(s1, dc);
			dst_base += dst_stride;
		} while (--height);
	} else {
		dc = CLIP255I(dc);
		dc = EXTEND_BYTE(dc);
		do {
			uint32_t s0, s1;
			s0 = *(uint32_t *)dst_base;
			s1 = *(uint32_t *)(dst_base + 4);
			*(uint32_t *)dst_base = ADD_CLIP(s0, dc);
			*(uint32_t *)(dst_base + 4) = ADD_CLIP(s1, dc);
			dst_base += dst_stride;
		} while (--height);
	}
}
#endif /* __RENESAS_VERSION__ */

#endif /* FAST_DECODE */
