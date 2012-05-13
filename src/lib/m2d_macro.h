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

#ifndef _M2D_MACRO_H_
#define _M2D_MACRO_H_

#include "m2types.h"

#if defined(__RENESAS_VERSION__) && (defined(_SH4ALDSP) || defined(_SH4A))
#include <umachine.h>
#define bswap32(a) end_cnvl(a)
#define MUL_EXTEND

#pragma inline_asm(read4_unalign)
static uint32_t read4_unalign(const uint32_t *src)
{
	MOVUA	@R4,R0
}

#elif defined(__GNUC__) && defined(__sh__) /* SuperH */

#define MUL_EXTEND
static inline uint32_t bswap32(uint32_t a) {
	__asm__(
		"swap.b %0,%0\n\t"
		"swap.w %0,%0\n\t"
		"swap.b %0,%0\n\t"
		: "=r" (a) : "0" (a) );
	return a;
}
static inline uint32_t read4_unalign(uint32_t *src)
{
	__asm__(
		"movua.l @%0,%0"
		: "=r" (src) : "0" (src) );
	return src;
}
#elif defined(_M_IX86) || defined(__i386__) || defined(__x86_64__)

	#if defined(__i386__) || defined(__x86_64__)
		static uint32_t bswap32(uint32_t a) {
			__asm__( "bswap %0" : "=r" (a) : "0" (a) );
			return a;
		}
	#else
		static uint32_t bswap32(uint32_t data) {
			__asm {
				mov	eax, data
				bswap	eax
			}
		}
	#endif
#define read4_unalign(s) *(uint32_t *)(s)

#else

	#define bswap32(a) ( ( (a) << 24 ) | ( ( (a) & 0xff00 ) << 8 ) | ( ( (a) >> 8 ) & 0xff00 ) | ( ( (a) >> 24 ) & 0xff ) )
static __inline uint32_t read4_unalign(const uint32_t *s) {
	const uint8_t *src = (const uint8_t *)s;
#ifdef _BIG_ENDIAN_
	uint32_t d;
	d = *src++;
	d = (d << 8) | *src++;
	d = (d << 8) | *src++;
	d = (d << 8) | *src;
	return d;
#else
	return (src[3] << 24) | (src[2] << 16) | (src[1] << 8) | src[0];
#endif
}
#endif

#if 1
extern "C" {
extern const uint8_t * const m2d_cliptable;
extern const uint8_t m2d_cliptable_h[];
}
#define CLIP255H(a) (m2d_cliptable_h[(a) & 0x3ff])
#define CLIP255C(a) (m2d_cliptable[(int)(a)])
#define CLIP255I(a) CLIP255C(a)
#else
#define CLIP255C(a) (!((a) & ~255) ? (a) : ~((unsigned)(a) >> 16))
#define CLIP255I(a) (!((a) & ~255) ? (a) : ((unsigned)~(a) >> 24))
#define CLIP255H(s) (!((s) & 0x300) ? (s) : -((((s) >> 9) & 1) ^ 1));
#endif

#define EXTEND_BYTE(a) ((a) = (a) * 0x01010101)

#endif /* _M2D_MACRO_H_ */
