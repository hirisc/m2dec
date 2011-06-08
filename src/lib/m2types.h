#ifndef _M2TYPES_H_
#define _M2TYPES_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __RENESAS_VERSION__

//typedef long intptr_t;
//typedef unsigned long uintptr_t;
typedef unsigned int uint32_t;
typedef signed int int32_t;

#elif defined(__GNUC__)

#include <stddef.h>
#include <stdint.h>

#elif defined(_M_IX86)

#undef __LIBM2DEC_API
#ifdef _WINDLL
#define __LIBM2DEC_API __declspec(dllexport)
#elif defined(USE_DLL)
#define __LIBM2DEC_API __declspec(dllimport)
#else
#define __LIBM2DEC_API
#endif
#include <stddef.h>
typedef unsigned int uint32_t;
typedef signed int int32_t;

#endif /* _M_IX86 */



#ifdef __RENESAS_VERSION__
#define __inline
#pragma inline(dec_bits_error, endofbuffer, endofbuffer_check, mix_bit, flush_bits, bit_width)
typedef int intptr_t;
typedef unsigned uintptr_t;
typedef signed char int8_t;
typedef unsigned short uint16_t;
typedef signed short int16_t;
typedef unsigned char uint8_t;

#elif defined(_M_IX86)

#include <stddef.h>
typedef unsigned int uint32_t;
typedef signed char int8_t;
typedef unsigned short uint16_t;
typedef signed short int16_t;
typedef signed int int32_t;
typedef unsigned char uint8_t;

#else

#include <stdint.h>

#endif

typedef unsigned char byte_t;

#ifndef __x86_64__
typedef unsigned long long uint64_t;
#endif
#if SIZEOF_INT_P > 4
#define ENCBIT64
typedef uint64_t cache_t;
#else
typedef unsigned int cache_t;
#endif

#if defined(__RENESAS_VERSION__) && (defined(_SH4ALDSP) || defined(_SH4A))
#include <umachine.h>
#define ABS(a) ((0 <= (a)) ? (a) : -(a))
#define bswap32(a) end_cnvl(a)
#define MUL_EXTEND

#ifdef UNIT_TEST
#pragma inline_asm(read4_unalign)
static uint32_t read4_unalign(const uint32_t *src)
{
	MOVUA	@R4,R0
}
#endif

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
#ifdef __cplusplus
}
#endif

#define CLIP255C(a) (!((a) & ~255) ? (a) : ~((unsigned)(a) >> 16))
#define CLIP255I(a) (!((a) & ~255) ? (a) : ((unsigned)~(a) >> 24))
#define EXTEND_BYTE(a) ((a) = (a) * 0x01010101)

#endif /* _M2TYPES_H_ */
