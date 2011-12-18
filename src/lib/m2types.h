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

#ifdef __cplusplus
}
#endif

#endif /* _M2TYPES_H_ */
