/** Bitstream reader class.
 *  Copyright 2007, 2008 Takayuki Minegishi
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

#ifndef _BIT_IO_H_
#define _BIT_IO_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <setjmp.h>
#ifndef _M_IX86
#include "config.h"
#endif
#include "types.h"

enum {
	CACHE_BITS = sizeof(cache_t) * 8
};

struct dec_bits_t {
	cache_t cache_;
	int cache_len_;
	const byte_t *buf_;
	void (*load_bytes)(struct dec_bits_t *, int bytes);
	const byte_t *buf_tail_;
	const byte_t *buf_head_;
	int (*error_func_)(void *);
	void *error_arg_;
	jmp_buf *jmp;
};

typedef struct dec_bits_t dec_bits;


int dec_bits_open(dec_bits *ths, void (*loadbytes_func)(dec_bits *, int bytes));
void dec_bits_close(dec_bits *ths);
void dec_bits_set_callback(dec_bits *ths, int (*error_func)(void *), void *error_arg, jmp_buf *jmp);

int dec_bits_set_data(dec_bits *ths, const byte_t *buf, size_t buf_len);
cache_t show_bits(dec_bits *ths, int bit_len);
cache_t show_onebit(dec_bits *ths);
cache_t get_bits(dec_bits *ths, int bit_len);
cache_t get_onebit(dec_bits *ths);
void skip_bits(dec_bits *ths, int bit_len);
int not_aligned_bits(dec_bits *ths);
void byte_align(dec_bits *ths);
void skip_bytes(dec_bits *ths, int byte_len);
const byte_t *dec_bits_current(dec_bits *ths);
const byte_t *dec_bits_tail(dec_bits *ths);
const byte_t *dec_bits_load_next(dec_bits *ths, int *data_size);

void dec_bits_dump(dec_bits *ths);

#ifdef __cplusplus
void dec_bits_cachefill(dec_bits *ths);
static inline cache_t get_onebit_inline(dec_bits *ths)
{
	int cache_len = ths->cache_len_;
	cache_t cache;
	if (cache_len <= 0) {
		dec_bits_cachefill(ths);
	}
	cache = ths->cache_;
	ths->cache_ = cache * 2;
	ths->cache_len_--;
	return ((intptr_t)cache < 0);
}
}
#endif

#endif /* _BIT_IO_H_ */
