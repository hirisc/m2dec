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

#include <string.h>
#include <assert.h>

#include "bitio.h"

enum {
	CACHE_BITS = sizeof(cache_t) * 8
};

#ifdef DEBUG_BITIO
static  void BIT_SANITY_CHECK(dec_bits *ths)
{
	int cache_len = ths->cache_len_;
	int loaded_bytes = (unsigned)(cache_len + 7) >> 3;
	const byte_t *buf = ths->buf_ - loaded_bytes;
	cache_t pattern = 0;
	while (loaded_bytes != 0) {
		pattern = (pattern << 8) | *buf++;
		loaded_bytes--;
	}
	pattern <<= CACHE_BITS - cache_len;
	assert(pattern == ths->cache_);
}
#define BIT_SANITY_CHECK(x) bit_sanity_check(x)
#else
#define BIT_SANITY_CHECK(x)
#endif

static int endofbuffer_check(dec_bits *ths);
static void load_bytes(dec_bits *ths, int read_bytes);

const byte_t *dec_bits_load_next(dec_bits *ths, int *data_size)
{
	const byte_t *data;
	if (ths->error_func_(ths->error_arg_) < 0) {
		return 0;
	}
	ths->cache_len_ = 0;
	data = dec_bits_current(ths);
	*data_size = (int)(dec_bits_tail(ths) - data);
	return data;
}

void dec_bits_cachefill(dec_bits *ths) {
	int cache_len;
	intptr_t read_bytes;
	const byte_t *buf;
	intptr_t available_bytes;

	cache_len = ths->cache_len_;
	buf = ths->buf_;
	read_bytes = (CACHE_BITS / 8) - ((unsigned)(cache_len + 7) >> 3);
	available_bytes = ths->buf_tail_ - buf;
	if (available_bytes < read_bytes) {
		if (0 < available_bytes) {
			ths->load_bytes(ths, (int)available_bytes);
		} else if (endofbuffer_check(ths) < 0) {
			return;
		} else {
			ths->load_bytes(ths, (int)read_bytes);
		}
	} else {
		ths->load_bytes(ths, (int)read_bytes);
	}
}

static void load_bytes(dec_bits *ths, int read_bytes)
{
	int cache_len;
	int shift_bits;
	const byte_t *buf;
	cache_t cache;

	cache_len = ths->cache_len_;
	ths->cache_len_ = read_bytes * 8 + cache_len;
	shift_bits = (sizeof(cache) - read_bytes) * 8 - cache_len;
	buf = ths->buf_;
	cache = 0;
	do {
		cache = (cache << 8) | *buf++;
	} while (--read_bytes);
	ths->cache_ = ths->cache_ | (cache << shift_bits);
	ths->buf_ = buf;
}

/** Check whether extra data are available
 */
static int endofbuffer_check(dec_bits *ths)
{
	cache_t cache_ = ths->cache_;
	int cache_len_ = ths->cache_len_;
	int ret = ths->error_func_(ths->error_arg_);
	if ((ret < 0) && (cache_len_ == 0)) {
		longjmp(ths->jmp, 1);
		/* NOTREACHED */
	}
	ths->cache_ = cache_;
	ths->cache_len_ = cache_len_;
	return ret;
}

/**Returns specified number of bits, while read position shall be unchanged.
 *Read bits are LSB-aligned.
 */
uint32_t show_bits(dec_bits *ths, int pat_len)
{
	int cache_len;
/*	assert((unsigned)pat_len <= CACHE_BITS - 8 + (not_aligned_bits(ths))); */
	cache_len = ths->cache_len_;
	BIT_SANITY_CHECK(ths);

	if (cache_len < pat_len) {
		dec_bits_cachefill(ths);
	}
	return (uint32_t)(ths->cache_ >> (CACHE_BITS - pat_len));
}

uint32_t show_onebit(dec_bits *ths)
{
	return (uint32_t)((intptr_t)ths->cache_ < 0);
}

/**Returns specified number of bits as well as read position shall be incremented.
 *Read bits are LSB-aligned.
 */
uint32_t get_bits(dec_bits *ths, int pat_len)
{
	int cache_len;
	cache_t cache;
	assert((unsigned)pat_len <= CACHE_BITS);
	cache_len = ths->cache_len_;

	if (cache_len < pat_len) {
		dec_bits_cachefill(ths);
	}
	cache = ths->cache_;
	ths->cache_ = cache << pat_len;
	ths->cache_len_ -= pat_len;
	return (uint32_t)(cache >> (CACHE_BITS - pat_len));
}

uint32_t get_onebit(dec_bits *ths)
{
	int cache_len = ths->cache_len_;
	cache_t cache;
	if (cache_len <= 0) {
		dec_bits_cachefill(ths);
	}
	cache = ths->cache_;
	ths->cache_ = cache * 2;
	ths->cache_len_--;
	return (uint32_t)((intptr_t)cache < 0);
}

/**Discards specified number of bits. Read position shall be incremented.
 */
void skip_bits(dec_bits *ths, int pat_len)
{
	int cache_len = ths->cache_len_ - pat_len;

	assert(0 < pat_len);
	if (0 <= cache_len) {
		ths->cache_ <<= pat_len;
		ths->cache_len_ = cache_len;
	} else {
		unsigned bits_shortage = (unsigned)(-cache_len);
		if (bits_shortage & ~7) {
			get_bits(ths, bits_shortage);
		} else {
			ths->cache_ = 0;
			ths->cache_len_ = 0;
			dec_bits_cachefill(ths);
			ths->cache_ <<= bits_shortage;
			ths->cache_len_ = CACHE_BITS - bits_shortage;
		}
	}
}

int not_aligned_bits(dec_bits *ths)
{
	return ths->cache_len_ & 7;
}


void byte_align(dec_bits *ths)
{
	int not_aligned = ths->cache_len_ & 7;
	if (not_aligned) {
		skip_bits(ths, not_aligned);
	}
}

/**Discards specified number of bits. Read position shall be incremented.
 */
void skip_bytes(dec_bits *ths, int bytes)
{
	const byte_t *tail, *current;
	int rest;

	assert(0 < bytes);
	tail = dec_bits_tail(ths);
	for (;;) {
		current = dec_bits_current(ths);
		rest = (int)(tail - current);
		if (bytes <= rest) {
			dec_bits_set_data(ths, current + bytes, (size_t)(rest - bytes), ths->id);
			break;
		}
		bytes -= rest;
		dec_bits_set_data(ths, current + rest, (size_t)bytes, ths->id);
		show_bits(ths, 8);
	}
}

/**Returns current read position with adjustment for caching.
 */
const byte_t *dec_bits_current(dec_bits *ths)
{
	return ths->buf_ - ((unsigned)(ths->cache_len_ + 7) >> 3);
}

const byte_t *dec_bits_tail(dec_bits *ths)
{
	return ths->buf_tail_;
}

/**Specify input bitstream.
 */
int dec_bits_set_data(dec_bits *ths, const byte_t *buf, size_t buf_len, void *id)
{
	if ((buf == 0) || ((int)buf_len < 0)) {
		return -1;
	}
	ths->cache_ = 0;
	ths->cache_len_ = 0;
	ths->buf_ = buf;
	ths->buf_tail_ = buf + buf_len;
	ths->buf_head_ = buf;
	ths->id = id;
	return 0;
}

/**Deoder bitstream initializer.
 */
int dec_bits_open(dec_bits *ths, void (*loadbytes_func)(dec_bits *, int bytes))
{
	ths->load_bytes = loadbytes_func ? loadbytes_func : load_bytes;
	return 0;
}

static int error_func_dummy(void *p) {
	return -1;
}

void dec_bits_set_callback(dec_bits *ths, int (*error_func)(void *), void *error_arg)
{
	assert(ths != 0);
	ths->error_func_ = error_func ? error_func : error_func_dummy;
	ths->error_arg_ = error_arg ? error_arg : ths->error_arg_;
}

/**Decoder bitstream finalizer.
 */
void dec_bits_close(dec_bits *ths) {
	assert(ths);
	memset(ths, 0, sizeof(*ths));
}

void dec_bits_tell_error(dec_bits *ths)
{
	longjmp(ths->jmp, 1);
	/* NOTREACHED */
}
