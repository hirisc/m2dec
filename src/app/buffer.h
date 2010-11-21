/** MPEG-2 buffer manager for decoder application example
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

#ifndef __BUFFER_H__
#define __BUFFER_H__

#include <assert.h>
#include <deque>
#include <list>
#include <algorithm>

#ifndef __RENESAS_VERSION__
using namespace std;
#endif

class BufferPort {
	bool is_input_;
public:
	BufferPort(bool is_input) : is_input_(is_input) {
	}
};

struct AllocateArray {
	int *operator()(size_t length) {
		return new int[length];
	}
};

struct DeleteArray {
	void operator()(int *ptr) {
		delete[] ptr;
	}
};

struct Buffer {
	size_t length_;
	char *buf_;
	char *begin_;
	char *end_;
	Buffer(int length) : length_(length) {
		begin_ = buf_ = new char[length];
		end_ = buf_ + length;
	}
	~Buffer() {
		delete[] buf_;
	}
};

class Ring {
	int total_;
	int readpos_;
	int writepos_;
	int inc(int pos) {
		++pos;
		return (pos < total_) ? pos : 0;
	}
public:
	Ring(int total) : total_(total), readpos_(0), writepos_(0) {}
	bool empty() {
		return readpos_ == writepos_;
	}
	bool full() {
		return readpos_ == inc(writepos_);
	}
	int get_empty() {
		int rpos = readpos_;
		int wpos = writepos_;
		int nextpos = inc(wpos);
		if (rpos == nextpos) {
			return -1;
		} else if (rpos == wpos) {
			return rpos;
		} else {
			return nextpos;
		}
	}
	int pop_filled() {
		if (empty()) {
			return -1;
		} else {
			int pos = readpos_;
			readpos_ = inc(pos);
			return pos;
		}
	}
	int push() {
		int wpos = writepos_;
		int nextpos = inc(wpos);
		if (readpos_ == nextpos) {
			return -1;
		} else {
			writepos_ = nextpos;
			return wpos;
		}
	}
};

class BufferManager {
	list<BufferPort> ports_;
	size_t number_;
	size_t length_;
	Buffer **buffers_;
	Ring *ring_;
public:
	BufferManager(size_t number, size_t length) : number_(number), length_(length) {
		assert(number && (number < length));
		buffers_ = new Buffer*[number];
		ring_ = new Ring(number);
		for (Buffer **p = buffers_; p != buffers_ + number; ++p) {
			*p = new Buffer(length);
		}
	}
	~BufferManager() {
		ports_.clear();
		for (Buffer **p = buffers_; p != buffers_ + number_; ++p) {
			delete *p;
		}
		delete ring_;
		delete[] buffers_;
	}

	Buffer *PopFilledBuffer() {
		int pos = ring_->pop_filled();
		return (0 <= pos) ? buffers_[pos] : 0;
	}

	Buffer *GetEmptyBuffer() {
		int pos = ring_->get_empty();
		return (0 <= pos) ? buffers_[pos] : 0;
	}

	int PushFilledBuffer(Buffer *buf) {
		int pos = ring_->push();
		if (0 <= pos) {
//			printf("ptr: %08x:%08x\n", buffers_[pos], buf);

			buffers_[pos] = buf;
			return 0;
		} else {
			return -1;
		}
	}
};

#endif /* __BUFFER_H__ */

