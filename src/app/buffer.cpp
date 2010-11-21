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

#include <memory>
#include <assert.h>
#include "buffer.h"

const int BUF_LEN = 65536;
const int BUF_NUM = 2;

int test_buffer()
{
#ifdef UNIT_TEST
	auto_ptr<BufferManager> bufman(new BufferManager(BUF_NUM, BUF_LEN));
	Buffer *buf_out = bufman.get()->PopFilledBuffer();
	assert(buf_out == 0);
	Buffer *buf_in = bufman.get()->GetEmptyBuffer();
	assert(buf_in != 0);
	int ret = bufman.get()->PushFilledBuffer(buf_in);
	assert(ret == 0);
	buf_out = bufman.get()->PopFilledBuffer();
	assert(buf_out != 0);
	buf_out = bufman.get()->PopFilledBuffer();
	assert(buf_out == 0);
#endif

	return 0;
}

