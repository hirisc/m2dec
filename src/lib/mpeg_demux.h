/** Program Stream (PS) DeMUXer.
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

#ifndef _MPEG_DEMUX_H_
#define _MPEG_DEMUX_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "bitio.h"

typedef struct {
	dec_bits *stream;
	int shortage;
	const byte_t *packet_head;
	int packet_len;
	dec_bits stream_i;
} pes_demuxer_t;

int mpeg_demux_init(pes_demuxer_t *dmx, int (*callback_func)(void *), void *arg);
const byte_t *mpeg_demux_get_video(pes_demuxer_t *dmx, int *packet_size_p);
int mpeg_demux_request_set_data(pes_demuxer_t *dmx, const byte_t **indata_p, int *size_p, void *arg);

#ifdef __cplusplus
}
#endif

#endif /* _MPEG_DEMUX_H_ */
