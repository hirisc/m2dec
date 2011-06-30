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

#include <assert.h>
#include <string.h>
#include "mpeg_demux.h"
#include "mpeg2.h"

static inline int skip_packet(pes_demuxer_t *dmx);
static int video_element_packet(pes_demuxer_t *dmx);

static int process_system_packet(pes_demuxer_t *dmx)
{
	int err;
	int code;

	dec_bits *bit = dmx->stream;
	code = get_bits(bit, 8);
	err = 0;
	switch (code) {
	case 0xb9:
		return -1; /* end of iso */
		/* NOT REACHED */
		break;
	case 0xba:
		skip_bytes(bit, 8);
		break;
	case 0xbd:
		err = skip_packet(dmx);
		break;
	case 0xc0:
		err = skip_packet(dmx);
		break;
	case 0xe0:
		err = video_element_packet(dmx);
		break;
	default:
		err = skip_packet(dmx);
		break;
	}
	return err;
}

static inline int skip_packet(pes_demuxer_t *dmx)
{
	int block_size = get_bits(dmx->stream, 2 * 8);
	skip_bytes(dmx->stream, block_size);
	return (0 < block_size) ? 0 : -1;
}

static int video_element_packet(pes_demuxer_t *dmx)
{
	dec_bits *bit = dmx->stream;
	int val = get_bits(bit, 3 * 8);
	const byte_t *packet_tail;

	packet_tail = dec_bits_current(bit) + ((unsigned)val >> 8) - 1;
	dmx->shortage = (int)(packet_tail - dec_bits_tail(bit));
	if ((val & 0xc0) == 0x80) {
		val = get_bits(bit, 2 * 8) & 255;
		if (val) {
			skip_bytes(bit, val);
		}
	} else {
		val &= 255; 
		while (val == 255) {
			val = get_bits(bit, 8);
		}
		if (val & 0xc0) {
			if (val & 0x80) {
				return -1;
			}
			val = get_bits(bit, 2 * 8) & 255;
		}
		if (0x30 <= val) {
			if (val & 0xc0) {
				return -1;
			}
			skip_bytes(bit, 9);
		} else if (val & 0x20) {
			skip_bytes(bit, 4);
		} else if (val != 0x0f) {
			return -1;
		}
	}
	dmx->packet_head = dec_bits_current(bit);
	if (0 < dmx->shortage) {
		dmx->packet_len = (int)(dec_bits_tail(bit) - dmx->packet_head);
	} else {
		dmx->packet_len = (int)(packet_tail - dmx->packet_head);
	}
	skip_bytes(bit, dmx->packet_len);
	return 1;
}

int mpeg_demux_init(pes_demuxer_t *dmx, int (*callback_func)(void *), void *arg)
{
	if (dmx == 0) {
		return -1;
	}
	memset((void *)dmx, 0, sizeof(*dmx));
	dec_bits_open(dmx->stream = &dmx->stream_i, 0);
	dec_bits_set_callback(dmx->stream, callback_func, arg);
	return 0;
}

const byte_t *mpeg_demux_get_video(pes_demuxer_t *dmx, int *packet_size_p)
{
	dec_bits *stream;
	int err;

	stream = dmx->stream;
	if (setjmp(stream->jmp) != 0) {
		return 0;
	}
	if (0 < dmx->shortage) {
		const byte_t *current;
		int rest, min;

		show_bits(stream, 8);
		current = dec_bits_current(stream);
		rest = (int)(dec_bits_tail(stream) - current);
		min = (rest < dmx->shortage) ? rest : dmx->shortage;
		skip_bytes(stream, min);
		*packet_size_p = min;
		dmx->shortage -= rest;
		return current;
	} else {
		err = 0;
		do {
			if (err < 0 || m2d_find_mpeg_data(stream) < 0) {
				*packet_size_p = 0;
				return 0;
			}
			err = process_system_packet(dmx);
		} while (err <= 0);
		*packet_size_p = dmx->packet_len;
		return dmx->packet_head;
	}
}

