/** Yet Another H.264 decoder
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

#include <assert.h>
#include <string.h>
#include <setjmp.h>
#include <stdlib.h>
#include <limits.h>
#ifdef _MSC_VER
#include <crtdbg.h>
#define VC_CHECK assert(_CrtCheckMemory());
#else
#define VC_CHECK
#endif
#include "m2d.h"

typedef struct {
	int id;
	dec_bits *stream;
	int (*header_callback)(void *arg, void *seq_id);
	void *header_callback_arg;
	dec_bits stream_i;
} h265d_context;

int h265d_init(h265d_context *h2d, int dpb_max, int (*header_callback)(void *, void *), void *arg) {
	if (!h2d) {
		return -1;
	}
	memset(h2d, 0, sizeof(*h2d));
	h2d->stream = &h2d->stream_i;
	dec_bits_open(h2d->stream, m2d_load_bytes_skip03);
	return 0;
}

dec_bits *h265d_stream_pos(h265d_context *h2d) {
	return h2d->stream;
}

int h265d_get_info(h265d_context *h2d, m2d_info_t *info) {
	int src_width;
	if (!h2d || !info) {
		return -1;
	}
/*	h264d_sps *sps = &h2d->sps_i[h2d->pps_i[h2d->slice_header->pic_parameter_set_id].seq_parameter_set_id];
	info->src_width = src_width = sps->pic_width;
	info->src_height = sps->pic_height;
	info->disp_width = sps->pic_width;
	info->disp_height = sps->pic_height;
	info->frame_num = sps->num_ref_frames + 1;
	for (int i = 0; i < 4; ++i) {
		info->crop[i] = sps->frame_crop[i];
	}
	info->additional_size = sizeof(prev_mb_t) * ((src_width >> 4) + 1)
		+ sizeof(uint32_t) * (src_width >> 2) * 2
		+ (sizeof(deblock_info_t) + (sizeof(h264d_col_mb_t) * 17)) * ((src_width * info->src_height) >> 8)
		+ sizeof(h264d_col_pic_t) * 17;*/
	return 0;
}

int h265d_set_frames(h265d_context *h2d, int num_frame, m2d_frame_t *frame, uint8_t *second_frame, int second_frame_size) {
//	if (!h2d || (num_frame < 3) || (NUM_ARRAY(mb->frame->frames) < (size_t)num_frame) || !frame || !second_frame) {
	if (!h2d || (num_frame < 3) || (16 < (size_t)num_frame) || !frame || !second_frame) {
		return -1;
	}
/*	h265d_mb_current *mb = &h2d->mb_current;
	frames_init(mb, num_frame, frame);
	h2d->slice_header->reorder[0].ref_frames = mb->frame->refs[0];
	h2d->slice_header->reorder[1].ref_frames = mb->frame->refs[1];
	return init_mb_buffer(mb, second_frame, second_frame_size); */
	return 0;
}

typedef enum {
	VPS_NAL = 32,
	SPS_NAL = 33
} h265d_nal_t;

typedef enum {
	WRONG_PARAM = -1,
	OUT_OF_RANGE = -2
} h265d_error_t;

static int dispatch_one_nal(h265d_context *h2d, int code_type) {
	int err = 0;
	dec_bits *st = h2d->stream;

	switch (static_cast<h265d_nal_t>((code_type >> 1) & 31)) {
	case VPS_NAL:
		break;
	case SPS_NAL:
		break;
	default:
		throw OUT_OF_RANGE;
		break;
	}
	return err;
}

int h265d_decode_picture(h265d_context *h2d) {
	if (!h2d) {
		return -1;
	}
	dec_bits* stream = h2d->stream;
/*	h2d->slice_header->first_mb_in_slice = UINT_MAX;*/
	int err = 0;
	int code_type = 0;
	try {
		do {
			if (0 <= (err = m2d_find_mpeg_data(stream))) {
				code_type = get_bits(stream, 8);
				err = dispatch_one_nal(h2d, code_type);
			} else {
				throw OUT_OF_RANGE;
			}
			VC_CHECK;
		} while (err == 0);// || (code_type == SPS_NAL && 0 < err));
	} catch (const h265d_error_t& e) {
		return static_cast<int>(e);
	}
	return err;
}

int h265d_peek_decoded_frame(h265d_context *h2d, m2d_frame_t *frame, int bypass_dpb)
{
/*	h265d_frame_info_t *frm;
	int frame_idx;

	if (!h2d || !frame) {
		return -1;
	}
	frm = h2d->mb_current.frame;
	if (!bypass_dpb) {
		if (frm->dpb.is_ready) {
			frame_idx = dpb_force_peek(&frm->dpb);
		} else {
			frame_idx = frm->dpb.output;
		}
	} else {
		frame_idx = dpb_force_peek(&frm->dpb);
	}
	if (frame_idx < 0) {
		return 0;
	}
	*frame = frm->frames[frame_idx];
	return 1;*/
	return 0;
}

int h265d_get_decoded_frame(h265d_context *h2d, m2d_frame_t *frame, int bypass_dpb)
{
/*	h265d_frame_info_t *frm;
	int frame_idx;

	if (!h2d || !frame) {
		return -1;
	}
	frm = h2d->mb_current.frame;
	if (!bypass_dpb) {
		if (frm->dpb.is_ready) {
			frame_idx = dpb_force_pop(&frm->dpb);
		} else {
			frame_idx = frm->dpb.output;
			frm->dpb.output = -1;
		}
	} else {
		frame_idx = dpb_force_pop(&frm->dpb);
	}
	dump_dpb(&frm->dpb);
	if (frame_idx < 0) {
		return 0;
	}
	*frame = frm->frames[frame_idx];
	return 1;*/
	return 0;
}

static const m2d_func_table_t h265d_func_ = {
	sizeof(h265d_context),
	(int (*)(void *, int, int (*)(void *, void *), void *))h265d_init,
	(dec_bits *(*)(void *))h265d_stream_pos,
	(int (*)(void *, m2d_info_t *))h265d_get_info,
	(int (*)(void *, int, m2d_frame_t *, uint8_t *, int))h265d_set_frames,
	(int (*)(void *))h265d_decode_picture,
	(int (*)(void *, m2d_frame_t *, int))h265d_peek_decoded_frame,
	(int (*)(void *, m2d_frame_t *, int))h265d_get_decoded_frame
};

extern "C" {
extern const m2d_func_table_t * const h265d_func;
}

const m2d_func_table_t * const h265d_func = &h265d_func_;
