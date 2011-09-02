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

#ifdef __RENESAS_VERSION__
extern "C" {
void exit(int);
}
#endif

#include <functional>
#include <algorithm>
#include "bitio.h"
#include "h264.h"
#include "h264vld.h"

#define MIN(a, b) ((a) <= (b) ? (a) : (b))
static inline int ABS(int a) {
	return (0 <= a) ? a : -a;
}

#define READ_UE_RANGE(dst, st, max) {uint32_t t = ue_golomb(st); (dst) = t; if ((max) < t) return -1;}
#define READ_SE_RANGE(dst, st, min, max) {int32_t t = se_golomb(st); (dst) = t; if (t < (min) || (max) < t) return -1;}

#define UNPACK(a, num) (((a) >> ((num) * 4)) & 15)
#define PACK(a, val, num) ((a) | ((val) << ((num) * 4)))

static inline cache_t get_bits32(dec_bits *ths, int bit_len)
{
	if (bit_len <= 24) {
		return get_bits(ths, bit_len);
	} else {
		int rest = bit_len - 24;
		return (get_bits(ths, 24) << rest) | get_bits(ths, rest);
	}
}

static uint32_t ue_golomb(dec_bits *str)
{
	int bits, rest;
	int i;

	if (get_onebit_inline(str)) {
		return 0;
	}
	rest = 0;
	i = 16;
	do {
		bits = get_bits(str, 2);
		switch (bits) {
		case 0:
			rest += 2;
			break;
		case 1:
			return get_bits32(str, rest + 2) + ((bits << 2) << rest) - 1;
			/* NOTREACHED */
			break;
		case 2:
			/* FALLTHROUGH */
		case 3:
			return (rest ? get_bits32(str, rest) : 0) + (bits << rest) - 1;
			/* NOTREACHED */
			break;
		}
	} while (--i);
	return 0;
}

static int32_t se_golomb(dec_bits *stream)
{
	int32_t ue = ue_golomb(stream);
	int32_t t = (ue + 1) >> 1;
	return (ue & 1) ? t : -t;
}

static const int8_t me_golomb_lut[2][48] = {
	{
		47, 31, 15, 0, 23, 27, 29, 30,
		7, 11, 13, 14, 39, 43, 45, 46,
		16, 3, 5, 10, 12, 19, 21, 26,
		28, 35, 37, 42, 44, 1, 2, 4,
		8, 17, 18, 20, 24, 6, 9, 22,
		25, 32, 33, 34, 36, 40, 38, 41
	},
	{
		0, 16, 1, 2, 4, 8, 32, 3,
		5, 10, 12, 15, 47, 7, 11, 13,
		14, 6, 9, 31, 35, 37, 42, 44,
		33, 34, 36, 40, 39, 43, 45, 46,
		17, 18, 20, 24, 19, 21, 26, 28,
		23, 27, 29, 30, 22, 25, 38, 41
	}
};

static inline int32_t me_golomb(dec_bits *stream, const int8_t *me_lut)
{
	uint32_t ue = ue_golomb(stream);
	return me_lut[(ue < 48) ? ue : 0];
}

static inline int32_t te_golomb(dec_bits *stream, int range)
{
	if (range == 1) {
		return get_onebit_inline(stream) ^ 1;
	} else {
		int32_t ue = ue_golomb(stream);
		return (ue <= range) ? ue : range;
	}
}

static uint32_t get_32bits(dec_bits *stream)
{
	uint32_t t = get_bits(stream, 16);
	return (t << 16) | get_bits(stream, 16);
}


static void hrd_parameters(hrd_parameters_t *hrd, dec_bits *stream)
{
	int max;

	hrd->cpb_cnt_minus1 = max = ue_golomb(stream);
	hrd->bit_rate_scale = get_bits(stream, 4);
	hrd->cpb_size_scale = get_bits(stream, 4);
	hrd->cbr_flag = 0;
	for (int i = 0; i <= max; ++i) {
		hrd->bit_rate_value_minus1[i] = ue_golomb(stream);
		hrd->cpb_size_value_minus1[i] = ue_golomb(stream);
		hrd->cbr_flag |= get_onebit(stream) << i;
	}
	hrd->initial_cpb_removal_delay_length_minus1 = get_bits(stream, 5);
	hrd->cpb_removal_delay_length_minus1 = get_bits(stream, 5);
	hrd->dpb_output_delay_length_minus1 = get_bits(stream, 5);
	hrd->time_offset_length = get_bits(stream, 5);
}

static int vui_parameters(vui_parameters_t *vui, dec_bits *stream)
{
	if ((vui->aspect_ratio_info_present_flag = get_onebit(stream)) != 0) {
		if ((vui->aspect_ratio_idc = get_bits(stream, 8)) == EXTENDED_SAR) {
			vui->sar_width = get_bits(stream, 16);
			vui->sar_height = get_bits(stream, 16);
		}
	}
	if ((vui->overscan_info_present_flag = get_onebit(stream)) != 0) {
		vui->overscan_appropriate_flag = get_onebit(stream);
	}
	if ((vui->video_signal_type_present_flag = get_onebit(stream)) != 0) {
		vui->video_format = get_bits(stream, 3);
		vui->video_full_range_flag = get_onebit(stream);
		if ((vui->colour_description_present_flag = get_onebit(stream)) != 0) {
			vui->colour_primaries = get_bits(stream, 8);
			vui->transfer_characteristics = get_bits(stream, 8);
			vui->matrix_coefficients = get_bits(stream, 8);
		}
	}
	if ((vui->chroma_loc_info_present_flag = get_onebit(stream)) != 0) {
		vui->chroma_sample_loc_type_top_field = ue_golomb(stream);
		vui->chroma_sample_loc_type_bottom_field = ue_golomb(stream);
	}
	if ((vui->timing_info_present_flag = get_onebit(stream)) != 0) {
		vui->num_units_in_tick = get_32bits(stream);
		vui->time_scale = get_32bits(stream);
		vui->fixed_frame_rate_flag = get_onebit(stream);
	}
	if ((vui->nal_hrd_parameters_present_flag = get_onebit(stream)) != 0) {
		hrd_parameters(&vui->nal_hrd_parameters, stream);
	}
	if ((vui->vcl_hrd_parameters_present_flag = get_onebit(stream)) != 0) {
		hrd_parameters(&vui->vcl_hrd_parameters, stream);
	}
	if (vui->nal_hrd_parameters_present_flag || vui->vcl_hrd_parameters_present_flag) {
		vui->low_delay_hrd_flag = get_onebit(stream);
	}
	vui->pic_struct_present_flag = get_onebit(stream);
	if ((vui->bitstream_restriction_flag = get_onebit(stream)) != 0) {
		vui->motion_vectors_over_pic_boundaries_flag = get_onebit(stream);
		vui->max_bytes_per_pic_denom = ue_golomb(stream);
		vui->max_bits_per_mb_denom = ue_golomb(stream);
		vui->log2_max_mv_length_horizontal = ue_golomb(stream);
		vui->log2_max_mv_length_vertical = ue_golomb(stream);
		vui->num_reorder_frames = ue_golomb(stream);
		vui->max_dec_frame_buffering = ue_golomb(stream);
	}
	return 0;
}

static void read_poc_type1_cycle(h264d_sps *sps, dec_bits *st, int max_cycles)
{
	int32_t *offset = sps->offset_for_ref_frame;
	int32_t delta = 0;
	for (int i = 0; i < max_cycles; ++i) {
		delta += se_golomb(st);
		offset[i] = delta;
	}
}

static inline int max_dpb_mbs(int profile_idc, int level_idc, unsigned constrained_set)
{
	int max_dpb;
	if (profile_idc == 100) {
		if (level_idc == 9) {
			level_idc = 10; /* 1b */
		}
	} else {
		if ((level_idc == 10) && (constrained_set & 16)) {
			level_idc = 10; /* 1b */
		}
	}
	switch (level_idc) {
	case 10:
		max_dpb = 396;
		break;
	case 11:
		max_dpb = 900;
		break;
	case 12:
	case 13:
	case 20:
		max_dpb = 2376;
		break;
	case 21:
		max_dpb = 4752;
		break;
	case 22:
	case 30:
		max_dpb = 8100;
		break;
	case 31:
		max_dpb = 18000;
		break;
	case 32:
		max_dpb = 20480;
		break;
	case 40:
	case 41:
		max_dpb = 32768;
		break;
	case 42:
		max_dpb = 34816;
		break;
	case 50:
		max_dpb = 110400;
		break;
	case 51:
		max_dpb = 184320;
		break;
	default:
		max_dpb = -1;
		break;
	}
	return max_dpb;
}

static int read_seq_parameter_set(h264d_sps *sps, dec_bits *stream)
{
	uint32_t sps_profile_idc, sps_constraint_set_flag;
	uint32_t sps_level_idc, sps_id;
	uint32_t tmp;

	sps_profile_idc = get_bits(stream, 8);
	sps_constraint_set_flag = get_bits(stream, 8);
	sps_level_idc = get_bits(stream, 8);
	READ_UE_RANGE(sps_id, stream, 31);
	sps += sps_id;
	sps->profile_idc = sps_profile_idc;
	sps->constraint_set_flag = sps_constraint_set_flag;
	sps->level_idc = sps_level_idc;
	if ((80 <= sps->profile_idc) && (sps->profile_idc <= 83)) {
	}
	READ_UE_RANGE(tmp, stream, 27);
	sps->log2_max_frame_num = tmp + 4;
	READ_UE_RANGE(sps->poc_type, stream, 2);
	if (sps->poc_type == 0) {
		READ_UE_RANGE(tmp, stream, 27);
		sps->log2_max_poc_lsb = tmp + 4;
	} else if (sps->poc_type == 1) {
		sps->delta_pic_order_always_zero_flag = get_onebit(stream);
		sps->offset_for_non_ref_pic = se_golomb(stream);
		sps->offset_for_top_to_bottom_field = se_golomb(stream);
		READ_UE_RANGE(sps->num_ref_frames_in_pic_order_cnt_cycle, stream, 255);
		read_poc_type1_cycle(sps, stream, sps->num_ref_frames_in_pic_order_cnt_cycle);
	}
	READ_UE_RANGE(sps->num_ref_frames, stream, 16);
	sps->gaps_in_frame_num_value_allowed_flag = get_onebit(stream);
	sps->pic_width = (ue_golomb(stream) + 1) * 16;
	sps->pic_height = (ue_golomb(stream) + 1) * 16;
	sps->max_dpb_in_mbs = max_dpb_mbs(sps->profile_idc, sps->level_idc, sps->constraint_set_flag);
	if ((sps->frame_mbs_only_flag = get_onebit(stream)) == 0) {
		sps->mb_adaptive_frame_field_flag = get_onebit(stream);
	}
	sps->direct_8x8_inference_flag = get_onebit(stream);
	if ((sps->frame_cropping_flag = get_onebit(stream)) != 0) {
		for (int i = 0; i < 4; ++i) {
			sps->frame_crop[i] = ue_golomb(stream) * 2;
		}
	} else {
		memset(sps->frame_crop, 0, sizeof(sps->frame_crop));
	}
	if ((sps->vui_parameters_present_flag = get_onebit(stream)) != 0) {
		int err = vui_parameters(&sps->vui, stream);
		if (err < 0) {
			return err;
		}
	}
	return sps_id;
}

static int get_sei_message_size(dec_bits *stream);
static int skip_sei_message(dec_bits *stream);

static int skip_sei(dec_bits *stream)
{
	do {
		skip_sei_message(stream);
		byte_align(stream);
	} while (show_bits(stream, 8) != 0);
	return 0;
}

static int skip_sei_message(dec_bits *stream)
{
	int d;
	d = get_sei_message_size(stream);
	skip_bits(stream, d * 8);
	d = get_sei_message_size(stream);
	skip_bits(stream, d * 8);
	return 0;
}

static int get_sei_message_size(dec_bits *stream)
{
	int c;
	int d = -255;
	do {
		c = get_bits(stream, 8);
		d += 255;
	} while (c == 0xff);
	return d + c;
}

static int read_pic_parameter_set(h264d_pps *pps, dec_bits *stream)
{
	uint8_t pps_id;
	int tmp;

	READ_UE_RANGE(pps_id, stream, 255);
	pps += pps_id;
	READ_UE_RANGE(pps->seq_parameter_set_id, stream, 31);
	pps->entropy_coding_mode_flag = get_onebit(stream);
	pps->pic_order_present_flag = get_onebit(stream);
	if (0 < (pps->num_slice_groups_minus1 = ue_golomb(stream))) {
		/* FMO not implemented */
		return -1;
	}
	READ_UE_RANGE(pps->num_ref_idx_l0_active_minus1, stream, 31);
	READ_UE_RANGE(pps->num_ref_idx_l1_active_minus1, stream, 31);
	pps->weighted_pred_flag = get_onebit(stream);
	pps->weighted_bipred_idc = get_bits(stream, 2);
	READ_SE_RANGE(tmp, stream, -26, 25);
	pps->pic_init_qp = tmp + 26;
	READ_SE_RANGE(tmp, stream, -26, 25);
	pps->pic_init_qs = tmp + 26;
	READ_SE_RANGE(pps->chroma_qp_index[0], stream, -12, 12);
	pps->deblocking_filter_control_present_flag = get_onebit(stream);
	pps->constrained_intra_pred_flag = get_onebit(stream);
	pps->redundant_pic_cnt_present_flag = get_onebit(stream);

	return 0;
}

void h264d_load_bytes_skip03(dec_bits *ths, int read_bytes)
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
		byte_t c = *buf++;
		if (c == 3) {
			if (buf[-2] == 0 && buf[-3] == 0) {
				c = *buf++;
			}
		}
		cache = (cache << 8) | c;
	} while (--read_bytes);
	ths->cache_ = ths->cache_ | (cache << shift_bits);
	ths->buf_ = buf;
}

static inline void dpb_init(h264d_dpb_t *dpb, int maxsize);

static int header_dummyfunc(void *arg, int seq_id) {return 0;}

int h264d_init(h264d_context *h2d, int dpb_max, int (*header_callback)(void *, int), void *arg)
{
	if (!h2d) {
		return -1;
	}
	memset(h2d, 0, sizeof(*h2d));
	h2d->stream = &h2d->stream_i;
	h2d->slice_header = &h2d->slice_header_i;
	h2d->mb_current.bdirect = &h2d->mb_current.bdirect_i;
	h2d->mb_current.frame = &h2d->mb_current.frame_i;
	h2d->mb_current.cabac = &h2d->mb_current.cabac_i;
	h2d->header_callback = header_callback ? header_callback : header_dummyfunc;
	h2d->header_callback_arg = arg;
	h2d->mb_current.num_ref_idx_lx_active_minus1[0] = &h2d->slice_header->num_ref_idx_lx_active_minus1[0];
	h2d->mb_current.num_ref_idx_lx_active_minus1[1] = &h2d->slice_header->num_ref_idx_lx_active_minus1[1];
	dpb_init(&h2d->mb_current.frame->dpb, dpb_max);
	dec_bits_open(h2d->stream, h264d_load_bytes_skip03);
	return 0;
}

static dec_bits *h264d_stream_pos(h264d_context *h2d)
{
	return &h2d->stream_i;
}

static void set_mb_size(h264d_mb_current *mb, int width, int height);

int h264d_read_header(h264d_context *h2d, const byte_t *data, size_t len)
{
	dec_bits *st;
	int nal_type;
	int err;
	int sps_id;

	st = h2d->stream;
	dec_bits_open(st, h264d_load_bytes_skip03);
	err = dec_bits_set_data(st, data, len);
	if (err < 0) {
		return err;
	}
	if (setjmp(st->jmp) != 0) {
		return 0;
	}
	do {
		err = m2d_find_mpeg_data(st);
		if (err < 0) {
			return err;
		}
		nal_type = get_bits(st, 8) & 31;
	} while (nal_type != SPS_NAL);
	sps_id = read_seq_parameter_set(h2d->sps_i, st);
	if (sps_id < 0) {
		return sps_id;
	}
	set_mb_size(&h2d->mb_current, h2d->sps_i[sps_id].pic_width, h2d->sps_i[sps_id].pic_height);
	return 0;
}

int h264d_get_info(h264d_context *h2d, m2d_info_t *info)
{
	int src_width;
	if (!h2d || !info) {
		return -1;
	}
	h264d_sps *sps = &h2d->sps_i[h2d->pps_i[h2d->slice_header->pic_parameter_set_id].seq_parameter_set_id];
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
		+ sizeof(h264d_col_pic_t) * 17;
	return 0;
}

static int init_mb_buffer(h264d_mb_current *mb, uint8_t *buffer, int len)
{
	uint8_t *src = buffer;
	mb->mb_base = (prev_mb_t *)src;
	src += sizeof(*mb->mb_base) * (mb->max_x + 1);
	mb->top4x4pred_base = (int32_t *)src;
	src += sizeof(*mb->top4x4pred_base) * mb->max_x;
	mb->top4x4coef_base = (int32_t *)src;
	src += sizeof(*mb->top4x4coef_base) * mb->max_x;
	mb->deblock_base = (deblock_info_t *)src;
	int mb_num = mb->max_x * mb->max_y;
	src += sizeof(*mb->deblock_base) * mb_num;
	for (int i = 0; i < 16; ++i) {
		mb->frame->refs[1][i].col = (h264d_col_pic_t *)src;
		src += sizeof(h264d_col_pic_t) + sizeof(h264d_col_mb_t) * (mb_num - 1);
	}
	mb->frame->curr_col = (h264d_col_pic_t *)src;
	src += sizeof(h264d_col_pic_t) + sizeof(h264d_col_mb_t) * (mb_num - 1);
	return (uintptr_t)(buffer + len) < (uintptr_t)src ? -1 : 0;
}

static void build_4x4offset_table(int dst[], int stride) {
	int offset = 0;
	for (int i = 0; i < 4; ++i) {
		dst[0] = offset;
		dst[1] = offset + 4;
		dst[2] = offset + stride * 4;
		dst[3] = offset + (stride + 1) * 4;
		dst += 4;
		offset += (i & 1) ? (stride - 1) * 8 : 8;
	}
}

static void set_mb_size(h264d_mb_current *mb, int width, int height)
{
	mb->max_x = width >> 4;
	mb->max_y = height >> 4;
	build_4x4offset_table(mb->offset4x4, width);
}

/**Invoked just before each slice_data.
 */
static void set_mb_pos(h264d_mb_current *mb, int mbpos)
{
	int x, y;
	div_t d;

	d = div(mbpos, mb->max_x);
	mb->y = y = d.quot;
	mb->x = x = d.rem;
	mb->luma = mb->frame->curr_luma + y * mb->max_x * 16 * 16 + x * 16;
	mb->chroma = mb->frame->curr_chroma + y * mb->max_x * 16 * 8 + x * 16;
	mb->firstline = mb->max_x;
	mb->left4x4pred = 0x22222222;
	mb->prev_qp_delta = 0;
	memset(mb->top4x4pred_base, 0x22, mb->max_x * sizeof(*mb->top4x4pred_base));
	mb->top4x4pred = mb->top4x4pred_base + x;
	mb->top4x4coef = mb->top4x4coef_base + x;
	mb->deblock_curr = mb->deblock_base + mbpos;
	mb->left4x4pred = 0;
	*mb->top4x4pred = 0;
	mb->top4x4inter = mb->mb_base + 1 + x;
	mb->left4x4inter = mb->mb_base;
	mb->col_curr = mb->frame->curr_col->col_mb + y * mb->max_x + x;
	mb->cbf = 0;
}

static inline uint16_t cbf_top(uint32_t cbf)
{
	return ((cbf >> 16) & 0x700) | ((cbf >> 14) & 0xc0) | ((cbf >> 12) & 0x3c) | ((cbf >> 10) & 3);
}

static inline uint16_t cbf_left(uint32_t cbf)
{
	return ((cbf >> 16) & 0x600) | ((cbf >> 15) & 0x100) | ((cbf >> 14) & 0x80) | ((cbf >> 13) & 0x40) | ((cbf >> 12) & 0x38) | ((cbf >> 11) & 4) | ((cbf >> 6) & 2) | ((cbf >> 5) & 1);
}

static int increment_mb_pos(h264d_mb_current *mb)
{
	int mb_type;
	int x;

	mb_type = mb->type;
	mb->top4x4inter->type = mb_type;
	mb->left4x4inter->type = mb_type;
	mb->top4x4inter->cbp = mb->cbp;
	mb->top4x4inter->cbf = cbf_top(mb->cbf);
	mb->left4x4inter->cbp = mb->cbp;
	mb->left4x4inter->cbf = cbf_left(mb->cbf);
	mb->top4x4inter->chroma_pred_mode = mb->chroma_pred_mode;
	mb->left4x4inter->chroma_pred_mode = mb->chroma_pred_mode;
	mb->cbf = 0;
	x = mb->x + 1;
	mb->top4x4pred++;
	mb->top4x4coef++;
	mb->top4x4inter++;
	mb->col_curr++;
	mb->deblock_curr++;
	mb->luma += 16;
	mb->chroma += 16;
	if (mb->max_x <= x) {
		int stride;
		int y = mb->y + 1;
		stride = mb->max_x * 16;
		x = 0;
		mb->y = y;
		if (mb->max_y <= y) {
			return -1;
		}
		mb->luma += stride * 15;
		mb->chroma += stride * 7;
		mb->top4x4pred = mb->top4x4pred_base;
		mb->top4x4coef = mb->top4x4coef_base;
		mb->top4x4inter = mb->mb_base + 1;
	}
	mb->x = x;
	mb->deblock_curr->idc = 0;
	if (0 <= mb->firstline) {
		mb->firstline--;
	}
	return 0;
}

static inline void frames_init(h264d_mb_current *mb, int num_frame, const m2d_frame_t *frame)
{
	h264d_frame_info_t *frm = mb->frame;
	frm->num = num_frame;
	std::copy(frame, frame + num_frame, frm->frames);
	memset(frm->lru, 0, sizeof(frm->lru));
}

#define NUM_ARRAY(x) (sizeof(x) / sizeof(x[0]))

int h264d_set_frames(h264d_context *h2d, int num_frame, m2d_frame_t *frame, uint8_t *second_frame, int second_frame_size)
{
	h264d_mb_current *mb;

	if (!h2d || (num_frame < 3) || (NUM_ARRAY(mb->frame->frames) < (size_t)num_frame) || !frame || !second_frame) {
		return -1;
	}
	mb = &h2d->mb_current;
	frames_init(mb, num_frame, frame);
	h2d->slice_header->reorder[0].ref_frames = mb->frame->refs[0];
	h2d->slice_header->reorder[1].ref_frames = mb->frame->refs[1];
	return init_mb_buffer(mb, second_frame, second_frame_size);
}

static int h2d_dispatch_one_nal(h264d_context *h2d, int code_type);

int h264d_decode_picture(h264d_context *h2d)
{
	dec_bits *stream;
	int code_type;
	int err;

	if (!h2d) {
		return -1;
	}
	stream = h2d->stream;
	if (setjmp(stream->jmp) != 0) {
		return -2;
	}
	h2d->slice_header->first_mb_in_slice = UINT_MAX;
	err = 0;
	code_type = 0;
	do {
		if (0 <= (err = m2d_find_mpeg_data(stream))) {
			code_type = get_bits(stream, 8);
			err = h2d_dispatch_one_nal(h2d, code_type);
		} else {
			err = -2;
			break;
		}
		VC_CHECK;
	} while (err == 0 || (code_type == SPS_NAL && 0 < err));
#ifdef DUMP_COEF
	print_coefs();
#endif
	return err;
}

static inline void dpb_init(h264d_dpb_t *dpb, int maxsize)
{
	dpb->size = 0;
	dpb->max = maxsize;
	dpb->output = -1;
}

//#define DUMP_DPB
static void dump_dpb(h264d_dpb_t *dpb)
{
#if defined(DUMP_DPB) && !defined(NDEBUG) && !defined(__RENESAS_VERSION__)
	printf("DPB length: %d\n", dpb->size);
	for (int i = 0; i < dpb->size; ++i) {
		printf("\t[%d] :\tpoc: %d idx: %d\n", i, dpb->data[i].poc, dpb->data[i].frame_idx);
	}
#endif
}

static inline void dpb_insert_non_idr(h264d_dpb_t *dpb, int poc, int frame_idx)
{
	int size = dpb->size;
	h264d_dpb_elem_t *d = dpb->data;
	h264d_dpb_elem_t *end = d + size;
	if (size < dpb->max) {
		dpb->size = size + 1;
		dpb->output = -1;
		while (d->poc < poc && d != end) {
			++d;
		}
		memmove(d + 1, d, (end - d) * sizeof(*d));
	} else {
		dpb->output = dpb->data[0].frame_idx;
		do {
			++d;
		} while (d->poc < poc && d != end);
		--d;
		memmove(dpb->data, dpb->data + 1, (d - dpb->data) * sizeof(*d));
	}
	d->poc = poc;
	d->frame_idx = frame_idx;
}

static inline void dpb_insert_idr(h264d_dpb_t *dpb, int poc, int frame_idx)
{
	int size = dpb->size;
	if (size < dpb->max) {
		dpb->size = size + 1;
	} else {
		size--;
		dpb->output = dpb->data[0].frame_idx;
		memmove(dpb->data, dpb->data + 1, size * sizeof(dpb->data[0]));
	}
	dpb->data[size].poc = poc;
	dpb->data[size].frame_idx = frame_idx;
}


static int dpb_force_pop(h264d_dpb_t *dpb)
{
	int size = dpb->size;
	int pop_idx = dpb->output;
	if (0 <= pop_idx) {
		dpb->output = -1;
		return pop_idx;
	} else if (size == 0) {
		return -1;
	}
	size -= 1;
	dpb->size = size;
	dpb->output = -1;
	pop_idx = dpb->data[0].frame_idx;
	memmove(dpb->data, dpb->data + 1, size * sizeof(dpb->data[0]));
	return pop_idx;
}

static inline int dpb_exist(const h264d_dpb_t *dpb, int frame_idx) {
	const h264d_dpb_elem_t *d = dpb->data;
	const h264d_dpb_elem_t *end = d + dpb->size;
	while (d != end) {
		if (d->frame_idx == frame_idx) {
			return 1;
		}
		++d;
	}
	return 0;
}

int h264d_get_decoded_frame(h264d_context *h2d, m2d_frame_t *frame, int bypass_dpb)
{
	h264d_frame_info_t *frm;
	int frame_idx;

	if (!h2d || !frame) {
		return -1;
	}
	frm = h2d->mb_current.frame;
	if (!bypass_dpb) {
		if (h2d->slice_header->marking.idr || h2d->slice_header->marking.mmco5) {
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
	return 1;
}

static int read_slice(h264d_context *h2d, dec_bits *st);

static int h2d_dispatch_one_nal(h264d_context *h2d, int code_type)
{
	int err;
	dec_bits *st = h2d->stream;

	switch(code_type & 31) {
	case SLICE_NONIDR_NAL:
	case SLICE_IDR_NAL:
		h2d->id = code_type;
		err = read_slice(h2d, st);
		break;
	case SEI_NAL:
		err = skip_sei(st);
		break;
	case SPS_NAL:
		err = read_seq_parameter_set(h2d->sps_i, st);
		if (0 <= err) {
			set_mb_size(&h2d->mb_current, h2d->sps_i[err].pic_width, h2d->sps_i[err].pic_height);
			h2d->header_callback(h2d->header_callback_arg, err);
		}
		break;
	case PPS_NAL:
		err = read_pic_parameter_set(h2d->pps_i, st);
		break;
	default:
		err = 0;
		break;
	}
	return err;
}

static int slice_header(h264d_context *h2d, dec_bits *st);
static int slice_data(h264d_context *h2d, dec_bits *st);

static int read_slice(h264d_context *h2d, dec_bits *st)
{
	int err;

	err = slice_header(h2d, st);
	if (err < 0) {
		return err;
	}
	return slice_data(h2d, st);
}

static int ref_pic_list_reordering(h264d_reorder_t *rdr, dec_bits *st, int num_ref_frames, int num_frames, int max_num_frames);
static int pred_weight_table(h264d_weight_table_t *tbl, dec_bits *st, int active_num);
static int dec_ref_pic_marking(int nal_id, h264d_marking_t *mrk, dec_bits *st);

static int slice_type_adjust(int slice_type)
{
	return (SI_SLICE < slice_type) ? slice_type - SI_SLICE - 1 : slice_type;
}

static inline void find_empty_frame(h264d_mb_current *mb)
{
	h264d_frame_info_t *frm = mb->frame;
	h264d_ref_frame_t *refs0 = frm->refs[0];
	h264d_ref_frame_t *refs1 = frm->refs[1];
	int8_t *lru = frm->lru;
	int max_idx = 0;
	int max_val = -1;
	int frm_num = frm->num;
	h264d_dpb_t *dpb;

	dpb = &frm->dpb;
	for (int i = 0; i < frm_num; ++i) {
		if (dpb_exist(dpb, i)) {
			lru[i] = 0;
		} else {
			lru[i] += 1;
		}
	}
	for (int i = 0; i < 16; ++i) {
		if (refs0[i].in_use) {
			lru[refs0[i].frame_idx] = 0;
		}
		if (refs1[i].in_use) {
			lru[refs1[i].frame_idx] = 0;
		}
	}
	for (int i = 0; i < frm_num; ++i) {
		int val = lru[i];
		if (max_val < val) {
			max_val = val;
			max_idx = i;
		}
	}
	lru[max_idx] = 0;
	frm->index = max_idx;
	frm->curr_luma = frm->frames[max_idx].luma;
	frm->curr_chroma = frm->frames[max_idx].chroma;
}

static void qp_matrix(int16_t *matrix, int scale, int shift)
{
	static const int8_t normAdjust[6][3] = {
		{10, 16, 13},
		{11, 18, 14},
		{13, 20, 16},
		{14, 23, 18},
		{16, 25, 20},
		{18, 29, 23}
	};
	/* before zigzag-scan, order of normAdjust v shall be:
	   0, 2, 2, 0, 1, 0, 2, 2, 2, 2, 1, 0, 1, 2, 2, 1
	 */
	int v0 = normAdjust[scale][0] << shift;
	int v1 = normAdjust[scale][1] << shift;
	int v2 = normAdjust[scale][2] << shift;
	/* write backword for SuperH. */
	matrix += 15;
	*matrix = v1;
	*--matrix = v2;
	*--matrix = v2;
	*--matrix = v1;
	*--matrix = v0;
	*--matrix = v1;
	*--matrix = v2;
	*--matrix = v2;
	*--matrix = v2;
	*--matrix = v2;
	*--matrix = v0;
	*--matrix = v1;
	*--matrix = v0;
	*--matrix = v2;
	*--matrix = v2;
	*--matrix = v0;
}

static void set_qp(h264d_mb_current *mb, int qpy)
{
	static const int8_t qpc_adjust[22] = {
		29, 30, 31, 32, 32, 33, 34, 34,
		35, 35, 36, 36, 37, 37, 37, 38,
		38, 38, 39, 39, 39, 39
	};
	int div;
	int mod;
	int qpc;

	if (qpy < 0) {
		qpy += 52;
	} else if (52 <= qpy) {
		qpy -= 52;
	}
	mb->qp = qpy;
	div = (unsigned)qpy / 6;
	mod = qpy - div * 6;
	qp_matrix(mb->qmaty, mod, div);
	qpc = qpy + mb->pps->chroma_qp_index[0];
	if (0 < qpc) {
		if (30 <= qpc) {
			if (51 < qpc) {
				qpc = 51;
			}
			qpc = qpc_adjust[qpc - 30];
		}
	} else {
		qpc = 0;
	}
	mb->qp_chroma = qpc;
	if (qpy == qpc) {
		mb->qmatc_p = mb->qmaty;
	} else {
		div = (unsigned)qpc / 6;
		mod = qpc - div * 6;
		mb->qmatc_p = mb->qmatc;
		qp_matrix(mb->qmatc, mod, div);
	}
}

static inline void calc_poc0(h264d_slice_header *hdr, int log2_max_lsb, int lsb)
{
	int prev_lsb, prev_msb;
	int msb;
	int max_lsb_2;
	if (hdr->first_mb_in_slice != 0) {
		return;
	}
	if (hdr->marking.idr || hdr->marking.mmco5) {
		prev_msb = 0;
		if (hdr->marking.mmco5 && hdr->field_pic_flag && hdr->bottom_field_flag) {
			prev_lsb = hdr->poc0.lsb;
		} else {
			prev_lsb = 0;
		}
	} else  {
		prev_lsb = hdr->poc0.lsb;
		prev_msb = hdr->poc0.msb;
	}
	hdr->poc0.lsb = lsb;
	max_lsb_2 = (1 << log2_max_lsb) >> 1;
	if ((lsb < prev_lsb) && (max_lsb_2 <= (prev_lsb - lsb))) {
		msb = prev_msb + max_lsb_2 * 2;
	} else if ((prev_lsb < lsb) && (max_lsb_2 < (lsb - prev_lsb))) {
		msb = prev_msb - max_lsb_2 * 2;
	} else {
		msb = prev_msb;
	}
	hdr->poc0.msb = msb;
	hdr->poc = msb + lsb;
	hdr->poc_bottom = hdr->poc + hdr->poc0.delta_pic_order_cnt_bottom;
}

static inline void calc_poc1(h264d_slice_header *hdr, const h264d_sps *sps, int nal_id)
{
	unsigned frame_num;
	int poc;

	if (hdr->first_mb_in_slice != 0) {
		return;
	}
	frame_num = hdr->frame_num;
	if (!hdr->marking.idr && !hdr->marking.mmco5) {
		if (frame_num < hdr->prev_frame_num) {
			hdr->poc1.num_offset += 1 << sps->log2_max_frame_num;
		}
	} else {
		hdr->poc1.num_offset = 0;
	}
	if (sps->num_ref_frames_in_pic_order_cnt_cycle) {
		frame_num += hdr->poc1.num_offset;
		if (frame_num != 0) {
			int cycle_cnt = 0;
			int cycle_sum = sps->offset_for_ref_frame[sps->num_ref_frames_in_pic_order_cnt_cycle - 1];
			frame_num--;
			if (frame_num != 0 && !(nal_id & 0x60)) {
				frame_num--;
			}
			while (cycle_sum <= (int)frame_num) {
				frame_num -= cycle_sum;
				cycle_cnt++;
			}
			poc = cycle_cnt * cycle_sum + sps->offset_for_ref_frame[frame_num & 255];
		} else {
			poc = sps->offset_for_ref_frame[0];
		}
		if ((nal_id & 0x60) == 0) {
			poc += sps->offset_for_non_ref_pic;
		}
	} else {
		poc = 0;
	}
	hdr->poc = poc = poc + hdr->poc1.delta_pic_order_cnt[0];
	hdr->poc_bottom = poc + sps->offset_for_top_to_bottom_field + hdr->poc1.delta_pic_order_cnt[1];
}

static inline void calc_poc2(h264d_slice_header *hdr, const h264d_sps *sps, int nal_id)
{
	int poc;

	if (hdr->first_mb_in_slice != 0) {
		return;
	}
	if (hdr->marking.idr || hdr->marking.mmco5) {
		hdr->poc2_prev_frameoffset = 0;
		poc = 0;
	} else {
		uint32_t frame_num = hdr->frame_num;
		if (frame_num < hdr->prev_frame_num) {
			hdr->poc2_prev_frameoffset += (1 << sps->log2_max_frame_num);
		}
		poc = (frame_num + hdr->poc2_prev_frameoffset) * 2 - ((nal_id & 0x60) == 0);
	}
	hdr->poc = poc;
	hdr->poc_bottom = poc;
}

static inline void ref_pic_init_p(h264d_slice_header *hdr, int max_frame_num, int num_ref_frames);
static inline void ref_pic_init_b(h264d_slice_header *hdr, int num_ref_frames);
static inline void set_dpb_max(h264d_dpb_t *dpb, const h264d_sps *sps)
{
	if (dpb->max < 0) {
		/* FIXME: dpb shall exists for each sps, so this method would be unnecessary. */
		int dpb_num = sps->max_dpb_in_mbs / ((uint32_t)(sps->pic_width * sps->pic_height) >> 8);
		dpb->max = 16 < dpb_num ? 16 : dpb_num;
	}
}

static inline int find_col_idx(const h264d_ref_frame_t *ref0, int len, int col_frameidx)
{
	int i;
	if (col_frameidx < 0) {
		return -1;
	}
	for (i = 0; i < len; ++i) {
		if (ref0[i].frame_idx == col_frameidx) {
			break;
		}
	}
	return (i < len) ? i : -1;
}

static inline int CLIP_P(int lower, int upper, int val)
{
	return (val < lower) ? lower : ((upper < val) ? upper : val);
}

static inline int col_scale(const h264d_ref_frame_t *ref0, int poc1, int curr_poc)
{
	int poc0 = ref0->poc;
	if (poc1 == poc0) {
		return 256;
	} else {
		int td = CLIP_P(-128, 127, poc1 - poc0);
		int tb = CLIP_P(-128, 127, curr_poc - poc0);
		int tx = (16384 + ABS(td / 2)) / td;
		return CLIP_P(-1024, 1023, (tb * tx + 32) >> 6);
	}
}

static void create_map_col_to_list0(int8_t *map_col_to_list0, int16_t *scale, const h264d_ref_frame_t *ref0, const h264d_ref_frame_t *ref1, int curr_poc, int len)
{
	int poc1 = ref1->poc;
	const int8_t *map = ref1[0].col->map_col_frameidx;
	for (int i = 0; i < len; ++i) {
		map_col_to_list0[i] = find_col_idx(ref0, len, map[i]);
		scale[i] = col_scale(&ref0[i], poc1, curr_poc);
	}
}

static void pred_direct8x8_temporal(h264d_mb_current *mb, int blk_idx, prev8x8_t *pblk);
static void b_skip_mb_temporal(h264d_mb_current *mb, int8_t *ref_idx, h264d_vector_set_t *mv);
static void pred_direct8x8_spatial(h264d_mb_current *mb, int blk_idx, prev8x8_t *pblk);
static void b_skip_mb_spatial(h264d_mb_current *mb, int8_t *ref_idx, h264d_vector_set_t *mv);

const h264d_bdirect_functions_t bdirect_functions[2] = {
	{
		pred_direct8x8_temporal,
		b_skip_mb_temporal
	},
	{
		pred_direct8x8_spatial,
		b_skip_mb_spatial
	}
};

static int slice_header(h264d_context *h2d, dec_bits *st)
{
	h264d_slice_header *hdr = h2d->slice_header;
	h264d_sps *sps;
	h264d_pps *pps;
	h264d_mb_current *mb;
	uint32_t prev_first_mb = hdr->first_mb_in_slice;
	int slice_type;

	mb = &h2d->mb_current;
	if ((hdr->first_mb_in_slice = ue_golomb(st)) <= prev_first_mb) {
		if (prev_first_mb != UINT_MAX) {
			return -2;
		}
		find_empty_frame(mb);
		memset(mb->deblock_base, 0, sizeof(*mb->deblock_base) * mb->max_x * mb->max_y);
	}
	READ_UE_RANGE(slice_type, st, 9);
	hdr->slice_type = slice_type_adjust(slice_type);
	if (3U <= (unsigned)hdr->slice_type) {
		return -1;
	}
	READ_UE_RANGE(hdr->pic_parameter_set_id, st, 255);
	pps = &h2d->pps_i[hdr->pic_parameter_set_id];
	sps = &h2d->sps_i[pps->seq_parameter_set_id];
	mb->pps = pps;
	mb->is_constrained_intra = pps->constrained_intra_pred_flag;

	if (hdr->first_mb_in_slice <= prev_first_mb) {
		m2d_frame_t *frm = &mb->frame->frames[mb->frame->index];
		frm->width = sps->pic_width;
		frm->height = sps->pic_height;
		memcpy(frm->crop, sps->frame_crop, sizeof(sps->frame_crop));
	}
	hdr->frame_num = get_bits(st, sps->log2_max_frame_num);
	if (!sps->frame_mbs_only_flag) {
		if ((hdr->field_pic_flag = get_onebit(st)) != 0) {
			hdr->bottom_field_flag = get_onebit(st);
		}
	}
	if ((h2d->id & 31) == 5) {
		hdr->marking.idr = 1;
		READ_UE_RANGE(hdr->idr_pic_id, st, 65535);
	} else {
		hdr->marking.idr = 0;
	}
	set_mb_size(mb, sps->pic_width, sps->pic_height);
	set_dpb_max(&mb->frame->dpb, sps);
	set_mb_pos(mb, hdr->first_mb_in_slice);
	if (sps->poc_type == 0) {
		uint32_t lsb = get_bits(st, sps->log2_max_poc_lsb);
		if (!hdr->field_pic_flag && pps->pic_order_present_flag) {
			hdr->poc0.delta_pic_order_cnt_bottom = se_golomb(st);
		} else {
			hdr->poc0.delta_pic_order_cnt_bottom = 0;
		}
		calc_poc0(hdr, sps->log2_max_poc_lsb, lsb);
	} else if (sps->poc_type == 1) {
		if (!sps->delta_pic_order_always_zero_flag) {
			hdr->poc1.delta_pic_order_cnt[0] = se_golomb(st);
			if (!hdr->field_pic_flag && pps->pic_order_present_flag) {
				hdr->poc1.delta_pic_order_cnt[1] = se_golomb(st);
			}
		} else {
			hdr->poc1.delta_pic_order_cnt[0] = 0;
			hdr->poc1.delta_pic_order_cnt[1] = 0;
		}
		calc_poc1(hdr, sps, h2d->id);
	} else {
		calc_poc2(hdr, sps, h2d->id);
	}
	if (pps->redundant_pic_cnt_present_flag) {
		hdr->redundant_pic_cnt = ue_golomb(st);
	}
	int max_frame_num = 1 << sps->log2_max_frame_num;
	switch (hdr->slice_type)
	{
	case B_SLICE:
		hdr->direct_spatial_mv_pred_flag = get_onebit(st);
		/* FALLTHROUGH */
	case P_SLICE:
		if ((hdr->num_ref_idx_active_override_flag = get_onebit(st)) != 0) {
			READ_UE_RANGE(hdr->num_ref_idx_lx_active_minus1[0], st, 31);
			if (hdr->slice_type == B_SLICE) {
				READ_UE_RANGE(hdr->num_ref_idx_lx_active_minus1[1], st, 31);
			}
		} else {
			hdr->num_ref_idx_lx_active_minus1[0] = pps->num_ref_idx_l0_active_minus1;
			hdr->num_ref_idx_lx_active_minus1[1] = pps->num_ref_idx_l1_active_minus1;
		}
		if (hdr->slice_type == P_SLICE) {
			ref_pic_init_p(hdr, max_frame_num, sps->num_ref_frames);
		} else {
			ref_pic_init_b(hdr, sps->num_ref_frames);
		}
		if (ref_pic_list_reordering(&hdr->reorder[0], st, sps->num_ref_frames, hdr->frame_num, max_frame_num)) {
			return -1;
		}
		if (hdr->slice_type == B_SLICE) {
			if (ref_pic_list_reordering(&hdr->reorder[1], st, sps->num_ref_frames, hdr->frame_num, max_frame_num)) {
				return -1;
			}
			mb->sub_mb_ref_map = sub_mb_ref_map_b;
			if (hdr->direct_spatial_mv_pred_flag == 0) {
				mb->bdirect->func = &bdirect_functions[0];
				create_map_col_to_list0(mb->bdirect->map_col_to_list0, mb->bdirect->scale, mb->frame->refs[0], mb->frame->refs[1], hdr->poc, sps->num_ref_frames);
			} else {
				mb->bdirect->func = &bdirect_functions[1];
			}
			if (pps->weighted_bipred_idc == 1) {
				pred_weight_table(&hdr->pred_weight_table[0], st, hdr->num_ref_idx_lx_active_minus1[0]);
				pred_weight_table(&hdr->pred_weight_table[1], st, hdr->num_ref_idx_lx_active_minus1[1]);
			}
		} else {
			mb->sub_mb_ref_map = sub_mb_ref_map_p;
			if (pps->weighted_pred_flag) {
				pred_weight_table(&hdr->pred_weight_table[0], st, hdr->num_ref_idx_lx_active_minus1[0]);
			}
		}
	}
	if (h2d->id & 0x60) {
		if (dec_ref_pic_marking(h2d->id & 31, &hdr->marking, st) < 0) {
			return -1;
		}
	} else {
		hdr->marking.mmco5 = 0;
	}
	if (pps->entropy_coding_mode_flag) {
		if ((hdr->slice_type != I_SLICE)
		    && (hdr->slice_type != SI_SLICE)) {
			READ_UE_RANGE(hdr->cabac_init_idc, st, 2);
		}
	}
	hdr->qp_delta = se_golomb(st);
	set_qp(&h2d->mb_current, pps->pic_init_qp + hdr->qp_delta);
	if ((hdr->slice_type == SP_SLICE)
	    || (hdr->slice_type == SI_SLICE)) {
		if ((hdr->slice_type == SP_SLICE)) {
			hdr->sp_for_switch_flag = get_onebit(st);
		}
		hdr->qs_delta = se_golomb(st);
	}
	deblock_info_t *firstmb = h2d->mb_current.deblock_base + hdr->first_mb_in_slice;
	if (pps->deblocking_filter_control_present_flag) {
		READ_UE_RANGE(hdr->disable_deblocking_filter_idc, st, 2);
		if (hdr->disable_deblocking_filter_idc != 1) {
			READ_SE_RANGE(hdr->slice_alpha_c0_offset_div2, st, -6, 6);
			READ_SE_RANGE(hdr->slice_beta_offset_div2, st, -6, 6);
			ENC_SLICEHDR(firstmb->slicehdr, hdr->slice_alpha_c0_offset_div2, hdr->slice_beta_offset_div2);
		} else {
			ENC_SLICEHDR(firstmb->slicehdr, 0, 0);
		}
	} else {
		hdr->disable_deblocking_filter_idc = 0;
		ENC_SLICEHDR(firstmb->slicehdr, 0, 0);
	}
	firstmb->idc = hdr->disable_deblocking_filter_idc + 1;
	return 0;
}

static inline int calc_short_term(int idc, int num, int frame_num, int max_frame_num)
{
	int no_wrap;
	assert(idc == 0 || idc == 1);
	if (idc == 0) {
		no_wrap = frame_num - num - 1;
		while (no_wrap < 0) {
			no_wrap += max_frame_num;
		}
	} else {
		no_wrap = frame_num + num + 1;
		while (max_frame_num <= no_wrap) {
			no_wrap -= max_frame_num;
		}
	}
	return no_wrap;
}

//#define DUMP_REF_LIST
static void dump_ref_list(h264d_ref_frame_t *refs, int num_ref_frames)
{
#if defined(DUMP_REF_LIST) && !defined(NDEBUG) && !defined(__RENESAS_VERSION__)
	printf("refs(%d) use\tnum\tpoc\n", num_ref_frames);
	for (int i = 0; i < 16; ++i) {
		if (!refs[i].in_use) {
			break;
		}
		printf("\t%d,\t%d,\t%d\n", refs[i].in_use, refs[i].num, refs[i].poc);
	}
#endif
}

static int ref_pic_list_reordering(h264d_reorder_t *rdr, dec_bits *st, int num_ref_frames, int frame_num, int max_frame_num)
{
	assert((unsigned)num_ref_frames <= 16);
	if ((rdr->ref_pic_list_reordering_flag = get_onebit(st)) != 0) {
		h264d_ref_frame_t *refs = rdr->ref_frames;
		int refIdxLx = -1;
		while (0 < num_ref_frames && !refs[num_ref_frames - 1].in_use) {
			--num_ref_frames;
		}
		while (++refIdxLx < 16) {
			h264d_ref_frame_t tmp_ref;
			int idc;
			uint32_t num;
			int curr_idx;
			int mode;

			READ_UE_RANGE(idc, st, 3);
			if (3 <= idc) {
				if (3 < idc) {
					return -1;
				}
				break;
			} else {
				num = ue_golomb(st);
			}
			if (idc < 2) {
				num = calc_short_term(idc, num, frame_num, max_frame_num);
				frame_num = num;
				mode = SHORT_TERM;
			} else {
				mode = LONG_TERM;
			}
			curr_idx = refIdxLx;
			while (!(refs[curr_idx].num == num && refs[curr_idx].in_use == mode) && curr_idx < num_ref_frames) {
				curr_idx++;
			}
			if (num_ref_frames <= curr_idx || refIdxLx == curr_idx) {
				continue;
			}
			tmp_ref = refs[curr_idx];
			memmove(&refs[refIdxLx + 1], &refs[refIdxLx], (curr_idx - refIdxLx) * sizeof(refs[0]));
			refs[refIdxLx] = tmp_ref;
		}
	}
	dump_ref_list(rdr->ref_frames, num_ref_frames);
	return 0;
}

static int pred_weight_table(h264d_weight_table_t *tbl, dec_bits *st, int active_num)
{
	unsigned log2_luma_denom;
	assert(0);
	READ_UE_RANGE(log2_luma_denom, st, 7);
	for (int i = 0; i < active_num; ++i) {
		int weight, offset;
		if (get_onebit(st)) {
		int weight, offset;
			READ_SE_RANGE(weight, st, -128, 127);
			READ_SE_RANGE(offset, st, -128, 127);
		} else {
			weight = 1;
			offset = 0;
		}
		tbl->luma_weight[i] = (weight << log2_luma_denom) + offset;
	}
	return 0;
}

static int dec_ref_pic_marking(int nal_unit_type, h264d_marking_t *mrk, dec_bits *st)
{
	uint32_t tmp = get_onebit(st);
	int op5_detect = 0;

	if (nal_unit_type == 5) {
		mrk->no_output_of_prior_pic_flag = tmp;
		mrk->long_term_reference_flag = get_onebit(st);
	} else {
		mrk->adaptive_ref_pic_marking_mode_flag = tmp;
		if (tmp) {
			h264d_mmco *mmco = mrk->mmco;
			int i = 16;
			do {
				READ_UE_RANGE(mmco->op, st, 6);
				if (mmco->op == 0) {
					break;
				} else if (mmco->op == 5) {
					op5_detect = 1;
				} else {
					tmp = ue_golomb(st);
					switch (mmco->op) {
					case 3:
						mmco->arg2 = ue_golomb(st);
						/* FALLTHROUGH */
					case 1:
					case 2:
					case 4:
					case 6:
						mmco->arg1 = tmp;
						break;
					}
				}
				mmco++;
			} while (--i);
		}
	}
	mrk->mmco5 = op5_detect;
	return 0;
}

static inline int get_nC(int na, int nb)
{
	int nc;
	if (na < 0) {
		if (nb < 0) {
			nc = 0;
		} else {
			nc = nb;
		}
	} else if (nb < 0) {
		nc = na;
	} else {
		nc = (na + nb + 1) >> 1;
	}
	return nc;
}

static inline int level_prefix(dec_bits *st)
{
	int val, len;

	const vlc_t *d = &level_prefix_bit8[show_bits(st, 8)];
	val = d->pattern;
	len = d->length;
	while (len < 0) {
		skip_bits(st, 8);
		d = &level_prefix_bit8[show_bits(st, 8)];
		val += d->pattern;
		len = d->length;
	}
	skip_bits(st, len);
	return val;
}

static int8_t total_zeros16(dec_bits *st, int total_coeff)
{
	const vlc_t *d;
	int zeros = 0;
	switch (total_coeff) {
	case 0:
		/* FALLTHROUGH */
	case 1:
		zeros = m2d_dec_vld_unary(st, total_zeros1_bit6, 6);
		break;
	case 2:
		d = &total_zeros2_bit6[show_bits(st, 6)];
		zeros = d->pattern;
		skip_bits(st, d->length);
		break;
	case 3:
		d = &total_zeros3_bit6[show_bits(st, 6)];
		zeros = d->pattern;
		skip_bits(st, d->length);
		break;
	case 4:
		d = &total_zeros4_bit5[show_bits(st, 5)];
		zeros = d->pattern;
		skip_bits(st, d->length);
		break;
	case 5:
		d = &total_zeros5_bit5[show_bits(st, 5)];
		zeros = d->pattern;
		skip_bits(st, d->length);
		break;
	case 6:
		d = &total_zeros6_bit6[show_bits(st, 6)];
		zeros = d->pattern;
		skip_bits(st, d->length);
		break;
	case 7:
		d = &total_zeros7_bit6[show_bits(st, 6)];
		zeros = d->pattern;
		skip_bits(st, d->length);
		break;
	case 8:
		d = &total_zeros8_bit6[show_bits(st, 6)];
		zeros = d->pattern;
		skip_bits(st, d->length);
		break;
	case 9:
		d = &total_zeros9_bit6[show_bits(st, 6)];
		zeros = d->pattern;
		skip_bits(st, d->length);
		break;
	case 10:
		d = &total_zeros10_bit5[show_bits(st, 5)];
		zeros = d->pattern;
		skip_bits(st, d->length);
		break;
	case 11:
		d = &total_zeros11_bit4[show_bits(st, 4)];
		zeros = d->pattern;
		skip_bits(st, d->length);
		break;
	case 12:
		d = &total_zeros12_bit4[show_bits(st, 4)];
		zeros = d->pattern;
		skip_bits(st, d->length);
		break;
	case 13:
		d = &total_zeros13_bit3[show_bits(st, 3)];
		zeros = d->pattern;
		skip_bits(st, d->length);
		break;
	case 14:
		d = &total_zeros14_bit2[show_bits(st, 2)];
		zeros = d->pattern;
		skip_bits(st, d->length);
		break;
	case 15:
		zeros = get_onebit_inline(st);
		break;
	}
	return zeros;
}

/**Read total_zeros for Chroma DC.
 */
static int8_t total_zeros4(dec_bits *st, int total_coeff)
{
	if (get_onebit(st)) {
		return 0;
	}
	if (total_coeff == 1) {
		if (get_onebit(st)) {
			return 1;
		} else {
			return 3 - get_onebit(st);
		}
	} else if (total_coeff == 2) {
		return 2 - get_onebit(st);
	} else {
		return 1;
	}
}

static int8_t run_before(dec_bits *st, int zeros_left)
{
	const vlc_t *d;
	int r;

	switch (zeros_left) {
	case 1:
		r = get_onebit(st) ^ 1;
		break;
	case 2:
		d = &run_before_2_bit2[show_bits(st, 2)];
		r = d->pattern;
		skip_bits(st, d->length);
		break;
	case 3:
		r = 3 - get_bits(st, 2);
		break;
	case 4:
		d = &run_before_4_bit3[show_bits(st, 3)];
		r = d->pattern;
		skip_bits(st, d->length);
		break;
	case 5:
		d = &run_before_5_bit3[show_bits(st, 3)];
		r = d->pattern;
		skip_bits(st, d->length);
		break;
	case 6:
		d = &run_before_6_bit3[show_bits(st, 3)];
		r = d->pattern;
		skip_bits(st, d->length);
		break;
	default:
		r = m2d_dec_vld_unary(st, run_before_7_bit3, 3);
		break;
	}
	return r;
}

static inline void coeff_writeback(int *coeff, int total_coeff, const int8_t *run, const int *level, const int16_t *qmat, uint32_t dc_mask)
{
	int idx;
	idx = (dc_mask >> 4) - 1;
	for (int i = total_coeff - 1; 0 <= i; --i) {
		idx = idx + 1 + run[i];
		idx &= 15; /* for bit error */
		coeff[idx] = level[i] * qmat[idx & dc_mask];
	}
	
}

struct residual_block_cavlc {
	int operator()(h264d_mb_current *mb, int na, int nb, dec_bits *st, int *coeff, int num_coeff, const int16_t *qmat, int avail, int pos4x4, int cat, uint32_t dc_mask) const {
		int level[16];
		int8_t run[16];
		int val;
		int trailing_ones;
		int total_coeff;
		int zeros_left;

		if (num_coeff <= 4) {
			val = m2d_dec_vld_unary(st, total_ones_nc_chroma_bit6, 6);
		} else {
			int nc = get_nC(na, nb);
			if (nc < 2) {
				val = m2d_dec_vld_unary(st, total_ones_nc02_bit6, 6);
			} else if (nc < 4) {
				val = m2d_dec_vld_unary(st, total_ones_nc24_bit6, 6);
			} else if (nc < 8) {
				val = m2d_dec_vld_unary(st, total_ones_nc48_bit6, 6);
			} else {
				val = m2d_dec_vld_unary(st, total_ones_nc8_bit6, 6);
			}
		}
		trailing_ones = val >> 5;
		total_coeff = val & 31;
		if (total_coeff == 0) {
			return 0;
		}
		memset(coeff + (dc_mask >> 4), 0, sizeof(*coeff) * num_coeff);
		int suffix_len = ((10 < total_coeff) && (trailing_ones < 3));
		for (int i = 0; i < total_coeff; ++i) {
			if (i < trailing_ones) {
				level[i] = 1 - get_onebit_inline(st) * 2;
			} else {
				int lvl_prefix = level_prefix(st);
				int lvl = lvl_prefix << suffix_len;
				if (0 < suffix_len || 14 <= lvl_prefix) {
					int size = suffix_len;
					if (lvl_prefix == 14 && !size) {
						size = 4;
					} else if (lvl_prefix == 15) {
						size = 12;
					}
					if (size) {
						lvl += get_bits(st, size);
					}
				}
				if (suffix_len == 0 && lvl_prefix == 15) {
					lvl += 15;
				}
				if (i == trailing_ones && trailing_ones < 3) {
					lvl += 2;
				}
				level[i] = lvl = ((lvl & 1) ? (-lvl - 1) : (lvl + 2)) >> 1;
				suffix_len = suffix_len ? suffix_len : 1;
				if (suffix_len < 6 && ((3 << (suffix_len - 1)) < ABS(lvl))) {
					suffix_len++;
				}
			}
		}
		if (total_coeff < num_coeff) {
			if (4 < num_coeff) {
				zeros_left = total_zeros16(st, total_coeff);
			} else {
				zeros_left = total_zeros4(st, total_coeff);
			}
		} else {
			zeros_left = 0;
		}
		for (int i = 0; i < total_coeff - 1; ++i) {
			int r;
			run[i] = r = (0 < zeros_left) ? run_before(st, zeros_left) : 0;
			zeros_left -= r;
		}
		run[total_coeff - 1] = zeros_left;
		coeff_writeback(coeff, total_coeff, run, level, qmat, dc_mask);
		return total_coeff <= 15 ? total_coeff : 15;
	}
};

static void intra_chroma_dc_transform(const int *src, int *dst);
static void ac4x4transform_dconly_chroma(uint8_t *dst, int dc, int stride);
static void ac4x4transform_acdc_chroma(uint8_t *dst, const int *coeff, int stride);

template <typename F0>
static inline int residual_chroma(h264d_mb_current *mb, uint32_t cbp, dec_bits *st, int avail, F0 ResidualBlock)
{
	int coeff[16];
	int dc[2][4];
	const int16_t *qmat;
	uint8_t *chroma;
	int stride;
	int *dcp;

	cbp >>= 4;
	if (!cbp) {
		mb->left4x4coef &= 0x0000ffff;
		*mb->top4x4coef &= 0x0000ffff;
		return 0;
	}
	qmat = mb->qmatc_p;
	for (int i = 0; i < 2; ++i) {
		if (ResidualBlock(mb, 0, 0, st, coeff, 4, qmat, avail, 16 + i, 3, 0)) {
			intra_chroma_dc_transform(coeff, dc[i]);
		} else {
			memset(dc[i], 0, sizeof(dc[0][0]) * 4);
		}
	}
	chroma = mb->chroma;
	stride = mb->max_x * 16;
	dcp = dc[0];
	if (cbp & 2) {
		int c0, c1, c2, c3;
		uint32_t left = mb->left4x4coef >> 16;
		uint32_t top = *mb->top4x4coef >> 16;
		for (int i = 0; i < 2; ++i) {
			int c0left, c2left;
			int c0top, c1top;
			if (avail & 1) {
				c0left = UNPACK(left, 0);
				c2left = UNPACK(left, 1);
			} else {
				c0left = c2left = -1;
			}
			if (avail & 2) {
				c0top = UNPACK(top, 0);
				c1top = UNPACK(top, 1);
			} else {
				c0top = c1top = -1;
			}
			if ((c0 = ResidualBlock(mb, c0left, c0top, st, coeff, 15, qmat, avail, 18 + i * 4, 4, 0x1f)) != 0) {
				coeff[0] = *dcp++;
				ac4x4transform_acdc_chroma(chroma, coeff, stride);
			} else {
				ac4x4transform_dconly_chroma(chroma, *dcp++, stride);
			}
			if ((c1 = ResidualBlock(mb, c0, c1top, st, coeff, 15, qmat, avail, 19 + i * 4, 4, 0x1f)) != 0) {
				coeff[0] = *dcp++;
				ac4x4transform_acdc_chroma(chroma + 8, coeff, stride);
			} else {
				ac4x4transform_dconly_chroma(chroma + 8, *dcp++, stride);
			}
			if ((c2 = ResidualBlock(mb, c2left, c0, st, coeff, 15, qmat, avail, 20 + i * 4, 4, 0x1f)) != 0) {
				coeff[0] = *dcp++;
				ac4x4transform_acdc_chroma(chroma + stride * 4, coeff, stride);
			} else {
				ac4x4transform_dconly_chroma(chroma + stride * 4, *dcp++, stride);
			}
			if ((c3 = ResidualBlock(mb, c2, c1, st, coeff, 15, qmat, avail, 21 + i * 4, 4, 0x1f)) != 0) {
				coeff[0] = *dcp++;
				ac4x4transform_acdc_chroma(chroma + stride * 4 + 8, coeff, stride);
			} else {
				ac4x4transform_dconly_chroma(chroma + stride * 4 + 8, *dcp++, stride);
			}
			left = ((left >> 8) & 0xff) | (c3 << 12) | (c1 << 8);
			top = ((top >> 8) & 0xff)| (c3 << 12) | (c2 << 8);
			chroma++;
		}
		mb->left4x4coef = (mb->left4x4coef & 0x0000ffff) | (left << 16);
		*mb->top4x4coef = (*mb->top4x4coef & 0x0000ffff) | (top << 16);
	} else {
		ac4x4transform_dconly_chroma(chroma, *dcp++, stride);
		ac4x4transform_dconly_chroma(chroma + 8, *dcp++, stride);
		ac4x4transform_dconly_chroma(chroma + stride * 4, *dcp++, stride);
		ac4x4transform_dconly_chroma(chroma + stride * 4 + 8, *dcp++, stride);
		ac4x4transform_dconly_chroma(chroma + 1, *dcp++, stride);
		ac4x4transform_dconly_chroma(chroma + 9, *dcp++, stride);
		ac4x4transform_dconly_chroma(chroma + stride * 4 + 1, *dcp++, stride);
		ac4x4transform_dconly_chroma(chroma + stride * 4 + 9, *dcp++, stride);
		mb->left4x4coef &= 0x0000ffff;
		*mb->top4x4coef &= 0x0000ffff;
	}
	VC_CHECK;
	return 0;
}

/* intra4x4 prediction */
template <int N>
static uint32_t sum_top(uint8_t *src, int stride)
{
	uint32_t dc = 0;
	int i = N / 4;
	src -= stride;
	do {
		dc += *src++;
		dc += *src++;
		dc += *src++;
		dc += *src++;
	} while (--i);
	return dc;
}

template <int N>
static uint32_t sum_left(uint8_t *src, int stride)
{
	uint32_t dc = 0;
	int i = N / 4;
	src--;
	do {
		dc += *src;
		src += stride;
		dc += *src;
		src += stride;
		dc += *src;
		src += stride;
		dc += *src;
		src += stride;
	} while (--i);
	return dc;
}

struct intra4x4pred_mode_cavlc {
	int operator()(int a, int b, dec_bits *st, h264d_cabac_t *cb) const {
		int pred = MIN(a, b);
		if (!get_onebit_inline(st)) {
			int rem = get_bits(st, 3);
			pred = (rem < pred) ? rem : rem + 1;
		}
		return pred;
	}
};

static int intra4x4pred_dc(uint8_t *dst, int stride, int avail)
{
	uint32_t dc;
	
	if (avail & 1) {
		if (avail & 2) {
			dc = (sum_left<4>(dst, stride) + sum_top<4>(dst, stride) + 4) >> 3;
		} else {
			dc = (sum_left<4>(dst, stride) + 2) >> 2;
		}
	} else if (avail & 2) {
		dc = (sum_top<4>(dst, stride) + 2) >> 2;
	} else {
		dc = 0x80;
	}
	dc = dc * 0x01010101U;
	*(uint32_t *)dst = dc;
	dst += stride;
	*(uint32_t *)dst = dc;
	dst += stride;
	*(uint32_t *)dst = dc;
	dst += stride;
	*(uint32_t *)dst = dc;
	return 0;
}

static int intra4x4pred_horiz(uint8_t *dst, int stride, int avail)
{
	if (!(avail & 1)) {
		return -1;
	}
	*(uint32_t *)dst = dst[-1] * 0x01010101U;
	dst = dst + stride;
	*(uint32_t *)dst = dst[-1] * 0x01010101U;
	dst = dst + stride;
	*(uint32_t *)dst = dst[-1] * 0x01010101U;
	dst = dst + stride;
	*(uint32_t *)dst = dst[-1] * 0x01010101U;
	return 0;
}

static int intra4x4pred_vert(uint8_t *dst, int stride, int avail)
{
	uint32_t *src;
	uint32_t t0;

	if (!(avail & 2)) {
		return -1;
	}
	src = (uint32_t *)(dst - stride);

	t0 = *src;
	*(uint32_t *)dst = t0;
	dst += stride;
	*(uint32_t *)dst = t0;
	dst += stride;
	*(uint32_t *)dst = t0;
	dst += stride;
	*(uint32_t *)dst = t0;
	return 0;
}

#define FIR3(a, b, c) (((a) + (b) * 2 + (c) + 2) >> 2)
#define FIR2(a, b) (((a) + (b) + 1) >> 1)


/**Intra 4x4 prediction Diagonal Down Left.
 */
static int intra4x4pred_ddl(uint8_t *dst, int stride, int avail)
{
	uint8_t *src;
	uint32_t t0, t1, t2, d0;
	src = dst - stride;
	t0 = *src++;
	t1 = *src++;
	t2 = *src++;
	d0 = FIR3(t0, t1, t2);
	t0 = *src++;
	if (avail & 4) {
		d0 = (d0 << 8) | FIR3(t1, t2, t0);
		t1 = *src++;
		d0 = (d0 << 8) | FIR3(t2, t0, t1);
		t2 = *src++;
		d0 = (d0 << 8) | FIR3(t0, t1, t2);
#ifndef WORDS_BIGENDIAN
		d0 = bswap32(d0);
#endif
		t0 = *src++;
		*(uint32_t *)dst = d0;
		dst += stride;
#ifdef WORDS_BIGENDIAN
		d0 = (d0 << 8) | (FIR3(t1, t2, t0);
#else
		d0 = (FIR3(t1, t2, t0) << 24) | (d0 >> 8);
#endif
		t1 = *src++;
		*(uint32_t *)dst = d0;
		dst += stride;
#ifdef WORDS_BIGENDIAN
		d0 = (d0 << 8) | FIR3(t1, t2, t0);
#else
		d0 = (FIR3(t2, t0, t1) << 24) | (d0 >> 8);
#endif
		*(uint32_t *)dst = d0;
		dst += stride;
#ifdef WORDS_BIGENDIAN
		d0 = (d0 << 8) | FIR3(t1, t2, t0);
#else
		d0 = (FIR3(t0, t1, t1) << 24) | (d0 >> 8);
#endif
	} else {
		d0 = (d0 << 8) | FIR3(t1, t2, t0);
		d0 = (d0 << 8) | FIR3(t2, t0, t0);
		d0 = (d0 << 8) | t0;
#ifndef WORDS_BIGENDIAN
		d0 = bswap32(d0);
		t0 = t0 << 24;
#endif
		*(uint32_t *)dst = d0;
		dst += stride;
#ifdef WORDS_BIGENDIAN
		d0 = (d0 << 8) | t0;
#else
		d0 = t0 | (d0 >> 8);
#endif
		*(uint32_t *)dst = d0;
		dst += stride;
#ifdef WORDS_BIGENDIAN
		d0 = (d0 << 8) | t0;
#else
		d0 = t0 | (d0 >> 8);
#endif
		*(uint32_t *)dst = d0;
		dst += stride;
#ifdef WORDS_BIGENDIAN
		d0 = (d0 << 8) | t0;
#else
		d0 = t0 | (d0 >> 8);
#endif
	}
	*(uint32_t *)dst = d0;
	return 0;
}

/** Intra 4x4 prediction Diagonal Down Right.
 */
static int intra4x4pred_ddr(uint8_t *dst, int stride, int avail)
{
	uint8_t *src;
	uint32_t t0, t1, t2, t3, d0;
	if ((avail & 3) != 3) {
		return -1;
	}
	src = dst - stride - 1;
	t0 = *src++;
	t1 = *src++;
	t2 = *src++;
	d0 = FIR3(t0, t1, t2);
	t3 = *src++;
	d0 = (d0 << 8) | FIR3(t1, t2, t3);
	d0 = (d0 << 8) | FIR3(t2, t3, *src);
	src = dst - 1;
	t3 =  *src;
	d0 = (FIR3(t3, t0, t1) << 24) | d0;
#ifndef WORDS_BIGENDIAN
	d0 = bswap32(d0);
#endif
	src += stride;
	*(uint32_t *)dst = d0;

	t2 = *src;
	dst += stride;
#ifdef WORDS_BIGENDIAN
	d0 = (d0 >> 8) | (FIR3(t2, t3, t0) << 24);
#else
	d0 = (d0 << 8) | FIR3(t2, t3, t0);
#endif
	src += stride;
	*(uint32_t *)dst = d0;

	t1 = *src;
	dst += stride;
#ifdef WORDS_BIGENDIAN
	d0 = (d0 >> 8) | (FIR3(t1, t2, t3) << 24);
#else
	d0 = (d0 << 8) | FIR3(t1, t2, t3);
#endif
	src += stride;
	*(uint32_t *)dst = d0;

	t0 = *src;
	dst += stride;
#ifdef WORDS_BIGENDIAN
	d0 = (d0 >> 8) | (FIR3(t0, t1, t2) << 24);
#else
	d0 = (d0 << 8) | FIR3(t0, t1, t2);
#endif
	*(uint32_t *)dst = d0;
	return 0;
}

/** Intra 4x4 prediction Vertical Right.
 */
static int intra4x4pred_vr(uint8_t *dst, int stride, int avail)
{
	uint8_t *src;
	uint32_t t0, t1, t2, t3, t4, d0, d1;
	if ((avail & 3) != 3) {
		return -1;
	}
	src = dst - stride - 1;
	t0 = *src++;
	t1 = *src++;
	d0 = FIR2(t0, t1);
	t2 = *src++;
	d0 = (FIR2(t1, t2) << 8) | d0;
	d1 = FIR3(t0, t1, t2) << 8;
	t3 = *src++;
	d0 = (FIR2(t2, t3) << 16) | d0;
	d1 = (FIR3(t1, t2, t3) << 16) | d1;

	t4 = *src++;
	d0 = (FIR2(t3, t4) << 24) | d0;
	d1 = (FIR3(t2, t3, t4) << 24) | d1;
#ifdef WORDS_BIGENDIAN
	d0 = bswap32(d0);
#endif
	src = dst - 1;
	*(uint32_t *)dst = d0;

	t4 = *src;
	dst += stride;
	d1 = d1 | FIR3(t4, t0, t1);
#ifdef WORDS_BIGENDIAN
	d0 = bswap32(d0);
#endif
	src += stride;
	*(uint32_t *)dst = d1;

	t3 = *src;
	dst += stride;
#ifdef WORDS_BIGENDIAN
	d0 = (FIR3(t3, t4, t0) << 24)| (d0 >> 8);
#else
	d0 = (d0 << 8) | FIR3(t3, t4, t0);
#endif
	src += stride;
	*(uint32_t *)dst = d0;

	t2 = *src;
	dst += stride;
#ifdef WORDS_BIGENDIAN
	d1 = (FIR3(t2, t3, t4) << 24)| (d1 >> 8);
#else
	d1 = (d1 << 8) | FIR3(t2, t3, t4);
#endif
	*(uint32_t *)dst = d1;
	return 0;
}

/** Intra 4x4 prediction Horizontal Down.
 */
static int intra4x4pred_hd(uint8_t *dst, int stride, int avail)
{
	uint8_t *src;
	uint32_t t0, t1, t2, d0;
	if ((avail & 3) != 3) {
		return -1;
	}
	src = dst - stride - 1;
	t0 = *src++;
	t1 = *src++;
	t2 = *src++;
	d0 = FIR3(t1, t2, *src);
	src = dst - 1;
	d0 = (d0 << 8) | FIR3(t0, t1, t2);
	t2 = *src;
	d0 = (d0 << 8) | FIR3(t1, t0, t2);
	src += stride;
	d0 = (d0 << 8) | FIR2(t0, t2);
#ifdef WORDS_BIGENDIAN
	d0 = bswap32(d0);
#endif
	*(uint32_t *)dst = d0;

	t1 = *src;
	dst += stride;
	src += stride;
#ifdef WORDS_BIGENDIAN
	d0 = (FIR3(t0, t2, t1) << 16) | (d0 >> 16);
	d0 = (FIR2(t2, t1) << 24) | d0;
#else
	d0 = (d0 << 8) | FIR3(t0, t2, t1);
	d0 = (d0 << 8) | FIR2(t2, t1);
#endif
	*(uint32_t *)dst = d0;

	t0 = *src;
	dst += stride;
	src += stride;
#ifdef WORDS_BIGENDIAN
	d0 = (FIR3(t2, t1, t0) << 16) | (d0 >> 16);
	d0 = (FIR2(t1, t0) << 24) | d0;
#else
	d0 = (d0 << 8) | FIR3(t2, t1, t0);
	d0 = (d0 << 8) | FIR2(t1, t0);
#endif
	*(uint32_t *)dst = d0;

	t2 = *src;
	dst += stride;
#ifdef WORDS_BIGENDIAN
	d0 = (FIR3(t1, t0, t2) << 16) | (d0 >> 16);
	d0 = (FIR2(t0, t2) << 24) | d0;
#else
	d0 = (d0 << 8) | FIR3(t1, t0, t2);
	d0 = (d0 << 8) | FIR2(t0, t2);
#endif
	*(uint32_t *)dst = d0;
	return 0;
}

/** Intra 4x4 prediction Vertical Left.
 */
static int intra4x4pred_vl(uint8_t *dst, int stride, int avail)
{
	uint8_t *src;
	uint32_t t0, t1, t2, d0, d1;

	src = dst - stride;
	t0 = *src++;
	t1 = *src++;
	t2 = *src++;
	d0 = FIR2(t0, t1);
#ifdef WORDS_BIGENDIAN
	d0 = (d0 << 8) | FIR2(t1, t2);
#else
	d0 = (FIR2(t1, t2) << 8) | d0;
#endif
	d1 = FIR3(t0, t1, t2);
	t0 = *src++;
#ifdef WORDS_BIGENDIAN
	d0 = (d0 << 8) | FIR2(t2, t0);
	d1 = (d1 << 8) | FIR3(t1, t2, t0);
#else
	d0 = (FIR2(t2, t0) << 16) | d0;
	d1 = (FIR3(t1, t2, t0) << 8) | d1;
#endif
	if (avail & 4) {
		t1 = *src++;
#ifdef WORDS_BIGENDIAN
		d0 = (d0 << 8) | FIR2(t1, t0);
#else
		d0 = (FIR2(t1, t0) << 24) | d0;
#endif
		*(uint32_t *)dst = d0;

#ifdef WORDS_BIGENDIAN
		d1 = (d1 << 8) | FIR3(t1, t0, t2) << 16);
#else
		d1 = (FIR3(t1, t0, t2) << 16) | d1;
#endif
		dst += stride;
		t2 = *src++;
#ifdef WORDS_BIGENDIAN
		d1 = (d1 << 8) | FIR3(t2, t1, t0);
#else
		d1 = (FIR3(t2, t1, t0) << 24) | d1;
#endif
		*(uint32_t *)dst = d1;

		dst += stride;
#ifdef WORDS_BIGENDIAN
		d0 = (d0 << 8) | FIR2(t2, t1);
#else
		d0 = (FIR2(t2, t1) << 24) | (d0 >> 8);
#endif
		*(uint32_t *)dst = d0;

		t0 = *src;
		dst += stride;
#ifdef WORDS_BIGENDIAN
		d1 = (d1 << 8) | FIR3(t1, t2, t0);
#else
		d1 = (FIR3(t1, t2, t0) << 24) | (d1 >> 8);
#endif
	} else {
		t1 = FIR3(t2, t0, t0);
#ifdef WORDS_BIGENDIAN
		d0 = (d0 << 8) | t0;
#else
		t0 <<= 24;
		d0 = t0 | d0;
#endif
		*(uint32_t *)dst = d0;
		dst += stride;
#ifdef WORDS_BIGENDIAN
		d1 = (d1 << 16) | (t1 << 8) | t0;
#else
		d1 = (t1 << 16) | d1 | t0;
#endif
		*(uint32_t *)dst = d1;
		dst += stride;
#ifdef WORDS_BIGENDIAN
		d0 = (d0 << 8) | t0;
#else
		d0 = (d0 >> 8) | t0;
#endif
		*(uint32_t *)dst = d0;
		dst += stride;
#ifdef WORDS_BIGENDIAN
		d1 = (d1 << 8) | t0;
#else
		d1 = (d1 >> 8) | t0;
#endif
	}
	*(uint32_t *)dst = d1;
	return 0;
}

/** Intra 4x4 prediction Horizontal Up.
 */
static int intra4x4pred_hu(uint8_t *dst, int stride, int avail)
{
	uint8_t *src;
	uint32_t t0, t1, t2, d0;
	if (!(avail & 1)) {
		return -1;
	}
	src = dst - 1;
	t0 = *src;
	src += stride;
	t1 = *src;
	src += stride;
	d0 = FIR2(t0, t1);
	t2 = *src;
	src += stride;
	d0 = (FIR3(t0, t1, t2) << 8) | d0;
	d0 = (FIR2(t1, t2) << 16) | d0;
	t0 = *src;
	d0 = (FIR3(t1, t2, t0) << 24) | d0;
#ifdef WORDS_BIGENDIAN
	d0 = bswap32(d0);
#endif
	*(uint32_t *)dst = d0;

	dst += stride;
#ifdef WORDS_BIGENDIAN
	d0 = (d0 << 8) | FIR2(t2, t0);
	d0 = (d0 << 8) | FIR3(t2, t0, t0);
#else
	d0 = (FIR2(t2, t0) << 16) | (d0 >> 16);
	d0 = (FIR3(t2, t0, t0) << 24) | d0;
#endif
	*(uint32_t *)dst = d0;

	dst += stride;
#ifdef WORDS_BIGENDIAN
	t1 = (t0 << 8) | d0;
	d0 = (d0 << 16) | t1;
#else
	t1 = (t0 << 24) | (t0 << 16);
	d0 = t1 | (d0 >> 16);
#endif
	*(uint32_t *)dst = d0;

	dst += stride;
#ifdef WORDS_BIGENDIAN
	d0 = (d0 << 16) | t1;
#else
	d0 = t1 | (d0 >> 16 );
#endif
	*(uint32_t *)dst = d0;
	return 0;
}

static int (* const intra4x4pred_func[9])(uint8_t *dst, int stride, int avail) = {
	intra4x4pred_vert,
	intra4x4pred_horiz,
	intra4x4pred_dc,
	intra4x4pred_ddl,
	intra4x4pred_ddr,
	intra4x4pred_vr,
	intra4x4pred_hd,
	intra4x4pred_vl,
	intra4x4pred_hu
};

template <typename F>
static int mb_pred_intra4x4(h264d_mb_current *mb, dec_bits *st, int avail, int8_t *pred4x4, F I4x4PredMode) {
	uint32_t left = mb->left4x4pred;
	uint32_t top = *mb->top4x4pred;
	h264d_cabac_t *cb = mb->cabac;
	pred4x4[0] = I4x4PredMode(avail & 2 ? UNPACK(left, 0) : 2, avail & 1 ? UNPACK(top, 0) : 2, st, cb);
	pred4x4[1] = I4x4PredMode(avail & 2 ? pred4x4[0] : 2, UNPACK(top, 1), st, cb);
	pred4x4[2] = I4x4PredMode(UNPACK(left, 1), avail & 1 ? pred4x4[0] : 2, st, cb);
	pred4x4[3] = I4x4PredMode(pred4x4[2], pred4x4[1], st, cb);
	pred4x4[4] = I4x4PredMode(avail & 2 ? pred4x4[1] : 2, UNPACK(top, 2), st, cb);
	pred4x4[5] = I4x4PredMode(avail & 2 ? pred4x4[4] : 2, UNPACK(top, 3), st, cb);
	pred4x4[6] = I4x4PredMode(pred4x4[3], pred4x4[4], st, cb);
	pred4x4[7] = I4x4PredMode(pred4x4[6], pred4x4[5], st, cb);

	pred4x4[8] = I4x4PredMode(UNPACK(left, 2), avail & 1 ? pred4x4[2] : 2, st, cb);
	pred4x4[9] = I4x4PredMode(pred4x4[8], pred4x4[3], st, cb);
	pred4x4[10] = I4x4PredMode(UNPACK(left, 3), avail & 1 ? pred4x4[8] : 2, st, cb);
	pred4x4[11] = I4x4PredMode(pred4x4[10], pred4x4[9], st, cb);
	pred4x4[12] = I4x4PredMode(pred4x4[9], pred4x4[6], st, cb);
	pred4x4[13] = I4x4PredMode(pred4x4[12], pred4x4[7], st, cb);
	pred4x4[14] = I4x4PredMode(pred4x4[11], pred4x4[12], st, cb);
	pred4x4[15] = I4x4PredMode(pred4x4[14], pred4x4[13], st, cb);

	mb->left4x4pred = (pred4x4[15] << 12) | (pred4x4[13] << 8) | (pred4x4[7] << 4)| pred4x4[5];
	*mb->top4x4pred = (pred4x4[15] << 12) | (pred4x4[14] << 8) | (pred4x4[11] << 4)| pred4x4[10];
	return 0;
}

static inline void fill_dc_if_unavailable(h264d_mb_current *mb, int avail)
{
	if (!(avail & 1)) {
		mb->left4x4pred = 0x22222222;
	}
	if (!(avail & 2)) {
		*mb->top4x4pred = 0x22222222;
	}
}

static int mb_intra_chroma_pred_dc(uint8_t *dst, int stride, int avail);
static int mb_intra_chroma_pred_horiz(uint8_t *dst, int stride, int avail);
static int mb_intra_chroma_pred_planer(uint8_t *dst, int stride, int avail);

template <int N>
static int mb_intra16xpred_vert(uint8_t *dst, int stride, int avail)
{
	uint32_t *src;
	uint32_t t0, t1, t2, t3;
	int i;

	if (!(avail & 2)) {
		return -1;
	}
	src = (uint32_t *)(dst - stride);
	t0 = *src++;
	t1 = *src++;
	t2 = *src++;
	t3 = *src++;
	i = N;
	do {
		*((uint32_t *)dst) = t0;
		*((uint32_t *)dst + 1) = t1;
		*((uint32_t *)dst + 2) = t2;
		*((uint32_t *)dst + 3) = t3;
		dst += stride;
	} while (--i);
	return 0;
}


static int (* const intra_chroma_pred[4])(uint8_t *dst, int stride, int avail) = {
	mb_intra_chroma_pred_dc,
	mb_intra_chroma_pred_horiz,
	mb_intra16xpred_vert<8>,
	mb_intra_chroma_pred_planer
};

static inline void ac4x4transform_maybe(uint8_t *dst, const int *coeff, int stride, int num_coeff);
static void mb_intra_save_info(h264d_mb_current *mb)
{
	mb->lefttop_ref[0] = mb->top4x4inter->ref[1][0];
	mb->lefttop_ref[1] = mb->top4x4inter->ref[1][1];
	mb->lefttop_mv[0].vector = mb->top4x4inter->mov[3].mv[0].vector;
	mb->lefttop_mv[1].vector = mb->top4x4inter->mov[3].mv[1].vector;
	mb->left4x4inter->direct8x8 = 0;
	mb->top4x4inter->direct8x8 = 0;
	memset(mb->left4x4inter->mov, 0, sizeof(mb->left4x4inter->mov));
	memset(mb->left4x4inter->mvd, 0, sizeof(mb->left4x4inter->mvd));
	memset(mb->top4x4inter->mov, 0, sizeof(mb->top4x4inter->mov));
	memset(mb->top4x4inter->mvd, 0, sizeof(mb->top4x4inter->mvd));
	memset(mb->left4x4inter->ref, -1, sizeof(mb->left4x4inter->ref));
	memset(mb->top4x4inter->ref, -1, sizeof(mb->top4x4inter->ref));
	mb->col_curr->type = COL_MB16x16;
	memset(mb->col_curr->ref, -1, sizeof(mb->col_curr->ref));
}

static inline void store_strength_intra(h264d_mb_current *mb) {
	deblock_info_t *deb = mb->deblock_curr;
	deb->qpy = mb->qp;
	deb->qpc = mb->qp_chroma;
	deb->str4_horiz = 1;
	deb->str4_vert = 1;
	deb->str_horiz = 0xffffffff;
	deb->str_vert = 0xffffffff;
}

template <typename F>
static inline void luma_intra4x4_with_residual(h264d_mb_current *mb, dec_bits *st, uint32_t cbp, int avail, int avail_intra, const int8_t *pr, int stride,
						    F ResidualBlock)
{
	int coeff[16];
	uint32_t top, left;
	int c0, c1, c2, c3, c4, c5;
	uint8_t *luma = mb->luma;
	const int *offset = mb->offset4x4;
	const int16_t *qmat = mb->qmaty;

	if (cbp & 1) {
		intra4x4pred_func[*pr++](luma, stride, avail_intra | (avail_intra & 2 ? 4 : 0));
		c0 = ResidualBlock(mb, avail & 1 ? UNPACK(mb->left4x4coef, 0) : -1, avail & 2 ? UNPACK(*mb->top4x4coef, 0) : -1, st, coeff, 16, qmat, avail_intra, 0, 2, 0xf);
		ac4x4transform_maybe(luma, coeff, stride, c0);
		intra4x4pred_func[*pr++](luma + 4, stride, avail_intra | (avail_intra & 2 ? 5 : 1));
		c1 = ResidualBlock(mb, c0, avail & 2 ? UNPACK(*mb->top4x4coef, 1) : -1, st, coeff, 16, qmat, avail_intra, 1, 2, 0xf);
		ac4x4transform_maybe(luma + 4, coeff, stride, c1);
		intra4x4pred_func[*pr++](luma + offset[2], stride, avail_intra | 6);
		c2 = ResidualBlock(mb, avail & 1 ? UNPACK(mb->left4x4coef, 1) : -1, c0, st, coeff, 16, qmat, avail_intra, 2, 2, 0xf);
		ac4x4transform_maybe(luma + offset[2], coeff, stride, c2);
		intra4x4pred_func[*pr++](luma + offset[3], stride, 3);
		c3 = ResidualBlock(mb, c2, c1, st, coeff, 16, qmat, avail_intra, 3, 2, 0xf);
		ac4x4transform_maybe(luma + offset[3], coeff, stride, c3);
	} else {
		intra4x4pred_func[*pr++](luma, stride, avail_intra | (avail_intra & 2 ? 4 : 0));
		intra4x4pred_func[*pr++](luma + 4, stride, avail_intra | (avail_intra & 2 ? 5 : 1));
		intra4x4pred_func[*pr++](luma + offset[2], stride, avail_intra | 6);
		intra4x4pred_func[*pr++](luma + offset[3], stride, 3);
		c0 = 0;
		c1 = 0;
		c2 = 0;
		c3 = 0;
	}
	if (cbp & 2) {
		intra4x4pred_func[*pr++](luma + offset[4], stride, avail_intra | (avail_intra & 2 ? 5 : 1));
		c0 = ResidualBlock(mb, c1, avail & 2 ? UNPACK(*mb->top4x4coef, 2) : -1, st, coeff, 16, qmat, avail_intra, 4, 2, 0xf);
		ac4x4transform_maybe(luma + offset[4], coeff, stride, c0);
		intra4x4pred_func[*pr++](luma + offset[5], stride, avail_intra | 1);
		c1 = ResidualBlock(mb, c0, avail & 2 ? UNPACK(*mb->top4x4coef, 3) : -1, st, coeff, 16, qmat, avail_intra, 5, 2, 0xf);
		left = PACK(0, c1, 0);
		ac4x4transform_maybe(luma + offset[5], coeff, stride, c1);
		intra4x4pred_func[*pr++](luma + offset[6], stride, 7);
		c4 = ResidualBlock(mb, c3, c0, st, coeff, 16, qmat, avail_intra, 6, 2, 0xf);
		ac4x4transform_maybe(luma + offset[6], coeff, stride, c4);
		intra4x4pred_func[*pr++](luma + offset[7], stride, 3);
		c5 = ResidualBlock(mb, c4, c1, st, coeff, 16, qmat, avail_intra, 7, 2, 0xf);
		left = PACK(left, c5, 1);
		ac4x4transform_maybe(luma + offset[7], coeff, stride, c5);
	} else {
		intra4x4pred_func[*pr++](luma + offset[4], stride, avail_intra | (avail_intra & 2 ? 5 : 1));
		intra4x4pred_func[*pr++](luma + offset[5], stride, avail_intra | 1);
		intra4x4pred_func[*pr++](luma + offset[6], stride, 7);
		intra4x4pred_func[*pr++](luma + offset[7], stride, 3);
		c0 = 0;
		c1 = 0;
		c4 = 0;
		c5 = 0;
		left = 0;
	}
	if (cbp & 4) {
		intra4x4pred_func[*pr++](luma + offset[8], stride, avail_intra | 6);
		c0 = ResidualBlock(mb, avail & 1 ? UNPACK(mb->left4x4coef, 2) : -1, c2, st, coeff, 16, qmat, avail_intra, 8, 2, 0xf);
		ac4x4transform_maybe(luma + offset[8], coeff, stride, c0);
		intra4x4pred_func[*pr++](luma + offset[9], stride, 7);
		c1 = ResidualBlock(mb, c0, c3, st, coeff, 16, qmat, avail_intra, 9, 2, 0xf);
		ac4x4transform_maybe(luma + offset[9], coeff, stride, c1);
		intra4x4pred_func[*pr++](luma + offset[10], stride, avail_intra | 6);
		c2 = ResidualBlock(mb, avail & 1 ? UNPACK(mb->left4x4coef, 3) : -1, c0, st, coeff, 16, qmat, avail_intra, 10, 2, 0xf);
		top = PACK(0, c2, 0);
		ac4x4transform_maybe(luma + offset[10], coeff, stride, c2);
		intra4x4pred_func[*pr++](luma + offset[11], stride, 3);
		c3 = ResidualBlock(mb, c2, c1, st, coeff, 16, qmat, avail_intra, 11, 2, 0xf);
		top = PACK(top, c3, 1);
		ac4x4transform_maybe(luma + offset[11], coeff, stride, c3);
	} else {
		intra4x4pred_func[*pr++](luma + offset[8], stride, avail_intra | 6);
		intra4x4pred_func[*pr++](luma + offset[9], stride, 7);
		intra4x4pred_func[*pr++](luma + offset[10], stride, avail_intra | 6);
		intra4x4pred_func[*pr++](luma + offset[11], stride, 3);
		c0 = 0;
		c1 = 0;
		c2 = 0;
		c3 = 0;
		top = 0;
	}
	if (cbp & 8) {
		intra4x4pred_func[*pr++](luma + offset[12], stride, 7);
		c0 = ResidualBlock(mb, c1, c4, st, coeff, 16, qmat, avail_intra, 12, 2, 0xf);
		ac4x4transform_maybe(luma + offset[12], coeff, stride, c0);
		intra4x4pred_func[*pr++](luma + offset[13], stride, 3);
		c1 = ResidualBlock(mb, c0, c5, st, coeff, 16, qmat, avail_intra, 13, 2, 0xf);
		left = PACK(left, c1, 2);
		ac4x4transform_maybe(luma + offset[13], coeff, stride, c1);
		intra4x4pred_func[*pr++](luma + offset[14], stride, 7);
		c2 = ResidualBlock(mb, c3, c0, st, coeff, 16, qmat, avail_intra, 14, 2, 0xf);
		top = PACK(top, c2, 2);
		ac4x4transform_maybe(luma + offset[14], coeff, stride, c2);
		intra4x4pred_func[*pr++](luma + offset[15], stride, 3);
		c3 = ResidualBlock(mb, c2, c1, st, coeff, 16, qmat, avail_intra, 15, 2, 0xf);
		ac4x4transform_maybe(luma + offset[15], coeff, stride, c3);
	} else {
		intra4x4pred_func[*pr++](luma + offset[12], stride, 7);
		intra4x4pred_func[*pr++](luma + offset[13], stride, 3);
		intra4x4pred_func[*pr++](luma + offset[14], stride, 7);
		intra4x4pred_func[*pr++](luma + offset[15], stride, 3);
		c3 = 0; 
	}
	mb->left4x4coef = (mb->left4x4coef & 0xffff0000) | PACK(left, c3, 3);
	*mb->top4x4coef = (*mb->top4x4coef & 0xffff0000) | PACK(top, c3, 3);
}

static inline void luma_intra4x4_pred(h264d_mb_current *mb, int avail_intra, const int8_t *pr, int stride)
{
	uint8_t *luma = mb->luma;
	const int *offset = mb->offset4x4;
	intra4x4pred_func[*pr++](luma, stride, avail_intra | (avail_intra & 2 ? 4 : 0));
	intra4x4pred_func[*pr++](luma + 4, stride, avail_intra | (avail_intra & 2 ? 5 : 1));
	intra4x4pred_func[*pr++](luma + offset[2], stride, avail_intra | 6);
	intra4x4pred_func[*pr++](luma + offset[3], stride, 3);
	intra4x4pred_func[*pr++](luma + offset[4], stride, avail_intra | (avail_intra & 2 ? 5 : 1));
	intra4x4pred_func[*pr++](luma + offset[5], stride, avail_intra | 1);
	intra4x4pred_func[*pr++](luma + offset[6], stride, 7);
	intra4x4pred_func[*pr++](luma + offset[7], stride, 3);
	intra4x4pred_func[*pr++](luma + offset[8], stride, avail_intra | 6);
	intra4x4pred_func[*pr++](luma + offset[9], stride, 7);
	intra4x4pred_func[*pr++](luma + offset[10], stride, avail_intra | 6);
	intra4x4pred_func[*pr++](luma + offset[11], stride, 3);
	intra4x4pred_func[*pr++](luma + offset[12], stride, 7);
	intra4x4pred_func[*pr++](luma + offset[13], stride, 3);
	intra4x4pred_func[*pr++](luma + offset[14], stride, 7);
	intra4x4pred_func[*pr](luma + offset[15], stride, 3);
	mb->left4x4coef &= 0xffff0000;
	*mb->top4x4coef &= 0xffff0000;
}

template <typename F0, typename F1, typename F2, typename F3, typename F4>
static inline int mb_intra4x4(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail,
				F0 Intra4x4PredMode,
				F1 IntraChromaPredMode,
				F2 CodedBlockPattern,
				F3 QpDelta,
				F4 ResidualBlock)
{
	int8_t pred4x4[16];
	uint32_t intra_chroma_pred_mode;
	int stride;
	uint32_t cbp;
	int avail_intra;

	avail_intra = avail;
	if (mb->is_constrained_intra) {
		avail_intra &= ~((MB_IPCM < mb->top4x4inter[1].type) * 4 | ((MB_IPCM < mb->top4x4inter->type) * 2) | (MB_IPCM < mb->left4x4inter->type));
	}
	fill_dc_if_unavailable(mb, avail_intra);
	mb_pred_intra4x4(mb, st, avail_intra, pred4x4, Intra4x4PredMode);
	VC_CHECK;
	intra_chroma_pred_mode = IntraChromaPredMode(mb, st, avail_intra);
	stride = mb->max_x * 16;
	intra_chroma_pred[intra_chroma_pred_mode](mb->chroma, stride, avail_intra);
	cbp = CodedBlockPattern(mb, st, avail);
	if (cbp) {
		int32_t qp_delta = QpDelta(mb, st, avail);
		if (qp_delta) {
			set_qp(mb, mb->qp + qp_delta);
		}
	} else {
		mb->prev_qp_delta = 0;
	}
	if (cbp & 15) {
		luma_intra4x4_with_residual(mb, st, cbp, avail, avail_intra, pred4x4, stride, ResidualBlock);
	} else {
		luma_intra4x4_pred(mb, avail, pred4x4, stride);
	}
	store_strength_intra(mb);
	mb_intra_save_info(mb);
	mb->cbp = cbp;
	VC_CHECK;
	return residual_chroma(mb, cbp, st, avail, ResidualBlock);
}

struct intra_chroma_pred_mode_cavlc {
	uint32_t operator()(h264d_mb_current *mb, dec_bits *st, int avail) const {
		uint32_t pred_mode = ue_golomb(st);
		pred_mode = pred_mode <= 3 ? pred_mode : 0;
		mb->chroma_pred_mode = pred_mode;
		return pred_mode;
	}
};

struct cbp_intra_cavlc {
	uint32_t operator()(h264d_mb_current *mb, dec_bits *st, int avail) const {
		return me_golomb(st, me_golomb_lut[0]);
	}
};

struct cbp_inter_cavlc {
	uint32_t operator()(h264d_mb_current *mb, dec_bits *st, int avail) const {
		return me_golomb(st, me_golomb_lut[1]);
	}
};

struct qp_delta_cavlc {
	int operator()(h264d_mb_current *mb, dec_bits *st, int avail) const {
		int delta = se_golomb(st);
		return (delta < -26) ? -26 : ((25 < delta) ? 25 : delta);
	}
};

static int mb_intra4x4_cavlc(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_intra4x4(mb, mbc, st, avail, intra4x4pred_mode_cavlc(), intra_chroma_pred_mode_cavlc(), cbp_intra_cavlc(), qp_delta_cavlc(), residual_block_cavlc());
} 

static int mb_intra16x16pred_dc(uint8_t *dst, int stride, int avail)
{
	uint32_t dc;
	int i;
	
	if (avail & 1) {
		if (avail & 2) {
			dc = (sum_left<16>(dst, stride) + sum_top<16>(dst, stride) + 16) >> 5;
		} else {
			dc = (sum_left<16>(dst, stride) + 8) >> 4;
		}
	} else if (avail & 2) {
		dc = (sum_top<16>(dst, stride) + 8) >> 4;
	} else {
		dc = 0x80;
	}
	dc = dc * 0x01010101U;
	i = 16;
	do {
		*(uint32_t *)dst = dc;
		*((uint32_t *)dst + 1) = dc;
		*((uint32_t *)dst + 2) = dc;
		*((uint32_t *)dst + 3) = dc;
		dst += stride;
	} while (--i);
	return 0;
}

static int mb_intra16x16pred_horiz(uint8_t *dst, int stride, int avail)
{
	int i;

	if (!(avail & 1)) {
		return -1;
	}
	i = 16;
	do {
		uint32_t t0 = dst[-1] * 0x01010101U;
		*(uint32_t *)dst = t0;
		*((uint32_t *)dst + 1) = t0;
		*((uint32_t *)dst + 2) = t0;
		*((uint32_t *)dst + 3) = t0;
		dst = dst + stride;
	} while (--i);
	return 0;
}

static int mb_intra16x16pred_planer(uint8_t *dst, int stride, int avail)
{
	const uint8_t *src, *src2;
	int t0, p0, h, v;
	int y;
	
	src = dst - stride;
	p0 = src[15];
	src -= 1;
	t0 = p0 - src[0];
	h = t0;
	t0 = t0 + src[15] - src[1];
	h += t0;
	t0 = t0 + src[14] - src[2];
	h += t0;
	t0 = t0 + src[13] - src[3];
	h += t0;
	t0 = t0 + src[12] - src[4];
	h += t0;
	t0 = t0 + src[11] - src[5];
	h += t0;
	t0 = t0 + src[10] - src[6];
	h += t0;
	t0 = t0 + src[9] - src[7];
	h += t0;
	h = ((h * 5) + 32) >> 6;

	src = dst - 1;
	src2 = src + (stride * 15);
	src = src - stride;
	t0 = src2[0];
	p0 = (p0 + t0) * 16;
	t0 = t0 - src[0];
	v = t0;
	src2 -= stride;
	src += stride;
	t0 = t0 + src2[0] - src[0];
	v += t0;
	src2 -= stride;
	src += stride;
	t0 = t0 + src2[0] - src[0];
	v += t0;
	src2 -= stride;
	src += stride;
	t0 = t0 + src2[0] - src[0];
	v += t0;
	src2 -= stride;
	src += stride;
	t0 = t0 + src2[0] - src[0];
	v += t0;
	src2 -= stride;
	src += stride;
	t0 = t0 + src2[0] - src[0];
	v += t0;
	src2 -= stride;
	src += stride;
	t0 = t0 + src2[0] - src[0];
	v += t0;
	src2 -= stride;
	src += stride;
	t0 = t0 + src2[0] - src[0];
	v += t0;
	v = ((v * 5) + 32) >> 6;

	dst += 16 + (stride * 15);
	p0 = p0 + ((h + v) * 8) + 16;
	y = 16;
	do {
		uint8_t *d = dst;
		int x = 16;
		t0 = p0;
		d = dst;
		do {
			int s = t0 >> 5;
			*--d = CLIP255C(s);
			t0 -= h;
		} while (--x);
		p0 -= v;
		dst -= stride;
	} while (--y);
	return 0;
}

/** Inverse 16x16 luma DC transformation with inverse zigzag scan.
 * Output is 4x4 block scan order.
 */
static void intra16x16_dc_transform(const int *src, int *dst)
{
	int c0, c1, c2, c3;
	int t0, t1;

	c0 = src[0] + src[1] + src[5] + src[6];
	c1 = src[2] + src[4] + src[7] + src[12];
	c2 = src[3] + src[8] + src[11] + src[13];
	c3 = src[9] + src[10] + src[14] + src[15];
	t0 = c0 + c1;
	t1 = c2 + c3;
	dst[0] = (t0 + t1 + 2) >> 2;
	dst[2] = (t0 - t1 + 2) >> 2;
	t0 = c0 - c1;
	t1 = c2 - c3;
	dst[8] = (t0 - t1 + 2) >> 2;
	dst[10] = (t0 + t1 + 2) >> 2;

	c0 = src[0] + src[1] - src[5] - src[6];
	c1 = src[2] + src[4] - src[7] - src[12];
	c2 = src[3] + src[8] - src[11] - src[13];
	c3 = src[9] + src[10] - src[14] - src[15];
	t0 = c0 + c1;
	t1 = c2 + c3;
	dst[1] = (t0 + t1 + 2) >> 2;
	dst[3] = (t0 - t1 + 2) >> 2;
	t0 = c0 - c1;
	t1 = c2 - c3;
	dst[9] = (t0 - t1 + 2) >> 2;
	dst[11] = (t0 + t1 + 2) >> 2;

	c0 = src[0] - src[1] - src[5] + src[6];
	c1 = src[2] - src[4] - src[7] + src[12];
	c2 = src[3] - src[8] - src[11] + src[13];
	c3 = src[9] - src[10] - src[14] + src[15];
	t0 = c0 + c1;
	t1 = c2 + c3;
	dst[4] = (t0 + t1 + 2) >> 2;
	dst[6] = (t0 - t1 + 2) >> 2;
	t0 = c0 - c1;
	t1 = c2 - c3;
	dst[12] = (t0 - t1 + 2) >> 2;
	dst[14] = (t0 + t1 + 2) >> 2;

	c0 = src[0] - src[1] + src[5] - src[6];
	c1 = src[2] - src[4] + src[7] - src[12];
	c2 = src[3] - src[8] + src[11] - src[13];
	c3 = src[9] - src[10] + src[14] - src[15];
	t0 = c0 + c1;
	t1 = c2 + c3;
	dst[5] = (t0 + t1 + 2) >> 2;
	dst[7] = (t0 - t1 + 2) >> 2;
	t0 = c0 - c1;
	t1 = c2 - c3;
	dst[13] = (t0 - t1 + 2) >> 2;
	dst[15] = (t0 + t1 + 2) >> 2;
}

struct AddSaturate {
	uint32_t operator()(uint32_t x, uint32_t y) const {
		uint32_t msk;
		msk = ((x & y) + (((x ^ y) >> 1) & 0x7f7f7f7f)) & ~0x7f7f7f7f;
		msk = (msk << 1) - (msk >> 7);
		return ((x + y) - msk) | msk;
	}
};

struct SubSaturate {
	uint32_t operator()(uint32_t x, uint32_t y) const {
		uint32_t msk;
		msk = ((~x & y) + (((~x ^ y) >> 1) & 0x7f7f7f7f)) & ~0x7f7f7f7f;
		msk = (msk << 1) - (msk >> 7);
		return (x | msk) - (y | msk);
	}
};

template<typename T>
static inline void ac4x4transform_dconly_part(uint8_t *dst, uint32_t dc, int stride, T saturate)
{
	int y = 4;
	dc = dc * 0x01010101;
	do {
		*(uint32_t *)dst = saturate(*(uint32_t *)dst, dc);
		dst += stride;
	} while (--y);
}

static void ac4x4transform_dconly(uint8_t *dst, int dc, int stride)
{
	dc = (dc + 32) >> 6;
	if (dc < 0) {
		ac4x4transform_dconly_part(dst, (uint32_t)(-dc), stride, SubSaturate());
	} else {
		ac4x4transform_dconly_part(dst, (uint32_t)dc, stride, AddSaturate());
	}
}

static void ac4x4transform_dconly_chroma(uint8_t *dst, int dc, int stride)
{
	int y;
	dc = (dc + 32) >> 6;
	y = 4;
	do {
		int t;
		t = dst[0] + dc;
		dst[0] = CLIP255C(t);
		t = dst[2] + dc;
		dst[2] = CLIP255C(t);
		t = dst[4] + dc;
		dst[4] = CLIP255C(t);
		t = dst[6] + dc;
		dst[6] = CLIP255C(t);
		dst += stride;
	} while (--y);
}

static inline void transform4x4_vert(int *dst, int d0, int d1, int d2, int d3)
{
	int t0 = d0 + d2;
	int t1 = d0 - d2;
	int t2 = (d1 >> 1) - d3;
	int t3 = d1 + (d3 >> 1);
	dst[0] = t0 + t3;
	dst[4] = t1 + t2;
	dst[8] = t1 - t2;
	dst[12] = t0 - t3;
}

template <int N>
static inline void transform4x4_horiz_loop(uint8_t *dst, const int *src, int stride)
{
	int x = 4;
	do {
		int e0, e1, e2, e3;
		int f0, f1, f2, f3;
		uint8_t *d;

		int t0;
		e0 = *src++;
		e1 = *src++;
		e2 = *src++;
		e3 = *src++;
		d = dst;
		f0 = e0 + e2;
		f1 = e0 - e2;
		f2 = (e1 >> 1) - e3;
		f3 = e1 + (e3 >> 1);
		t0 = *d + ((f0 + f3) >> 6);
		*d = CLIP255C(t0);
		d += stride;
		t0 = *d + ((f1 + f2) >> 6);
		*d = CLIP255C(t0);
		d += stride;
		t0 = *d + ((f1 - f2) >> 6);
		*d = CLIP255C(t0);
		d += stride;
		t0 = *d + ((f0 - f3) >> 6);
		*d = CLIP255C(t0);
		dst += N;
	} while (--x);
}

/** Read coefficients in inverse zigzag order and then reconstruct.
 */
static void ac4x4transform_acdc(uint8_t *dst, const int *coeff, int stride)
{
	int tmp[16];
	transform4x4_vert(tmp, coeff[0] + 32, coeff[1], coeff[5], coeff[6]);
	transform4x4_vert(tmp + 1, coeff[2], coeff[4], coeff[7], coeff[12]);
	transform4x4_vert(tmp + 2, coeff[3], coeff[8], coeff[11], coeff[13]);
	transform4x4_vert(tmp + 3, coeff[9], coeff[10], coeff[14], coeff[15]);
	transform4x4_horiz_loop<1>(dst, tmp, stride);
}

static inline void ac4x4transform_maybe(uint8_t *dst, const int *coeff, int stride, int num_coeff)
{
	if (num_coeff) {
		ac4x4transform_acdc(dst, coeff, stride);
	}
}

static inline void ac4x4transform(uint8_t *dst, int *coeff, int stride, int num_coeff, int dc)
{
	if (num_coeff) {
		coeff[0] = dc;
		ac4x4transform_acdc(dst, coeff, stride);
	} else {
		ac4x4transform_dconly(dst, dc, stride);
	}
}

/** Read coefficients in inverse zigzag order and then reconstruct, chroma part.
 */
static void ac4x4transform_acdc_chroma(uint8_t *dst, const int *coeff, int stride)
{
	int tmp[16];
	transform4x4_vert(tmp, coeff[0] + 32, coeff[1], coeff[5], coeff[6]);
	transform4x4_vert(tmp + 1, coeff[2], coeff[4], coeff[7], coeff[12]);
	transform4x4_vert(tmp + 2, coeff[3], coeff[8], coeff[11], coeff[13]);
	transform4x4_vert(tmp + 3, coeff[9], coeff[10], coeff[14], coeff[15]);
	transform4x4_horiz_loop<2>(dst, tmp, stride);
}

/** Inverse 8x8 chroma DC transformation with inverse zigzag scan.
 * Output is 4x4 block scan order.
 */
static void intra_chroma_dc_transform(const int *src, int *dst)
{
	int c0, c1, c2, c3;
	int t0, t1;

	c0 = src[0];
	c1 = src[1];
	c2 = src[2];
	c3 = src[3];
	t0 = c0 + c1;
	t1 = c2 + c3;
	dst[0] = (t0 + t1) >> 1;
	dst[2] = (t0 - t1) >> 1;
	t0 = c0 - c1;
	t1 = c2 - c3;
	dst[1] = (t0 + t1) >> 1;
	dst[3] = (t0 - t1) >> 1;
}

template <typename F0, typename F1, typename F2>
static int mb_intra16x16_dconly(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail,
				F0 IntraChromaPredMode,
				F1 QpDelta,
				F2 ResidualBlock)
{
	int coeff[16];
	int dc[16];
	uint8_t *luma;
	int stride;
	int avail_intra;
	uint32_t intra_chroma_pred_mode;
	int32_t qp_delta;
	const int *offset;

	luma = mb->luma;
	stride = mb->max_x * 16;
	avail_intra = avail;
	if (mb->is_constrained_intra) {
		avail_intra &= ~((MB_IPCM < mb->top4x4inter[1].type) * 4 | ((MB_IPCM < mb->top4x4inter->type) * 2) | (MB_IPCM < mb->left4x4inter->type));
	}
	mbc->mb_pred(luma, stride, avail_intra);
	intra_chroma_pred_mode = IntraChromaPredMode(mb, st, avail_intra);
	intra_chroma_pred[intra_chroma_pred_mode](mb->chroma, stride, avail_intra);
	qp_delta = QpDelta(mb, st, avail);
	if (qp_delta) {
		set_qp(mb, mb->qp + qp_delta);
	}
	if (ResidualBlock(mb, avail & 1 ? UNPACK(mb->left4x4coef, 0) : -1, avail & 2 ? UNPACK(*mb->top4x4coef, 0) : -1, st, coeff, 16, mb->qmaty, avail_intra, 26, 0, 0)) {
		intra16x16_dc_transform(coeff, dc);
		offset = mb->offset4x4;
		for (int i = 0; i < 16; ++i) {
			ac4x4transform_dconly(luma + *offset++, dc[i], stride);
		}
	}
	mb->left4x4coef &= 0xffff0000;
	*mb->top4x4coef &= 0xffff0000;
	mb->left4x4pred = 0x22222222;
	*mb->top4x4pred = 0x22222222;
	store_strength_intra(mb);
	mb_intra_save_info(mb);
	mb->cbp = mbc->cbp;
	return residual_chroma(mb, mbc->cbp, st, avail, ResidualBlock);
}

static int mb_intra16x16_dconly_cavlc(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_intra16x16_dconly(mb, mbc, st, avail, intra_chroma_pred_mode_cavlc(), qp_delta_cavlc(), residual_block_cavlc());
}

template <typename F0, typename F1, typename F2>
static int mb_intra16x16_acdc(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail,
				F0 IntraChromaPredMode,
				F1 QpDelta,
				F2 ResidualBlock)
{
	int dc[16];
	int coeff[16];
	uint8_t *luma;
	int stride;
	int avail_intra;
	uint32_t intra_chroma_pred_mode;
	int32_t qp_delta;
	const int16_t *qmat;
	uint32_t top, left;
	int c0, c1, c2, c3, c4, c5;
	int na, nb;
	const int *offset;
	int *dcp;

	luma = mb->luma;
	stride = mb->max_x * 16;
	avail_intra = avail;
	if (mb->is_constrained_intra) {
		avail_intra &= ~((MB_IPCM < mb->top4x4inter[1].type) * 4 | ((MB_IPCM < mb->top4x4inter->type) * 2) | (MB_IPCM < mb->left4x4inter->type));
	}
	mbc->mb_pred(luma, stride, avail_intra);
	intra_chroma_pred_mode = IntraChromaPredMode(mb, st, avail_intra);
	intra_chroma_pred[intra_chroma_pred_mode](mb->chroma, stride, avail_intra);
	qp_delta = QpDelta(mb, st, avail);
	if (qp_delta) {
		set_qp(mb, mb->qp + qp_delta);
	}

	na = avail & 1 ? UNPACK(mb->left4x4coef, 0) : -1;
	nb = avail & 2 ? UNPACK(*mb->top4x4coef, 0) : -1;
	qmat = mb->qmaty;
	if (ResidualBlock(mb, na, nb, st, coeff, 16, qmat, avail_intra, 26, 0, 0)) {
		intra16x16_dc_transform(coeff, dc);
	} else {
		memset(dc, 0, sizeof(dc));
	}
	offset = mb->offset4x4;
	dcp = dc;
	c0 = ResidualBlock(mb, na, nb, st, coeff, 15, qmat, avail_intra, 0, 1, 0x1f);
	ac4x4transform(luma + *offset++, coeff, stride, c0, *dcp++);
	c1 = ResidualBlock(mb, c0, avail & 2 ? UNPACK(*mb->top4x4coef, 1) : -1, st, coeff, 15, qmat, avail_intra, 1, 1, 0x1f);
	ac4x4transform(luma + *offset++, coeff, stride, c1, *dcp++);
	c2 = ResidualBlock(mb, avail & 1 ? UNPACK(mb->left4x4coef, 1) : -1, c0, st, coeff, 15, qmat, avail_intra, 2, 1, 0x1f);
	ac4x4transform(luma + *offset++, coeff, stride, c2, *dcp++);
	c3 = ResidualBlock(mb, c2, c1, st, coeff, 15, qmat, avail_intra, 3, 1, 0x1f);
	ac4x4transform(luma + *offset++, coeff, stride, c3, *dcp++);

	c0 = ResidualBlock(mb, c1, avail & 2 ? UNPACK(*mb->top4x4coef, 2) : -1, st, coeff, 15, qmat, avail_intra, 4, 1, 0x1f);
	ac4x4transform(luma + *offset++, coeff, stride, c0, *dcp++);
	c1 = ResidualBlock(mb, c0, avail & 2 ? UNPACK(*mb->top4x4coef, 3) : -1, st, coeff, 15, qmat, avail_intra, 5, 1, 0x1f);
	left = mb->left4x4coef & 0xffff0000;
	left = PACK(left, c1, 0);
	ac4x4transform(luma + *offset++, coeff, stride, c1, *dcp++);
	c4 = ResidualBlock(mb, c3, c0, st, coeff, 15, qmat, avail_intra, 6, 1, 0x1f);
	ac4x4transform(luma + *offset++, coeff, stride, c4, *dcp++);
	c5 = ResidualBlock(mb, c4, c1, st, coeff, 15, qmat, avail_intra, 7, 1, 0x1f);
	left = PACK(left, c5, 1);
	ac4x4transform(luma + *offset++, coeff, stride, c5, *dcp++);

	c0 = ResidualBlock(mb, avail & 1 ? UNPACK(mb->left4x4coef, 2) : -1, c2, st, coeff, 15, qmat, avail_intra, 8, 1, 0x1f);
	ac4x4transform(luma + *offset++, coeff, stride, c0, *dcp++);
	c1 = ResidualBlock(mb, c0, c3, st, coeff, 15, qmat, avail_intra, 9, 1, 0x1f);
	ac4x4transform(luma + *offset++, coeff, stride, c1, *dcp++);
	c2 = ResidualBlock(mb, avail & 1 ? UNPACK(mb->left4x4coef, 3) : -1, c0, st, coeff, 15, qmat, avail_intra, 10, 1, 0x1f);
	top = *mb->top4x4coef & 0xffff0000;
	top = PACK(top, c2, 0);
	ac4x4transform(luma + *offset++, coeff, stride, c2, *dcp++);
	c3 = ResidualBlock(mb, c2, c1, st, coeff, 15, qmat, avail_intra, 11, 1, 0x1f);
	top = PACK(top, c3, 1);
	ac4x4transform(luma + *offset++, coeff, stride, c3, *dcp++);

	c0 = ResidualBlock(mb, c1, c4, st, coeff, 15, qmat, avail_intra, 12, 1, 0x1f);
	ac4x4transform(luma + *offset++, coeff, stride, c0, *dcp++);
	c1 = ResidualBlock(mb, c0, c5, st, coeff, 15, qmat, avail_intra, 13, 1, 0x1f);
	left = PACK(left, c1, 2);
	ac4x4transform(luma + *offset++, coeff, stride, c1, *dcp++);
	c2 = ResidualBlock(mb, c3, c0, st, coeff, 15, qmat, avail_intra, 14, 1, 0x1f);
	top = PACK(top, c2, 2);
	ac4x4transform(luma + *offset++, coeff, stride, c2, *dcp++);
	c3 = ResidualBlock(mb, c2, c1, st, coeff, 15, qmat, avail_intra, 15, 1, 0x1f);
	ac4x4transform(luma + *offset++, coeff, stride, c3, *dcp++);

	mb->left4x4coef = PACK(left, c3, 3);
	*mb->top4x4coef = PACK(top, c3, 3);
	mb->left4x4pred = 0x22222222;
	*mb->top4x4pred = 0x22222222;
	store_strength_intra(mb);
	mb_intra_save_info(mb);
	mb->cbp = mbc->cbp;
	return residual_chroma(mb, mbc->cbp, st, avail, ResidualBlock);
}

static int mb_intra16x16_acdc_cavlc(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_intra16x16_acdc(mb, mbc, st, avail, intra_chroma_pred_mode_cavlc(), qp_delta_cavlc(), residual_block_cavlc());
}

/** Sum of top of 4x4 block, NV12 chroma part.
 */
static inline uint32_t sum_top_chroma(const uint8_t *src, int stride) {
	src -= stride;
	return src[0] + src[2] + src[4] + src[6];
}

static inline void fill_4x4_chroma(uint8_t *dst, uint32_t dc, int stride)
{
	int y;
	dc = dc * 0x00010001U;
#ifdef WORDS_BIGENDIAN
	dc = bswap32(dc);
#endif
	y = 4;
	do {
		*(uint32_t *)dst = dc;
		*((uint32_t *)dst + 1) = dc;
		dst += stride;
	} while (--y);
}

static int mb_intra_chroma_pred_dc(uint8_t *dst, int stride, int avail)
{
	uint32_t dc0, dc1, dc2, dc3;

	if (avail & 1) {
		if (avail & 2) {
			uint32_t left0, left1, top0, top1;

			left0 = sum_left<4>(dst - 1, stride);
			left1 = sum_left<4>(dst, stride);
			top0 = sum_top_chroma(dst, stride);
			top1 = sum_top_chroma(dst + 1, stride);
			dc0 = ((left0 + top0 + 4) >> 3) | (((left1 + top1 + 4) >> 3) << 8);

			top0 = sum_top_chroma(dst + 8, stride);
			top1 = sum_top_chroma(dst + 9, stride);
			dc1 = ((top0 + 2) >> 2) | (((top1 + 2) >> 2) << 8);

			left0 = sum_left<4>(dst + 4 * stride - 1, stride);
			left1 = sum_left<4>(dst + 4 * stride, stride);
			dc2 = ((left0 + 2) >> 2) | (((left1 + 2) >> 2) << 8);
			dc3 = ((left0 + top0 + 4) >> 3) | (((left1 + top1 + 4) >> 3) << 8);
		} else {
			dc1 = dc0 = (((sum_left<4>(dst, stride) + 2) >> 2) << 8) | ((sum_left<4>(dst - 1, stride) + 2) >> 2);
			dc3 = dc2 = (((sum_left<4>(dst + 4 * stride, stride) + 2) >> 2) << 8) | ((sum_left<4>(dst + 4 * stride - 1, stride) + 2) >> 2);
		}
	} else if (avail & 2) {
		int l0, l1, t0, t1;
		l0 = sum_top_chroma(dst + 1, stride);
		l1 = sum_top_chroma(dst + 0, stride);
		dc2 = dc0 = (((l0 + 2) >> 2) << 8) | ((l1 + 2) >> 2);
		t0 = sum_top_chroma(dst + 9, stride);
		t1 = sum_top_chroma(dst + 8, stride);
		dc3 = dc1 = (((t0 + 2) >> 2) << 8) | ((t1 + 2) >> 2);
	} else {
		dc0 = dc1 = dc2 = dc3 = 0x8080;
	}
	fill_4x4_chroma(dst, dc0, stride);
	fill_4x4_chroma(dst + 8, dc1, stride);
	fill_4x4_chroma(dst + 4 * stride, dc2, stride);
	fill_4x4_chroma(dst + 4 * (stride + 2), dc3, stride);
	return 0;
}

static int mb_intra_chroma_pred_horiz(uint8_t *dst, int stride, int avail)
{
	int i;

	if (!(avail & 1)) {
		return -1;
	}
	i = 8;
	do {
		uint32_t t0 = *((uint16_t *)dst - 1) * 0x00010001U;
		*(uint32_t *)dst = t0;
		*((uint32_t *)dst + 1) = t0;
		*((uint32_t *)dst + 2) = t0;
		*((uint32_t *)dst + 3) = t0;
		dst = dst + stride;
	} while (--i);
	return 0;
}

static int mb_intra_chroma_pred_planer(uint8_t *dst, int stride, int avail)
{
	const uint8_t *src, *src2;
	int a0, a1, h0, h1, v0, v1;
	int y;

	src = dst - stride + 14;
	a0 = src[0];
	a1 = src[1];
	src = src - 16;
	h0 = ((int)src[10] - src[6]) + ((int)src[12] - src[4]) * 2 + ((int)src[14] - src[2]) * 3 + (a0 - src[0]) * 4;
	h1 = ((int)src[11] - src[7]) + ((int)src[13] - src[5]) * 2 + ((int)src[15] - src[3]) * 3 + (a1 - src[1]) * 4;
	h0 = (h0 * 17 + 16) >> 5;
	h1 = (h1 * 17 + 16) >> 5;

	src = dst + (stride * 7) - 2;
	a0 = (a0 + src[0]) * 16;
	a1 = (a1 + src[1]) * 16;

	src = dst + stride * 4 - 2;
	src2 = src - stride * 2;
	v0 = (int)src[0] - src2[0];
	v1 = (int)src[1] - src2[1];
	src += stride;
	src2 -= stride;
	v0 += ((int)src[0] - src2[0]) * 2;
	v1 += ((int)src[1] - src2[1]) * 2;
	src += stride;
	src2 -= stride;
	v0 += ((int)src[0] - src2[0]) * 3;
	v1 += ((int)src[1] - src2[1]) * 3;
	src += stride;
	src2 -= stride;
	v0 += ((int)src[0] - src2[0]) * 4;
	v1 += ((int)src[1] - src2[1]) * 4;
	v0 = (v0 * 17 + 16) >> 5;
	v1 = (v1 * 17 + 16) >> 5;

	a0 = a0 - ((h0 + v0) * 3) + 16;
	a1 = a1 - ((h1 + v1) * 3) + 16;
	y = 8;
	do {
		int at0 = a0;
		int at1 = a1;
		int x = 8;
		uint8_t *d = dst;
		do {
			int t0, t1;
			t0 = at0 >> 5;
			d[0] = CLIP255C(t0);
			at0 += h0;
			t1 = at1 >> 5;
			d[1] = CLIP255C(t1);
			at1 += h1;
			d += 2;
		} while (--x);
		a0 += v0;
		a1 += v1;
		dst += stride;
	} while (--y);
	return 0;
}

template <int STEP>
static inline void intrapcm_block(uint8_t *dst, int stride, dec_bits *st)
{
	int y = 16 / STEP;
	do {
		for (int x = 0; x < 16; x += (STEP * 2)) {
			dst[x] = get_bits(st, 8);
			dst[x + STEP] = get_bits(st, 8);
		}
		dst += stride;
	} while (--y);
}

static inline void intrapcm_luma(uint8_t *dst, int stride, dec_bits *st)
{
#ifdef WORDS_BIGENDIAN
	int y = 16;
	do {
		*(uint32_t *)dst = get_bits(st, 32);
		*(uint32_t *)(dst + 4) = get_bits(st, 32);
		*(uint32_t *)(dst + 8) = get_bits(st, 32);
		*(uint32_t *)(dst + 12) = get_bits(st, 32);
		dst += stride;
	} while (--y);
#else
	intrapcm_block<1>(dst, stride, st);
#endif
}

static int mb_intrapcm(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	int stride;

	stride = mb->max_x * 16;
	byte_align(st);
	intrapcm_luma(mb->luma, stride, st);
	intrapcm_block<2>(mb->chroma, stride, st);
	intrapcm_block<2>(mb->chroma + 1, stride, st);
	mb->left4x4coef = 0xffffffff;
	*mb->top4x4coef = 0xffffffff;
	mb->left4x4pred = 0x22222222;
	*mb->top4x4pred = 0x22222222;
	mb->deblock_curr->qpy = 0;
	mb->deblock_curr->qpc = 0;
	mb->prev_qp_delta = 0;
	mb->left4x4inter->cbp = 0x3f;
	mb->left4x4inter->cbf = 0xffff;
	mb->top4x4inter->cbp = 0x3f;
	mb->top4x4inter->cbf = 0xffff;
	mb_intra_save_info(mb);
	return 0;
}


static void copy_inter16x_align(const uint8_t *src, uint8_t *dst, int height, int src_stride, int stride)
{
	do {
		((uint32_t *)dst)[0] = ((uint32_t *)src)[0];
		((uint32_t *)dst)[1] = ((uint32_t *)src)[1];
		((uint32_t *)dst)[2] = ((uint32_t *)src)[2];
		((uint32_t *)dst)[3] = ((uint32_t *)src)[3];
		VC_CHECK;
		dst += stride;
		src += src_stride;
	} while (--height);
}

static void copy_inter8x_align(const uint8_t *src, uint8_t *dst, int height, int src_stride, int stride)
{
	do {
		((uint32_t *)dst)[0] = ((uint32_t *)src)[0];
		((uint32_t *)dst)[1] = ((uint32_t *)src)[1];
		dst += stride;
		src += src_stride;
	} while (--height);
}

static void copy_inter4x_align(const uint8_t *src, uint8_t *dst, int height, int src_stride, int stride)
{
	do {
		((uint32_t *)dst)[0] = ((uint32_t *)src)[0];
		dst += stride;
		src += src_stride;
	} while (--height);
}

static void (* const copy_inter_align[3])(const uint8_t *src, uint8_t *dst, int height, int src_stride, int stride) = {
	copy_inter4x_align,
	copy_inter8x_align,
	copy_inter16x_align
};

static inline void copy_inter(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	width = (unsigned)width >> 2;
	if (((intptr_t)src & 3) == 0) {
		copy_inter_align[width >> 1](src, dst, height, src_stride, stride);
		VC_CHECK;
	} else {
		do {
			const uint8_t *s = src;
			uint8_t *d = dst;
			int x = width;
			do {
				d[0] = *s++;
				d[1] = *s++;
				d[2] = *s++;
				d[3] = *s++;
				d += 4;
			} while (--x);
			dst += stride;
			src += src_stride;
		} while (--height);
	}
}

static inline int inter_pred_mvoffset_luma(int mvint_x, int mvint_y, int stride)
{
	return mvint_y * stride + mvint_x;
}

typedef struct {
	const uint8_t *src_luma;
	const uint8_t *src_chroma;
	uint8_t *dst_chroma;
	int16_t pos_x, pos_y;
} mb_pred_t;

static inline void filter_chroma_horiz(const uint8_t *src, uint8_t *dst, int width, int height, int frac, int src_stride, int dst_stride)
{
	int c0 = (8 - frac) * 8;
	int c1 = frac * 8;
	width >>= 1;
	do {
		int x = width;
		const uint8_t *s = src;
		uint8_t *d = dst;
		int s0 = *s++;
		int s1 = *s++;
		do {
			int s2 = *s++;
			int s3 = *s++;
			d[0] = (s2 * c1 + s0 * c0 + 32) >> 6;
			d[1] = (s3 * c1 + s1 * c0 + 32) >> 6;
			s0 = s2;
			s1 = s3;
			d += 2;
		} while (--x);
		src += src_stride;
		dst += dst_stride;
	} while (--height);
}

static inline void filter_chroma_vert(const uint8_t *src, uint8_t *dst, int width, int height, int frac, int src_stride, int dst_stride)
{
	int c0 = (8 - frac) * 8;
	int c1 = frac * 8;

	do {
		const uint8_t *s = src++;
		uint8_t *d = dst++;
		int t0 = *s;
		s += src_stride;
		int y = height;
		do {
			int t1 = *s;
			*d = (t0 * c0 + t1 * c1 + 32) >> 6;
			t0 = t1;
			s += src_stride;
			d += dst_stride;
		} while (--y);
	} while (--width);
}

static inline void filter_chroma_vert_horiz(const uint8_t *src, uint8_t *dst, int width, int height, int fracx, int fracy, int src_stride, int stride)
{
	int c0, c1, c2, c3;

	c0 = (8 - fracx) * (8 - fracy);
	c3 = fracx * fracy;
	c1 = fracx * (8 - fracy);
	c2 = (8 - fracx) * fracy;

	do {
		const uint8_t *src0, *src1;
		uint8_t *d = dst;
		int x = width;
		src0 = src;
		src += src_stride;
		src1 = src;
		do {
			int t = src0[2] * c1 + src1[2] * c3 + 32;
			*d++ = (t + *src0++ * c0 + *src1++ * c2) >> 6;
		} while (--x);
		dst += stride;
	} while (--height);
}

static inline void extend_left_chroma(uint8_t *dst, int left, int width, int height)
{
	do {
		int c = *(int16_t *)(dst + left);
		int x = left >> 1;
		int16_t *d = (int16_t *)dst;
		do {
			*d++ = c;
		} while (--x);
		dst += width;
	} while (--height);
}

static inline void fill_left_top_chroma(const uint8_t *src, uint8_t *buf, int left, int top, int width, int height, int stride)
{
	uint8_t *dst;
	int y;

	src += top * stride + left;
	left = (width == left) ? left - 2 : left;
	dst = buf + left;
	for (y = 0; y < top; ++y) {
		memcpy(dst, src, width - left);
		dst += width;
	}
	for (; y < height; ++y) {
		memcpy(dst, src, width - left);
		src += stride;
		dst += width;
	}
	if (left) {
		extend_left_chroma(buf, left, width, height);
	}
}

static inline void fill_left_bottom_chroma(const uint8_t *src, uint8_t *buf, int left, int bottom, int width, int height, int stride)
{
	uint8_t *dst;
	int y;

	src += left;
	dst = buf + left;
	for (y = 0; y < height - bottom; ++y) {
		memcpy(dst, src, width - left);
		src += stride;
		dst += width;
	}
	src = dst - width;
	for (; y < height; ++y) {
		memcpy(dst, src, width - left);
		dst += width;
	}
	if (left) {
		extend_left_chroma(buf, left, width, height);
	}
}

static inline void extend_right_chroma(uint8_t *dst, int right, int width, int height)
{
	dst = dst + width - right;
	do {
		int c = ((int16_t *)dst)[-1];
		int x = right >> 1;
		int16_t *d = (int16_t *)dst;
		do {
			*d++ = c;
		} while (--x);
		dst += width;
	} while (--height);
}

static inline void fill_right_top_chroma(const uint8_t *src, uint8_t *buf, int right, int top, int width, int height, int stride)
{
	uint8_t *dst;
	int y;

	src = src + top * stride;
	dst = buf;
	for (y = 0; y < top; ++y) {
		memcpy(dst, src, width - right);
		dst += width;
	}
	for (; y < height; ++y) {
		memcpy(dst, src, width - right);
		src += stride;
		dst += width;
	}
	if (right) {
		extend_right_chroma(buf, right, width, height);
	}
}

static inline void fill_right_bottom_chroma(const uint8_t *src, uint8_t *buf, int right, int bottom, int width, int height, int stride)
{
	uint8_t *dst;
	int y;

	dst = buf;
	for (y = 0; y < height - bottom; ++y) {
		memcpy(dst, src, width - right);
		src += stride;
		dst += width;
	}
	src -= stride;
	for (; y < height; ++y) {
		memcpy(dst, src, width - right);
		dst += width;
	}
	if (right) {
		extend_right_chroma(buf, right, width, height);
	}
}

static inline void fill_rect_umv_chroma(const uint8_t *src, uint8_t *buf, int width, int height, int stride, int vert_size, int posx, int posy)
{
	int left, right, top, bottom;

	left = -posx;
	top = -posy;
	if (0 < left) {
		if (0 < top) {
			/* left, top */
			fill_left_top_chroma(src, buf, left, top, width, height, stride);
		} else {
			bottom = posy - vert_size + height;
			if (0 < bottom) {
				/* left, bottom */
				fill_left_bottom_chroma(src, buf, left, bottom, width, height, stride);
			} else {
				/* left */
				fill_left_top_chroma(src, buf, left, 0, width, height, stride);
			}
		}
	} else {
		right = posx - stride + width;
		if (0 < top) {
			if (0 < right) {
				/* top, right */
				fill_right_top_chroma(src, buf, right, top, width, height, stride);
			} else {
				/* top */
				fill_left_top_chroma(src, buf, 0, top, width, height, stride);
			}
		} else {
			bottom = posy - vert_size + height;
			if (0 < right) {
				if (0 < bottom) {
					/* right, bottom */
					fill_right_bottom_chroma(src, buf, right, bottom, width, height, stride);
				} else {
					/* right */
					fill_right_top_chroma(src, buf, right, 0, width, height, stride);
				}
			} else {
				if (0 < bottom) {
					/* bottom */
					fill_right_bottom_chroma(src, buf, 0, bottom, width, height, stride);
				} else {
					/* in range */
					fill_right_top_chroma(src, buf, 0, 0, width, height, stride);
				}
			}
		}
	}
}


typedef struct {
	int mvx, mvy;
} chroma_filter_info_t;

static inline void chroma_inter_umv(const uint8_t *src, uint8_t *dst, int posx, int posy, int width, int height, int src_stride, int vert_size, int dst_stride, chroma_filter_info_t *filter)
{
	uint32_t buf[18 * 9 / sizeof(uint32_t) + 1];

	if (posx < -width) {
		src += -width - posx;
		posx = -width;
	} else if (src_stride - 2 < posx) {
		src -= posx - src_stride + 2;
		posx = src_stride - 2;
	}
	if (posy < -height) {
		src -= (height + posy) * src_stride;
		posy = -height;
	} else if (vert_size - 1 < posy) {
		src -= (posy - vert_size + 1) * src_stride;
		posy = vert_size - 1;
	}

	if (filter) {
		width += 2;
		height += 1;
	}
	fill_rect_umv_chroma(src, (uint8_t *)buf, width, height, src_stride, vert_size, posx, posy);
	if (filter) {
		width -= 2;
		height -= 1;
		filter_chroma_vert_horiz((const uint8_t *)buf, dst, width, height, filter->mvx, filter->mvy, width + 2, dst_stride);
	} else {
		copy_inter_align[(unsigned)width >> 3]((const uint8_t *)buf, dst, height, width, dst_stride);
	}
}

static void inter_pred_chroma(const mb_pred_t *pred, int mvx, int mvy, int width, int height, int src_stride, int vert_stride, int dst_stride, uint8_t *dst)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	mvx &= 7;
	mvy &= 7;
	const uint8_t *src_chroma = pred->src_chroma + inter_pred_mvoffset_luma(posx, posy, src_stride);
	if (mvx || mvy) {
		if ((unsigned)posx <= (unsigned)(src_stride - width - 2) && (unsigned)posy <= (unsigned)(vert_stride - height - 1)) {
			if (mvy) {
				if (mvx) {
					filter_chroma_vert_horiz(src_chroma, dst, width, height, mvx, mvy, src_stride, dst_stride);
				} else {
					filter_chroma_vert(src_chroma, dst, width, height, mvy, src_stride, dst_stride);
				}
			} else {
				filter_chroma_horiz(src_chroma, dst, width, height, mvx, src_stride, dst_stride);
			}
		} else {
			/* UMV */
			chroma_filter_info_t f = {
				mvx, mvy
			};
			chroma_inter_umv(src_chroma, dst, posx, posy, width, height, src_stride, vert_stride, dst_stride, &f);
		}

	} else {
		if ((unsigned)posx <= (unsigned)(src_stride - width) && (unsigned)posy <= (unsigned)(vert_stride - height)) {
			copy_inter(src_chroma, dst, width, height, src_stride, dst_stride);
		} else {
			chroma_inter_umv(src_chroma, dst, posx, posy, width, height, src_stride, vert_stride, dst_stride, 0);
		}
	}
}

static inline uint32_t AVERAGE2(uint32_t s1, uint32_t s2)
{
	uint32_t x = s1 ^ s2;
	return (s1 & s2) + ((x & ~0x01010101) >> 1) + (x & 0x01010101);
}

static inline void add_bidir(const uint8_t *src, uint8_t *dst, int width, int height, int stride)
{
	int x_len = (uint32_t)width >> 2;
	do {
		int x = x_len;
		const uint32_t *s = (const uint32_t *)src;
		uint32_t *d = (uint32_t *)dst;
		do {
			*d = AVERAGE2(*s++, *d);
			d++;
		} while (--x);
		src += width;
		dst += stride;
	} while (--height);
}

static void inter_pred_chroma_bidir(const mb_pred_t *pred, int mvx, int mvy, int width, int height, int src_stride, int vert_stride)
{
	uint32_t tmp[(16 * 8) / sizeof(uint32_t)];
	inter_pred_chroma(pred, mvx, mvy, width, height, src_stride, vert_stride, width, (uint8_t *)tmp);
	add_bidir((const uint8_t *)tmp, pred->dst_chroma, width, height, src_stride);
}

static inline void inter_pred_luma_filter02_core(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	width >>= 1;
	do {
		int c0, c1, c2, c3, c4;
		const uint8_t *s = src;
		uint8_t *d = dst;
		int x = width;

		c0 = *s++;
		c1 = *s++;
		c2 = *s++;
		c3 = *s++;
		c4 = *s++;
		do {
			int t, c5, c6;
			c5 = *s++;
			t = (((c2 + c3) * 4 - c1 - c4) * 5 + c0 + c5 + 16) >> 5;
			d[0] = CLIP255C(t);
			c6 = *s++;
			t = (((c3 + c4) * 4 - c2 - c5) * 5 + c1 + c6 + 16) >> 5;
			d[1] = CLIP255C(t);
			c0 = c2;
			c1 = c3;
			c2 = c4;
			c3 = c5;
			c4 = c6;
			d += 2;
		} while (--x);
		src += src_stride;
		dst += stride;
	} while (--height);
}

static inline void inter_pred_luma_filter20_core(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	height >>= 1;
	do {
		int c0, c1, c2, c3, c4;
		const uint8_t *s = src++;
		uint8_t *d = dst++;
		int y = height;

		c0 = *s;
		s += src_stride;
		c1 = *s;
		s += src_stride;
		c2 = *s;
		s += src_stride;
		c3 = *s;
		s += src_stride;
		c4 = *s;
		s += src_stride;
		do {
			int t, c5, c6;
			c5 = *s;
			t = (((c2 + c3) * 4 - c1 - c4) * 5 + c0 + c5 + 16) >> 5;
			s += src_stride;
			d[0] = CLIP255C(t);
			c6 = *s;
			d += stride;
			t = (((c3 + c4) * 4 - c2 - c5) * 5 + c1 + c6 + 16) >> 5;
			s += src_stride;
			d[0] = CLIP255C(t);
			c0 = c2;
			c1 = c3;
			c2 = c4;
			c3 = c5;
			c4 = c6;
			d += stride;
		} while (--y);
	} while (--width);
}

static inline void filter_1_3_v_post(const uint8_t *src0, const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	height >>= 1;
	do {
		int c0, c1, c2, c3, c4;
		const uint8_t *s;
		const uint8_t *s0;
		uint8_t *d;
		int y = height;

		s = src++;
		s0 = src0++;
		d = dst++;
		c0 = *s;
		s += src_stride;
		c1 = *s;
		s += src_stride;
		c2 = *s;
		s += src_stride;
		c3 = *s;
		s += src_stride;
		c4 = *s;
		s += src_stride;
		do {
			int t, c5, c6;
			c5 = *s;
			t = (((c2 + c3) * 4 - c1 - c4) * 5 + c0 + c5 + 16) >> 5;
			*d = (CLIP255I(t) + *s0 + 1) >> 1;
			s += src_stride;
			s0 += stride;
			d += stride;
			c6 = *s;
			t = (((c3 + c4) * 4 - c2 - c5) * 5 + c1 + c6 + 16) >> 5;
			*d = (CLIP255I(t) + *s0 + 1) >> 1;
			s += src_stride;
			s0 += stride;
			d += stride;
			c0 = c2;
			c1 = c3;
			c2 = c4;
			c3 = c5;
			c4 = c6;
		} while (--y);
	} while (--width);
}

template <typename F>
static inline void inter_pred_luma_filter22_horiz(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride, F Pred)
{
	int buf[16 * 21];
	const uint8_t *ss = src;
	int *dd = buf;
	int y = height + 5;

	width >>= 1;
	do {
		int c0, c1, c2, c3, c4;
		const uint8_t *s = ss;
		int x = width;

		c0 = *s++;
		c1 = *s++;
		c2 = *s++;
		c3 = *s++;
		c4 = *s++;
		do {
			int c5, c6;
			c5 = *s++;
			dd[0] = ((c2 + c3) * 4 - c1 - c4) * 5 + c0 + c5;
			c6 = *s++;
			dd[1] = ((c3 + c4) * 4 - c2 - c5) * 5 + c1 + c6;
			c0 = c2;
			c1 = c3;
			c2 = c4;
			c3 = c5;
			c4 = c6;
			dd += 2;
		} while (--x);
		ss += src_stride;
	} while (--y);

	width = width * 2;
	dd = buf;
	y = width;
	height >>= 1;
	do {
		int c0, c1, c2, c3, c4;
		int yy = height;
		int *d = dd++;
		uint8_t *dest = dst++;

		c0 = *d;
		d += width;
		c1 = *d;
		d += width;
		c2 = *d;
		d += width;
		c3 = *d;
		d += width;
		c4 = *d;
		d += width;
		do {
			int c5, c6;
			c5 = *d;
			d += width;
			*dest = Pred(c0, c1, c2, c3, c4, c5);
			c6 = *d;
			dest += stride;
			d += width;
			*dest = Pred(c1, c2, c3, c4, c5, c6);
			c0 = c2;
			c1 = c3;
			c2 = c4;
			c3 = c5;
			c4 = c6;
			dest += stride;
		} while (--yy);
	} while (--y);
}

template <typename F>
static inline void inter_pred_luma_filter22_vert(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride, F Pred)
{
	int buf[21 * 16];
	const uint8_t *ss = src;
	int *dd = buf;
	int i = width + 5;

	height >>= 1;
	do {
		int c0, c1, c2, c3, c4;
		const uint8_t *s = ss++;
		int *d = dd++;
		int y = height;

		c0 = *s;
		s += src_stride;
		c1 = *s;
		s += src_stride;
		c2 = *s;
		s += src_stride;
		c3 = *s;
		s += src_stride;
		c4 = *s;
		s += src_stride;
		do {
			int c5, c6;
			c5 = *s;
			s += src_stride;
			*d = ((c2 + c3) * 4 - c1 - c4) * 5 + c0 + c5;
			d += width + 5;
			c6 = *s;
			*d = ((c3 + c4) * 4 - c2 - c5) * 5 + c1 + c6;
			s += src_stride;
			c0 = c2;
			c1 = c3;
			c2 = c4;
			c3 = c5;
			c4 = c6;
			d += width + 5;
		} while (--y);
	} while (--i);

	height = height * 2;
	dd = buf;
	i = height;
	width >>= 1;
	do {
		int c0, c1, c2, c3, c4;
		int x = width;
		uint8_t *dest = dst;

		c0 = *dd++;
		c1 = *dd++;
		c2 = *dd++;
		c3 = *dd++;
		c4 = *dd++;
		do {
			int c5, c6;
			c5 = *dd++;
			dest[0] = Pred(c0, c1, c2, c3, c4, c5);
			c6 = *dd++;
			dest[1] = Pred(c1, c2, c3, c4, c5, c6);
			c0 = c2;
			c1 = c3;
			c2 = c4;
			c3 = c5;
			c4 = c6;
			dest += 2;
		} while (--x);
		dst += stride;
	} while (--i);
}

struct PPred12 {
	int operator()(int c0, int c1, int c2, int c3, int c4, int c5) const {
		int t = (((c2 + c3) * 4 - c1 - c4) * 5 + c0 + c5 + 512) >> 10;
		int c = (c2 + 16) >> 5;
		return (CLIP255I(t) + CLIP255I(c) + 1) >> 1;
	}
};

struct PPred22 {
	int operator()(int c0, int c1, int c2, int c3, int c4, int c5) const {
		int t = (((c2 + c3) * 4 - c1 - c4) * 5 + c0 + c5 + 512) >> 10;
		return CLIP255C(t);
	}
};

struct PPred32 {
	int operator()(int c0, int c1, int c2, int c3, int c4, int c5) const {
		int t = (((c2 + c3) * 4 - c1 - c4) * 5 + c0 + c5 + 512) >> 10;
		int c = (c3 + 16) >> 5;
		return (CLIP255I(t) + CLIP255I(c) + 1) >> 1;
	}
};

static inline void inter_pred_luma_filter_add(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	stride -= width;
	src_stride -= width;
	width >>= 1;
	do {
		int x = width;
		do {
			dst[0] = (dst[0] + *src++ + 1) >> 1;
			dst[1] = (dst[1] + *src++ + 1) >> 1;
			dst += 2;
		} while (--x);
		src += src_stride;
		dst += stride;
	} while (--height);
}

static void inter_pred_luma_filter00(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	src = src + 2 + src_stride * 2;
	do {
		memcpy(dst, src, width);
		src += src_stride;
		dst += stride;
	} while (--height);
}

static void inter_pred_luma_filter01(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	inter_pred_luma_filter02_core(src + src_stride * 2, dst, width, height, src_stride, stride);
	inter_pred_luma_filter_add(src + src_stride * 2 + 2, dst, width, height, src_stride, stride);
}

static void inter_pred_luma_filter02(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	inter_pred_luma_filter02_core(src + src_stride * 2, dst, width, height, src_stride, stride);
}

static void inter_pred_luma_filter03(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	inter_pred_luma_filter02_core(src + src_stride * 2, dst, width, height, src_stride, stride);
	inter_pred_luma_filter_add(src + src_stride * 2 + 3, dst, width, height, src_stride, stride);
}

static void inter_pred_luma_filter10(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	inter_pred_luma_filter20_core(src + 2, dst, width, height, src_stride, stride);
	inter_pred_luma_filter_add(src + src_stride * 2 + 2, dst, width, height, src_stride, stride);
}

static void inter_pred_luma_filter11(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	inter_pred_luma_filter02_core(src + src_stride * 2, dst, width, height, src_stride, stride);
	filter_1_3_v_post(dst, src + 2, dst, width, height, src_stride, stride);
}

static void inter_pred_luma_filter12(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	inter_pred_luma_filter22_horiz(src, dst, width, height, src_stride, stride, PPred12());
}

static void inter_pred_luma_filter13(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	inter_pred_luma_filter02_core(src + src_stride * 2, dst, width, height, src_stride, stride);
	filter_1_3_v_post(dst, src + 3, dst, width, height, src_stride, stride);
}

static void inter_pred_luma_filter20(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	inter_pred_luma_filter20_core(src + 2, dst, width, height, src_stride, stride);
}

static void inter_pred_luma_filter21(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	inter_pred_luma_filter22_vert(src, dst, width, height, src_stride, stride, PPred12());
}

static inline void inter_pred_luma_filter22(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	inter_pred_luma_filter22_horiz(src, dst, width, height, src_stride, stride, PPred22());
}

static void inter_pred_luma_filter23(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	inter_pred_luma_filter22_vert(src, dst, width, height, src_stride, stride, PPred32());
}

static void inter_pred_luma_filter30(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	inter_pred_luma_filter20_core(src + 2, dst, width, height, src_stride, stride);
	inter_pred_luma_filter_add(src + src_stride * 3 + 2, dst, width, height, src_stride, stride);
}

static void inter_pred_luma_filter31(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	inter_pred_luma_filter02_core(src + src_stride * 3, dst, width, height, src_stride, stride);
	filter_1_3_v_post(dst, src + 2, dst, width, height, src_stride, stride);
}

static void inter_pred_luma_filter32(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	inter_pred_luma_filter22_horiz(src, dst, width, height, src_stride, stride, PPred32());
}

static void inter_pred_luma_filter33(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	inter_pred_luma_filter02_core(src + src_stride * 3, dst, width, height, src_stride, stride);
	filter_1_3_v_post(dst, src + 3, dst, width, height, src_stride, stride);
}

static void (* const inter_pred_luma_filter[4][4])(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int dst_stride) = {
	{
		inter_pred_luma_filter00,
		inter_pred_luma_filter01,
		inter_pred_luma_filter02,
		inter_pred_luma_filter03
	},
	{
		inter_pred_luma_filter10,
		inter_pred_luma_filter11,
		inter_pred_luma_filter12,
		inter_pred_luma_filter13
	},
	{
		inter_pred_luma_filter20,
		inter_pred_luma_filter21,
		inter_pred_luma_filter22,
		inter_pred_luma_filter23
	},
	{
		inter_pred_luma_filter30,
		inter_pred_luma_filter31,
		inter_pred_luma_filter32,
		inter_pred_luma_filter33
	},
};

static inline void extend_left_luma(uint8_t *dst, int left, int width, int height)
{
	do {
		int c = dst[left];
		memset(dst, c, left);
		dst += width;
	} while (--height);
}

static inline void fill_left_top(const uint8_t *src, uint8_t *buf, int left, int top, int width, int height, int stride)
{
	uint8_t *dst;
	int y;

	assert(left <= width);
	src += top * stride + left;
	dst = buf + left;
	for (y = 0; y < top; ++y) {
		memcpy(dst, src, width - left);
		dst += width;
	}
	for (; y < height; ++y) {
		memcpy(dst, src, width - left);
		src += stride;
		dst += width;
	}
	if (left) {
		extend_left_luma(buf, left, width, height);
	}
}

static inline void fill_left_bottom(const uint8_t *src, uint8_t *buf, int left, int bottom, int width, int height, int stride)
{
	uint8_t *dst;
	int y;

	assert(left <= width);
	src += left;
	dst = buf + left;
	for (y = 0; y < height - bottom; ++y) {
		memcpy(dst, src, width - left);
		src += stride;
		dst += width;
	}
	src = dst - width;
	for (; y < height; ++y) {
		memcpy(dst, src, width - left);
		dst += width;
	}
	if (left) {
		extend_left_luma(buf, left, width, height);
	}
}

static inline void extend_right_luma(uint8_t *dst, int right, int width, int height)
{
	dst = dst + width - right;
	do {
		int c = ((int8_t *)dst)[-1];
		memset(dst, c, right);
		dst += width;
	} while (--height);
}

static inline void fill_right_top(const uint8_t *src, uint8_t *buf, int right, int top, int width, int height, int stride)
{
	uint8_t *dst;
	int y;

	assert(right <= width);
	src = src + top * stride;
	dst = buf;
	for (y = 0; y < top; ++y) {
		memcpy(dst, src, width - right);
		dst += width;
	}
	for (; y < height; ++y) {
		memcpy(dst, src, width - right);
		src += stride;
		dst += width;
	}
	if (right) {
		extend_right_luma(buf, right, width, height);
	}
}

static inline void fill_right_bottom(const uint8_t *src, uint8_t *buf, int right, int bottom, int width, int height, int stride)
{
	uint8_t *dst;
	int y;

	assert(right <= width);
	dst = buf;
	for (y = 0; y < height - bottom; ++y) {
		memcpy(dst, src, width - right);
		src += stride;
		dst += width;
	}
	src -= stride;
	for (; y < height; ++y) {
		memcpy(dst, src, width - right);
		dst += width;
	}
	if (right) {
		extend_right_luma(buf, right, width, height);
	}
}

static inline void fill_rect_umv_luma(const uint8_t *src, uint8_t *buf, int width, int height, int stride, int vert_size, int posx, int posy)
{
	int left, right, top, bottom;

	left = -posx;
	top = -posy;
	if (0 < left) {
		if (0 < top) {
			/* left, top */
			fill_left_top(src, buf, left, top, width, height, stride);
		} else {
			bottom = posy - vert_size + height;
			if (0 < bottom) {
				/* left, bottom */
				fill_left_bottom(src, buf, left, bottom, width, height, stride);
			} else {
				/* left */
				fill_left_top(src, buf, left, 0, width, height, stride);
			}
		}
	} else {
		right = posx - stride + width;
		if (0 < top) {
			if (0 < right) {
				/* top, right */
				fill_right_top(src, buf, right, top, width, height, stride);
			} else {
				/* top */
				fill_left_top(src, buf, 0, top, width, height, stride);
			}
		} else {
			bottom = posy - vert_size + height;
			if (0 < right) {
				if (0 < bottom) {
					/* right, bottom */
					fill_right_bottom(src, buf, right, bottom, width, height, stride);
				} else {
					/* right */
					fill_right_top(src, buf, right, 0, width, height, stride);
				}
			} else {
				if (0 < bottom) {
					/* bottom */
					fill_right_bottom(src, buf, 0, bottom, width, height, stride);
				} else {
					/* in range */
					fill_right_top(src, buf, 0, 0, width, height, stride);
				}
			}
		}
	}
}

static void inter_pred_luma_umv(const mb_pred_t *pred, int width, int height, int src_stride, int vert_size, int dst_stride, int fracx, int fracy, uint8_t *dst)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	const uint8_t *src;
	uint8_t buf[22 * 22];

	width += 6;
	height += 6;
	src = pred->src_luma;
	/* rewind if beyond boundary */
	if (posx < 3 - width) {
		src += 3 - width - posx;
		posx = 3 - width;
	} else if (src_stride - 1 < posx - 2) {
		src -= posx - src_stride - 1;
		posx = src_stride + 1;
	}
	if (posy < 3 - height) {
		src += (3 - height - posy) * src_stride;
		posy = 3 - height;
	} else if (vert_size - 1 < posy - 2) {
		src -= (posy - vert_size - 1) * src_stride;
		posy = vert_size + 1;
	}

	posx = posx - 2;
	posy = posy - 2;
	fill_rect_umv_luma(src, buf, width, height, src_stride, vert_size, posx, posy);
	width -= 6;
	height -= 6;
	inter_pred_luma_filter[fracy][fracx](buf, dst, width, height, width + 6, dst_stride);
}

static void inter_pred_luma_frac00(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, int dst_stride, uint8_t *dst)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if ((unsigned)posx <= (unsigned)(stride - width) && (unsigned)posy <= (unsigned)(vert_stride - height)) {
		const uint8_t *src_luma = pred->src_luma + inter_pred_mvoffset_luma(2, 2, stride);
		copy_inter(src_luma, dst, width, height, stride, dst_stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, dst_stride, 0, 0, dst);
	}
}

static void inter_pred_luma_frac01(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, int dst_stride, uint8_t *dst)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && (unsigned)posy < (unsigned)(vert_stride - height)) {
		inter_pred_luma_filter01(pred->src_luma, dst, width, height, stride, dst_stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, dst_stride, 1, 0, dst);
	}
}

static void inter_pred_luma_frac02(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, int dst_stride, uint8_t *dst) 
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && (unsigned)posy < (unsigned)(vert_stride - height)) {
		inter_pred_luma_filter02(pred->src_luma, dst, width, height, stride, dst_stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, dst_stride, 2, 0, dst);
	}
}

static void inter_pred_luma_frac03(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, int dst_stride, uint8_t *dst)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && (unsigned)posy < (unsigned)(vert_stride - height)) {
		inter_pred_luma_filter03(pred->src_luma, dst, width, height, stride, dst_stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, dst_stride, 3, 0, dst);
	}
}

static void inter_pred_luma_frac10(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, int dst_stride, uint8_t *dst)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if ((unsigned)posx < (unsigned)(stride - width) && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter10(pred->src_luma, dst, width, height, stride, dst_stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, dst_stride, 0, 1, dst);
	}
}

static void inter_pred_luma_frac11(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, int dst_stride, uint8_t *dst)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter11(pred->src_luma, dst, width, height, stride, dst_stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, dst_stride, 1, 1, dst);
	}
}

static void inter_pred_luma_frac12(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, int dst_stride, uint8_t *dst)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter12(pred->src_luma, dst, width, height, stride, dst_stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, dst_stride, 2, 1, dst);
	}
}

static void inter_pred_luma_frac13(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, int dst_stride, uint8_t *dst)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter13(pred->src_luma, dst, width, height, stride, dst_stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, dst_stride, 3, 1, dst);
	}
}

static void inter_pred_luma_frac20(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, int dst_stride, uint8_t *dst)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if ((unsigned)posx < (unsigned)(stride - width) && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter20(pred->src_luma, dst, width, height, stride, dst_stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, dst_stride, 0, 2, dst);
	}
}

static void inter_pred_luma_frac21(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, int dst_stride, uint8_t *dst)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter21(pred->src_luma, dst, width, height, stride, dst_stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, dst_stride, 1, 2, dst);
	}
}

static void inter_pred_luma_frac22(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, int dst_stride, uint8_t *dst)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter22(pred->src_luma, dst, width, height, stride, dst_stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, dst_stride, 2, 2, dst);
	}
}

static void inter_pred_luma_frac23(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, int dst_stride, uint8_t *dst)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter23(pred->src_luma, dst, width, height, stride, dst_stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, dst_stride, 3, 2, dst);
	}
}

static void inter_pred_luma_frac30(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, int dst_stride, uint8_t *dst)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if ((unsigned)posx < (unsigned)(stride - width) && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter30(pred->src_luma, dst, width, height, stride, dst_stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, dst_stride, 0, 3, dst);
	}
}

static void inter_pred_luma_frac31(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, int dst_stride, uint8_t *dst)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter31(pred->src_luma, dst, width, height, stride, dst_stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, dst_stride, 1, 3, dst);
	}
}

static void inter_pred_luma_frac32(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, int dst_stride, uint8_t *dst)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter32(pred->src_luma, dst, width, height, stride, dst_stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, dst_stride, 2, 3, dst);
	}
}

static void inter_pred_luma_frac33(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, int dst_stride, uint8_t *dst)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter33(pred->src_luma, dst, width, height, stride, dst_stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, dst_stride, 3, 3, dst);
	}
}

template <typename F>
static inline void inter_pred_luma_bidir_latter(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, uint8_t *dst,
						F InterPred)
{
	uint32_t tmp[(16 * 16) / sizeof(uint32_t)];

	InterPred(pred, width, height, stride, vert_stride, width, (uint8_t *)tmp);
	add_bidir((const uint8_t *)tmp, dst, width, height, stride);
}

static void inter_pred_luma_bidir_latter00(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, uint8_t *dst)
{
	inter_pred_luma_bidir_latter(pred, width, height, stride, vert_stride, dst, inter_pred_luma_frac00);
}

static void inter_pred_luma_bidir_latter01(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, uint8_t *dst)
{
	inter_pred_luma_bidir_latter(pred, width, height, stride, vert_stride, dst, inter_pred_luma_frac01);
}

static void inter_pred_luma_bidir_latter02(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, uint8_t *dst)
{
	inter_pred_luma_bidir_latter(pred, width, height, stride, vert_stride, dst, inter_pred_luma_frac02);
}

static void inter_pred_luma_bidir_latter03(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, uint8_t *dst)
{
	inter_pred_luma_bidir_latter(pred, width, height, stride, vert_stride, dst, inter_pred_luma_frac03);
}

static void inter_pred_luma_bidir_latter10(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, uint8_t *dst)
{
	inter_pred_luma_bidir_latter(pred, width, height, stride, vert_stride, dst, inter_pred_luma_frac10);
}

static void inter_pred_luma_bidir_latter11(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, uint8_t *dst)
{
	inter_pred_luma_bidir_latter(pred, width, height, stride, vert_stride, dst, inter_pred_luma_frac11);
}

static void inter_pred_luma_bidir_latter12(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, uint8_t *dst)
{
	inter_pred_luma_bidir_latter(pred, width, height, stride, vert_stride, dst, inter_pred_luma_frac12);
}

static void inter_pred_luma_bidir_latter13(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, uint8_t *dst)
{
	inter_pred_luma_bidir_latter(pred, width, height, stride, vert_stride, dst, inter_pred_luma_frac13);
}

static void inter_pred_luma_bidir_latter20(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, uint8_t *dst)
{
	inter_pred_luma_bidir_latter(pred, width, height, stride, vert_stride, dst, inter_pred_luma_frac20);
}

static void inter_pred_luma_bidir_latter21(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, uint8_t *dst)
{
	inter_pred_luma_bidir_latter(pred, width, height, stride, vert_stride, dst, inter_pred_luma_frac21);
}

static void inter_pred_luma_bidir_latter22(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, uint8_t *dst)
{
	inter_pred_luma_bidir_latter(pred, width, height, stride, vert_stride, dst, inter_pred_luma_frac22);
}

static void inter_pred_luma_bidir_latter23(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, uint8_t *dst)
{
	inter_pred_luma_bidir_latter(pred, width, height, stride, vert_stride, dst, inter_pred_luma_frac23);
}

static void inter_pred_luma_bidir_latter30(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, uint8_t *dst)
{
	inter_pred_luma_bidir_latter(pred, width, height, stride, vert_stride, dst, inter_pred_luma_frac30);
}

static void inter_pred_luma_bidir_latter31(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, uint8_t *dst)
{
	inter_pred_luma_bidir_latter(pred, width, height, stride, vert_stride, dst, inter_pred_luma_frac31);
}

static void inter_pred_luma_bidir_latter32(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, uint8_t *dst)
{
	inter_pred_luma_bidir_latter(pred, width, height, stride, vert_stride, dst, inter_pred_luma_frac32);
}

static void inter_pred_luma_bidir_latter33(const mb_pred_t *pred, int width, int height, int stride, int vert_stride, uint8_t *dst)
{
	inter_pred_luma_bidir_latter(pred, width, height, stride, vert_stride, dst, inter_pred_luma_frac33);
}

static void (* const inter_pred_luma[4][4])(const mb_pred_t *pred, int width, int height, int src_stride, int vert_stride, int dst_stride, uint8_t *dst) = {
	{
		inter_pred_luma_frac00,
		inter_pred_luma_frac01,
		inter_pred_luma_frac02,
		inter_pred_luma_frac03
	},
	{
		inter_pred_luma_frac10,
		inter_pred_luma_frac11,
		inter_pred_luma_frac12,
		inter_pred_luma_frac13
	},
	{
		inter_pred_luma_frac20,
		inter_pred_luma_frac21,
		inter_pred_luma_frac22,
		inter_pred_luma_frac23
	},
	{
		inter_pred_luma_frac30,
		inter_pred_luma_frac31,
		inter_pred_luma_frac32,
		inter_pred_luma_frac33
	}
};

static void (* const inter_pred_luma_bidir[4][4])(const mb_pred_t *pred, int width, int height, int src_stride, int vert_stride, uint8_t *dst) = {
	{
		inter_pred_luma_bidir_latter00,
		inter_pred_luma_bidir_latter01,
		inter_pred_luma_bidir_latter02,
		inter_pred_luma_bidir_latter03
	},
	{
		inter_pred_luma_bidir_latter10,
		inter_pred_luma_bidir_latter11,
		inter_pred_luma_bidir_latter12,
		inter_pred_luma_bidir_latter13
	},
	{
		inter_pred_luma_bidir_latter20,
		inter_pred_luma_bidir_latter21,
		inter_pred_luma_bidir_latter22,
		inter_pred_luma_bidir_latter23
	},
	{
		inter_pred_luma_bidir_latter30,
		inter_pred_luma_bidir_latter31,
		inter_pred_luma_bidir_latter32,
		inter_pred_luma_bidir_latter33
	}
};

template <typename F0, typename F1>
static uint32_t residual_luma_inter(h264d_mb_current *mb, uint32_t cbp, dec_bits *st, int avail,
				    F0 QpDelta,
				    F1 ResidualBlock)
{
	int coeff[16];
	int32_t qp_delta;
	const int16_t *qmat;
	const int *offset;
	uint32_t top, left;
	int c0, c1, c2, c3, c4, c5;
	uint32_t str_map; /* deblocking filter strength 3, vertical */
	uint8_t *luma;
	int stride;

	qp_delta = QpDelta(mb, st, avail);
	if (qp_delta) {
		set_qp(mb, mb->qp + qp_delta);
	}
	qmat = mb->qmaty;
	offset = mb->offset4x4;
	luma = mb->luma;
	stride = mb->max_x * 16;
	str_map = 0;
	if (cbp & 1) {
		if ((c0 = ResidualBlock(mb, avail & 1 ? UNPACK(mb->left4x4coef, 0) : -1, avail & 2 ? UNPACK(*mb->top4x4coef, 0) : -1, st, coeff, 16, qmat, avail, 0, 2, 0xf)) != 0) {
			ac4x4transform_acdc(luma + offset[0], coeff, stride);
			str_map = 0x2;
		}
		if ((c1 = ResidualBlock(mb, c0, avail & 2 ? UNPACK(*mb->top4x4coef, 1) : -1, st, coeff, 16, qmat, avail, 1, 2, 0xf)) != 0) {
			ac4x4transform_acdc(luma + offset[1], coeff, stride);
			str_map |= 0x8;
		}
		if ((c2 = ResidualBlock(mb, avail & 1 ? UNPACK(mb->left4x4coef, 1) : -1, c0, st, coeff, 16, qmat, avail, 2, 2, 0xf)) != 0) {
			ac4x4transform_acdc(luma + offset[2], coeff, stride);
			str_map |= 0x200;
		}
		if ((c3 = ResidualBlock(mb, c2, c1, st, coeff, 16, qmat, avail, 3, 2, 0xf)) != 0) {
			ac4x4transform_acdc(luma + offset[3], coeff, stride);
			str_map |= 0x800;
		}
	} else {
		c0 = 0;
		c1 = 0;
		c2 = 0;
		c3 = 0;
	}
	if (cbp & 2) {
		if ((c0 = ResidualBlock(mb, c1, avail & 2 ? UNPACK(*mb->top4x4coef, 2) : -1, st, coeff, 16, qmat, avail, 4, 2, 0xf)) != 0) {
			ac4x4transform_acdc(luma + offset[4], coeff, stride);
			str_map |= 0x20;
		}
		if ((c1 = ResidualBlock(mb, c0, avail & 2 ? UNPACK(*mb->top4x4coef, 3) : -1, st, coeff, 16, qmat, avail, 5, 2, 0xf)) != 0) {
			left = PACK(0, c1, 0);
			str_map |= 0x80;
			ac4x4transform_acdc(luma + offset[5], coeff, stride);
		} else {
			left = 0;
		}
		if ((c4 = ResidualBlock(mb, c3, c0, st, coeff, 16, qmat, avail, 6, 2, 0xf)) != 0) {
			ac4x4transform_acdc(luma + offset[6], coeff, stride);
			str_map |= 0x2000;
		}
		if ((c5 = ResidualBlock(mb, c4, c1, st, coeff, 16, qmat, avail, 7, 2, 0xf)) != 0) {
			left = PACK(left, c5, 1);
			str_map |= 0x8000;
			ac4x4transform_acdc(luma + offset[7], coeff, stride);
		}
	} else {
		c0 = 0;
		c1 = 0;
		c4 = 0;
		c5 = 0;
		left = 0;
	}
	if (cbp & 4) {
		if ((c0 = ResidualBlock(mb, avail & 1 ? UNPACK(mb->left4x4coef, 2) : -1, c2, st, coeff, 16, qmat, avail, 8, 2, 0xf)) != 0) {
			ac4x4transform_acdc(luma + offset[8], coeff, stride);
			str_map |= 0x20000;
		}
		if ((c1 = ResidualBlock(mb, c0, c3, st, coeff, 16, qmat, avail, 9, 2, 0xf)) != 0) {
			ac4x4transform_acdc(luma + offset[9], coeff, stride);
			str_map |= 0x80000;
		}
		if ((c2 = ResidualBlock(mb, avail & 1 ? UNPACK(mb->left4x4coef, 3) : -1, c0, st, coeff, 16, qmat, avail, 10, 2, 0xf)) != 0) {
			top = PACK(0, c2, 0);
			str_map |= 0x2000000;
			ac4x4transform_acdc(luma + offset[10], coeff, stride);
		} else {
			top = 0;
		}
		if ((c3 = ResidualBlock(mb, c2, c1, st, coeff, 16, qmat, avail, 11, 2, 0xf)) != 0) {
			top = PACK(top, c3, 1);
			str_map |= 0x8000000;
			ac4x4transform_acdc(luma + offset[11], coeff, stride);
		}
	} else {
		c0 = 0;
		c1 = 0;
		c2 = 0;
		c3 = 0;
		top = 0;
	}
	if (cbp & 8) {
		if ((c0 = ResidualBlock(mb, c1, c4, st, coeff, 16, qmat, avail, 12, 2, 0xf)) != 0) {
			ac4x4transform_acdc(luma + offset[12], coeff, stride);
			str_map |= 0x200000;
		}
		if ((c1 = ResidualBlock(mb, c0, c5, st, coeff, 16, qmat, avail, 13, 2, 0xf)) != 0) {
			left = PACK(left, c1, 2);
			str_map |= 0x800000;
			ac4x4transform_acdc(luma + offset[13], coeff, stride);
		}
		if ((c2 = ResidualBlock(mb, c3, c0, st, coeff, 16, qmat, avail, 14, 2, 0xf)) != 0) {
			top = PACK(top, c2, 2);
			str_map |= 0x20000000;
			ac4x4transform_acdc(luma + offset[14], coeff, stride);
		}
		if ((c3 = ResidualBlock(mb, c2, c1, st, coeff, 16, qmat, avail, 15, 2, 0xf)) != 0) {
			str_map |= 0x80000000;
			ac4x4transform_acdc(luma + offset[15], coeff, stride);
		}
	} else {
		c3 = 0; 
	}
	mb->left4x4coef = (mb->left4x4coef & 0xffff0000) | PACK(left, c3, 3);
	*mb->top4x4coef = (*mb->top4x4coef & 0xffff0000) | PACK(top, c3, 3);
	return str_map;
}

static inline int MEDIAN(int a, int b, int c)
{
	return (a <= b) ? ((b <= c) ? b : (a <= c ? c : a)) : (a <= c ? a : (b <= c ? c : b));
}

static const int16_t zero_mv[2 * 8] = {
	0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0
};

static const int8_t non_ref[4] = {
	-1, -1, -1, -1
};

static const h264d_vector_t zero_mov[2] = {
	{0U}, {0U}
};

static inline void determine_pmv(const int16_t *mva, const int16_t *mvb, const int16_t *mvc, int16_t pmv[], int avail, int idx_map)
{
	static const int not_one_hot = 0xe9;
	int pmvx, pmvy;
	if (((avail & 7) == 1) || (idx_map == 1)) {
		pmvx = mva[0];
		pmvy = mva[1];
	} else if (not_one_hot & (1 << idx_map)) {
		pmvx = MEDIAN(mva[0], mvb[0], mvc[0]);
		pmvy = MEDIAN(mva[1], mvb[1], mvc[1]);
	} else if (idx_map == 2) {
		pmvx = mvb[0];
		pmvy = mvb[1];
	} else {
		pmvx = mvc[0];
		pmvy = mvc[1];
	}
	pmv[0] = pmvx;
	pmv[1] = pmvy;
}


static inline void calc_mv16x16(h264d_mb_current *mb, int16_t *pmv, const int16_t *&mvd_a, const int16_t *&mvd_b, int lx, int ref_idx, int avail)
{
	prev_mb_t *pmb;
	const int16_t *mva, *mvb, *mvc;
	int idx_map;

	if (avail & 1) {
		pmb = mb->left4x4inter;
		idx_map = (ref_idx == pmb->ref[0][lx]);
		mva = pmb->mov[0].mv[lx].v;
		mvd_a = pmb->mvd[0].mv[lx].v;
	} else {
		idx_map = 0;
		mvd_a = mva = zero_mv;
	}
	if (avail & 2) {
		pmb = mb->top4x4inter;
		idx_map |= (ref_idx == pmb->ref[0][lx]) * 2;
		mvb = pmb->mov[0].mv[lx].v;
		mvd_b = pmb->mvd[0].mv[lx].v;
	} else {
		mvd_b = mvb = zero_mv;
	}
	if (avail & 4) {
		pmb = mb->top4x4inter + 1;
		idx_map |= (ref_idx == pmb->ref[0][lx]) * 4;
		mvc = pmb->mov[0].mv[lx].v;
	} else if (avail & 8) {
		idx_map |= (ref_idx == mb->lefttop_ref[lx]) * 4;
		mvc = mb->lefttop_mv[lx].v;
	} else {
		mvc = zero_mv;
	}
	determine_pmv(mva, mvb, mvc, pmv, avail, idx_map);
}

static inline void inter_pred8x8(h264d_mb_current *mb, const int16_t mv[], int width, int height, int frame_idx, int offsetx, int offsety, int bidir)
{
	m2d_frame_t *frms = &(mb->frame->frames[frame_idx]);
	int stride, vert_size;
	int posx, posy;
	int mvx, mvy;
	mb_pred_t pred;

	stride = mb->max_x * 16;
	vert_size = mb->max_y * 16;
	pred.src_luma = frms->luma;
	pred.src_chroma = frms->chroma;
	pred.dst_chroma = mb->chroma + (offsety >> 1) * stride + offsetx;
	mvx = mv[0];
	mvy = mv[1];
	pred.pos_x = posx = mb->x * 16 + (mvx >> 2) + offsetx;
	pred.pos_y = posy = mb->y * 16 + (mvy >> 2) + offsety;
	pred.src_luma = pred.src_luma + inter_pred_mvoffset_luma(posx - 2, posy - 2, stride);
	if (bidir == 0) {
		inter_pred_luma[mvy & 3][mvx & 3](&pred, width, height, stride, vert_size, stride, mb->luma + offsety * stride + offsetx);
	} else {
		inter_pred_luma_bidir[mvy & 3][mvx & 3](&pred, width, height, stride, vert_size, mb->luma + offsety * stride + offsetx);
	}
	pred.pos_x = mb->x * 16 + (mvx >> 3) * 2 + offsetx;
	pred.pos_y = mb->y * 8 + (mvy >> 3) + (offsety >> 1);
	height >>= 1;
	vert_size >>= 1;
	if (bidir == 0) {
		inter_pred_chroma(&pred, mvx, mvy, width, height, stride, vert_size, stride, pred.dst_chroma);
	} else {
		inter_pred_chroma_bidir(&pred, mvx, mvy, width, height, stride, vert_size);
	}
}

static inline void direct_zero_pred_block(const uint8_t *src0, const uint8_t *src1, uint8_t *dst, int width, int height, int stride)
{
	int x_len = (uint32_t)width >> 2;
	do {
		int x = x_len;
		const uint32_t *s0 = (const uint32_t *)src0;
		const uint32_t *s1 = (const uint32_t *)src1;
		uint32_t *d = (uint32_t *)dst;
		do {
			*d++ = AVERAGE2(*s0++, *s1++);
		} while (--x);
		src0 += stride;
		src1 += stride;
		dst += stride;
	} while (--height);
}

static inline void direct_zero_pred(h264d_mb_current *mb, int width, int height, int xoffset, int yoffset)
{
	m2d_frame_t *frms0 = &(mb->frame->frames[mb->frame->refs[0][0].frame_idx]);
	m2d_frame_t *frms1 = &(mb->frame->frames[mb->frame->refs[1][0].frame_idx]);
	int stride = mb->max_x * 16;
	int offset;

	offset = (mb->y * 16 + yoffset) * stride + mb->x * 16 + xoffset;
	direct_zero_pred_block(frms0->luma + offset, frms1->luma + offset, mb->luma + yoffset * stride + xoffset, width, height, stride);
	offset = (((mb->y * 16 + yoffset) * stride) >> 1) + mb->x * 16 + xoffset;
	direct_zero_pred_block(frms0->chroma + offset, frms1->chroma + offset, mb->chroma + ((yoffset * stride) >> 1) + xoffset, width, height >> 1, stride);
}

static inline uint32_t spread(uint32_t a)
{
	return ((a & 0xc0) << 18) | ((a & 0x30) << 12) | ((a & 0xc) << 6) | (a & 3);
}

static inline uint32_t transposition(uint32_t a)
{
	return (spread((a >> 24) & 0xff) << 6) | (spread((a >> 16) & 0xff) << 4) | (spread((a >> 8) & 0xff) << 2) | spread(a & 0xff);
}

static inline uint32_t expand_coef_str(uint32_t str) {
	return str | (str << 8);
}

static inline uint32_t str_previous_coef(uint32_t map, uint32_t prev4x4)
{
	if (prev4x4) {
		for (int i = 0; i < 4; ++i) {
			if ((prev4x4 & 0xf) != 0) {
				map |= 2 << (i * 2);
			}
			prev4x4 >>= 4;
		}
	}
	return map;
}

static inline int DIF_SQUARE(int a, int b) {
	int t = a - b;
	return t * t;
}

static inline bool DIF_ABS_LARGER_THAN4(int a, int b) {
	return 16 <= DIF_SQUARE(a, b);
}

static inline bool ABS_LARGER_THAN4(uint32_t x)
{
	return 16 <= x * x;
}

template <int MV_STEP>
static inline uint32_t str_mv_calc16x16_bidir_both(uint32_t str, int offset, const h264d_vector_set_t *mvxy, const prev_mb_t *prev)
{
	uint32_t mask = 2 << (offset * 2);
	for (int j = 0; j < 2; ++j) {
		if (!(str & mask)) {
			const int16_t *mv_prev = &(prev->mov[j + offset].mv[0].v[0]);
			int prev0x = *mv_prev++;
			int prev0y = *mv_prev++;
			int prev1x = *mv_prev++;
			int prev1y = *mv_prev;
			if ((DIF_ABS_LARGER_THAN4(mvxy->mv[0].v[0], prev0x)
				|| DIF_ABS_LARGER_THAN4(mvxy->mv[0].v[1], prev0y)
				|| DIF_ABS_LARGER_THAN4(mvxy->mv[1].v[0], prev1x)
				|| DIF_ABS_LARGER_THAN4(mvxy->mv[1].v[1], prev1y))
				&& (DIF_ABS_LARGER_THAN4(mvxy->mv[0].v[0], prev1x)
				|| DIF_ABS_LARGER_THAN4(mvxy->mv[0].v[1], prev1y)
				|| DIF_ABS_LARGER_THAN4(mvxy->mv[1].v[0], prev0x)
				|| DIF_ABS_LARGER_THAN4(mvxy->mv[1].v[1], prev0y))) {
				str = str | (mask >> 1);
			}
		}
		mask <<= 2;
		mvxy += MV_STEP;
	}
	return str;
}

template <int MV_STEP>
static inline uint32_t str_mv_calc16x16_bidir_one(uint32_t str, int ref0, int prev_ref0, int offset, const h264d_vector_set_t *mvxy, const prev_mb_t *prev)
{
	int lx0 = (ref0 != prev_ref0);
	int lx1 = lx0 ^ 1;
	uint32_t mask = 2 << (offset * 2);
	for (int j = 0; j < 2; ++j) {
		if (!(str & mask)) {
			const int16_t *mv_prev = &(prev->mov[j + offset].mv[0].v[0]);
			if (DIF_ABS_LARGER_THAN4(mvxy->mv[lx0].v[0], *mv_prev++)
				|| DIF_ABS_LARGER_THAN4(mvxy->mv[lx0].v[1], *mv_prev++)
				|| DIF_ABS_LARGER_THAN4(mvxy->mv[lx1].v[0], *mv_prev++)
				|| DIF_ABS_LARGER_THAN4(mvxy->mv[lx1].v[1], *mv_prev)) {
				str = str | (mask >> 1);
			}
		}
		mask <<= 2;
		mvxy += MV_STEP;
	}
	return str;
}

template <int MV_STEP>
static inline uint32_t str_mv_calc16x16_bidir(uint32_t str, int ref0, int ref1, int prev_ref0, int offset, const h264d_vector_set_t *mvxy, const prev_mb_t *prev)
{
	if (ref0 == ref1) {
		return str_mv_calc16x16_bidir_both<MV_STEP>(str, offset, mvxy, prev);
	} else {
		return str_mv_calc16x16_bidir_one<MV_STEP>(str, ref0, prev_ref0, offset, mvxy, prev);
	}
}

template <int MV_STEP>
static inline uint32_t str_mv_calc16x16_onedir(uint32_t str, int ref0, int ref1, int prev_ref0, int offset, const h264d_vector_set_t *mvxy, const prev_mb_t *prev)
{
	int lx_curr, lx_prev;
	if (0 <= ref0) {
		lx_curr = 0;
		lx_prev = (ref0 != prev_ref0);
	} else {
		lx_curr = 1;
		lx_prev = (ref1 != prev_ref0);
	}
	uint32_t mask = 2 << (offset * 2);
	for (int j = 0; j < 2; ++j) {
		if ((str & mask) == 0) {
			if (DIF_ABS_LARGER_THAN4(mvxy->mv[lx_curr].v[0], prev->mov[j + offset].mv[lx_prev].v[0])
				|| DIF_ABS_LARGER_THAN4(mvxy->mv[lx_curr].v[1], prev->mov[j + offset].mv[lx_prev].v[1])) {
				str |= mask >> 1;
			}
		}
		mask <<= 2;
		mvxy += MV_STEP;
	}
	return str;
}

static inline int frame_idx_of_ref(const h264d_mb_current *mb, int ref_idx, int lx) {
	return (0 <= ref_idx) ? mb->frame->refs[lx][ref_idx].frame_idx : -1;
}

template <int MV_STEP>
static inline uint32_t str_mv_calc16x16_mv(uint32_t str, int ref0, int ref1, int prev_ref0, int offset, const h264d_vector_set_t *mvxy, const prev_mb_t *prev)
{
	if ((0 <= ref0) && (0 <= ref1)) {
		return str_mv_calc16x16_bidir<MV_STEP>(str, ref0, ref1, prev_ref0, offset, mvxy, prev);
	} else {
		return str_mv_calc16x16_onedir<MV_STEP>(str, ref0, ref1, prev_ref0, offset, mvxy, prev);
	}
}

static inline uint32_t str_mv_calc16x16(const h264d_mb_current *mb, uint32_t str, const h264d_vector_set_t *mvxy, const int8_t ref_idx[], const prev_mb_t *prev)
{
	int ref0 = frame_idx_of_ref(mb, *ref_idx++, 0);
	int ref1 = frame_idx_of_ref(mb, *ref_idx, 1);
	uint32_t mask = 0xa;
	for (int i = 0; i < 2; ++i) {
		if ((str & mask) != mask) {
			int prev0 = frame_idx_of_ref(mb, prev->ref[i][0], 0);
			int prev1 = frame_idx_of_ref(mb, prev->ref[i][1], 1);
			if (((prev0 != ref0) || (prev1 != ref1)) && ((prev1 != ref0) || (prev0 != ref1))) {
				uint32_t m = mask >> 1;
				str = str | (((str >> 1) ^ m) & m);
			} else {
				str = str_mv_calc16x16_mv<0>(str, ref0, ref1, prev0, i * 2, mvxy, prev);
			}
		}
		mask <<= 4;
	}
	return str;
}

static uint32_t store_str_inter16xedge(const h264d_mb_current *mb, const prev_mb_t *prev4x4inter, int8_t& str4, const h264d_vector_set_t *mv, const int8_t ref_idx[], uint32_t str, uint32_t coeff4x4)
{
	if (prev4x4inter->type <= MB_IPCM) {
		str4 = 1;
		str |= 0xaa;
	} else {
		str = str_previous_coef(str, coeff4x4);
		str = str_mv_calc16x16(mb, str, mv, ref_idx, prev4x4inter);
	}
	return str;
}

static void store_info_inter16x16(h264d_mb_current *mb, const h264d_vector_set_t mv[], const int8_t ref_idx[], uint32_t str_vert, uint32_t str_horiz, uint32_t left4x4, uint32_t top4x4)
{
	deblock_info_t *deb = mb->deblock_curr;
	deb->qpy = mb->qp;
	deb->qpc = mb->qp_chroma;
	if (mb->y != 0) {
		str_vert = store_str_inter16xedge(mb, mb->top4x4inter, deb->str4_vert, mv, ref_idx, str_vert, top4x4);
	}
	deb->str_vert = str_vert;
	if (mb->x != 0) {
		str_horiz = store_str_inter16xedge(mb, mb->left4x4inter, deb->str4_horiz, mv, ref_idx, str_horiz, left4x4);
	}
	deb->str_horiz = str_horiz;
	*mb->top4x4pred = 0x22222222;
	mb->left4x4pred = 0x22222222;

	mb->left4x4inter->direct8x8 = 0;
	mb->top4x4inter->direct8x8 = 0;
	for (int i = 0; i < 2; ++i) {
		mb->lefttop_ref[i] = mb->top4x4inter->ref[1][i];
		mb->lefttop_mv[i].vector = mb->top4x4inter->mov[3].mv[i].vector;
		for (int j = 0; j < 2; ++j) {
			mb->top4x4inter->ref[j][i] = ref_idx[i];
			mb->left4x4inter->ref[j][i] = ref_idx[i];
		}
	}
	for (int i = 0; i < 4; ++i) {
		mb->left4x4inter->mov[i] = mv[0];
		mb->left4x4inter->mvd[i] = mv[1];
		mb->top4x4inter->mov[i] = mv[0];
		mb->top4x4inter->mvd[i] = mv[1];
	}
	int refcol;
	uint32_t mvcol;
	if (0 <= ref_idx[0]) {
		refcol = ref_idx[0];
		mvcol = mv->mv[0].vector;
	} else {
		refcol = ref_idx[1];
		mvcol = mv->mv[1].vector;
	}
	mb->col_curr->type = COL_MB16x16;
	memset(mb->col_curr->ref, refcol, sizeof(mb->col_curr->ref));
	h264d_vector_t *mvdst = mb->col_curr->mv;
	for (int i = 0; i < 16; ++i) {
		mvdst[i].vector = mvcol;
	}
}

template <typename F0 ,typename F1, typename F2, typename F3, typename F4>
static int mb_inter16x16(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail,
			 F0 RefIdx16x16,
			 F1 MvdXY,
			 F2 CodedBlockPattern,
			 F3 QpDelta,
			 F4 ResidualBlock)
{
	const int16_t *mvd_a;
	const int16_t *mvd_b;
	h264d_vector_set_t mv[2];
	uint32_t predmap;
	uint32_t cbp;
	uint32_t str_vert, str_horiz;
	int bidir;
	int8_t ref_idx[2];

	predmap = mbc->cbp;
	for (int lx = 0; lx < 2; ++lx) {
		ref_idx[lx] = (predmap & (1 << lx)) ? RefIdx16x16(mb, st, lx, avail) : -1;
	}
	memset(&mv, 0, sizeof(mv));
	bidir = 0;
	for (int lx = 0; lx < 2; ++lx) {
		if (predmap & (1 << lx)) {
			calc_mv16x16(mb, mv[0].mv[lx].v, mvd_a, mvd_b, lx, ref_idx[lx], avail);
			MvdXY(mb, st, mv[1].mv[lx].v, mvd_a, mvd_b);
			mv[0].mv[lx].v[0] += mv[1].mv[lx].v[0];
			mv[0].mv[lx].v[1] += mv[1].mv[lx].v[1];
			inter_pred8x8(mb, mv[0].mv[lx].v, 16, 16, mb->frame->refs[lx][ref_idx[lx]].frame_idx, 0, 0, bidir++);
		}
	}
	uint32_t left4x4 = mb->left4x4coef;
	uint32_t top4x4 = *mb->top4x4coef;
	mb->cbp = cbp = CodedBlockPattern(mb, st, avail);
	if (cbp) {
		str_vert = residual_luma_inter(mb, cbp, st, avail, QpDelta, ResidualBlock);
		str_horiz = expand_coef_str(transposition(str_vert));
		str_vert = expand_coef_str(str_vert);
	} else {
		mb->prev_qp_delta = 0;
		mb->left4x4coef = 0;
		*mb->top4x4coef = 0;
		str_vert = 0;
		str_horiz = 0;
	}
	store_info_inter16x16(mb, &mv[0], ref_idx, str_vert, str_horiz, left4x4, top4x4);
	return residual_chroma(mb, cbp, st, avail, ResidualBlock);
}

static inline void calc_mv16x8top(h264d_mb_current *mb, int16_t pmv[], const int16_t *&mvd_a, const int16_t *&mvd_b, int lx, int ref_idx, int avail)
{
	prev_mb_t *pmb;
	const int16_t *mva, *mvb, *mvc;
	int idx_map;

	if (avail & 2) {
		pmb = mb->top4x4inter;
		mvd_b = pmb->mvd[0].mv[lx].v;
		if (ref_idx == pmb->ref[0][lx]) {
			pmv[0] = pmb->mov[0].mv[lx].v[0];
			pmv[1] = pmb->mov[0].mv[lx].v[1];
			mvd_a = (avail & 1) ? mb->left4x4inter->mvd[0].mv[lx].v : zero_mv;
			return;
		}
		mvb = pmb->mov[0].mv[lx].v;
	} else {
		mvd_b = mvb = zero_mv;
	}
	if (avail & 1) {
		pmb = mb->left4x4inter;
		idx_map = (ref_idx == pmb->ref[0][lx]);
		mva = pmb->mov[0].mv[lx].v;
		mvd_a = pmb->mvd[0].mv[lx].v;
	} else {
		mvd_a = mva = zero_mv;
		idx_map = 0;
	}
	if (avail & 4) {
		pmb = mb->top4x4inter + 1;
		idx_map |= (ref_idx == pmb->ref[0][lx]) * 4;
		mvc = pmb->mov[0].mv[lx].v;
	} else if (avail & 8) {
		idx_map |= (ref_idx == mb->lefttop_ref[lx]) * 4;
		mvc = mb->lefttop_mv[lx].v;
	} else {
		mvc = zero_mv;
	}
	determine_pmv(mva, mvb, mvc, pmv, avail, idx_map);
}

static inline void calc_mv16x8bottom(h264d_mb_current *mb, int16_t mv[], const int16_t *&mvd_a, const int16_t *&mvd_b, int lx, int ref_idx, int avail, int prev_ref, const h264d_vector_set_t *prev_mv)
{
	prev_mb_t *pmb;
	const int16_t *mva, *mvb, *mvc;
	int idx_map;

	if (avail & 1) {
		pmb = mb->left4x4inter;
		mvd_a = pmb->mvd[2].mv[lx].v;
		if (ref_idx == pmb->ref[1][lx]) {
			mv[0] = pmb->mov[2].mv[lx].v[0];
			mv[1] = pmb->mov[2].mv[lx].v[1];
			mvd_b = prev_mv[2].mv[lx].v;
			return;
		}
		idx_map = (ref_idx == pmb->ref[0][lx]) * 4;
		mva = pmb->mov[2].mv[lx].v;
		mvc = pmb->mov[1].mv[lx].v;
	} else {
		idx_map = 0;
		mvd_a = mva = zero_mv;
		mvc = zero_mv;
	}
	/* upper block: always exists */
	mvb = prev_mv->mv[lx].v;
	mvd_b = prev_mv[2].mv[lx].v;
	idx_map |= (ref_idx == prev_ref) * 2;
	avail |= 2;
	determine_pmv(mva, mvb, mvc, mv, avail, idx_map);
}

template <int IS_HORIZ, int MV_WIDTH, int MV_STEP2>
static inline uint32_t str_mv_calc16x8_left(const h264d_mb_current *mb, uint32_t str, const int8_t ref_idx[], const h264d_vector_set_t *mv, const prev_mb_t *prev)
{
	for (int i = 0; i < 2; ++i) {
		uint32_t mask = 0xa << (i * 4);
		if ((str & mask) != mask) {
			int prev_ref0 = frame_idx_of_ref(mb, prev->ref[i][0], 0);
			int prev_ref1 = frame_idx_of_ref(mb, prev->ref[i][1], 1);
			int ref0 = frame_idx_of_ref(mb, ref_idx[0], 0);
			int ref1 = frame_idx_of_ref(mb, ref_idx[1], 1);
			if (((prev_ref0 != ref0) || (prev_ref1 != ref1))
			    && ((prev_ref1 != ref0) || (prev_ref0 != ref1))) {
				uint32_t m = mask >> 1;
				str |= (((str >> 1) ^ m) & m);
			} else {
				str = str_mv_calc16x16_mv<MV_WIDTH * (MV_STEP2 / 2)>(str, ref0, ref1, prev_ref0, i * 2, &mv[0], prev);
			}
		}
		ref_idx += (MV_WIDTH == 1) ? 2 : 4;
		mv += MV_WIDTH * MV_STEP2;
	}
	return str;
}


static inline bool is_str_mv_calc16x8_center_bidir(int top_ref0, int bot_ref0, const h264d_vector_set_t *mv)
{
	const int16_t *mv_top0, *mv_top1;
	const int16_t *mv_bot0, *mv_bot1;
	if (top_ref0 == bot_ref0) {
		mv_top0 = mv[0].mv[0].v;
		mv_top1 = mv[0].mv[1].v;
	} else {
		mv_top1 = mv[0].mv[0].v;
		mv_top0 = mv[0].mv[1].v;
	}
	mv_bot0 = mv[1].mv[0].v;
	mv_bot1 = mv[1].mv[1].v;
	return (DIF_ABS_LARGER_THAN4(*mv_top0++, *mv_bot0++)
		|| DIF_ABS_LARGER_THAN4(*mv_top1++, *mv_bot1++)
		|| DIF_ABS_LARGER_THAN4(*mv_top0, *mv_bot0)
		|| DIF_ABS_LARGER_THAN4(*mv_top1, *mv_bot1));
}

static inline bool is_str_mv_calc16x8_center_onedir(int top_ref0, int bot_ref0, const h264d_vector_set_t *mv)
{
	const int16_t *top_mv = mv[0].mv[top_ref0 < 0].v;
	const int16_t *bot_mv = mv[1].mv[bot_ref0 < 0].v;
	return (DIF_ABS_LARGER_THAN4(*top_mv++, *bot_mv++)
		|| DIF_ABS_LARGER_THAN4(*top_mv, *bot_mv));
}

static inline uint32_t str_mv_calc16x8_vert(const h264d_mb_current *mb, uint32_t str, const int8_t ref_idx[], const h264d_vector_set_t *mv)
{
	if ((str & 0xaa0000) == 0xaa0000) {
		return str;
	}
	int top_ref0 = frame_idx_of_ref(mb, *ref_idx++, 0);
	int top_ref1 = frame_idx_of_ref(mb, *ref_idx++, 1);
	int bot_ref0 = frame_idx_of_ref(mb, *ref_idx++, 0);
	int bot_ref1 = frame_idx_of_ref(mb, *ref_idx, 1);
	if ((((top_ref0 != bot_ref0) || (top_ref1 != bot_ref1)) && ((top_ref1 != bot_ref0) || (top_ref0 != bot_ref1)))
		|| (((0 <= top_ref0) && (0 <= top_ref1)) ? is_str_mv_calc16x8_center_bidir : is_str_mv_calc16x8_center_onedir)(top_ref0, bot_ref0, mv)) {
		uint32_t mask = 0x550000;
		str |= (((str >> 1) ^ mask) & mask);
	}
	return str;
}

static inline void store_col16x8(h264d_col_mb_t *col, const int8_t *ref_idx, const h264d_vector_set_t *mv)
{
	h264d_vector_t *mvdst = col->mv;
	int8_t *refdst = col->ref;
	col->type = COL_MB16x8;
	for (int y = 0; y < 2; ++y) {
		int refcol;
		uint32_t mvcol;
		if (0 <= ref_idx[0]) {
			refcol = ref_idx[0];
			mvcol = mv[y].mv[0].vector;
		} else {
			refcol = ref_idx[1];
			mvcol = mv[y].mv[1].vector;
		}
		refdst[0] = refcol;
		refdst[1] = refcol;
		for (int i = 0; i < 16 / 2; ++i) {
			mvdst[i].vector = mvcol; 
		}
		ref_idx += 2;
		refdst += 2;
		mvdst += 8;
	}
}

template <int IS_HORIZ, int MV_WIDTH, int MV_STEP2>
static uint32_t store_str_inter8xedge(const h264d_mb_current *mb, const prev_mb_t *prev4x4inter, int8_t& str4,  const h264d_vector_set_t mv[], const int8_t ref_idx[], uint32_t str, uint32_t coeff4x4)
{
	if (prev4x4inter->type <= MB_IPCM) {
		str4 = 1;
		str |= 0xaa;
	} else {
		str = str_previous_coef(str, coeff4x4);
		str = str_mv_calc16x8_left<IS_HORIZ, MV_WIDTH, MV_STEP2>(mb, str, ref_idx, mv, prev4x4inter);
	}
	return str;
}

static void store_info_inter16x8(h264d_mb_current *mb, const h264d_vector_set_t mv[], const int8_t ref_idx[], uint32_t str_vert, uint32_t str_horiz, uint32_t left4x4, uint32_t top4x4)
{
	deblock_info_t *deb = mb->deblock_curr;
	deb->qpy = mb->qp;
	deb->qpc = mb->qp_chroma;
	if (mb->y != 0) {
		str_vert = store_str_inter16xedge(mb, mb->top4x4inter, deb->str4_vert, mv, ref_idx, str_vert, top4x4);
	}
	deb->str_vert = str_mv_calc16x8_vert(mb, str_vert, ref_idx, mv);
	if (mb->x != 0) {
		str_horiz = store_str_inter8xedge<1, 1, 1>(mb, mb->left4x4inter, deb->str4_horiz, mv, ref_idx, str_horiz, left4x4);
	}
	deb->str_horiz = str_horiz;

	mb->left4x4pred = 0x22222222;
	*mb->top4x4pred = 0x22222222;
	mb->lefttop_ref[0] = mb->top4x4inter->ref[1][0];
	mb->lefttop_ref[1] = mb->top4x4inter->ref[1][1];
	mb->lefttop_mv[0].vector = mb->top4x4inter->mov[3].mv[0].vector;
	mb->lefttop_mv[1].vector = mb->top4x4inter->mov[3].mv[1].vector;
	mb->left4x4inter->direct8x8 = 0;
	mb->top4x4inter->direct8x8 = 0;
	for (int i = 0; i < 4; ++i) {
		mb->top4x4inter->mov[i] = mv[1];
		mb->top4x4inter->mvd[i] = mv[3];
	}
	for (int i = 0; i < 2; ++i) {
		mb->top4x4inter->ref[i][0] = ref_idx[2];
		mb->top4x4inter->ref[i][1] = ref_idx[3];
		mb->left4x4inter->mov[i] = mv[0];
		mb->left4x4inter->mvd[i] = mv[2];
		mb->left4x4inter->mov[2 + i] = mv[1];
		mb->left4x4inter->mvd[2 + i] = mv[3];
	}
	mb->left4x4inter->ref[0][0] = ref_idx[0];
	mb->left4x4inter->ref[0][1] = ref_idx[1];
	mb->left4x4inter->ref[1][0] = ref_idx[2];
	mb->left4x4inter->ref[1][1] = ref_idx[3];
	store_col16x8(mb->col_curr, ref_idx, mv);
}

template <typename F0 ,typename F1, typename F2, typename F3, typename F4>
static int mb_inter16x8(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail,
			 F0 RefIdx16x8,
			 F1 MvdXY,
			 F2 CodedBlockPattern,
			 F3 QpDelta,
			 F4 ResidualBlock)
{
	uint32_t left4x4, top4x4;
	uint32_t refmap;
	uint32_t cbp;
	uint32_t str_vert, str_horiz;
	const int16_t *mvd_a, *mvd_b;
	int bidir0, bidir1;
	h264d_vector_set_t mv[2][2];
	int8_t ref_idx[4];

	refmap = mbc->cbp;
	RefIdx16x8(mb, st, ref_idx, refmap, avail);
	memset(mv, 0, sizeof(mv));
	bidir0 = 0;
	bidir1 = 0;
	for (int lx = 0; lx < 2; ++lx) {
		if (refmap & 1) {
			calc_mv16x8top(mb, mv[0][0].mv[lx].v, mvd_a, mvd_b, lx, ref_idx[lx], avail);
			MvdXY(mb, st, mv[1][0].mv[lx].v, mvd_a, mvd_b);
			mv[0][0].mv[lx].v[0] += mv[1][0].mv[lx].v[0];
			mv[0][0].mv[lx].v[1] += mv[1][0].mv[lx].v[1];
			inter_pred8x8(mb, mv[0][0].mv[lx].v, 16, 8, mb->frame->refs[lx][ref_idx[lx]].frame_idx, 0, 0, bidir0);
			bidir0++;
		}
		if (refmap & 2) {
			calc_mv16x8bottom(mb, mv[0][1].mv[lx].v, mvd_a, mvd_b, lx, ref_idx[lx + 2], avail, ref_idx[lx], &mv[0][0]);
			MvdXY(mb, st, mv[1][1].mv[lx].v, mvd_a, mvd_b);
			mv[0][1].mv[lx].v[0] += mv[1][1].mv[lx].v[0];
			mv[0][1].mv[lx].v[1] += mv[1][1].mv[lx].v[1];
			inter_pred8x8(mb, mv[0][1].mv[lx].v, 16, 8, mb->frame->refs[lx][ref_idx[lx + 2]].frame_idx, 0, 8, bidir1);
			bidir1++;
		}
		refmap >>= 2;
	}
	left4x4 = mb->left4x4coef;
	top4x4 = *mb->top4x4coef;
	mb->cbp = cbp = CodedBlockPattern(mb, st, avail);
	if (cbp) {
		str_vert = residual_luma_inter(mb, cbp, st, avail, QpDelta, ResidualBlock);
		str_horiz = expand_coef_str(transposition(str_vert));
		str_vert = expand_coef_str(str_vert);
	} else {
		mb->prev_qp_delta = 0;
		mb->left4x4coef = 0;
		*mb->top4x4coef = 0;
		str_horiz = 0;
		str_vert = 0;
	}
	store_info_inter16x8(mb, &mv[0][0], ref_idx, str_vert, str_horiz, left4x4, top4x4);
	return residual_chroma(mb, cbp, st, avail, ResidualBlock);
}

static inline void calc_mv8x16left(h264d_mb_current *mb, int16_t pmv[], const int16_t *&mvd_a, const int16_t *&mvd_b, int lx, int ref_idx, int avail)
{
	prev_mb_t *pmb;
	const int16_t *mva, *mvb, *mvc;
	int idx_map;

	if (avail & 1) {
		pmb = mb->left4x4inter;
		mvd_a = pmb->mvd[0].mv[lx].v;
		if (ref_idx == pmb->ref[0][lx]) {
			pmv[0] = pmb->mov[0].mv[lx].v[0];
			pmv[1] = pmb->mov[0].mv[lx].v[1];
			mvd_b = (avail & 2) ? mb->top4x4inter->mvd[0].mv[lx].v : zero_mv;
			return;
		}
		mva = pmb->mov[0].mv[lx].v;
	} else {
		mvd_a = mva = zero_mv;
	}
	idx_map = 0;
	if (avail & 2) {
		pmb = mb->top4x4inter;
		idx_map |= (ref_idx == pmb->ref[0][lx]) * 2;
		idx_map |= (ref_idx == pmb->ref[1][lx]) * 4;
		avail |= 4;
		mvb = pmb->mov[0].mv[lx].v;
		mvd_b = pmb->mvd[0].mv[lx].v;
		mvc = pmb->mov[2].mv[lx].v;
	} else {
		mvd_b = mvb = zero_mv;
		avail &= ~4;
		if (avail & 8) {
			idx_map |= (ref_idx == mb->lefttop_ref[lx]) * 4;
			mvc = mb->lefttop_mv[lx].v;
		} else {
			mvc = zero_mv;
		}
	}
	determine_pmv(mva, mvb, mvc, pmv, avail, idx_map);
}

static inline void calc_mv8x16right(h264d_mb_current *mb, int16_t pmv[], const int16_t *&mvd_a, const int16_t *&mvd_b, int lx, int ref_idx, int avail, int prev_ref, const h264d_vector_set_t *prev_mv)
{
	prev_mb_t *pmb;
	const int16_t *mva, *mvb, *mvc;
	int idx_map;

	idx_map = 0;
	if (avail & 4) {
		pmb = mb->top4x4inter + 1;
		if (ref_idx == pmb->ref[0][lx]) {
			pmv[0] = pmb->mov[0].mv[lx].v[0];
			pmv[1] = pmb->mov[0].mv[lx].v[1];
			mvd_a = prev_mv[2].mv[lx].v;
			mvd_b = (avail & 2) ? mb->top4x4inter->mvd[2].mv[lx].v : zero_mv;
			return;
		}
		mvc = pmb->mov[0].mv[lx].v;
	} else if (avail & 2) {
		pmb = mb->top4x4inter;
		idx_map = (ref_idx == pmb->ref[0][lx]) * 4;
		mvd_b = pmb->mvd[2].mv[lx].v;
		if (idx_map) {
			pmv[0] = pmb->mov[1].mv[lx].v[0];
			pmv[1] = pmb->mov[1].mv[lx].v[1];
			mvd_a = prev_mv[2].mv[lx].v;
			return;
		} else {
			mvc = pmb->mov[1].mv[lx].v;
		}
	} else {
		mvc = zero_mv;
	}
	/* left block are always available */
	idx_map |= (ref_idx == prev_ref);
	mva = prev_mv->mv[lx].v;
	mvd_a = prev_mv[2].mv[lx].v;
	avail |= 1;
	if (avail & 2) {
		pmb = mb->top4x4inter;
		idx_map |= (ref_idx == pmb->ref[1][lx]) * 2;
		mvb = pmb->mov[2].mv[lx].v;
		mvd_b = pmb->mvd[2].mv[lx].v;
	} else {
		mvd_b = mvb = zero_mv;
	}
	determine_pmv(mva, mvb, mvc, pmv, avail, idx_map);
}

static inline void store_col8x16(h264d_col_mb_t *col, const int8_t *ref_idx, const h264d_vector_set_t *mv)
{
	h264d_vector_t *mvdst = col->mv;
	int8_t *refdst = col->ref;
	col->type = COL_MB8x16;
	for (int x = 0; x < 2; ++x) {
		int refcol;
		uint32_t mvcol;
		if (0 <= ref_idx[0]) {
			refcol = ref_idx[0];
			mvcol = mv[x].mv[0].vector;
		} else {
			refcol = ref_idx[1];
			mvcol = mv[x].mv[1].vector;
		}
		refdst[0] = refcol;
		refdst[2] = refcol;
		uint32_t *dst = &mvdst[0].vector;
		int i = 4;
		do {
			dst[0] = mvcol;
			dst[1] = mvcol;
			dst += 4;
		} while (--i);
		ref_idx += 2;
		refdst += 1;
		mvdst += 2;
	}
}

static void store_info_inter8x16(h264d_mb_current *mb, const h264d_vector_set_t mv[], const int8_t ref_idx[], uint32_t str_vert, uint32_t str_horiz, uint32_t left4x4, uint32_t top4x4)
{
	deblock_info_t *deb = mb->deblock_curr;
	deb->qpy = mb->qp;
	deb->qpc = mb->qp_chroma;
	if (mb->y != 0) {
		str_vert = store_str_inter8xedge<0, 1, 1>(mb, mb->top4x4inter, deb->str4_vert, mv, ref_idx, str_vert, top4x4);
	}
	deb->str_vert = str_vert;
	if (mb->x != 0) {
		str_horiz = store_str_inter16xedge(mb, mb->left4x4inter, deb->str4_horiz, mv, ref_idx, str_horiz, left4x4);
	}
	deb->str_horiz = str_mv_calc16x8_vert(mb, str_horiz, ref_idx, mv);

	mb->left4x4pred = 0x22222222;
	*mb->top4x4pred = 0x22222222;
	mb->left4x4inter->direct8x8 = 0;
	mb->top4x4inter->direct8x8 = 0;
	for (int i = 0; i < 2; ++i) {
		mb->lefttop_ref[i] = mb->top4x4inter->ref[1][i];
		mb->top4x4inter->ref[i][0] = ref_idx[i * 2];
		mb->top4x4inter->ref[i][1] = ref_idx[i * 2 + 1];
		mb->left4x4inter->ref[i][0] = ref_idx[2];
		mb->left4x4inter->ref[i][1] = ref_idx[3];
		mb->lefttop_mv[i].vector = mb->top4x4inter->mov[3].mv[i].vector;
		mb->top4x4inter->mov[i] = mv[0];
		mb->top4x4inter->mvd[i] = mv[2];
		mb->top4x4inter->mov[i + 2] = mv[1];
		mb->top4x4inter->mvd[i + 2] = mv[3];
	}
	for (int i = 0; i < 4; ++i) {
		mb->left4x4inter->mov[i] = mv[1];
		mb->left4x4inter->mvd[i] = mv[3];
	}
	store_col8x16(mb->col_curr, ref_idx, mv);
}

template <typename F0 ,typename F1, typename F2, typename F3, typename F4>
static int mb_inter8x16(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail,
			 F0 RefIdx8x16,
			 F1 MvdXY,
			 F2 CodedBlockPattern,
			 F3 QpDelta,
			 F4 ResidualBlock)
{
	h264d_vector_set_t mv[2][2];
	const int16_t *mvd_a, *mvd_b;
	int bidir0, bidir1;
	uint32_t refmap;
	uint32_t cbp;
	uint32_t top4x4, left4x4;
	uint32_t str_vert, str_horiz;
	int8_t ref_idx[4];

	refmap = mbc->cbp;
	RefIdx8x16(mb, st, ref_idx, refmap, avail);
	memset(mv, 0, sizeof(mv));
	bidir0 = 0;
	bidir1 = 0;
	for (int lx = 0; lx < 2; ++lx) {
		if (refmap & 1) {
			calc_mv8x16left(mb, mv[0][0].mv[lx].v, mvd_a, mvd_b, lx, ref_idx[lx], avail);
			MvdXY(mb, st, mv[1][0].mv[lx].v, mvd_a, mvd_b);
			mv[0][0].mv[lx].v[0] += mv[1][0].mv[lx].v[0];
			mv[0][0].mv[lx].v[1] += mv[1][0].mv[lx].v[1];
			inter_pred8x8(mb, mv[0][0].mv[lx].v, 8, 16, mb->frame->refs[lx][ref_idx[lx]].frame_idx, 0, 0, bidir0);
			bidir0++;
		}
		if (refmap & 2) {
			calc_mv8x16right(mb, mv[0][1].mv[lx].v, mvd_a, mvd_b, lx, ref_idx[lx + 2], avail, ref_idx[lx], &mv[0][0]);
			MvdXY(mb, st, mv[1][1].mv[lx].v, mvd_a, mvd_b);
			mv[0][1].mv[lx].v[0] += mv[1][1].mv[lx].v[0];
			mv[0][1].mv[lx].v[1] += mv[1][1].mv[lx].v[1];
			inter_pred8x8(mb, mv[0][1].mv[lx].v, 8, 16, mb->frame->refs[lx][ref_idx[lx + 2]].frame_idx, 8, 0, bidir1);
			bidir1++;
		}
		refmap >>= 2;
	}
	left4x4 = mb->left4x4coef;
	top4x4 = *mb->top4x4coef;
	mb->cbp = cbp = CodedBlockPattern(mb, st, avail);
	if (cbp) {
		str_vert = residual_luma_inter(mb, cbp, st, avail, QpDelta, ResidualBlock);
		str_horiz = expand_coef_str(transposition(str_vert));
		str_vert = expand_coef_str(str_vert);
	} else {
		mb->prev_qp_delta = 0;
		mb->left4x4coef = 0;
		*mb->top4x4coef = 0;
		str_vert = 0;
		str_horiz = 0;
	}
	store_info_inter8x16(mb, &mv[0][0], ref_idx, str_vert, str_horiz, left4x4, top4x4);
	return residual_chroma(mb, cbp, st, avail, ResidualBlock);
}

static inline void calc_mv8x8_sub8x8(h264d_mb_current *mb, int16_t *pmv, const int16_t *&mvd_a, const int16_t *&mvd_b, int avail, int lx, int ref_idx, int blk_idx, prev8x8_t *pblk)
{
	const int16_t *mva;
	const int16_t *mvb;
	const int16_t *mvc;
	prev_mb_t *pmb;
	int idx_map;

	if (blk_idx & 1) {
		idx_map = (ref_idx == pblk[blk_idx - 1].ref[lx]);
		mva = pblk[blk_idx - 1].mv[1][lx].v;
		mvd_a = pblk[blk_idx - 1].mvd[1][lx].v;
		avail |= 1;
	} else if (avail & 1) {
		pmb = mb->left4x4inter;
		idx_map = (ref_idx == pmb->ref[blk_idx >> 1][lx]);
		mva = pmb->mov[blk_idx].mv[lx].v;
		mvd_a = pmb->mvd[blk_idx].mv[lx].v;
	} else {
		idx_map = 0;
		mvd_a = mva = zero_mv;
	}

	if (blk_idx & 2) {
		idx_map |= (ref_idx == pblk[blk_idx - 2].ref[lx]) * 2;
		mvb = pblk[blk_idx - 2].mv[2][lx].v;
		mvd_b = pblk[blk_idx - 2].mvd[2][lx].v;
		avail |= 2;
	} else if (avail & 2) {
		pmb = mb->top4x4inter;
		idx_map |= (ref_idx == pmb->ref[blk_idx][lx]) * 2;
		mvb = pmb->mov[blk_idx * 2].mv[lx].v;
		mvd_b = pmb->mvd[blk_idx * 2].mv[lx].v;
	} else {
		mvd_b = mvb = zero_mv;
	}

	switch (blk_idx) {
	case 0:
		if (avail & 2) {
			pmb = mb->top4x4inter;
			idx_map |= (ref_idx == pmb->ref[1][lx]) * 4;
			mvc = pmb->mov[2].mv[lx].v;
			avail |= 4;
		} else if (avail & 8) {
			idx_map |= (ref_idx == mb->lefttop_ref[lx]) * 4;
			mvc = mb->lefttop_mv[lx].v;
			avail |= 4;
		} else {
			avail &= ~4;
			mvc = zero_mv;
		}
		break;
	case 1:
		if (avail & 4) {
			pmb = mb->top4x4inter + 1;
			idx_map |= (ref_idx == pmb->ref[0][lx]) * 4;
			mvc = pmb->mov[0].mv[lx].v;
		} else if (avail & 2) {
			pmb = mb->top4x4inter;
			idx_map |= (ref_idx == pmb->ref[0][lx]) * 4;
			mvc = pmb->mov[1].mv[lx].v;
		} else {
			mvc = zero_mv;
		}
		break;
	case 2:
		idx_map |= (ref_idx == pblk[1].ref[lx]) * 4;
		mvc = pblk[1].mv[2][lx].v;
		avail |= 4;
		break;
	case 3:
		idx_map |= (ref_idx == pblk[0].ref[lx]) * 4;
		mvc = pblk[0].mv[3][lx].v;
		avail |= 4;
		break;
	}
	determine_pmv(mva, mvb, mvc, pmv, avail, idx_map);
}

static inline void calc_mv8x8_sub8x4(h264d_mb_current *mb, int16_t *pmv, const int16_t *&mvd_a, const int16_t *&mvd_b, int avail, int lx, int ref_idx, int blk_idx, prev8x8_t *pblk, int y)
{
	const int16_t *mva;
	const int16_t *mvb;
	const int16_t *mvc;
	prev_mb_t *pmb;
	int idx_map;

	if (blk_idx & 1) {
		idx_map = (ref_idx == pblk[blk_idx - 1].ref[lx]);
		mva = pblk[blk_idx - 1].mv[y * 2 + 1][lx].v;
		mvd_a = pblk[blk_idx - 1].mvd[y * 2 + 1][lx].v;
		avail |= 1;
	} else if (avail & 1) {
		pmb = mb->left4x4inter;
		idx_map = (ref_idx == pmb->ref[blk_idx >> 1][lx]);
		mva = pmb->mov[(blk_idx & 2) + y].mv[lx].v;
		mvd_a = pmb->mvd[(blk_idx & 2) + y].mv[lx].v;
	} else {
		idx_map = 0;
		mvd_a = mva = zero_mv;
	}

	if (y != 0) {
		idx_map |= 2;
		mvb = pblk[blk_idx].mv[0][lx].v;
		mvd_b = pblk[blk_idx].mvd[0][lx].v;
		avail |= 2;
	} else if (blk_idx & 2) {
		idx_map |= (ref_idx == pblk[blk_idx - 2].ref[lx]) * 2;
		mvb = pblk[blk_idx - 2].mv[2][lx].v;
		mvd_b = pblk[blk_idx - 2].mvd[2][lx].v;
		avail |= 2;
	} else if (avail & 2) {
		pmb = mb->top4x4inter;
		idx_map |= (ref_idx == pmb->ref[blk_idx & 1][lx]) * 2;
		mvb = pmb->mov[blk_idx * 2].mv[lx].v;
		mvd_b = pmb->mvd[blk_idx * 2].mv[lx].v;
	} else {
		mvd_b = mvb = zero_mv;
	}

	switch (blk_idx) {
	case 0:
		if (y == 0) {
			if (avail & 2) {
				pmb = mb->top4x4inter;
				idx_map |= (ref_idx == pmb->ref[1][lx]) * 4;
				avail |= 4;
				mvc = pmb->mov[2].mv[lx].v;
			} else if (avail & 8) {
				idx_map |= (ref_idx == mb->lefttop_ref[lx]) * 4;
				avail |= 4;
				mvc = mb->lefttop_mv[lx].v;
			} else {
				avail &= ~4;
				mvc = zero_mv;
			}
		} else if (avail & 1) {
			pmb = mb->left4x4inter;
			idx_map |= (ref_idx == pmb->ref[0][lx]) * 4;
			mvc = pmb->mov[0].mv[lx].v;
			avail |= 4;
		} else {
			avail &= ~4;
			mvc = zero_mv;
		}
		break;
	case 1:
		if (y == 0) {
			if (avail & 4) {
				pmb = mb->top4x4inter + 1;
				idx_map |= (ref_idx == pmb->ref[0][lx]) * 4;
				mvc = pmb->mov[0].mv[lx].v;
				avail |= 4;
			} else if (avail & 2) {
				pmb = mb->top4x4inter;
				idx_map |= (ref_idx == pmb->ref[0][lx]) * 4;
				mvc = pmb->mov[1].mv[lx].v;
				avail |= 4;
			} else {
				mvc = zero_mv;
			}
		} else {
			idx_map |= (ref_idx == pblk[0].ref[lx]) * 4;
			mvc = pblk[0].mv[1][lx].v;
			avail |= 4;
		}
		break;
	case 2:
		if (y == 0) {
			idx_map |= (ref_idx == pblk[1].ref[lx]) * 4;
			mvc = pblk[1].mv[2][lx].v;
			avail |= 4;
		} else if (avail & 1) {
			pmb = mb->left4x4inter;
			idx_map |= (ref_idx == pmb->ref[1][lx]) * 4;
			mvc = pmb->mov[2].mv[lx].v;
			avail |= 4;
		} else {
			avail &= ~4;
			mvc = zero_mv;
		}
		break;
	case 3:
		idx_map |= (ref_idx == pblk[y * 2].ref[lx]) * 4;
		mvc = pblk[y * 2].mv[3 - y * 2][lx].v;
		avail |= 4;
		break;
	}
	determine_pmv(mva, mvb, mvc, pmv, avail, idx_map);
}

static inline void calc_mv8x8_sub4x8(h264d_mb_current *mb, int16_t *pmv, const int16_t *&mvd_a, const int16_t *&mvd_b, int avail, int lx, int ref_idx, int blk_idx, prev8x8_t *pblk, int x)
{
	const int16_t *mva;
	const int16_t *mvb;
	const int16_t *mvc;
	prev_mb_t *pmb;
	int idx_map;

	if (x != 0) {
		idx_map = 1;
		mva = pblk[blk_idx].mv[0][lx].v;
		mvd_a = pblk[blk_idx].mvd[0][lx].v;
		avail |= 1;
	} else if (blk_idx & 1) {
		idx_map = (ref_idx == pblk[blk_idx - 1].ref[lx]);
		mva = pblk[blk_idx - 1].mv[1][lx].v;
		mvd_a = pblk[blk_idx - 1].mvd[1][lx].v;
		avail |= 1;
	} else if (avail & 1) {
		pmb = mb->left4x4inter;
		idx_map = (ref_idx == pmb->ref[blk_idx >> 1][lx]);
		mva = pmb->mov[blk_idx].mv[lx].v; /* blk_idx shall be 0 or 2 */
		mvd_a = pmb->mvd[blk_idx].mv[lx].v;
	} else {
		idx_map = 0;
		mvd_a = mva = zero_mv;
	}

	if (blk_idx & 2) {
		idx_map |= (ref_idx == pblk[blk_idx - 2].ref[lx]) * 2;
		mvb = pblk[blk_idx - 2].mv[2 + x][lx].v;
		mvd_b = pblk[blk_idx - 2].mvd[2 + x][lx].v;
		avail |= 2;
	} else if (avail & 2) {
		pmb = mb->top4x4inter;
		idx_map |= (ref_idx == pmb->ref[blk_idx & 1][lx]) * 2;
		mvb = pmb->mov[blk_idx * 2 + x].mv[lx].v;
		mvd_b = pmb->mvd[blk_idx * 2 + x].mv[lx].v;
	} else {
		mvd_b = mvb = zero_mv;
	}

	switch (blk_idx) {
	case 0:
		if (avail & 2) {
			pmb = mb->top4x4inter;
			idx_map |= (ref_idx == pmb->ref[x][lx]) * 4;
			mvc = pmb->mov[x + 1].mv[lx].v;
			avail |= 4;
		} else {
			avail &= ~4;
			if (x == 0 && (avail & 8)) {
				idx_map |= (ref_idx == mb->lefttop_ref[lx]) * 4;
				mvc = mb->lefttop_mv[lx].v;
			} else {
				mvc = zero_mv;
			}
		}
		break;
	case 1:
		if (x == 0) {
			if (avail & 2) {
				pmb = mb->top4x4inter;
				idx_map |= (ref_idx == pmb->ref[1][lx]) * 4;
				mvc = pmb->mov[3].mv[lx].v;
				avail |= 4;
			} else {
				avail &= ~4;
				mvc = zero_mv;
			}
		} else {
			if (avail & 4) {
				pmb = mb->top4x4inter + 1;
				idx_map |= (ref_idx == pmb->ref[0][lx]) * 4;
				mvc = pmb->mov[0].mv[lx].v;
			} else if (avail & 2) {
				pmb = mb->top4x4inter;
				idx_map |= (ref_idx == pmb->ref[1][lx]) * 4;
				if (0 <= pmb->ref[1][lx]) {
					mvc = pmb->mov[2].mv[lx].v;
				} else {
					mvc = zero_mv;
				}
			} else {
				mvc = zero_mv;
			}
		}
		break;
	case 2:
		avail |= 4;
		idx_map |= (ref_idx == pblk[x].ref[lx]) * 4;
		mvc = pblk[x].mv[3 - x][lx].v;
		break;
	case 3:
		avail |= 4;
		idx_map |= (ref_idx == pblk[1].ref[lx]) * 4;
		mvc = pblk[1].mv[3 - x][lx].v;
		break;
	}
	determine_pmv(mva, mvb, mvc, pmv, avail, idx_map);
}

static inline void calc_mv8x8_sub4x4(h264d_mb_current *mb, int16_t *pmv, const int16_t *&mvd_a, const int16_t *&mvd_b, int avail, int lx, int ref_idx, int blk_idx, prev8x8_t *pblk, int xy)
{
	const int16_t *mva;
	const int16_t *mvb;
	const int16_t *mvc;
	prev_mb_t *pmb;
	int idx_map;

	if (xy & 1) {
		idx_map = 1;
		mva = pblk[blk_idx].mv[xy - 1][lx].v;
		mvd_a = pblk[blk_idx].mvd[xy - 1][lx].v;
		avail |= 1;
	} else if (blk_idx & 1) {
		idx_map = (ref_idx == pblk[blk_idx - 1].ref[lx]);
		mva = pblk[blk_idx - 1].mv[xy + 1][lx].v; /* xy shall be 0 or 2 */
		mvd_a = pblk[blk_idx - 1].mvd[xy + 1][lx].v;
		avail |= 1;
	} else if (avail & 1) {
		pmb = mb->left4x4inter;
		idx_map = (ref_idx == pmb->ref[blk_idx >> 1][lx]);
		mva = pmb->mov[blk_idx + (xy >> 1)].mv[lx].v; /* blk_idx shall be 0 or 2 */
		mvd_a = pmb->mvd[blk_idx + (xy >> 1)].mv[lx].v;
	} else {
		idx_map = 0;
		mvd_a = mva = zero_mv;
	}

	if (xy & 2) {
		idx_map |= 2;
		mvb = pblk[blk_idx].mv[xy - 2][lx].v;
		mvd_b = pblk[blk_idx].mvd[xy - 2][lx].v;
		avail |= 2;
	} else if (blk_idx & 2) {
		idx_map |= (ref_idx == pblk[blk_idx - 2].ref[lx]) * 2;
		mvb = pblk[blk_idx - 2].mv[2 + (xy & 1)][lx].v;
		mvd_b = pblk[blk_idx - 2].mvd[2 + (xy & 1)][lx].v;
		avail |= 2;
	} else if (avail & 2) {
		pmb = mb->top4x4inter;
		idx_map |= (ref_idx == pmb->ref[blk_idx & 1][lx]) * 2;
		mvb = pmb->mov[blk_idx * 2 + (xy & 1)].mv[lx].v;
		mvd_b = pmb->mvd[blk_idx * 2 + (xy & 1)].mv[lx].v;
	} else {
		mvd_b = mvb = zero_mv;
	}

	switch (blk_idx) {
	case 0:
		switch (xy) {
		case 0:
			if (avail & 2) {
				pmb = mb->top4x4inter;
				idx_map |= (ref_idx == pmb->ref[0][lx]) * 4;
				avail |= 4;
				mvc = pmb->mov[1].mv[lx].v;
			} else if (avail & 8) {
				avail &= ~4;
				idx_map |= (ref_idx == mb->lefttop_ref[lx]) * 4;
				mvc = mb->lefttop_mv[lx].v;
			} else {
				avail &= ~4;
				mvc = zero_mv;
			}
			break;
		case 1:
			if (avail & 2) {
				pmb = mb->top4x4inter;
				idx_map |= (ref_idx == pmb->ref[1][lx]) * 4;
				avail |= 4;
				mvc = pmb->mov[2].mv[lx].v;
			} else {
				avail &= ~4;
				mvc = zero_mv;
			}
			break;
		case 2:
			idx_map |= 4;
			avail |= 4;
			mvc = pblk[blk_idx].mv[1][lx].v;
			break;
		case 3:
			idx_map |= 4;
			avail |= 4;
			mvc = pblk[blk_idx].mv[0][lx].v;
			break;
		}
		break;
	case 1:
		switch (xy) {
		case 0:
			if (avail & 2) {
				pmb = mb->top4x4inter;
			idx_map |= (ref_idx == pmb->ref[1][lx]) * 4;
				mvc = pmb->mov[3].mv[lx].v;
				avail |= 4;
			} else {
				avail &= ~4;
				mvc = zero_mv;
			}
			break;
		case 1:
			if (avail & 4) {
				pmb = mb->top4x4inter + 1;
				idx_map |= (ref_idx == pmb->ref[0][lx]) * 4;
				mvc = pmb->mov[0].mv[lx].v;
			} else if (avail & 2) {
				pmb = mb->top4x4inter;
				idx_map |= (ref_idx == pmb->ref[1][lx]) * 4;
				mvc = pmb->mov[2].mv[lx].v;
				avail |= 4;
			} else {
				mvc = zero_mv;
			}
			break;
		case 2:
		case 3:
			idx_map |= 4;
			avail |= 4;
			mvc = pblk[blk_idx].mv[3 - xy][lx].v;
			break;
		}
		break;
	case 2:
		avail |= 4;
		switch (xy) {
		case 0:
		case 1:
			idx_map |= (ref_idx == pblk[xy].ref[lx]) * 4;
			mvc = pblk[xy].mv[3 - xy][lx].v;
			break;
		case 2:
		case 3:
			idx_map |= 4;
			mvc = pblk[2].mv[3 - xy][lx].v;
			break;
		}
		break;
	case 3:
		avail |= 4;
		switch (xy) {
		case 0:
		case 1:
			idx_map |= (ref_idx == pblk[1].ref[lx]) * 4;
			mvc = pblk[1].mv[3 - xy][lx].v;
			break;
		case 2:
		case 3:
			idx_map |= 4;
			mvc = pblk[3].mv[3 - xy][lx].v;
			break;
		}
		break;
	}
	determine_pmv(mva, mvb, mvc, pmv, avail, idx_map);
}

static inline void b_skip_ref_mv(int8_t *ref_idx, int16_t *mv, const int8_t *ref_a, const int8_t *ref_b, const int8_t *ref_c, const h264d_vector_t *mv_a, const h264d_vector_t *mv_b, const h264d_vector_t *mv_c)
{
	for (int lx = 0; lx < 2; ++lx) {
		int ra = *ref_a++;
		int rb = *ref_b++;
		int rc = *ref_c++;
		int ref = MIN((unsigned)ra, (unsigned)rb);
		ref = MIN((unsigned)ref, (unsigned)rc);
		if (ref < 0) {
			*(uint32_t *)mv = 0U;
		} else if ((ra == ref) && (rb != ref) && (rc != ref)) {
			*(uint32_t *)mv = mv_a->vector;
		} else if ((ra != ref) && (rb == ref) && (rc != ref)) {
			*(uint32_t *)mv = mv_b->vector;
		} else if ((ra != ref) && (rb != ref) && (rc == ref)) {
			*(uint32_t *)mv = mv_c->vector;
		} else {
			mv[0] = MEDIAN(mv_a->v[0], mv_b->v[0], mv_c->v[0]);
			mv[1] = MEDIAN(mv_a->v[1], mv_b->v[1], mv_c->v[1]);
		}
		*ref_idx++ = ref;
		mv_a++;
		mv_b++;
		mv_c++;
		mv += 2;
	}
}

static void b_direct_ref_mv_calc(h264d_mb_current *mb, int avail, int8_t *ref_idx, int16_t *mv)
{
	const int8_t *ref_a;
	const int8_t *ref_b;
	const int8_t *ref_c;
	const h264d_vector_t *mv_a;
	const h264d_vector_t *mv_b;
	const h264d_vector_t *mv_c;

	if (avail & 1) {
		ref_a = mb->left4x4inter->ref[0];
		mv_a = mb->left4x4inter->mov[0].mv;
	} else {
		ref_a = non_ref;
		mv_a = zero_mov;
	}
	if (avail & 2) {
		ref_b = mb->top4x4inter->ref[0];
		mv_b = mb->top4x4inter->mov[0].mv;
	} else {
		ref_b = non_ref;
		mv_b = zero_mov;
	}
	if (avail & 4) {
		ref_c = mb->top4x4inter[1].ref[0];
		mv_c = mb->top4x4inter[1].mov[0].mv;
	} else if (avail & 8) {
		ref_c = mb->lefttop_ref;
		mv_c = mb->lefttop_mv;
	} else {
		ref_c = non_ref;
		mv_c = zero_mov;
	}
	b_skip_ref_mv(ref_idx, mv, ref_a, ref_b, ref_c, mv_a, mv_b, mv_c);
}

static inline void direct_mv_pred(h264d_mb_current *mb, const int8_t *ref_idx, const h264d_vector_t *mv, int xsize, int ysize, int xoffset, int yoffset)
{
	int bidir = 0;
	for (int lx = 0; lx < 2; ++lx) {
		int ref = *ref_idx++;
		if (0 <= ref) {
			inter_pred8x8(mb, mv->v, xsize, ysize, mb->frame->refs[lx][ref].frame_idx, xoffset, yoffset, bidir++);
		}
		mv++;
	}
}

/*
 4x4:   0x00000001
 4x8:   0x00000101
 8x4:   0x00000005
 8x8:   0x00000505
 8x16:  0x05050505
 16x8:  0x00005555
 16x16: 0x55555555
*/
#define COLBITS(X, Y) ((((X) == 4) ? 1 : (((X) == 8) ? 5 : 0x55)) * (((Y) == 4) ? 1 : (((Y) == 8) ? 0x101 : 0x01010101)))

template <int N, int X, int Y>
struct pred_direct_col_block_bidir {
	void operator()(h264d_mb_current *mb, const h264d_vector_t *mvcol, h264d_vector_t *mv_curr, const int8_t *ref_idx, int xofs, int yofs) const {
		if ((mv_curr[0].vector != 0U || mv_curr[1].vector != 0U) && ((unsigned)(mvcol->v[0] + 1) <= 2U) && ((unsigned)(mvcol->v[1] + 1) <= 2U)) {
			mv_curr[0].vector = 0U;
			mv_curr[1].vector = 0U;
			direct_zero_pred(mb, X, Y, xofs, yofs);
		} else {
			direct_mv_pred(mb, ref_idx, mv_curr, X, Y, xofs, yofs);
		}
	}
};

template <int LX, int N, int X, int Y>
struct pred_direct_col_block_onedir {
	void operator()(h264d_mb_current *mb, const h264d_vector_t *mvcol, h264d_vector_t *mv_curr, const int8_t *ref_idx, int xofs, int yofs) const {
		if ((mv_curr[LX].vector != 0U) && ((unsigned)(mvcol->v[0] + 1) <= 2U) && ((unsigned)(mvcol->v[1] + 1) <= 2U)) {
			mv_curr[LX].vector = 0U;
		}
		direct_mv_pred(mb, ref_idx, mv_curr, X, Y, xofs, yofs);
	}
};

static void pred_direct8x8_block8x8_nocol(h264d_mb_current *mb,  const h264d_vector_t *mvcol, const int8_t *ref_idx, h264d_vector_t *mv, int blk_idx)
{
	direct_mv_pred(mb, ref_idx, mv, 8, 8, (blk_idx & 1) * 8, (blk_idx & 2) * 4);
}

template <int N, int BLOCK, typename F0>
static inline void pred_direct_block(h264d_mb_current *mb, const h264d_vector_t *mvcol, const int8_t *ref_idx, h264d_vector_t *mv, int blk_idx,
				    F0 PredDirectBlock)
{
	int xoffset = (blk_idx & 1) * 8;
	int yoffset = (blk_idx & 2) * 4;
	for (int i = 0; i < 64 / (BLOCK * BLOCK); ++i) {
		int xofs = xoffset + (i & 1) * 4;
		int yofs = yoffset + (i & 2) * 2;
		h264d_vector_t *mv_curr = &mv[(i & 2) * (N / 4) + (i & 1) * 2];
		PredDirectBlock(mb, mvcol, mv_curr, ref_idx, xofs, yofs);
		mvcol += (i & 1) ? 3 : 1;
	}
}

static void pred_direct8x8_block8x8_l0(h264d_mb_current *mb,  const h264d_vector_t *mvcol, const int8_t *ref_idx, h264d_vector_t *mv, int blk_idx)
{
	pred_direct_block<8, 8>(mb, mvcol, ref_idx, mv, blk_idx, pred_direct_col_block_onedir<0, 8, 8, 8>());
}

static void pred_direct8x8_block8x8_l1(h264d_mb_current *mb,  const h264d_vector_t *mvcol, const int8_t *ref_idx, h264d_vector_t *mv, int blk_idx)
{
	pred_direct_block<8, 8>(mb, mvcol, ref_idx, mv, blk_idx, pred_direct_col_block_onedir<1, 8, 8, 8>());
}

static void pred_direct8x8_block8x8_bidir(h264d_mb_current *mb,  const h264d_vector_t *mvcol, const int8_t *ref_idx, h264d_vector_t *mv, int blk_idx)
{
	pred_direct_block<8, 8>(mb, mvcol, ref_idx, mv, blk_idx, pred_direct_col_block_bidir<8, 8, 8>());
}

static void pred_direct8x8_block4x4_l0(h264d_mb_current *mb,  const h264d_vector_t *mvcol, const int8_t *ref_idx, h264d_vector_t *mv, int blk_idx)
{
	pred_direct_block<8, 4>(mb, mvcol, ref_idx, mv, blk_idx, pred_direct_col_block_onedir<0, 8, 4, 4>());
}

static void pred_direct8x8_block4x4_l1(h264d_mb_current *mb,  const h264d_vector_t *mvcol, const int8_t *ref_idx, h264d_vector_t *mv, int blk_idx)
{
	pred_direct_block<8, 4>(mb, mvcol, ref_idx, mv, blk_idx, pred_direct_col_block_onedir<1, 8, 4, 4>());
}

static void pred_direct8x8_block4x4_bidir(h264d_mb_current *mb,  const h264d_vector_t *mvcol, const int8_t *ref_idx, h264d_vector_t *mv, int blk_idx)
{
	pred_direct_block<8, 4>(mb, mvcol, ref_idx, mv, blk_idx, pred_direct_col_block_bidir<8, 4, 4>());
}

static void pred_direct8x8_spatial(h264d_mb_current *mb, int blk_idx, prev8x8_t *pblk)
{
	static void (* const pred_direct8x8_block_lut[4][4])(h264d_mb_current *mb, const h264d_vector_t *mvcol, const int8_t *ref_idx, h264d_vector_t *mv, int blk_idx) = {
		{
			pred_direct8x8_block8x8_nocol,
			pred_direct8x8_block8x8_l0,
			pred_direct8x8_block8x8_l1,
			pred_direct8x8_block8x8_bidir
		},
		{
			pred_direct8x8_block8x8_nocol,
			pred_direct8x8_block8x8_l0,
			pred_direct8x8_block8x8_l1,
			pred_direct8x8_block8x8_bidir
		},
		{
			pred_direct8x8_block8x8_nocol,
			pred_direct8x8_block8x8_l0,
			pred_direct8x8_block8x8_l1,
			pred_direct8x8_block8x8_bidir
		},
		{
			/* COL_MB8x8 */
			pred_direct8x8_block8x8_nocol,
			pred_direct8x8_block4x4_l0,
			pred_direct8x8_block4x4_l1,
			pred_direct8x8_block4x4_bidir
		}
	};
	int xoffset = (blk_idx & 1) * 8;
	int yoffset = (blk_idx & 2) * 4;
	pblk += blk_idx;
	int8_t *ref_idx = pblk->ref;
	if ((0 <= ref_idx[0]) || (0 <= ref_idx[1])) {
		const h264d_ref_frame_t *colpic = &(mb->frame->refs[1][0]);
		const h264d_col_mb_t *col_mb = &colpic->col->col_mb[mb->y * mb->max_x + mb->x];
		if ((colpic->in_use == SHORT_TERM) && (col_mb->ref[blk_idx] == 0)) {
			const h264d_vector_t *mvcol = &col_mb->mv[((unsigned)xoffset >> 2) + yoffset];
			int refs = (ref_idx[0] == 0) + (ref_idx[1] == 0) * 2;
			pred_direct8x8_block_lut[col_mb->type][refs](mb, mvcol, pblk->ref, pblk->mv[0], blk_idx);
		} else {
			direct_mv_pred(mb, ref_idx, &pblk->mv[0][0], 8, 8, xoffset, yoffset);
		}
	} else {
		pblk->ref[0] = 0;
		pblk->ref[1] = 0;
		memset(pblk->mv, 0, sizeof(pblk->mv));
		direct_zero_pred(mb, 8, 8, xoffset, yoffset);
	}
}

static void pred_direct8x8(h264d_mb_current *mb, int blk_idx, prev8x8_t *pblk)
{
	/* do nothing */
}

static void sub_mb8x8_direct(h264d_mb_current *mb, dec_bits *st, int avail, int blk_idx, prev8x8_t *pblk, int lx)
{
	if (lx == 0) {
		mb->bdirect->func->direct8x8(mb, blk_idx, pblk);
	}
}

template<typename F0>
static void sub_mb8x8_mv(h264d_mb_current *mb, dec_bits *st, int avail, int blk_idx, prev8x8_t *pblk, int lx,
			 F0 MvdXY)
{
	prev8x8_t *p = pblk + blk_idx;
	int idx = p->ref[lx];
	if (0 <= idx) {
		const int16_t *mvd_a;
		const int16_t *mvd_b;
		h264d_vector_t mv;
		h264d_vector_t mvd;

		calc_mv8x8_sub8x8(mb, mv.v, mvd_a, mvd_b, avail, lx, idx, blk_idx, pblk);
		MvdXY(mb, st, mvd.v, mvd_a, mvd_b);
		mv.v[0] += mvd.v[0];
		mv.v[1] += mvd.v[1];
		for (int i = 0; i < 4; ++i) {
			p->mv[i][lx].vector = mv.vector;
			p->mvd[i][lx].vector = mvd.vector;
		}
	}
}

template<typename F0>
static void sub_mb8x4_mv(h264d_mb_current *mb, dec_bits *st, int avail, int blk_idx, prev8x8_t *pblk, int lx,
			 F0 MvdXY)
{
	prev8x8_t *p = pblk + blk_idx;
	int idx = p->ref[lx];
	if (0 <= idx) {
		for (int y = 0; y < 2; ++y) {
			const int16_t *mvd_a;
			const int16_t *mvd_b;
			h264d_vector_t mv;
			h264d_vector_t mvd;

			calc_mv8x8_sub8x4(mb, mv.v, mvd_a, mvd_b, avail, lx, idx, blk_idx, pblk, y);
			MvdXY(mb, st, mvd.v, mvd_a, mvd_b);
			mv.v[0] += mvd.v[0];
			mv.v[1] += mvd.v[1];
			p->mv[y * 2][lx].vector = mv.vector;
			p->mvd[y * 2][lx].vector = mvd.vector;
			p->mv[y * 2 + 1][lx].vector = mv.vector;
			p->mvd[y * 2 + 1][lx].vector = mvd.vector;
		}
	}
}

template <typename F0>
static void sub_mb4x8_mv(h264d_mb_current *mb, dec_bits *st, int avail, int blk_idx, prev8x8_t *pblk, int lx,
			 F0 MvdXY)
{
	prev8x8_t *p = pblk + blk_idx;
	int idx = p->ref[lx];
	if (0 <= idx) {
		for (int x = 0; x < 2; ++x) {
			const int16_t *mvd_a;
			const int16_t *mvd_b;
			h264d_vector_t mv;
			h264d_vector_t mvd;

			calc_mv8x8_sub4x8(mb, mv.v, mvd_a, mvd_b, avail, lx, idx, blk_idx, pblk, x);
			MvdXY(mb, st, mvd.v, mvd_a, mvd_b);
			mv.v[0] += mvd.v[0];
			mv.v[1] += mvd.v[1];
			p->mv[x][lx].vector = mv.vector;
			p->mvd[x][lx].vector = mvd.vector;
			p->mv[x + 2][lx].vector = mv.vector;
			p->mvd[x + 2][lx].vector = mvd.vector;
		}
	}
}

template <typename F0>
static void sub_mb4x4_mv(h264d_mb_current *mb, dec_bits *st, int avail, int blk_idx, prev8x8_t *pblk, int lx,
			 F0 MvdXY)
{
	prev8x8_t *p = pblk + blk_idx;
	int idx = p->ref[lx];
	if (0 <= idx) {
		for (int xy = 0; xy < 4; ++xy) {
			const int16_t *mvd_a;
			const int16_t *mvd_b;
			h264d_vector_t mv;
			h264d_vector_t mvd;

			calc_mv8x8_sub4x4(mb, mv.v, mvd_a, mvd_b, avail, lx, idx, blk_idx, pblk, xy);
			MvdXY(mb, st, mvd.v, mvd_a, mvd_b);
			mv.v[0] += mvd.v[0];
			mv.v[1] += mvd.v[1];
			p->mv[xy][lx].vector = mv.vector;
			p->mvd[xy][lx].vector = mvd.vector;
		}
	}
}

struct mvd_xy_cavlc {
	void operator()(h264d_mb_current *mb, dec_bits *st, int16_t mv[], const int16_t mva[], const int16_t mvb[]) const {
	       mv[0] = se_golomb(st);
	       mv[1] = se_golomb(st);
	}
};

static void sub_mb8x8_mv_cavlc(h264d_mb_current *mb, dec_bits *st, int avail, int blk_idx, prev8x8_t *pblk, int lx)
{
	sub_mb8x8_mv(mb, st, avail, blk_idx, pblk, lx, mvd_xy_cavlc());
}

static void sub_mb8x4_mv_cavlc(h264d_mb_current *mb, dec_bits *st, int avail, int blk_idx, prev8x8_t *pblk, int lx)
{
	sub_mb8x4_mv(mb, st, avail, blk_idx, pblk, lx, mvd_xy_cavlc());
}

static void sub_mb4x8_mv_cavlc(h264d_mb_current *mb, dec_bits *st, int avail, int blk_idx, prev8x8_t *pblk, int lx)
{
	sub_mb4x8_mv(mb, st, avail, blk_idx, pblk, lx, mvd_xy_cavlc());
}

static void sub_mb4x4_mv_cavlc(h264d_mb_current *mb, dec_bits *st, int avail, int blk_idx, prev8x8_t *pblk, int lx)
{
	sub_mb4x4_mv(mb, st, avail, blk_idx, pblk, lx, mvd_xy_cavlc());
}

static void (* const sub_mb_p_cavlc[4])(h264d_mb_current *mb, dec_bits *st, int avail, int blk_idx, prev8x8_t *pblk, int lx) = {
	sub_mb8x8_mv_cavlc,
	sub_mb8x4_mv_cavlc,
	sub_mb4x8_mv_cavlc,
	sub_mb4x4_mv_cavlc,
};

static void (* const sub_mb_b_cavlc[13])(h264d_mb_current *mb, dec_bits *st, int avail, int blk_idx, prev8x8_t *pblk, int lx) = {
	sub_mb8x8_direct,
	sub_mb8x8_mv_cavlc,
	sub_mb8x8_mv_cavlc,
	sub_mb8x8_mv_cavlc,
	sub_mb8x4_mv_cavlc,
	sub_mb4x8_mv_cavlc,
	sub_mb8x4_mv_cavlc,
	sub_mb4x8_mv_cavlc,
	sub_mb8x4_mv_cavlc,
	sub_mb4x8_mv_cavlc,
	sub_mb4x4_mv_cavlc,
	sub_mb4x4_mv_cavlc,
	sub_mb4x4_mv_cavlc
};

struct sub_mbs_p_cavlc {
	void operator()(h264d_mb_current *mb, dec_bits *st, int avail, int8_t sub_mb_type[], prev8x8_t curr_blk[], int lx) const {
		if (lx == 0) {
			for (int i = 0; i < 4; ++i) {
				sub_mb_p_cavlc[sub_mb_type[i]](mb, st, avail, i, curr_blk, lx);
			}
		}
	}
};

struct sub_mbs_b_cavlc {
	void operator()(h264d_mb_current *mb, dec_bits *st, int avail, int8_t sub_mb_type[], prev8x8_t curr_blk[], int lx) const {
		for (int i = 0; i < 4; ++i) {
			sub_mb_b_cavlc[sub_mb_type[i]](mb, st, avail, i, curr_blk, lx);
		}
	}
};

static void sub_mb8x8_dec(h264d_mb_current *mb, int blk_idx, prev8x8_t *pblk)
{
	int bidir = 0;
	prev8x8_t *p = pblk + blk_idx;
	for (int lx = 0; lx < 2; ++lx) {
		int idx = p->ref[lx];
		if (0 <= idx) {
			inter_pred8x8(mb, p->mv[0][lx].v, 8, 8, mb->frame->refs[lx][idx].frame_idx, (blk_idx & 1) * 8, (blk_idx & 2) * 4, bidir++);
		}
	}
}

static void sub_mb8x4_dec(h264d_mb_current *mb, int blk_idx, prev8x8_t *pblk)

{
	prev8x8_t *p = pblk + blk_idx;
	int bidir = 0;
	for (int lx = 0; lx < 2; ++lx) {
		int idx = p->ref[lx];
		if (0 <= idx) {
			int frm_idx = mb->frame->refs[lx][idx].frame_idx;
			for (int y = 0; y < 2; ++y) {
				inter_pred8x8(mb, p->mv[y * 2][lx].v, 8, 4, frm_idx, (blk_idx & 1) * 8, ((blk_idx & 2) + y) * 4, bidir);
			}
			bidir++;
		}
	}
}

static void sub_mb4x8_dec(h264d_mb_current *mb, int blk_idx, prev8x8_t *pblk)
{
	prev8x8_t *p = pblk + blk_idx;
	int bidir = 0;
	for (int lx = 0; lx < 2; ++lx) {
		int idx = p->ref[lx];
		if (0 <= idx) {
			int frm_idx = mb->frame->refs[lx][idx].frame_idx;
			for (int x = 0; x < 2; ++x) {
				inter_pred8x8(mb, p->mv[x][lx].v, 4, 8, frm_idx, (blk_idx & 1) * 8 + x * 4, (blk_idx & 2) * 4, bidir);
			}
			bidir++;
		}
	}
}

static void sub_mb4x4_dec(h264d_mb_current *mb, int blk_idx, prev8x8_t *pblk)
{
	prev8x8_t *p = pblk + blk_idx;
	int bidir = 0;
	for (int lx = 0; lx < 2; ++lx) {
		int idx = p->ref[lx];
		if (0 <= idx) {
			int frm_idx = mb->frame->refs[lx][idx].frame_idx;
			for (int xy = 0; xy < 4; ++xy) {
				inter_pred8x8(mb, p->mv[xy][lx].v, 4, 4, frm_idx, (blk_idx & 1) * 8 + (xy & 1) * 4, (blk_idx & 2) * 4 + (xy & 2) * 2, bidir);
			}
			bidir++;
		}
	}
}

static void (* const sub_mb_dec_p[4])(h264d_mb_current *mb, int blk_idx, prev8x8_t *pblk) = {
	sub_mb8x8_dec,
	sub_mb8x4_dec,
	sub_mb4x8_dec,
	sub_mb4x4_dec
};

static void (* const sub_mb_dec_b[13])(h264d_mb_current *mb, int blk_idx, prev8x8_t *pblk) = {
	pred_direct8x8,
	sub_mb8x8_dec,
	sub_mb8x8_dec,
	sub_mb8x8_dec,
	sub_mb8x4_dec,
	sub_mb4x8_dec,
	sub_mb8x4_dec,
	sub_mb4x8_dec,
	sub_mb8x4_dec,
	sub_mb4x8_dec,
	sub_mb4x4_dec,
	sub_mb4x4_dec,
	sub_mb4x4_dec
};

struct sub_mbs_dec_p {
	void operator()(h264d_mb_current *mb, const int8_t sub_mb_type[], prev8x8_t curr_blk[], int avail) {
		for (int i = 0; i < 4; ++i) {
			sub_mb_dec_p[sub_mb_type[i]](mb, i, curr_blk);
		}
	}
};

struct sub_mbs_dec_b {
	void operator()(h264d_mb_current *mb, const int8_t sub_mb_type[], prev8x8_t curr_blk[], int avail) {
		for (int i = 0; i < 4; ++i) {
			sub_mb_dec_b[sub_mb_type[i]](mb, i, curr_blk);
		}
	}
};

template <int N>
static inline uint32_t str_mv_calc8x8_edge_bidir(uint32_t str, int ref0, int prev_ref0, int offset, const prev8x8_t *p, const prev_mb_t *prev)
{
	int lx = (ref0 != prev_ref0);
	for (int j = 0; j < 2; ++j) {
		if (((str & (2 << ((j + offset) * 2))) == 0)
			&& (DIF_ABS_LARGER_THAN4(p->mv[j * N][lx].v[0], prev->mov[j + offset].mv[0].v[0])
			|| DIF_ABS_LARGER_THAN4(p->mv[j * N][lx].v[1], prev->mov[j + offset].mv[0].v[1])
			|| DIF_ABS_LARGER_THAN4(p->mv[j * N][lx ^ 1].v[0], prev->mov[j + offset].mv[1].v[0])
			|| DIF_ABS_LARGER_THAN4(p->mv[j * N][lx ^ 1].v[1], prev->mov[j + offset].mv[1].v[1]))) {
			str = str | (1 << ((j + offset) * 2));
		}
	}
	return str;
}

template <int N>
static inline uint32_t str_mv_calc8x8_edge_onedir(uint32_t str, int ref0, int ref1, int prev_ref0, int offset, const prev8x8_t *p, const prev_mb_t *prev)
{
	int lx_s, lx_d;
	if (0 <= ref0) {
		lx_s = 0;
		lx_d = (ref0 != prev_ref0);
	} else {
		lx_s = 1;
		lx_d = (ref1 != prev_ref0);
	}
	for (int j = 0; j < 2; ++j) {
		if (((str & (2 << ((j + offset) * 2))) == 0)
			&& (DIF_ABS_LARGER_THAN4(p->mv[j * N][lx_s].v[0], prev->mov[j + offset].mv[lx_d].v[0])
			|| DIF_ABS_LARGER_THAN4(p->mv[j * N][lx_s].v[1], prev->mov[j + offset].mv[lx_d].v[1]))) {
			str = str | (1 << ((j + offset) * 2));
		}
	}
	return str;
}

template <int N>
static inline uint32_t str_mv_calc8x8_mv_edge(uint32_t str, int ref0, int ref1, int prev_ref0, int offset, const prev8x8_t *p, const prev_mb_t *prev)
{
	if ((0 <= ref0) && (0 <= ref1)) {
		return str_mv_calc8x8_edge_bidir<N>(str, ref0, prev_ref0, offset, p, prev);
	} else {
		return str_mv_calc8x8_edge_onedir<N>(str, ref0, ref1, prev_ref0, offset, p, prev);
	}
}

template <int N>
static uint32_t str_mv_calc8x8_edge(const h264d_mb_current *mb, uint32_t str, const prev8x8_t *p, const prev_mb_t *prev)
{
	for (int i = 0; i < 2; ++i) {
		uint32_t mask = 0xa << (i * 4);
		if ((str & mask) != mask) {
			int prev_ref0 = frame_idx_of_ref(mb, prev->ref[i][0], 0);
			int prev_ref1 = frame_idx_of_ref(mb, prev->ref[i][1], 1);
			int ref0 = frame_idx_of_ref(mb, p[i * N].ref[0], 0);
			int ref1 = frame_idx_of_ref(mb, p[i * N].ref[1], 1);
			if (((prev_ref0 != ref0) || (prev_ref1 != ref1)) && ((prev_ref1 != ref0) || (prev_ref0 != ref1))) {
				mask >>= 1;
				str |= (((str >> 1) ^ mask) & mask);
			} else {
				str = str_mv_calc8x8_mv_edge<N>(str, ref0, ref1, prev_ref0, i * 2, &p[i * N], prev);
			}
		}
	}
	return str;
}

template <int N>
static inline uint32_t str_mv_calc8x8_mid_bidir(uint32_t str, int offset, const prev8x8_t *p)
{
	for (int j = 0; j < 2; ++j) {
		if ((str & (2 << ((j + offset) * 2))) == 0) {
			if (DIF_ABS_LARGER_THAN4(p->mv[j * N][0].v[0], p->mv[j * N + (3 - N)][0].v[0])
				|| DIF_ABS_LARGER_THAN4(p->mv[j * N][0].v[1], p->mv[j * N + (3 - N)][0].v[1])
				|| DIF_ABS_LARGER_THAN4(p->mv[j * N][1].v[0], p->mv[j * N + (3 - N)][1].v[0])
				|| DIF_ABS_LARGER_THAN4(p->mv[j * N][1].v[1], p->mv[j * N + (3 - N)][1].v[1])) {
				str = str | (1 << ((j + offset) * 2));
			}
		}
	}
	return str;
}

template <int N>
static inline uint32_t str_mv_calc8x8_mid_onedir(uint32_t str, int lx, int offset, const prev8x8_t *p)
{
	for (int j = 0; j < 2; ++j) {
		if ((str & (2 << ((j + offset) * 2))) == 0) {
			if (DIF_ABS_LARGER_THAN4(p->mv[j * N][lx].v[0], p->mv[j * N + (3 - N)][lx].v[0])
				|| DIF_ABS_LARGER_THAN4(p->mv[j * N][lx].v[1], p->mv[j * N + (3 - N)][lx].v[1])) {
				str = str | (1 << ((j + offset) * 2));
			}
		}
	}
	return str;
}

template <int N>
static inline uint32_t str_mv_calc8x8_mv_mid(uint32_t str, int ref0, int ref1, int offset, const prev8x8_t *p)
{
	if ((0 <= ref0) && (0 <= ref1)) {
		return str_mv_calc8x8_mid_bidir<N>(str, offset, p);
	} else {
		return str_mv_calc8x8_mid_onedir<N>(str, (0 <= ref1), offset, p);
	}
}

template <int N>
static inline uint32_t str_mv_calc8x8_half_bidir(uint32_t str, int ref0, int prev_ref0, int offset, const prev8x8_t *p)
{
	int lx = (ref0 != prev_ref0);
	for (int j = 0; j < 2; ++j) {
		if ((str & (2 << ((j + offset) * 2))) == 0) {
			const int16_t *mv0 = p[0].mv[j * N + (3 - N)][0].v;
			const int16_t *mv1a = p[3 - N].mv[j * N][lx].v;
			const int16_t *mv1b = p[3 - N].mv[j * N][lx ^ 1].v;
			if (DIF_ABS_LARGER_THAN4(*mv0++, *mv1a++)
				|| DIF_ABS_LARGER_THAN4(*mv0++, *mv1a)
				|| DIF_ABS_LARGER_THAN4(*mv0++, *mv1b++)
				|| DIF_ABS_LARGER_THAN4(*mv0, *mv1b)) {
				str = str | (1 << ((j + offset) * 2));
			}
		}
	}
	return str;
}

template <int N>
static inline uint32_t str_mv_calc8x8_half_onedir(uint32_t str, int ref0, int ref1, int prev_ref0, int offset, const prev8x8_t *p)
{
	int lx_d, lx_s;
	if (0 <= ref0) {
		lx_d = 0;
		lx_s = (ref0 != prev_ref0);
	} else {
		lx_d = 1;
		lx_s = (ref1 != prev_ref0);
	}
	for (int j = 0; j < 2; ++j) {
		if ((str & (2 << ((j + offset) * 2))) == 0) {
			const int16_t *mv0 = p[0].mv[j * N + (3 - N)][lx_s].v;
			const int16_t *mv1 = p[3 - N].mv[j * N][lx_d].v;
			if (DIF_ABS_LARGER_THAN4(*mv0++, *mv1++) || DIF_ABS_LARGER_THAN4(*mv0, *mv1)) {
				str = str | (1 << ((j + offset) * 2));
			}
		}
	}
	return str;
}

template <int N>
static inline uint32_t str_mv_calc8x8_half_mv(uint32_t str, int ref0, int ref1, int prev_ref0, int offset, const prev8x8_t *p)
{
	if ((0 <= ref0) && (0 <= ref1)) {
		return str_mv_calc8x8_half_bidir<N>(str, ref0, prev_ref0, offset, p);
	} else {
		return str_mv_calc8x8_half_onedir<N>(str, ref0, ref1, prev_ref0, offset, p);
	}
}

template <int N>
static inline uint32_t str_mv_calc8x8_half(const h264d_mb_current *mb, uint32_t str, int offset, const prev8x8_t *p)
{
	int prev_ref0 = frame_idx_of_ref(mb, p[0].ref[0], 0);
	int prev_ref1 = frame_idx_of_ref(mb, p[0].ref[1], 1);
	int ref0 = frame_idx_of_ref(mb, p[3 - N].ref[0], 0);
	int ref1 = frame_idx_of_ref(mb, p[3 - N].ref[1], 1);
	if (((prev_ref0 != ref0) || (prev_ref1 != ref1)) && ((prev_ref1 != ref0) || (prev_ref0 != ref1))) {
		uint32_t mask = 5 << (offset * 2);
		str |= (((str >> 1) ^ mask) & mask);
	} else {
		str = str_mv_calc8x8_half_mv<N>(str, ref0, ref1, prev_ref0, offset, p);
	}
	return str;
}

template <int N>
static uint32_t str_mv_calc8x8_inner(const h264d_mb_current *mb, uint32_t str, const prev8x8_t *p)
{
	for (int i = 0; i < 2; ++i) {
		uint32_t mask = 0xa00 << (i * 4);
		if ((str & mask) != mask) {
			str = str_mv_calc8x8_mv_mid<N>(str, p[i * N].ref[0], p[i * N].ref[1], i * 2 + 4, p + i * N);
		}
	}
	for (int i = 0; i < 2; ++i) {
		uint32_t mask = 0xa0000 << (i * 4);
		if ((str & mask) != mask) {
			str = str_mv_calc8x8_half<N>(mb, str, i * 2 + 8, p + i * N);
		}
	}
	for (int i = 0; i < 2; ++i) {
		uint32_t mask = 0xa000000 << (i * 4);
		if ((str & mask) != mask) {
			str = str_mv_calc8x8_mv_mid<N>(str, p[i * N + (3 - N)].ref[0], p[i * N + (3 - N)].ref[1], i * 2 + 12, p + i * N + (3 - N));
		}
	}
	return str;
}

static inline void store_col8x8(h264d_col_mb_t *col_mb, const prev8x8_t *curr_blk)
{
	int8_t *refdst = col_mb->ref;
	h264d_vector_t *mvdst = col_mb->mv;

	col_mb->type = COL_MB8x8;
	for (int blk = 0; blk < 4; ++blk) {
		int refcol = curr_blk[blk].ref[0];
		const h264d_vector_t *mvcol;

		if (0 <= refcol) {
			mvcol = &curr_blk[blk].mv[0][0];
		} else {
			mvcol = &curr_blk[blk].mv[0][1];
			refcol = curr_blk[blk].ref[1];
		}
		refdst[blk] = refcol;
		mvdst[0].vector = mvcol[0].vector;
		mvdst[1].vector = mvcol[2].vector;
		mvdst[4].vector = mvcol[4].vector;
		mvdst[5].vector = mvcol[6].vector;
		mvdst = mvdst + ((blk & 1) * 4) + 2;
	}
}

struct store_direct8x8_info_b {
	void operator()(h264d_mb_current *mb, const int8_t *sub_mb_type) {
		mb->left4x4inter->direct8x8 = ((sub_mb_type[3] == 0) * 2) | (sub_mb_type[1] == 0);
		mb->top4x4inter->direct8x8 = ((sub_mb_type[3] == 0) * 2) | (sub_mb_type[2] == 0);
	}
};

struct store_direct8x8_info_p {
	void operator()(h264d_mb_current *mb, const int8_t *sub_by_type) {
		mb->left4x4inter->direct8x8 = 0;
		mb->top4x4inter->direct8x8 = 0;
	}
};

template <typename F0, typename F1, typename F2, typename F3, typename F4, typename F5, typename F6, typename F7>
static int mb_inter8x8(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail,
		       F0 SubMbTypes,
		       F1 RefIdx8x8,
		       F2 SubMbsMv,
		       F3 SubMbsDec,
		       F4 CodedBlockPattern,
		       F5 QpDelta,
		       F6 StoreDirect8x8Info,
		       F7 ResidualBlock)
{
	prev8x8_t curr_blk[4];
	int8_t sub_mb_type[4];
	uint32_t cbp;
	uint32_t str_vert, str_horiz;
	uint32_t left4x4, top4x4;

	memset(curr_blk, 0, sizeof(curr_blk));
	for (int i = 0; i < 4; ++i) {
		curr_blk[i].ref[0] = -1;
		curr_blk[i].ref[1] = -1;
	}
	if (SubMbTypes(mb, st, sub_mb_type, curr_blk, avail) < 0) {
		return -1;
	}
	for (int lx = 0; lx < 2; ++lx) {
		RefIdx8x8(mb, st, sub_mb_type, curr_blk, avail, lx);
	}
	for (int lx = 0; lx < 2; ++lx) {
		SubMbsMv(mb, st, avail, sub_mb_type, curr_blk, lx);
	}
	SubMbsDec(mb, sub_mb_type, curr_blk, avail);
	mb->cbp = cbp = CodedBlockPattern(mb, st, avail);
	left4x4 = mb->left4x4coef;
	top4x4 = *mb->top4x4coef;
	if (cbp) {
		str_vert = residual_luma_inter(mb, cbp, st, avail, QpDelta, ResidualBlock);
		str_horiz = expand_coef_str(transposition(str_vert));
		str_vert = expand_coef_str(str_vert);
	} else {
		mb->prev_qp_delta = 0;
		mb->left4x4coef = 0;
		*mb->top4x4coef = 0;
		str_vert = 0;
		str_horiz = 0;
	}
	deblock_info_t *deb = mb->deblock_curr;
	deb->qpy = mb->qp;
	deb->qpc = mb->qp_chroma;

	if (mb->y != 0) {
		if (mb->top4x4inter->type <= MB_IPCM) {
			deb->str4_vert = 1;
			str_vert |= 0xaa;
		} else {
			str_vert = str_previous_coef(str_vert, top4x4);
			str_vert = str_mv_calc8x8_edge<1>(mb, str_vert, curr_blk, mb->top4x4inter);
		}
	}
	deb->str_vert = str_mv_calc8x8_inner<1>(mb, str_vert, curr_blk);
	if (mb->x != 0) {
		if (mb->left4x4inter->type <= MB_IPCM) {
			deb->str4_horiz = 1;
			str_horiz |= 0xaa;
		} else {
			str_horiz = str_previous_coef(str_horiz, left4x4);
			str_horiz = str_mv_calc8x8_edge<2>(mb, str_horiz, curr_blk, mb->left4x4inter);
		}
	}
	deb->str_horiz = str_mv_calc8x8_inner<2>(mb, str_horiz, curr_blk);
	mb->left4x4pred = 0x22222222;
	*mb->top4x4pred = 0x22222222;
	for (int i = 0; i < 2; ++i) {
		mb->lefttop_mv[i].vector = mb->top4x4inter->mov[3].mv[i].vector;
		mb->top4x4inter->mov[0].mv[i].vector = curr_blk[2].mv[2][i].vector;
		mb->top4x4inter->mov[1].mv[i].vector = curr_blk[2].mv[3][i].vector;
		mb->top4x4inter->mov[2].mv[i].vector = curr_blk[3].mv[2][i].vector;
		mb->top4x4inter->mov[3].mv[i].vector = curr_blk[3].mv[3][i].vector;
		mb->top4x4inter->mvd[0].mv[i].vector = curr_blk[2].mvd[2][i].vector;
		mb->top4x4inter->mvd[1].mv[i].vector = curr_blk[2].mvd[3][i].vector;
		mb->top4x4inter->mvd[2].mv[i].vector = curr_blk[3].mvd[2][i].vector;
		mb->top4x4inter->mvd[3].mv[i].vector = curr_blk[3].mvd[3][i].vector;
		mb->lefttop_ref[i] = mb->top4x4inter->ref[1][i];
		mb->top4x4inter->ref[0][i] = curr_blk[2].ref[i];
		mb->top4x4inter->ref[1][i] = curr_blk[3].ref[i];
		mb->left4x4inter->ref[0][i] = curr_blk[1].ref[i];
		mb->left4x4inter->ref[1][i] = curr_blk[3].ref[i];
	}
	for (int i = 0; i < 4; ++i) {
		const prev8x8_t *p = &curr_blk[(i & 2) + 1];
		int idx = (i & 1) * 2 + 1;
		for (int j = 0; j < 2; ++j) {
			mb->left4x4inter->mov[i].mv[j].vector = p->mv[idx][j].vector;
			mb->left4x4inter->mvd[i].mv[j].vector = p->mvd[idx][j].vector;
		}
	}
	StoreDirect8x8Info(mb, sub_mb_type);
	store_col8x8(mb->col_curr, curr_blk);
	return residual_chroma(mb, cbp, st, avail, ResidualBlock);
}

template <int N, int IS_HORIZ>
static inline uint32_t str_mv_calc8x8_inner_onedir_center(const h264d_vector_set_t *mv, int top_lx, int bot_lx, uint32_t str, int shift)
{
	uint32_t mask = (N == 8 ? 0x000a0000U : 0x00020000U) << shift;
	const int16_t *mv_top = mv[IS_HORIZ ? 4 / N : 16 / N].mv[top_lx].v;
	const int16_t *mv_bot = mv[IS_HORIZ ? 8 / N : 16 * 2 / N].mv[bot_lx].v;
	uint32_t bits = 0U;
	for (int x = 0; x < 8 / N; ++x) {
		if ((str & mask) != mask) {
			if (DIF_ABS_LARGER_THAN4(mv_top[0], mv_bot[0])
				|| DIF_ABS_LARGER_THAN4(mv_top[1], mv_bot[1])) {
				bits |= (mask >> 1);
			}
		}
		mv_top += (IS_HORIZ) ? 16 * 4 / N : 4;
		mv_bot += (IS_HORIZ) ? 16 * 4 / N : 4;
		mask <<= 2;
	}
	return bits;
}

template <int N, int IS_HORIZ>
static inline uint32_t str_mv_calc4x4_inner_onedir(const h264d_vector_set_t *mv, int lx, int shift, uint32_t str)
{
	uint32_t bits = 0U;
	if (N <= 4) {
		uint32_t mask = 0x000200U << shift;
		const int16_t *mv_top = mv[0].mv[lx].v;
		const int16_t *mv_bot = mv[IS_HORIZ ? 4 / N : 16 / N].mv[lx].v;
		for (int x = 0; x < 2; ++x) {
			if ((str & mask) != mask) {
				if (DIF_ABS_LARGER_THAN4(mv_top[0], mv_bot[0])
					|| DIF_ABS_LARGER_THAN4(mv_top[1], mv_bot[1])) {
					bits |= (mask >> 1);
				}
			}
			mv_top += (IS_HORIZ) ? 16 * 4 / N : 4;
			mv_bot += (IS_HORIZ) ? 16 * 4 / N : 4;
			mask <<= 2;
		}
	}
	return bits;
}

template <int N, int IS_HORIZ>
static inline uint32_t str_mv_calc8x8_inner_onedir(const h264d_vector_set_t *mv, int top_lx, int bot_lx, uint32_t str, int shift)
{
	return str_mv_calc4x4_inner_onedir<N, IS_HORIZ>(mv + (IS_HORIZ ? 2 : 8), bot_lx, 16 + shift, str) | str_mv_calc8x8_inner_onedir_center<N, IS_HORIZ>(mv, top_lx, bot_lx, str, shift) | str_mv_calc4x4_inner_onedir<N, IS_HORIZ>(mv, top_lx, shift, str); 
}


template <int N, int IS_HORIZ>
static inline uint32_t str_mv_calc8x8_inner_bidir_center(const h264d_vector_set_t *mv, int lx, uint32_t str, int shift)
{
	uint32_t mask = (N == 8 ? 0x000a0000U : 0x00020000U) << shift;
	const int16_t *mv_top = mv[IS_HORIZ ? 4 / N : 16 / N].mv[0].v;
	const int16_t *mv_bot0 = mv[IS_HORIZ ? 8 / N : 16 * 2 / N].mv[lx].v;
	const int16_t *mv_bot1 = mv[IS_HORIZ ? 8 / N : 16 * 2 / N].mv[lx ^ 1].v;
	uint32_t bits = 0U;
	for (int x = 0; x < 8 / N; ++x) {
		if ((str & mask) != mask) {
			if (DIF_ABS_LARGER_THAN4(mv_top[0], mv_bot0[0])
				|| DIF_ABS_LARGER_THAN4(mv_top[1], mv_bot0[1])
				|| DIF_ABS_LARGER_THAN4(mv_top[2], mv_bot1[0])
				|| DIF_ABS_LARGER_THAN4(mv_top[3], mv_bot1[1])) {
				bits |= (mask >> 1);
			}
		}
		mv_top += (IS_HORIZ) ? 16 * 4 / N : 4;
		mv_bot0 += (IS_HORIZ) ? 16 * 4 / N : 4;
		mv_bot1 += (IS_HORIZ) ? 16 * 4 / N : 4;
		mask <<= 2;
	}
	return bits;
}

template <int N, int IS_HORIZ>
static inline uint32_t str_mv_calc4x4_inner_bidir(const h264d_vector_set_t *mv, int shift, uint32_t str)
{
	uint32_t bits = 0U;
	if (N <= 4) {
		uint32_t mask = 0x000200U << shift;
		const int16_t *mv_top = mv[0].mv[0].v;
		const int16_t *mv_bot = mv[IS_HORIZ ? 4 / N : 16 / N].mv[0].v;
		for (int x = 0; x < 2; ++x) {
			if ((str & mask) != mask) {
				if (DIF_ABS_LARGER_THAN4(mv_top[0], mv_bot[0])
					|| DIF_ABS_LARGER_THAN4(mv_top[1], mv_bot[1])
					|| DIF_ABS_LARGER_THAN4(mv_top[2], mv_bot[2])
					|| DIF_ABS_LARGER_THAN4(mv_top[3], mv_bot[3])) {
					bits |= (mask >> 1);
				}
			}
			mv_top += (IS_HORIZ) ? 16 * 4 / N : 4;
			mv_bot += (IS_HORIZ) ? 16 * 4 / N : 4;
			mask <<= 2;
		}
	}
	return bits;
}

template <int N, int IS_HORIZ>
static inline uint32_t str_mv_calc8x8_inner_bidir(const h264d_vector_set_t *mv, int lx, uint32_t str, int shift)
{
	return str_mv_calc4x4_inner_bidir<N, IS_HORIZ>(mv + (IS_HORIZ ? 2 : 8), 16 + shift, str) | str_mv_calc8x8_inner_bidir_center<N, IS_HORIZ>(mv, lx, str, shift) | str_mv_calc4x4_inner_bidir<N, IS_HORIZ>(mv, shift, str); 
}

template <int N, int IS_HORIZ>
static inline uint32_t str_mv_calc8x8_inner(const h264d_mb_current *mb, uint32_t str, const int8_t ref_idx[], const h264d_vector_set_t *mv)
{
	uint32_t mask = 0U;
	for (int x = 0; x < 2; ++x) {
		int top_ref0 = frame_idx_of_ref(mb, ref_idx[0], 0);
		int top_ref1 = frame_idx_of_ref(mb, ref_idx[1], 1);
		int bot_ref0 = frame_idx_of_ref(mb, ref_idx[IS_HORIZ ? 2 : 4], 0);
		int bot_ref1 = frame_idx_of_ref(mb, ref_idx[IS_HORIZ ? 3 : 5], 1);
		int shift = x * 4;
		uint32_t str_bits;
		if ((top_ref0 != bot_ref0 || top_ref1 != bot_ref1) && (top_ref1 != bot_ref0 || top_ref0 != bot_ref1)) {
			if ((0 <= top_ref0) && (0 <= top_ref1)) {
				str_bits = str_mv_calc4x4_inner_bidir<N, IS_HORIZ>(mv + (IS_HORIZ ? 2 : 8), 16 + shift, str) | str_mv_calc4x4_inner_bidir<N, IS_HORIZ>(mv, shift, str) | (0x00050000U << shift);
			} else {
				str_bits = str_mv_calc4x4_inner_onedir<N, IS_HORIZ>(mv + (IS_HORIZ ? 2 : 8), 0 <= bot_ref1, 16 + shift, str) | str_mv_calc4x4_inner_onedir<N, IS_HORIZ>(mv, 0 <= top_ref1, shift, str) | (0x00050000U << shift);
			}
		} else {
			if ((0 <= top_ref0) && (0 <= top_ref1)) {
				str_bits = str_mv_calc8x8_inner_bidir<N, IS_HORIZ>(mv, top_ref0 != bot_ref0, str, shift);
			} else {
				str_bits = str_mv_calc8x8_inner_onedir<N, IS_HORIZ>(mv, top_ref0 < 0, bot_ref0 < 0, str, shift);
			}
		}
		mask |= str_bits;
		ref_idx += IS_HORIZ ? 4 : 2;
		mv += IS_HORIZ ? (64 * 2 / (N * N)) : (8 / N);
	}
	return str | (((str >> 1) ^ mask) & mask);
}

template <int N>
static void store_info_inter8x8(h264d_mb_current *mb, const h264d_vector_set_t mv[], const int8_t ref_idx[], uint32_t str_vert, uint32_t str_horiz, uint32_t left4x4, uint32_t top4x4)
{
	deblock_info_t *deb = mb->deblock_curr;
	deb->qpy = mb->qp;
	deb->qpc = mb->qp_chroma;
	if (mb->y != 0) {
		str_vert = store_str_inter8xedge<0, 1, 8 / N>(mb, mb->top4x4inter, deb->str4_vert, mv, ref_idx, str_vert, top4x4);
	}
	deb->str_vert = str_mv_calc8x8_inner<N, 0>(mb, str_vert, ref_idx, mv);
	if (mb->x != 0) {
		str_horiz = store_str_inter8xedge<1, 4, 8 / N>(mb, mb->left4x4inter, deb->str4_horiz, mv, ref_idx, str_horiz, left4x4);
	}
	deb->str_horiz = str_mv_calc8x8_inner<N, 1>(mb, str_horiz, ref_idx, mv);

	mb->left4x4pred = 0x22222222;
	*mb->top4x4pred = 0x22222222;
	for (int i = 0; i < 2; ++i) {
		mb->lefttop_ref[i] = mb->top4x4inter->ref[1][i];
		mb->lefttop_mv[i].vector = mb->top4x4inter->mov[3].mv[i].vector;
		mb->top4x4inter->ref[i][0] = ref_idx[i * 2 + 4];
		mb->top4x4inter->ref[i][1] = ref_idx[i * 2 + 5];
		mb->left4x4inter->ref[i][0] = ref_idx[i * 4 + 2];
		mb->left4x4inter->ref[i][1] = ref_idx[i * 4 + 3];
	}
	for (int i = 0; i < 4; ++i) {
		mb->top4x4inter->mov[i] = mv[i + 12];
		mb->left4x4inter->mov[i] = mv[i * 4 + 3];
	}
	memset(mb->top4x4inter->mvd, 0, sizeof(mb->top4x4inter->mvd));
	memset(mb->left4x4inter->mvd, 0, sizeof(mb->left4x4inter->mvd));

	h264d_col_mb_t *col_mb = mb->col_curr;
	int8_t *refdst = col_mb->ref;
	h264d_vector_t *mvdst = col_mb->mv;

	col_mb->type = COL_MB8x8;
	for (int blk = 0; blk < 4; ++blk) {
		int refcol = ref_idx[blk * 2];
		int lx;

		if (0 <= refcol) {
			lx = 0;
		} else {
			lx = 1;
			refcol = ref_idx[blk * 2 + 1];
		}
		refdst[blk] = refcol;
		mvdst[0].vector = mv[0].mv[lx].vector;
		mvdst[1].vector = mv[1].mv[lx].vector;
		mvdst[4].vector = mv[4].mv[lx].vector;
		mvdst[5].vector = mv[5].mv[lx].vector;
		if (blk & 1) {
			mvdst += 8;
			mv += 8;
		} else {
			mvdst += 2;
			mv += 2;
		}
	}
}

static void (* const store_info_inter[2][4])(h264d_mb_current *mb, const h264d_vector_set_t mv[], const int8_t ref_idx[], uint32_t str_vert, uint32_t str_horiz, uint32_t left4x4, uint32_t top4x4) = {
	{
		store_info_inter16x16,
		store_info_inter16x8,
		store_info_inter8x16,
		store_info_inter8x8<4>
	},
	{
		store_info_inter16x16,
		store_info_inter16x8,
		store_info_inter8x16,
		store_info_inter8x8<8>
	},
};

template <typename F0, typename F1, typename F2>
static int mb_bdirect16x16(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail,
			   F0 CodedBlockPattern,
			   F1 QpDelta,
			   F2 ResidualBlock)
{
	h264d_vector_set_t mv[16 * 2];
	int8_t ref_idx[2 * 4];
	uint32_t str_vert, str_horiz;
	uint32_t left4x4, top4x4;
	uint32_t cbp;

	mb->bdirect->func->direct16x16(mb, ref_idx, mv);
	left4x4 = mb->left4x4coef;
	top4x4 = *mb->top4x4coef;
	mb->cbp = cbp = CodedBlockPattern(mb, st, avail);
	if (cbp) {
		str_vert = residual_luma_inter(mb, cbp, st, avail, QpDelta, ResidualBlock);
		str_horiz = expand_coef_str(transposition(str_vert));
		str_vert = expand_coef_str(str_vert);
	} else {
		mb->prev_qp_delta = 0;
		mb->left4x4coef = 0;
		*mb->top4x4coef = 0;
		str_vert = 0;
		str_horiz = 0;
	}
	const h264d_ref_frame_t *colpic = &(mb->frame->refs[1][0]);
	const h264d_col_mb_t *col_mb = &colpic->col->col_mb[mb->y * mb->max_x + mb->x];
	store_info_inter[0][col_mb->type](mb, mv, ref_idx, str_vert, str_horiz, left4x4, top4x4);
	mb->left4x4inter->direct8x8 = 3;
	mb->top4x4inter->direct8x8 = 3;
	return residual_chroma(mb, cbp, st, avail, ResidualBlock);
}

struct sub_mb_type_p_cavlc {
	int operator()(h264d_mb_current *mb, dec_bits *st, int8_t *sub_mb_type, prev8x8_t *curr_blk, int avail) const {
		for (int i = 0; i < 4; ++i) {
			READ_UE_RANGE(sub_mb_type[i], st, 3);
		}
		return 0;
	}
};


static inline void fill_direct8x8_mv(prev8x8_t *pblk, const prev8x8_t *src)
{
	uint32_t mov0 = src->mv[0][0].vector;
	uint32_t mov1 = src->mv[0][1].vector;
	pblk->ref[0] = src->ref[0];
	pblk->ref[1] = src->ref[1];
	for (int i = 0; i < 4; ++i) {
		pblk->mv[i][0].vector = mov0;
		pblk->mv[i][1].vector = mov1;
	}
}

template <typename F>
static inline int sub_mb_type_b_base(h264d_mb_current *mb, dec_bits *st, int8_t sub_mb_type[], prev8x8_t curr_blk[], int avail,
				      F SubMbTypeB)
{
	const prev8x8_t *ref_blk;
	int b_direct_cnt = 0;

	for (int i = 0; i < 4; ++i) {
		int type = SubMbTypeB(mb, st);
		if (type < 0) {
			return -1;
		}
		sub_mb_type[i] = type;
		if (type == 0) {
			if (b_direct_cnt++ == 0) {
				b_direct_ref_mv_calc(mb, avail, curr_blk[i].ref, curr_blk[i].mv[0][0].v);
				ref_blk = &curr_blk[i];
			}
			fill_direct8x8_mv(&curr_blk[i], ref_blk);
		}
	}
	return 0;
}


struct sub_mb_type_b_cavlc {
	int operator()(h264d_mb_current *mb, dec_bits *st) const {
		int sub_mb_type;
		READ_UE_RANGE(sub_mb_type, st, 12);
		return sub_mb_type;
	}
};

struct sub_mb_types_b_cavlc {
	int operator()(h264d_mb_current *mb, dec_bits *st, int8_t *sub_mb_type, prev8x8_t *curr_blk, int avail) const {
		return sub_mb_type_b_base(mb, st, sub_mb_type, curr_blk, avail, sub_mb_type_b_cavlc());
	}
};

struct ref_idx16x16_cavlc {
	 int operator()(h264d_mb_current *mb, dec_bits *st, int lx, int avail) const {
		int t = *(mb->num_ref_idx_lx_active_minus1[lx]);
		return t ? te_golomb(st, t) : 0;
	}
};

struct ref_idx16x8_cavlc {
	void operator()(h264d_mb_current *mb, dec_bits *st, int8_t *ref_idx, uint32_t blk_map, int avail) const {
		int8_t * const *num = mb->num_ref_idx_lx_active_minus1;
		for (int lx = 0; lx < 2; ++lx) {
			int t = *(num[0]);
			ref_idx[0] = (blk_map & 1) ? (t ? te_golomb(st, t) : 0) : -1;
			ref_idx[2] = (blk_map & 2) ? (t ? te_golomb(st, t) : 0) : -1;
			blk_map >>= 2;
			ref_idx++;
			num++;
		}
	}
};

struct ref_idx8x8_cavlc {
	void operator()(h264d_mb_current *mb, dec_bits *st, const int8_t *sub_mb_type, prev8x8_t *pblk, int avail, int lx) const {
		int t = (mb->type != MB_P8x8REF0) ? *(mb->num_ref_idx_lx_active_minus1[lx]) : 0;
		int dir = 1 << lx;
		const int8_t *sub_mb_ref_map = mb->sub_mb_ref_map;
		for (int i = 0; i < 4; ++i) {
			int sub_dir = sub_mb_ref_map[*sub_mb_type++];
			if (0 <= sub_dir) {
				pblk[i].ref[lx] = (dir & sub_dir) ? (t ? te_golomb(st, t) : 0) : -1;
			}
		}
	}
};

static int mb_inter16x16_cavlc(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_inter16x16(mb, mbc, st, avail, ref_idx16x16_cavlc(), mvd_xy_cavlc(), cbp_inter_cavlc(), qp_delta_cavlc(), residual_block_cavlc());
}

static int mb_inter16x8_cavlc(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_inter16x8(mb, mbc, st, avail, ref_idx16x8_cavlc(), mvd_xy_cavlc(), cbp_inter_cavlc(), qp_delta_cavlc(), residual_block_cavlc());
}

static int mb_inter8x16_cavlc(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_inter8x16(mb, mbc, st, avail, ref_idx16x8_cavlc(), mvd_xy_cavlc(), cbp_inter_cavlc(), qp_delta_cavlc(), residual_block_cavlc());
}

static int mb_inter8x8p_cavlc(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_inter8x8(mb, mbc, st, avail, sub_mb_type_p_cavlc(), ref_idx8x8_cavlc(), sub_mbs_p_cavlc(), sub_mbs_dec_p(), cbp_inter_cavlc(), qp_delta_cavlc(), store_direct8x8_info_p(), residual_block_cavlc());
}

static int mb_inter8x8b_cavlc(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_inter8x8(mb, mbc, st, avail, sub_mb_types_b_cavlc(), ref_idx8x8_cavlc(), sub_mbs_b_cavlc(), sub_mbs_dec_b(), cbp_inter_cavlc(), qp_delta_cavlc(), store_direct8x8_info_b(), residual_block_cavlc());
}

static int mb_bdirect16x16_cavlc(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_bdirect16x16(mb, mbc, st, avail, cbp_inter_cavlc(), qp_delta_cavlc(), residual_block_cavlc());
}

static const mb_code mb_decode[] = {
	{mb_intra4x4_cavlc, 0, 0},
	{mb_intra16x16_dconly_cavlc, mb_intra16xpred_vert<16>, 0},
	{mb_intra16x16_dconly_cavlc, mb_intra16x16pred_horiz, 0},
	{mb_intra16x16_dconly_cavlc, mb_intra16x16pred_dc, 0},
	{mb_intra16x16_dconly_cavlc, mb_intra16x16pred_planer, 0},
	{mb_intra16x16_dconly_cavlc, mb_intra16xpred_vert<16>, 0x10},
	{mb_intra16x16_dconly_cavlc, mb_intra16x16pred_horiz, 0x10},
	{mb_intra16x16_dconly_cavlc, mb_intra16x16pred_dc, 0x10},
	{mb_intra16x16_dconly_cavlc, mb_intra16x16pred_planer, 0x10},
	{mb_intra16x16_dconly_cavlc, mb_intra16xpred_vert<16>, 0x20},
	{mb_intra16x16_dconly_cavlc, mb_intra16x16pred_horiz, 0x20},
	{mb_intra16x16_dconly_cavlc, mb_intra16x16pred_dc, 0x20},
	{mb_intra16x16_dconly_cavlc, mb_intra16x16pred_planer, 0x20},
	{mb_intra16x16_acdc_cavlc, mb_intra16xpred_vert<16>, 0x0f},
	{mb_intra16x16_acdc_cavlc, mb_intra16x16pred_horiz, 0x0f},
	{mb_intra16x16_acdc_cavlc, mb_intra16x16pred_dc, 0x0f},
	{mb_intra16x16_acdc_cavlc, mb_intra16x16pred_planer, 0x0f},
	{mb_intra16x16_acdc_cavlc, mb_intra16xpred_vert<16>, 0x1f},
	{mb_intra16x16_acdc_cavlc, mb_intra16x16pred_horiz, 0x1f},
	{mb_intra16x16_acdc_cavlc, mb_intra16x16pred_dc, 0x1f},
	{mb_intra16x16_acdc_cavlc, mb_intra16x16pred_planer, 0x1f},
	{mb_intra16x16_acdc_cavlc, mb_intra16xpred_vert<16>, 0x2f},
	{mb_intra16x16_acdc_cavlc, mb_intra16x16pred_horiz, 0x2f},
	{mb_intra16x16_acdc_cavlc, mb_intra16x16pred_dc, 0x2f},
	{mb_intra16x16_acdc_cavlc, mb_intra16x16pred_planer, 0x2f},
	{mb_intrapcm, 0, 0},
	{mb_inter16x16_cavlc, 0, 1},
	{mb_inter16x8_cavlc, 0, 3},
	{mb_inter8x16_cavlc, 0, 3},
	{mb_inter8x8p_cavlc, 0, 0xf},
	{mb_inter8x8p_cavlc, 0, 0xf},
	{mb_bdirect16x16_cavlc, 0, 0},
	{mb_inter16x16_cavlc, 0, 1}, {mb_inter16x16_cavlc, 0, 2}, {mb_inter16x16_cavlc, 0, 3},
	{mb_inter16x8_cavlc, 0, 0x3}, {mb_inter8x16_cavlc, 0, 0x3},
	{mb_inter16x8_cavlc, 0, 0xc}, {mb_inter8x16_cavlc, 0, 0xc},
	{mb_inter16x8_cavlc, 0, 0x9}, {mb_inter8x16_cavlc, 0, 0x9},
	{mb_inter16x8_cavlc, 0, 0x6}, {mb_inter8x16_cavlc, 0, 0x6},
	{mb_inter16x8_cavlc, 0, 0xb}, {mb_inter8x16_cavlc, 0, 0xb},
	{mb_inter16x8_cavlc, 0, 0xe}, {mb_inter8x16_cavlc, 0, 0xe},
	{mb_inter16x8_cavlc, 0, 0x7}, {mb_inter8x16_cavlc, 0, 0x7},
	{mb_inter16x8_cavlc, 0, 0xd}, {mb_inter8x16_cavlc, 0, 0xd},
	{mb_inter16x8_cavlc, 0, 0xf}, {mb_inter8x16_cavlc, 0, 0xf},
	{mb_inter8x8b_cavlc, 0, 0}
};

/** Convert MB type number into unified order:
 * Intra < Inter < Bidirectional
 */
static int adjust_mb_type(int mb_type, int slice_type)
{
	if (slice_type == P_SLICE) {
		if (mb_type < 30) {
			mb_type -= 5;
			return mb_type < 0 ? mb_type + MB_BDIRECT16x16 : mb_type;
		} else {
			return -1;
		}
	} else if (slice_type == B_SLICE) {
		mb_type -= 23;
		return mb_type < 0 ? mb_type + 23 + MB_BDIRECT16x16 : mb_type;
	} else if ((slice_type == I_SLICE) && (mb_type <= 25)) {
		return  mb_type;
	} else {
		return -1;
	}
}

static inline int get_availability(h264d_mb_current *mb)
{
	int mbx, max_x, firstline;

	mbx = mb->x;
	max_x = mb->max_x;
	firstline = mb->firstline;
	return ((mbx != 0 && firstline < 0) * 8) /* bit3: top left */
		| ((mbx != max_x - 1 && firstline <= 1) * 4) /* bit2: top right */
		| ((firstline <= 0) * 2) /* bit1: top */
		| (mbx != 0 && firstline != max_x); /* bit0: left */
}

static inline int macroblock_layer_cabac(h264d_mb_current *mb, int slice_type, dec_bits *st);

static inline int macroblock_layer(h264d_mb_current *mb, int slice_type, dec_bits *st)
{
	const mb_code *mbc;
	int mbtype;
	int avail;

	READ_UE_RANGE(mbtype, st, 48);
	if ((mb->type = mbtype = adjust_mb_type(mbtype, slice_type)) < 0) {
		return -1;
	}
	mbc = &mb_decode[mbtype];
	avail = get_availability(mb);
	mbc->mb_dec(mb, mbc, st, avail);
	VC_CHECK;
	return 0;
}

static void calc_mv_pskip(h264d_mb_current *mb, int16_t mv[], int avail)
{
	prev_mb_t *pmb;
	int16_t pmv[2];
	const int16_t *mvd_a, *mvd_b; /* ignored */

	mv[0] = 0;
	mv[1] = 0;
	if ((avail & 3) != 3) {
		return;
	}
	pmb = mb->left4x4inter;
	if (pmb->ref[0][0] == 0 && pmb->mov[0].mv[0].vector == 0U) {
		return;
	}
	pmb = mb->top4x4inter;
	if (pmb->ref[0][0] == 0 && pmb->mov[0].mv[0].vector == 0U) {
		return;
	}
	calc_mv16x16(mb, pmv, mvd_a, mvd_b, 0, 0, avail);
	mv[0] = pmv[0];
	mv[1] = pmv[1];
}

static void p_skip_mb(h264d_mb_current *mb, int8_t *ref_idx, h264d_vector_set_t *mv)
{
	calc_mv_pskip(mb, mv->mv->v, get_availability(mb));
	memset(&mv[1], 0, sizeof(mv[1]));
	inter_pred8x8(mb, mv->mv->v, 16, 16, mb->frame->refs[0][0].frame_idx, 0, 0, 0);
}

template <int BLOCK>
static inline void fill_bskip_mv(h264d_vector_t *mv)
{
	uint32_t *mv32 = &(mv->vector);
	uint32_t d0 = mv32[0];
	uint32_t d1 = mv32[1];
	int i = (16 * 4 / BLOCK) - 1;
	do {
		mv32 += 2;
		mv32[0] = d0;
		mv32[1] = d1;
	} while (--i);
}

static void direct_mv_pred_nocol(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv)
{
	direct_mv_pred(mb, ref_idx, mv, 16, 16, 0, 0);
	((h264d_col_mb_t *)col_mb)->type = COL_MB16x16;
	memset(&mv[2], 0, sizeof(mv[2]) * 2);
}

template <typename F>
static inline void pred_direct16x16_col_base16x16(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv,
						      F PredDirectCol)
{
	if (col_mb->ref[0] == 0) {
		const h264d_vector_t *mvcol = &col_mb->mv[0];
		pred_direct_block<16, 8>(mb, mvcol, ref_idx, mv, 0, PredDirectCol);
	} else {
		direct_mv_pred(mb, ref_idx, mv, 16, 16, 0, 0);
	}
	memset(&mv[2], 0, sizeof(mv[2]) * 2);
}

template <typename F>
static inline void pred_direct16x16_col_base16x8(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv,
						     F PredDirectCol)
{
	memcpy(&mv[2], &mv[0], sizeof(mv[0]) * 2);
	for (int y = 0; y < 2; ++y) {
		if (col_mb->ref[y * 2] == 0) {
			const h264d_vector_t *mvcol = &col_mb->mv[y * 8];
			pred_direct_block<16, 8>(mb, mvcol, ref_idx, mv, y * 2, PredDirectCol);
		} else {
			direct_mv_pred(mb, ref_idx, mv, 16, 8, 0, y * 8);
		}
		mv += 2;
	}
	memset(&mv[0], 0, sizeof(mv[0]) * 4);
}

template <typename F>
static inline void pred_direct16x16_col_base8x16(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv,
						     F PredDirectCol)
{
	memcpy(&mv[2], &mv[0], sizeof(mv[0]) * 2);
	for (int x = 0; x < 2; ++x) {
		if (col_mb->ref[x] == 0) {
			const h264d_vector_t *mvcol = &col_mb->mv[x * 2];
			pred_direct_block<16, 8>(mb, mvcol, ref_idx, mv, x, PredDirectCol);
		} else {
			direct_mv_pred(mb, ref_idx, mv, 8, 16, x * 8, 0);
		}
		mv += 2;
	}
	memset(&mv[0], 0, sizeof(mv[0]) * 4);
}

template <int BLOCK, typename F>
static inline void pred_direct16x16_col_base8x8(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv,
					    F PredDirectCol)
{
	fill_bskip_mv<BLOCK>(mv);
	for (int blk8x8 = 0; blk8x8 < 4; ++blk8x8) {
		int yoffset = (blk8x8 & 2) * 4;
		if (col_mb->ref[blk8x8] == 0) {
			const h264d_vector_t *mvcol = &col_mb->mv[(blk8x8 & 1) * 2 + yoffset];
			pred_direct_block<16, BLOCK>(mb, mvcol, ref_idx, mv, blk8x8, PredDirectCol);
		} else {
			direct_mv_pred(mb, ref_idx, mv, 8, 8, (blk8x8 & 1) * 8, yoffset);
		}
		mv += (blk8x8 & 1) ? 16 * 4 / BLOCK - 16 / BLOCK : 16 / BLOCK;
	}
	memset(&mv[0], 0, sizeof(mv[0]) * 16 * 8 / BLOCK);
}

static void pred_direct16x16_col_ref1_16x16(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv)
{
	pred_direct16x16_col_base16x16(mb, col_mb, ref_idx, mv, pred_direct_col_block_onedir<0, 16, 16, 16>());
}

static void pred_direct16x16_col_ref2_16x16(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv)
{
	pred_direct16x16_col_base16x16(mb, col_mb, ref_idx, mv, pred_direct_col_block_onedir<1, 16, 16, 16>());
}

static void pred_direct16x16_col_ref3_16x16(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv)
{
	pred_direct16x16_col_base16x16(mb, col_mb, ref_idx, mv, pred_direct_col_block_bidir<16, 16, 16>());
}

static void pred_direct16x16_col_ref1_16x8(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv)
{
	pred_direct16x16_col_base16x8(mb, col_mb, ref_idx, mv, pred_direct_col_block_onedir<0, 16, 16, 8>());
}

static void pred_direct16x16_col_ref2_16x8(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv)
{
	pred_direct16x16_col_base16x8(mb, col_mb, ref_idx, mv, pred_direct_col_block_onedir<1, 16, 16, 8>());
}

static void pred_direct16x16_col_ref3_16x8(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv)
{
	pred_direct16x16_col_base16x8(mb, col_mb, ref_idx, mv, pred_direct_col_block_bidir<16, 16, 8>());
}

static void pred_direct16x16_col_ref1_8x16(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv)
{
	pred_direct16x16_col_base8x16(mb, col_mb, ref_idx, mv, pred_direct_col_block_onedir<0, 16, 8, 16>());
}

static void pred_direct16x16_col_ref2_8x16(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv)
{
	pred_direct16x16_col_base8x16(mb, col_mb, ref_idx, mv, pred_direct_col_block_onedir<1, 16, 8, 16>());
}

static void pred_direct16x16_col_ref3_8x16(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv)
{
	pred_direct16x16_col_base8x16(mb, col_mb, ref_idx, mv, pred_direct_col_block_bidir<16, 8, 16>());
}

static void pred_direct16x16_col_ref1_8x8(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv)
{
	pred_direct16x16_col_base8x8<8>(mb, col_mb, ref_idx, mv, pred_direct_col_block_onedir<0, 16, 8, 8>());
}

static void pred_direct16x16_col_ref2_8x8(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv)
{
	pred_direct16x16_col_base8x8<8>(mb, col_mb, ref_idx, mv, pred_direct_col_block_onedir<1, 16, 8, 8>());
}

static void pred_direct16x16_col_ref3_8x8(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv)
{
	pred_direct16x16_col_base8x8<8>(mb, col_mb, ref_idx, mv, pred_direct_col_block_bidir<16, 8, 8>());
}

static void pred_direct16x16_col_ref1_4x4(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv)
{
	pred_direct16x16_col_base8x8<4>(mb, col_mb, ref_idx, mv, pred_direct_col_block_onedir<0, 16, 4, 4>());
}

static void pred_direct16x16_col_ref2_4x4(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv)
{
	pred_direct16x16_col_base8x8<4>(mb, col_mb, ref_idx, mv, pred_direct_col_block_onedir<1, 16, 4, 4>());
}

static void pred_direct16x16_col_ref3_4x4(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv)
{
	pred_direct16x16_col_base8x8<4>(mb, col_mb, ref_idx, mv, pred_direct_col_block_bidir<16, 4, 4>());
}

static void (* const pred_direct16x16_col[5][4])(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, const int8_t *ref_idx, h264d_vector_t *mv) = {
	{
		direct_mv_pred_nocol,
		pred_direct16x16_col_ref1_16x16,
		pred_direct16x16_col_ref2_16x16,
		pred_direct16x16_col_ref3_16x16
	},
	{
		direct_mv_pred_nocol,
		pred_direct16x16_col_ref1_16x8,
		pred_direct16x16_col_ref2_16x8,
		pred_direct16x16_col_ref3_16x8
	},
	{
		direct_mv_pred_nocol,
		pred_direct16x16_col_ref1_8x16,
		pred_direct16x16_col_ref2_8x16,
		pred_direct16x16_col_ref3_8x16
	},
	{
		direct_mv_pred_nocol,
		pred_direct16x16_col_ref1_4x4,
		pred_direct16x16_col_ref2_4x4,
		pred_direct16x16_col_ref3_4x4
	},
	{
		direct_mv_pred_nocol,
		pred_direct16x16_col_ref1_8x8,
		pred_direct16x16_col_ref2_8x8,
		pred_direct16x16_col_ref3_8x8
	}
};

static void pred_direct16x16(h264d_mb_current *mb, int8_t *ref_idx, h264d_vector_t *mv)
{
	h264d_ref_frame_t *colpic = &(mb->frame->refs[1][0]);
	h264d_col_mb_t *col_mb = &colpic->col->col_mb[mb->y * mb->max_x + mb->x];
	if ((0 <= ref_idx[0]) || (0 <= ref_idx[1])) {
		if (colpic->in_use == SHORT_TERM) {
			int refs = (ref_idx[0] == 0) + (ref_idx[1] == 0) * 2;
			pred_direct16x16_col[col_mb->type][refs](mb, col_mb, ref_idx, mv);
		} else {
			col_mb->type = COL_MB16x16;
			memset(&mv[2], 0, sizeof(mv[2]) * 2);
			direct_mv_pred(mb, ref_idx, mv, 16, 16, 0, 0);
		}
	} else {
		ref_idx[0] = 0;
		ref_idx[1] = 0;
		col_mb->type = COL_MB16x16;
		memset(&mv[2], 0, sizeof(mv[2]) * 2);
		direct_zero_pred(mb, 16, 16, 0, 0);
	}
}

static void b_skip_mb_spatial(h264d_mb_current *mb, int8_t *ref_idx, h264d_vector_set_t *mv)
{
	int avail = get_availability(mb);
	b_direct_ref_mv_calc(mb, avail, ref_idx, mv->mv->v);
	for (int i = 1; i < 4; ++i) {
		ref_idx[i * 2] = ref_idx[0];
		ref_idx[i * 2 + 1] = ref_idx[1];
	}
	pred_direct16x16(mb, ref_idx, mv->mv);
}

template <int N, int BLOCK>
static inline void fill_temporal_mv(h264d_vector_t *mv)
{
	if (4 < BLOCK) { 
		uint32_t xy0 = mv[0].vector;
		uint32_t xy1 = mv[1].vector;
		for (int y = 0; y < BLOCK / 4; ++y) {
			for (int x = 0; x < BLOCK / 4; ++x) {
				mv[x * 2].vector = xy0;
				mv[x * 2 + 1].vector = xy1;
			}
			mv += N / 2;
		}
	}
}

struct tempral_vector_zero {
	void operator()(const h264d_vector_t *mvcol, int scale, h264d_vector_t *tmv) const {
		tmv[0].vector = mvcol->vector;
		tmv[1].vector = 0U;
	}
};

struct tempral_vector_nonzero {
	static inline void temporal_vector(int mvcol, int scale, int16_t& mv0, int16_t& mv1) {
		int t = (mvcol * scale + 128) >> 8;
		mv0 = t;
		mv1 = t - mvcol;
	}
	void operator()(const h264d_vector_t *mvcol, int scale, h264d_vector_t *tmv) const {
		temporal_vector(mvcol->v[0], scale, tmv[0].v[0], tmv[1].v[0]);
		temporal_vector(mvcol->v[1], scale, tmv[0].v[1], tmv[1].v[1]);
	}
};

template <int N, int BLOCK, int X, int Y, typename F>
static inline void temporal_direct_block_base(h264d_mb_current *mb, const h264d_vector_t *mvcol, int8_t *ref_idx, h264d_vector_t *mv, int blk_idx, int scale,
						 F TemporalVector)
{
	int xoffset = (blk_idx & 1) * 8;
	int yoffset = (blk_idx & 2) * 4;
	for (int i = 0; i < 64 / (BLOCK * BLOCK); ++i) {
		h264d_vector_t *tmv = &mv[N == 8 ? i * 2 : ((i & 2) * 4 + (i & 1) * 2)];
		TemporalVector(mvcol, scale, tmv);
		direct_mv_pred(mb, ref_idx, tmv, X, Y, xoffset + (i & 1) * 4, yoffset + (i & 2) * 2);
		mvcol += (i & 1) ? 3 : 1;
	}
}

template <int N, int BLOCK, int X, int Y>
static inline void temporal_direct_block(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, int8_t *ref_idx, h264d_vector_t *mv, int blk_idx)
{
	int map_idx = col_mb->ref[blk_idx];
	int ref = (0 <= map_idx) ? mb->bdirect->map_col_to_list0[map_idx] : 0;
	ref_idx[0] = ref;
	ref_idx[1] = 0;
	if ((0 <= map_idx) && (mb->frame->refs[0][ref].in_use != LONG_TERM)) {
		const h264d_vector_t *mvcol = &col_mb->mv[(blk_idx & 2) * 4 + (blk_idx & 1) * 2];
		temporal_direct_block_base<N, BLOCK, X, Y>(mb, mvcol, ref_idx, mv, blk_idx, mb->bdirect->scale[ref], tempral_vector_nonzero());
	} else {
		temporal_direct_block_base<N, BLOCK, X, Y>(mb, zero_mov, ref_idx, mv, blk_idx, 0, tempral_vector_zero());
	}
}

static void pred_direct8x8_temporal(h264d_mb_current *mb, int blk_idx, prev8x8_t *pblk)
{
	const h264d_ref_frame_t *colpic = &(mb->frame->refs[1][0]);
	const h264d_col_mb_t *col_mb = &colpic->col->col_mb[mb->y * mb->max_x + mb->x];
	pblk += blk_idx;
	if (col_mb->type == COL_MB8x8) {
		temporal_direct_block<8, 4, 4, 4>(mb, col_mb, pblk->ref, pblk->mv[0], blk_idx);
	} else {
		temporal_direct_block<8, 8, 8, 8>(mb, col_mb, pblk->ref, pblk->mv[0], blk_idx);
		memcpy(pblk->mv[1], pblk->mv[0], sizeof(pblk->mv[0]));
		memcpy(pblk->mv[2], pblk->mv[0], sizeof(pblk->mv[0]));
		memcpy(pblk->mv[3], pblk->mv[0], sizeof(pblk->mv[0]));
	}
}

static void temporal_direct16x16_block8x8_16x16(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, int8_t *ref_idx, h264d_vector_set_t *mv)
{
	temporal_direct_block<16, 8, 16, 16>(mb, col_mb, &ref_idx[0], &mv[0].mv[0], 0);
	memset(&mv[1], 0, sizeof(mv[1]));
}

static void temporal_direct16x16_block8x8_16x8(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, int8_t *ref_idx, h264d_vector_set_t *mv)
{
	for (int y = 0; y < 2; ++y) {
		temporal_direct_block<16, 8, 16, 8>(mb, col_mb, &ref_idx[y * 2], &mv[y].mv[0], y * 2);
	}
	memset(&mv[2], 0, sizeof(mv[2]) * 2);
}

static void temporal_direct16x16_block8x8_8x16(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, int8_t *ref_idx, h264d_vector_set_t *mv)
{
	for (int x = 0; x < 2; ++x) {
		temporal_direct_block<16, 8, 8, 16>(mb, col_mb, &ref_idx[x * 2], &mv[x].mv[0], x);
	}
	memset(&mv[2], 0, sizeof(mv[2]) * 2);
}

template <int N>
void temporal_direct16x16_blockNxN_8x8(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, int8_t *ref_idx, h264d_vector_set_t *mv)
{
	for (int blk_idx = 0; blk_idx < 4; ++blk_idx) {
		temporal_direct_block<16, N, N, N>(mb, col_mb, &ref_idx[blk_idx * 2], &mv[(blk_idx & 2) * 64 / (N * N) + (blk_idx & 1) * 8 / N].mv[0], blk_idx);
	}
//	memset(&mv[16 * 4 / N], 0, sizeof(mv[16 * 4 / N]) * 16 * 4 / N);
}

static void (* const temporal_direct16x16[2][4])(h264d_mb_current *mb, const h264d_col_mb_t *col_mb, int8_t *ref_idx, h264d_vector_set_t *mv) = {
	{
		temporal_direct16x16_block8x8_16x16,
		temporal_direct16x16_block8x8_16x8,
		temporal_direct16x16_block8x8_8x16,
		temporal_direct16x16_blockNxN_8x8<4>
	},
	{
		temporal_direct16x16_block8x8_16x16,
		temporal_direct16x16_block8x8_16x8,
		temporal_direct16x16_block8x8_8x16,
		temporal_direct16x16_blockNxN_8x8<8>
	}
};

static void b_skip_mb_temporal(h264d_mb_current *mb, int8_t *ref_idx, h264d_vector_set_t *mv)
{
	const h264d_ref_frame_t *colpic = &(mb->frame->refs[1][0]);
	const h264d_col_mb_t *col_mb = &colpic->col->col_mb[mb->y * mb->max_x + mb->x];
	temporal_direct16x16[0][col_mb->type](mb, col_mb, ref_idx, mv);
}

static int skip_mbs(h264d_mb_current *mb, uint32_t skip_mb_num, int slice_type)
{
	uint32_t max_mb_run = mb->max_x * mb->max_y - (mb->y * mb->max_x + mb->x);
	uint32_t left4x4, top4x4;
	int8_t *ref_idx;
	void (*skip_mb)(h264d_mb_current *mb, int8_t *ref_idx, h264d_vector_set_t *mv);
	int8_t ref_idx_b[2 * 4];
	static const int8_t ref_idx_p[2] = {
		0, -1
	};

	skip_mb_num = skip_mb_num < max_mb_run ? skip_mb_num : max_mb_run;
	mb->left4x4pred = 0x22222222;
	left4x4 = mb->left4x4coef;
	mb->left4x4coef = 0;
	mb->cbp = 0;
	mb->cbf = 0;
	if (slice_type == P_SLICE) {
		ref_idx = (int8_t *)ref_idx_p;
		skip_mb = p_skip_mb;
	} else {
		ref_idx = ref_idx_b;
		skip_mb = mb->bdirect->func->direct16x16;
	}
	do {
		h264d_vector_set_t mv[16 * 2];
		col_mbtype_t col_mb_type;

		skip_mb(mb, ref_idx, mv);
		*mb->top4x4pred = 0x22222222;
		top4x4 = *mb->top4x4coef;
		*mb->top4x4coef = 0;
		if (slice_type == B_SLICE) {
			const h264d_ref_frame_t *colpic = &(mb->frame->refs[1][0]);
			const h264d_col_mb_t *col_mb = &colpic->col->col_mb[mb->y * mb->max_x + mb->x];
			col_mb_type = col_mb->type;
		} else {
			col_mb_type = COL_MB16x16;
		}
		store_info_inter[0][col_mb_type](mb, mv, ref_idx, 0U, 0U, left4x4, top4x4);
		left4x4 = 0;
		mb->prev_qp_delta = 0;
		mb->type = MB_PSKIP;
		mb->left4x4inter->type = MB_PSKIP;
		mb->left4x4inter->mb_skip = 1;
		mb->left4x4inter->direct8x8 = 3;
		mb->top4x4inter->type = MB_PSKIP;
		mb->top4x4inter->direct8x8 = 3;
		mb->top4x4inter->mb_skip = 1;
		if (increment_mb_pos(mb) < 0) {
			return -1;
		}
	} while (--skip_mb_num);
	return 0;
}

static int more_rbsp_data(dec_bits *st)
{
	int bits;

	bits = not_aligned_bits(st);
	if (bits == 0) {
		bits = 8;
	}
	if (show_bits(st, bits) == (1U << (bits - 1))) {
		return (1 < (show_bits(st, bits + 24) & 0xffffff));
	} else {
		return 1;
	}
}

static int post_process(h264d_context *h2d, h264d_mb_current *mb);
static void init_cabac_context(h264d_cabac_t *cabac, int slice_qp, int idc);
static void init_cabac_engine(h264d_cabac_t *cb, dec_bits *st);
static inline int cabac_decode_terminate(h264d_cabac_t *cb, dec_bits *st);
static int mb_skip_cabac(h264d_mb_current *mb, dec_bits *st, int slice_type);

static int slice_data(h264d_context *h2d, dec_bits *st)
{
	h264d_slice_header *hdr = h2d->slice_header;
	h264d_pps *pps = &h2d->pps_i[hdr->pic_parameter_set_id];
	h264d_mb_current *mb = &h2d->mb_current;
	int is_ae = pps->entropy_coding_mode_flag;

	if (is_ae) {
		init_cabac_context(mb->cabac, mb->qp, (hdr->slice_type == I_SLICE) ? 0 : hdr->cabac_init_idc + 1);
		byte_align(st);
		init_cabac_engine(mb->cabac, st);
	}
	do {
		uint32_t skip_num;
		if ((hdr->slice_type != I_SLICE)
		    && (hdr->slice_type != SI_SLICE)) {
			skip_num = is_ae ? mb_skip_cabac(mb, st, hdr->slice_type) : ue_golomb(st);
			if (skip_num) {
				if (skip_mbs(mb, skip_num, hdr->slice_type) < 0) {
					break;
				}
				if (is_ae) {
					continue;
				}
			}
			if (!is_ae && !more_rbsp_data(st)) {
				break;
			}
		}
		if (is_ae) {
			macroblock_layer_cabac(mb, hdr->slice_type, st);
		} else {
			macroblock_layer(mb, hdr->slice_type, st);
		}
		mb->left4x4inter->mb_skip = 0;
		mb->top4x4inter->mb_skip = 0;
		if (increment_mb_pos(mb) < 0) {
			break;
		}
	} while (is_ae ? !cabac_decode_terminate(mb->cabac, st) : more_rbsp_data(st));
	return post_process(h2d, mb);
}

/**Pair of alpha, square of beta.
 */
static const int16_t deblock_ab2[52 - 16][2] = {
	{4, 4}, {4, 4}, {5, 4}, {6, 9},
	{7, 9}, {8, 9}, {9, 9}, {10, 16},
	{12, 16}, {13, 16}, {15, 36}, {17, 36},
	{20, 49}, {22, 49}, {25, 64}, {28, 64},
	{32, 81}, {36, 81}, {40, 100}, {45, 100},
	{50, 121}, {56, 121}, {63, 144}, {71, 144},
	{80, 169}, {90, 169}, {101, 196}, {113, 196},
	{127, 225}, {144, 225}, {162, 256}, {182, 256},
	{203, 289}, {226, 289}, {255, 324}, {255, 324}
};

static const int8_t str1_3_tc0[52 - 16][3] = {
	{0, 0, 0}, {0, 0, 1}, {0, 0, 1}, {0, 0, 1},
	{0, 0, 1}, {0, 1, 1}, {0, 1, 1}, {1, 1, 1},
	{1, 1, 1}, {1, 1, 1}, {1, 1, 1}, {1, 1, 2},
	{1, 1, 2}, {1, 1, 2}, {1, 1, 2}, {1, 2, 3},
	{1, 2, 3}, {2, 2, 3}, {2, 2, 4}, {2, 3, 4},
	{2, 3, 4}, {3, 3, 5}, {3, 4, 6}, {3, 4, 6},
	{4, 5, 7}, {4, 5, 8}, {4, 6, 9}, {5, 7, 10},
	{6, 8, 11}, {6, 8, 13}, {7, 10, 14}, {8, 11, 16},
	{9, 12, 18}, {10, 13, 20}, {11, 15, 23}, {13, 17, 25}
};


#define CLIP(a, x) ((unsigned)(x) <= (a) ? (x) : (a))

static inline void AlphaBetaTc0(int &alpha, int &beta2, const int8_t **tc0, int q, int alpha_offset, int beta_offset) {
	int a = q + alpha_offset;
	int b = q + beta_offset;
	if (a < 16 || b < 16) {
		alpha = 0;
	} else {
		a = a <= 51 ? a : 51;
		b = b <= 51 ? b : 51;
		alpha = deblock_ab2[a - 16][0];
		beta2 = deblock_ab2[b - 16][1];
		*tc0 = str1_3_tc0[a - 16];
	}
}

static inline int SQUARE(int x) {
	return x * x;
}

template<int N>
struct Strength4h {
	void operator()(uint8_t *dst, int q0, int q1, int p0, int p1, int a, int b2, int tc0) const {
		int q2 = dst[0];
		int p2 = dst[5 * N];
		int aq = SQUARE(q0 - q2);
		int ap = SQUARE(p0 - p2);
		int t;

		t = SQUARE(p0 - q0) < SQUARE((a >> 2) + 2);
		if (N == 1 && ap < b2 && t) {
			dst[3 * N] = ((p0 + p1 + q0) * 2 + q1 + p2 + 4) >> 3;
			dst[4 * N] = (q0 + p0 + p1 + p2 + 2) >> 2;
			dst[5 * N] = (dst[6 * N] * 2 + p2 * 3 + p1 + p0 + q0 + 4) >> 3;
		} else {
			dst[3 * N] = (p1 * 2 + p0 + q1 + 2) >> 2;
		}
		if (N == 1 && aq < b2 && t) {
			dst[2 * N] = ((q0 + q1 + p0) * 2 + p1 + q2 + 4) >> 3;
			dst[N] = (p0 + q0 + q1 + q2 + 2) >> 2;
			dst[0] = (dst[-N] * 2 + q2 * 3 + q1 + q0 + p0 + 4) >> 3;
		} else {
			dst[2 * N] = (q1 * 2 + q0 + p1 + 2) >> 2;
		}
	}
};

static inline int CLIP3(int x, int c) {
	return (x < -c) ? -c : (c < x) ? c : x;
}

template<int N>
struct Strength1_3h {
	void operator()(uint8_t *dst, int q0, int q1, int p0, int p1, int a, int b2, int tc0) const {
		int q2 = dst[0 * N];
		int p2 = dst[5 * N];
		int aq_smaller = SQUARE(q2 - q0) < b2;
		int ap_smaller = SQUARE(p2 - p0) < b2;
		if (N == 1 && tc0) {
			if (aq_smaller) {
				int t = (q2 + ((p0 + q0 + 1) >> 1) - (q1 * 2)) >> 1;
				if (t) {
					dst[1 * N] = CLIP3(t, tc0) + q1;
				}
			}
			if (ap_smaller) {
				int t = (p2 + ((q0 + p0 + 1) >> 1) - (p1 * 2)) >> 1;
				if (t) {
					dst[4 * N] = CLIP3(t, tc0) + p1;
				}
			}
		}
		tc0 = tc0 + ((N == 1) ? aq_smaller + ap_smaller : 1);
		if (tc0) {
			int delta = (((p0 - q0) * 4) + q1 - p1 + 4) >> 3;
			if (delta) {
				delta = CLIP3(delta, tc0);
				q0 = q0 + delta;
				p0 = p0 - delta;
				dst[2 * N] = CLIP255C(q0);
				dst[3 * N] = CLIP255C(p0);
			}
		}
	}
};

template <typename T, int N, int LEN>
static inline void deblock_luma_horiz(int a, int b2, uint8_t *luma, int tc0, int stride, T strength)
{
	uint8_t *dst = luma - 3 * N;
	assert(a);
	int a2 = SQUARE(a);
	for (int y = 0; y < LEN; ++y) {
		int q0, q1, p0, p1;
		int t;
		q0 = dst[2 * N];
		p0 = dst[3 * N];
		t = q0 - p0;
		if (SQUARE(t) < a2) {
			q1 = dst[1 * N];
			t = q1 - q0;
			if (SQUARE(t) < b2) {
				p1 = dst[4 * N];
				t = p0 - p1;
				if (SQUARE(t) < b2) {
					strength(dst, q0, q1, p0, p1, a, b2, tc0);
				}
			}
		}
		dst += stride;
	}
}

static inline void deblock_luma_horiz_str4(int a, int b2, uint8_t *luma, int stride) {
	deblock_luma_horiz<Strength4h<1>, 1, 16>(a, b2, luma, 0, stride, Strength4h<1>());
}

static inline void deblock_luma_horiz_str1_3(int a, int b2, const int8_t *tc0, uint8_t *luma, int str, int stride)
{
	for (int i = 0; i < 4; ++i) {
		int strength_1 = str & 3;
		if (strength_1) {
			deblock_luma_horiz<Strength1_3h<1>, 1, 4>(a, b2, luma, tc0[strength_1 - 1], stride, Strength1_3h<1>());
		}
		luma += stride * 4;
		str = (unsigned)str >> 2;
	}
}

static inline void deblock_chroma_horiz_str4(int a, int b2, uint8_t *luma, int stride) {
	deblock_luma_horiz<Strength4h<2>, 2, 8>(a, b2, luma, 0, stride, Strength4h<2>());
}

static inline void deblock_chroma_horiz_str1_3(int a, int b2, const int8_t *tc0, uint8_t *chroma, int str, int stride)
{
	for (int i = 0; i < 4; ++i) {
		int strength_1 = str & 3;
		if (strength_1) {
			deblock_luma_horiz<Strength1_3h<2>, 2, 2>(a, b2, chroma, tc0[strength_1 - 1], stride, Strength1_3h<2>());
		}
		chroma += stride * 2;
		str = (unsigned)str >> 2;
	}
}


template <int N>
struct Strength4v {
	void operator()(uint8_t *dst, int q0, int q1, int p0, int p1, int a, int b2, int tc0, int stride) const {
		int q2 = dst[0];
		int p2 = dst[stride * 5];
		int aq = SQUARE(q0 - q2);
		int ap = SQUARE(p0 - p2);
		int t;

		t = SQUARE(p0 - q0) < SQUARE((a >> 2) + 2);
		if (N == 1 && ap < b2 && t) {
			dst[stride * 3] = ((p0 + p1 + q0) * 2 + q1 + p2 + 4) >> 3;
			dst[stride * 4] = (q0 + p0 + p1 + p2 + 2) >> 2;
			dst[stride * 5] = (dst[stride * 6] * 2 + p2 * 3 + p1 + p0 + q0 + 4) >> 3;
		} else {
			dst[stride * 3] = (p1 * 2 + p0 + q1 + 2) >> 2;
		}
		if (N == 1 && aq < b2 && t) {
			dst[stride * 2] = ((q0 + q1 + p0) * 2 + p1 + q2 + 4) >> 3;
			dst[stride * 1] = (p0 + q0 + q1 + q2 + 2) >> 2;
			dst[0] = (dst[-stride] * 2 + q2 * 3 + q1 + q0 + p0 + 4) >> 3;
		} else {
			dst[stride * 2] = (q1 * 2 + q0 + p1 + 2) >> 2;
		}
	}
};

template <int N>
struct Strength1_3v {
	void operator()(uint8_t *dst, int q0, int q1, int p0, int p1, int a, int b2, int tc0, int stride) const {
		int q2 = dst[0];
		int p2 = dst[stride * 5];
		int aq_smaller = SQUARE(q2 - q0) < b2;
		int ap_smaller = SQUARE(p2 - p0) < b2;
		if (N == 1 && tc0) {
			if (aq_smaller) {
				int t = (q2 + ((p0 + q0 + 1) >> 1) - (q1 * 2)) >> 1;
				if (t) {
					dst[stride * 1] = CLIP3(t, tc0) + q1;
				}
			}
			if (ap_smaller) {
				int t = (p2 + ((q0 + p0 + 1) >> 1) - (p1 * 2)) >> 1;
				if (t) {
					dst[stride * 4] = CLIP3(t, tc0) + p1;
				}
			}
		}
		tc0 = tc0 + ((N == 1) ? aq_smaller + ap_smaller : 1);
		if (tc0) {
			int delta = (((p0 - q0) * 4) + q1 - p1 + 4) >> 3;
			if (delta) {
				delta = CLIP3(delta, tc0);
				q0 = q0 + delta;
				p0 = p0 - delta;
				dst[stride * 2] = CLIP255C(q0);
				dst[stride * 3] = CLIP255C(p0);
			}
		}
	}
};

template <typename T, int LEN>
static inline void deblock_luma_vert(int a, int b2, uint8_t *luma, int tc0, int stride, T strength)
{
	uint8_t *dst = luma - stride * 3;
	assert(a);
	for (int x = 0; x < LEN; ++x) {
		int q0, q1, p0, p1;
		int t;
		q0 = dst[stride * 2];
		p0 = dst[stride * 3];
		t = q0 - p0;
		if (SQUARE(t) < SQUARE(a)) {
			q1 = dst[stride];
			t = q1 - q0;
			if (SQUARE(t) < b2) {
				p1 = dst[stride * 4];
				t = p0 - p1;
				if (SQUARE(t) < b2) {
					strength(dst, q0, q1, p0, p1, a, b2, tc0, stride);
				}
			}
		}
		dst += 1;
	}
}

static inline void deblock_luma_vert_str4(int a, int b2, uint8_t *luma, int stride) {
	deblock_luma_vert<Strength4v<1>, 16>(a, b2, luma, 0, stride, Strength4v<1>());
}

static inline void deblock_luma_vert_str1_3(int a, int b2, const int8_t *tc0, uint8_t *luma, int str, int stride)
{
	for (int i = 0; i < 4; ++i) {
		int strength_1 = str & 3;
		if (strength_1) {
			deblock_luma_vert<Strength1_3v<1>, 4>(a, b2, luma, tc0[strength_1 - 1], stride, Strength1_3v<1>());
		}
		luma += 4;
		str = (unsigned)str >> 2;
	}
}

static inline void deblock_chroma_vert_str4(int a, int b2, uint8_t *chroma, int stride) {
	deblock_luma_vert<Strength4v<2>, 16>(a, b2, chroma, 0, stride, Strength4v<2>());
}

static inline void deblock_chroma_vert_str1_3(int a, int b2, const int8_t *tc0, uint8_t *chroma, int str, int stride)
{
	for (int i = 0; i < 4; ++i) {
		int strength_1 = str & 3;
		if (strength_1) {
			deblock_luma_vert<Strength1_3v<2>, 4>(a, b2, chroma, tc0[strength_1 - 1], stride, Strength1_3v<2>());
		}
		chroma += 4;
		str >>= 2;
	}
}

static inline void deblock_pb(h264d_mb_current *mb)
{
	int qp;
	int a, b2;
	const int8_t *tc0;
	int alpha_offset, beta_offset;
	int max_x = mb->max_x;
	int max_y = mb->max_y;
	int stride = max_x * 16;
	uint8_t *luma = mb->frame->curr_luma;
	uint8_t *chroma = mb->frame->curr_chroma;
	deblock_info_t *curr = mb->deblock_base;
	int idc;

	for (int y = 0; y < max_y; ++y) {
		for (int x = 0; x < max_x; ++x) {
			uint32_t str;
			if (curr->idc) {
				idc = curr->idc - 1;
				DEC_SLICEHDR(curr->slicehdr, alpha_offset, beta_offset);
			}
			if (idc == 1) {
				curr++;
				luma += 16;
				chroma += 16;
				continue;
			}
			str = curr->str_horiz;
			if ((x != 0) && (!idc || mb->firstline != max_x) && (str & 255)) {
				/* alpha, beta of MB left edge */
				qp = (curr->qpy + (curr - 1)->qpy + 1) >> 1;
				AlphaBetaTc0(a, b2, &tc0, qp, alpha_offset, beta_offset);
				if (a) {
					if (curr->str4_horiz) {
						deblock_luma_horiz_str4(a, b2, luma, stride);
					} else {
						if (str & 255) {
							deblock_luma_horiz_str1_3(a, b2, tc0, luma, str, stride);
						}
					}
				}
				qp = (curr->qpc + (curr - 1)->qpc + 1) >> 1;
				AlphaBetaTc0(a, b2, &tc0, qp, alpha_offset, beta_offset);
				if (a) {
					if (curr->str4_horiz) {
						deblock_chroma_horiz_str4(a, b2, chroma, stride);
						deblock_chroma_horiz_str4(a, b2, chroma + 1, stride);
					} else {
						if (str & 255) {
							deblock_chroma_horiz_str1_3(a, b2, tc0, chroma, str, stride);
							deblock_chroma_horiz_str1_3(a, b2, tc0, chroma + 1, str, stride);
						}
					}
				}
			}
			if (str & ~255) {
				AlphaBetaTc0(a, b2, &tc0, curr->qpy, alpha_offset, beta_offset);
				if (a) {
					uint32_t str_tmp = str >> 8;
					if (str_tmp & 255) {
						deblock_luma_horiz_str1_3(a, b2, tc0, luma + 4, str_tmp, stride);
					}
					str_tmp >>= 8;
					if (str_tmp & 255) {
						deblock_luma_horiz_str1_3(a, b2, tc0, luma + 8, str_tmp, stride);
					}
					str_tmp >>= 8;
					if (str_tmp & 255) {
						deblock_luma_horiz_str1_3(a, b2, tc0, luma + 12, str_tmp, stride);
					}
				}
				if (curr->qpy != curr->qpc) {
					AlphaBetaTc0(a, b2, &tc0, curr->qpc, alpha_offset, beta_offset);
				}
				if (a) {
					str >>= 16;
					if (str & 0xff) {
						deblock_chroma_horiz_str1_3(a, b2, tc0, chroma + 8, str, stride);
						deblock_chroma_horiz_str1_3(a, b2, tc0, chroma + 8 + 1, str, stride);
					}
				}
			}
			str = curr->str_vert;
			if ((y != 0) && (!idc || mb->firstline < 0) && (str & 255)) {
				/* top edge of MB */
				qp = (curr->qpy + (curr - max_x)->qpy + 1) >> 1;
				AlphaBetaTc0(a, b2, &tc0, qp, alpha_offset, beta_offset);
				if (a) {
					if (curr->str4_vert) {
						deblock_luma_vert_str4(a, b2, luma, stride);
					} else {
						if (str & 255) {
							deblock_luma_vert_str1_3(a, b2, tc0, luma, str, stride);
						}
					}
				}
				qp = (curr->qpc + (curr - max_x)->qpc + 1) >> 1;
				AlphaBetaTc0(a, b2, &tc0, qp, alpha_offset, beta_offset);
				if (a) {
					if (curr->str4_vert) {
						deblock_chroma_vert_str4(a, b2, chroma, stride);
					} else {
						if (str & 255) {
							deblock_chroma_vert_str1_3(a, b2, tc0, chroma, str, stride);
						}
					}
				}
			}
			if (str & ~255) {
				AlphaBetaTc0(a, b2, &tc0, curr->qpy, alpha_offset, beta_offset);
				if (a) {
					uint32_t str_tmp = str >> 8;
					if (str_tmp & 255) {
						deblock_luma_vert_str1_3(a, b2, tc0, luma + stride * 4, str_tmp, stride);
					}
					str_tmp >>= 8;
					if (str_tmp & 255) {
						deblock_luma_vert_str1_3(a, b2, tc0, luma + stride * 8, str_tmp, stride);
					}
					str_tmp >>= 8;
					if (str_tmp & 255) {
						deblock_luma_vert_str1_3(a, b2, tc0, luma + stride * 12, str_tmp, stride);
					}
				}
				if (curr->qpy != curr->qpc) {
					AlphaBetaTc0(a, b2, &tc0, curr->qpc, alpha_offset, beta_offset);
				}
				if (a) {
					str >>= 16;
					if (str & 255) {
						deblock_chroma_vert_str1_3(a, b2, tc0, chroma + stride * 4, str, stride);
					}
				}
			}
			curr++;
			luma += 16;
			chroma += 16;
		}
		luma += stride * 15;
		chroma += stride * 7;
	}
}

static inline h264d_ref_frame_t *marking_sliding_window(h264d_ref_frame_t *refs, int frame_ptr, int frame_num, int max_frame_num, int num_ref_frames, int poc)
{
	int min_frm_num = INT_MAX;
	int min_idx = 0;
	int empty_idx = -1;
	int num_long = 0;
	int num_short = 0;

	for (int i = 0; i < 16; ++i) {
		int in_use = refs[i].in_use;
		if (in_use == NOT_IN_USE) {
			if (empty_idx < 0) {
				empty_idx = i;
			}
		} else if (in_use == SHORT_TERM) {
			int num = refs[i].num;
			if (frame_num < num) {
				num -= max_frame_num;
			}
			if (num < min_frm_num) {
				min_frm_num = num;
				min_idx = i;
			}
			num_short++;
		} else {
			num_long++;
		}
	}
	if (num_short + num_long < num_ref_frames) {
		refs +=	(0 <= empty_idx) ? empty_idx : num_ref_frames - 1;
	} else {
		refs += min_idx;
	}
	refs->in_use = SHORT_TERM;
	refs->frame_idx = frame_ptr;
	refs->num = frame_num;
	refs->poc = poc;
	return refs;
}

static void mmco_discard(h264d_ref_frame_t *refs, int in_use, uint32_t target_num)
{
	int i = 16;
	do {
		if (refs->num == target_num && refs->in_use == in_use) {
			refs->in_use = NOT_IN_USE;
			break;
		}
		refs++;
	} while (--i);
}

static void mmco_op1(const h264d_mmco *mmco, h264d_ref_frame_t *refs, int frame_ptr, int frame_num, int max_frame_num, int num_ref_frames, int poc)
{
	int num = frame_num - mmco->arg1 - 1;
	while (num < 0) {
		num += max_frame_num;
	}
	mmco_discard(refs, SHORT_TERM, num);
}

static void mmco_op2(const h264d_mmco *mmco, h264d_ref_frame_t *refs, int frame_ptr, int frame_num, int max_frame_num, int num_ref_frames, int poc)
{
	mmco_discard(refs, LONG_TERM, mmco->arg1);
}

static void mmco_op3(const h264d_mmco *mmco, h264d_ref_frame_t *refs, int frame_ptr, int frame_num, int max_frame_num, int num_ref_frames, int poc)
{
	uint32_t long_num = mmco->arg2;
	uint32_t target_num = frame_num - mmco->arg1 - 1;
	int i = 16;

	while ((int)target_num < 0) {
		target_num += max_frame_num;
	}
	do {
		int in_use = refs->in_use;
		if ((in_use == LONG_TERM) && (refs->num == long_num)) {
			refs->in_use = NOT_IN_USE;
		} else if ((in_use == SHORT_TERM) && (refs->num == target_num)) {
			refs->in_use = LONG_TERM;
			refs->num = long_num;
		}
		refs++;
	} while (--i);
}

static void mmco_op4(const h264d_mmco *mmco, h264d_ref_frame_t *refs, int frame_ptr, int frame_num, int max_frame_num, int num_ref_frames, int poc)
{
	int i = 16;
	uint32_t max_long_term_idx_plus1 = mmco->arg1;
	do {
		if (refs->in_use == LONG_TERM && max_long_term_idx_plus1 <= refs->num) {
			refs->in_use = NOT_IN_USE;
		}
		refs++;
	} while (--i);
}

static void mmco_op5(const h264d_mmco *mmco, h264d_ref_frame_t *refs, int frame_ptr, int frame_num, int max_frame_num, int num_ref_frames, int poc)
{
	int i = 16;
	do {
		refs->in_use = NOT_IN_USE;
		refs++;
	} while (--i);
}

static void mmco_op6(const h264d_mmco *mmco, h264d_ref_frame_t *refs, int frame_ptr, int frame_num, int max_frame_num, int num_ref_frames, int poc)
{
	h264d_ref_frame_t *ref = marking_sliding_window(refs, frame_ptr, frame_num, max_frame_num, num_ref_frames, poc);
	ref->in_use = LONG_TERM;
	ref->num = mmco->arg1;
}

static void (* const mmco_ops[6])(const h264d_mmco *mmco, h264d_ref_frame_t *refs, int frame_ptr, int frame_num, int max_frame_num, int num_ref_frames, int poc) = {
	mmco_op1, mmco_op2, mmco_op3,
	mmco_op4, mmco_op5, mmco_op6,
};

static inline int marking_mmco(h264d_marking_t *mrk, h264d_ref_frame_t *refs, int frame_ptr, int frame_num, int max_frame_num, int num_ref_frames, int poc)
{
	const h264d_mmco *mmco = mrk->mmco;
	int op5_detect = 0;
	int op6_detect = 0;
	int i = 16;
	do {
		int op = mmco->op;
		if (op == 0) {
			break;
		} else if (5 <= op) {
			if (op == 5) {
				op5_detect = 1;
			} else {
				op6_detect = 1;
			}
		}
		mmco_ops[op - 1](mmco, refs, frame_ptr, frame_num, max_frame_num, num_ref_frames, poc);
		mmco++;
	} while (--i);
	if (!op6_detect) {
		if (op5_detect) {
			frame_num = poc = 0;
		}
		marking_sliding_window(refs, frame_ptr, frame_num, max_frame_num, num_ref_frames, poc);
	}
	return op5_detect;
}

static inline void gap_mbs(h264d_slice_header *hdr, h264d_mb_current *mb, h264d_ref_frame_t *refs, int max_frame_num, int num_ref_frames)
{
	int frame_num = hdr->frame_num;
	int prev_frame_num = hdr->prev_frame_num;
	int gap = frame_num - prev_frame_num;
	while (gap < 0) {
		gap += max_frame_num;
	}
	if (0 < --gap) {
		int poc = hdr->poc;
		if (16 < gap) {
			gap = 16;
			prev_frame_num = frame_num - 17;
		}
		do {
			if (max_frame_num <= ++prev_frame_num) {
				prev_frame_num -= max_frame_num;
			}
			marking_sliding_window(refs, mb->frame->index, prev_frame_num, max_frame_num, num_ref_frames, poc);
		} while (--gap);
	}
}

static inline void post_ref_pic_marking(h264d_slice_header *hdr, int nal_unit_type, int max_frame_num, int num_ref_frames, h264d_mb_current *mb, int lx)
{
	h264d_ref_frame_t *refs = hdr->reorder[lx].ref_frames;
	h264d_marking_t *mrk = &hdr->marking;
	int frame_num = hdr->frame_num;
	int poc = hdr->poc;

	if (nal_unit_type == SLICE_IDR_NAL) {
		refs[0].in_use = mrk->long_term_reference_flag ? LONG_TERM : SHORT_TERM;
		refs[0].frame_idx = mb->frame->index;
		refs[0].num = frame_num;
		refs[0].poc = poc;
		for (int i = 1; i < 16; ++i) {
			refs[i].in_use = NOT_IN_USE;
		}
	} else {
		if (!hdr->marking.idr && !hdr->marking.mmco5) {
			gap_mbs(hdr, mb, refs, max_frame_num, num_ref_frames);
		}
		if (mrk->adaptive_ref_pic_marking_mode_flag) {
			if (marking_mmco(mrk, refs, mb->frame->index, frame_num, max_frame_num, num_ref_frames, poc)) {
				hdr->frame_num = 0;
			}
		} else {
			marking_sliding_window(refs, mb->frame->index, frame_num, max_frame_num, num_ref_frames, poc);
		}
	}
}

static inline void insert_dpb(h264d_dpb_t *dpb, int poc, int frame_idx, int is_idr)
{
	if (is_idr) {
		dpb_insert_idr(dpb, poc, frame_idx);
	} else {
		dpb_insert_non_idr(dpb, poc, frame_idx);
	}
}

struct frame_num_descent_p {
	static int unwrap(int s, int frame_num, int max_frame_num) {
		return (frame_num < s) ? s - max_frame_num : s;
	}
	bool operator()(const int& l, const int& r, int frame_num, int max_frame_num) const {
		return unwrap(l, frame_num, max_frame_num) > unwrap(r, frame_num, max_frame_num);
	}
};

struct poc_order_b_l0 {
	bool operator()(const int& l, const int& r, int curr_poc, int na) const {
		if (l < curr_poc) {
			return (curr_poc < r) || (l > r);
		} else {
			return (curr_poc < r) && (l < r);
		}
	}
};

struct poc_order_b_l1 {
	bool operator()(const int& l, const int& r, int curr_poc, int na) const {
		if (l > curr_poc) {
			return (curr_poc > r) || (l < r);
		} else {
			return (curr_poc > r) && (l > r);
		}
	}
};

struct get_frame_num {
	int operator()(const h264d_ref_frame_t&ref) const {
		return ref.num;
	}
};

struct get_poc {
	int operator()(const h264d_ref_frame_t&ref) const {
		return ref.poc;
	}
};

template <typename F0, typename F1>
static inline bool ref_list_order(const h264d_ref_frame_t& lhs, const h264d_ref_frame_t& rhs, int curr_num, int max_num,
				 F0 GetNum,
				 F1 LessShortTerm)
{
	int l_use = lhs.in_use;
	int r_use = rhs.in_use;
	if (l_use == SHORT_TERM) {
		if (r_use == SHORT_TERM) {
			return LessShortTerm(GetNum(lhs), GetNum(rhs), curr_num, max_num);
		} else {
			return true;
		}
	} else if (l_use == LONG_TERM) {
		if (r_use == SHORT_TERM) {
			return false;
		} else if (r_use == LONG_TERM) {
			return GetNum(lhs) < GetNum(rhs);
		} else {
			return true;
		}
	} else {
		return false;
	}
}

struct ref_list_less_p {
	ref_list_less_p(int curr_num, int max_num) : curr_num_(curr_num), max_num_(max_num) {}
	bool operator()(const h264d_ref_frame_t& lhs, const h264d_ref_frame_t& rhs) const {
		return ref_list_order(lhs, rhs, curr_num_, max_num_, get_frame_num(), frame_num_descent_p());
	}
private:
	int curr_num_;
	int max_num_;
};

struct ref_list_order_b_ref0 {
	ref_list_order_b_ref0(int curr_poc) : curr_poc_(curr_poc) {}
	bool operator()(const h264d_ref_frame_t& lhs, const h264d_ref_frame_t& rhs) const {
		return ref_list_order(lhs, rhs, curr_poc_, 0, get_poc(), poc_order_b_l0());
	}
private:
	int curr_poc_;
};

struct ref_list_order_b_ref1 {
	ref_list_order_b_ref1(int curr_poc) : curr_poc_(curr_poc) {}
	bool operator()(const h264d_ref_frame_t& lhs, const h264d_ref_frame_t& rhs) const {
		return ref_list_order(lhs, rhs, curr_poc_, 0, get_poc(), poc_order_b_l1());
	}
private:
	int curr_poc_;
};

static inline void ref_pic_init_p(h264d_slice_header *hdr, int max_frame_num, int num_ref_frames)
{
	h264d_ref_frame_t *ref = hdr->reorder[0].ref_frames;

	std::sort(ref, ref + num_ref_frames, ref_list_less_p(hdr->frame_num, max_frame_num));
	for (int i = num_ref_frames; i < 16; ++i) {
		ref[i].in_use = NOT_IN_USE;
	}
}

static inline bool is_same_list(const h264d_ref_frame_t *a, const h264d_ref_frame_t *b, int num_elem)
{
	return !memcmp(a, b, sizeof(*a) * num_elem); /* FIXME */
}

static inline void ref_pic_init_b(h264d_slice_header *hdr, int num_ref_frames)
{
	h264d_ref_frame_t *ref0 = hdr->reorder[0].ref_frames;
	h264d_ref_frame_t *ref1 = hdr->reorder[1].ref_frames;

	std::sort(ref0, ref0 + num_ref_frames, ref_list_order_b_ref0(hdr->poc));
	std::sort(ref1, ref1 + num_ref_frames, ref_list_order_b_ref1(hdr->poc));
	if ((1 < num_ref_frames) && is_same_list(ref0, ref1, num_ref_frames)) {
		std::swap(ref1[0], ref1[1]);
	}
	for (int i = num_ref_frames; i < 16; ++i) {
		ref0[i].in_use = NOT_IN_USE;
		ref1[i].in_use = NOT_IN_USE;
	}
}

static inline void record_map_col_ref_frameidx(int8_t *map, const h264d_ref_frame_t *refs1, int num_ref_frames)
{
	for (int i = 0; i < num_ref_frames; ++i) {
		map[i] = (int8_t)refs1[i].frame_idx;
	}
	memset(&map[num_ref_frames], refs1[0].frame_idx, 16 - num_ref_frames);
}

static inline h264d_ref_frame_t *find_l1_curr_pic(h264d_ref_frame_t *refs, int poc)
{
	h264d_ref_frame_t *refs_found = 0;
	for (int i = 0; i < 16; ++i) {
		if (refs->in_use) {
			if (refs->poc == poc) {
				return refs;
			}
			if (!refs_found) {
				refs_found = refs;
			}
		}
		++refs;
	}
	return refs_found ? refs_found : refs - 16;
}

static int post_process(h264d_context *h2d, h264d_mb_current *mb)
{
	h264d_slice_header *hdr;
	int nal_id;
	int is_filled = (mb->y >= mb->max_y);

	hdr = h2d->slice_header;
	if (is_filled) {
		h264d_frame_info_t *frame;
		deblock_pb(mb);
		h264d_sps *sps = &h2d->sps_i[h2d->pps_i[hdr->pic_parameter_set_id].seq_parameter_set_id];
		int max_frame_num = 1 << sps->log2_max_frame_num;
		int num_ref_frames = sps->num_ref_frames;
		nal_id = h2d->id;
		frame = h2d->mb_current.frame;
		if (nal_id & 0x60) {
			post_ref_pic_marking(hdr, nal_id & 31, max_frame_num, num_ref_frames, mb, 0);
			post_ref_pic_marking(hdr, nal_id & 31, max_frame_num, num_ref_frames, mb, 1);
			record_map_col_ref_frameidx(mb->frame->curr_col->map_col_frameidx, mb->frame->refs[1], num_ref_frames);
			std::swap(mb->frame->curr_col, find_l1_curr_pic(mb->frame->refs[1], hdr->marking.mmco5 ? 0 : hdr->poc)->col);
			insert_dpb(&frame->dpb, hdr->poc, frame->index, hdr->marking.idr | hdr->marking.mmco5);
		} else {
			dpb_insert_non_idr(&frame->dpb, hdr->poc, frame->index);
		}
		hdr->prev_frame_num = hdr->frame_num;
		hdr->first_mb_in_slice = mb->max_x * mb->max_x;
	}
	return is_filled;
}

static void init_cabac_engine(h264d_cabac_t *cb, dec_bits *st)
{
	cb->range = 0x1fe;
	cb->offset = get_bits(st, 9);
}

typedef struct {
	int8_t m, n;
} ctx_idx_mn_t;

static const ctx_idx_mn_t ctx_idx_mn_IPB[4][460] = {
	{
		/* 0 to 10 */
		{20, -15}, {2, 54}, {3, 74}, {20, -15},	{2, 54},
		{3, 74}, {-28, 127}, {-23, 104}, {-6, 53}, {-1, 54},
		{7, 51},
		/* 11 to 59, void */
		{0, 0}, {0, 0}, {0, 0}, {0, 0},
		{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
		{0, 0}, {0, 0}, {0, 0}, {0, 0},	{0, 0},
		{0, 0}, {0, 0},	{0, 0}, {0, 0}, {0, 0},
		{0, 0},	{0, 0},	{0, 0}, {0, 0}, {0, 0},
		{0, 0},	{0, 0}, {0, 0}, {0, 0},	{0, 0},
		{0, 0}, {0, 0}, {0, 0},	{0, 0}, {0, 0},
		{0, 0}, {0, 0},	{0, 0},	{0, 0}, {0, 0},
		{0, 0},	{0, 0},	{0, 0}, {0, 0},	{0, 0},
		{0, 0}, {0, 0},	{0, 0}, {0, 0}, {0, 0},
		/* 60 to 69 */
		{0, 41}, {0, 63}, {0, 63}, {0, 63}, {-9, 83},
		{4, 86}, {0, 97}, {-7, 72}, {13, 41}, {3, 62},
		/* 70 to 104 */
		{0, 11}, {1, 55}, {0, 69}, {-17, 127}, {-13, 102},
		{0, 82}, {-7, 74}, {-21, 107}, {-27, 127}, {-31, 127},
		{-24, 127}, {-18, 95}, {-27, 127}, {-21, 114}, {-30, 127},
		{-17, 123}, {-12, 115}, {-16, 122}, {-11, 115}, {-12, 63},
		{-2, 68}, {-15, 84}, {-13, 104}, {-3, 70}, {-8, 93},
		{-10, 90}, {-30, 127}, {-1, 74}, {-6, 97}, {-7, 91},
		{-20, 127}, {-4, 56}, {-5, 82}, {-7, 76}, {-22, 125},
		/* 105 to 135 */
		{-7, 93}, {-11, 87}, {-3, 77}, {-5, 71}, {-4, 63},
		{-4, 68}, {-12, 84}, {-7, 62}, {-7, 65}, {8, 61},
		{5, 56}, {-2, 66}, {1, 64}, {0, 61}, {-2, 78},
		{1, 50}, {7, 52}, {10, 35}, {0, 44}, {11, 38},
		{1, 45}, {0, 46}, {5, 44}, {31, 17}, {1, 51},
		{7, 50}, {28, 19}, {16, 33}, {14, 62}, {-13, 108},
		{-15, 100},
		/* 136 to 165 */
		{-13, 101}, {-13, 91}, {-12, 94}, {-10, 88}, {-16, 84},
		{-10, 86}, {-7, 83}, {-13, 87}, {-19, 94}, {1, 70},
		{0, 72}, {-5, 74}, {18, 59}, {-8, 102}, {-15, 100},
		{0, 95}, {-4, 75}, {2, 72}, {-11, 75}, {-3, 71},
		{15, 46}, {-13, 69}, {0, 62}, {0, 65}, {21, 37},
		{-15, 72}, {9, 57}, {16, 54}, {0, 62}, {12, 72},
		/* 166 to 196 */
		{24, 0}, {15, 9}, {8, 25}, {13, 18}, {15, 9},
		{13, 19}, {10, 37}, {12, 18}, {6, 29}, {20, 33},
		{15, 30}, {4, 45}, {1, 58}, {0, 62}, {7, 61},
		{12, 38}, {11, 45}, {15, 39}, {11, 42}, {13, 44},
		{16, 45}, {12, 41}, {10, 49}, {30, 34}, {18, 42},
		{10, 55}, {17, 51}, {17, 46}, {0, 89}, {26, -19},
		{22, -17},
		/* 197 to 226 */
		{26, -17}, {30, -25}, {28, -20}, {33, -23}, {37, -27},
		{33, -23}, {40, -28}, {38, -17}, {33, -11}, {40, -15},
		{41, -6}, {38, 1}, {41, 17}, {30, -6}, {27, 3},
		{26, 22}, {37, -16}, {35, -4}, {38, -8}, {38, -3},
		{37, 3}, {38, 5}, {42, 0}, {35, 16}, {39, 22},
		{14, 48}, {27, 37}, {21, 60}, {12, 68}, {2, 97},
		/* 227 to 251 */
		{-3, 71}, {-6, 42}, {-5, 50}, {-3, 54}, {-2, 62},
		{0, 58}, {1, 63}, {-2, 72}, {-1, 74}, {-9, 91},
		{-5, 67}, {-5, 27}, {-3, 39}, {-2, 44}, {0, 46},
		{-16, 64}, {-8, 68}, {-10, 78}, {-6, 77}, {-10, 86},
		{-12, 92}, {-15, 55}, {-10, 60}, {-6, 62}, {-4, 65},
		/* 252 to 275 */
		{-12, 73}, {-8, 76}, {-7, 80}, {-9, 88}, {-17, 110},
		{-11, 97}, {-20, 84}, {-11, 79}, {-6, 73}, {-4, 74},
		{-13, 86}, {-13, 96}, {-11, 97}, {-19, 117}, {-8, 78},
		{-5, 33}, {-4, 48}, {-2, 53}, {-3, 62}, {-13, 71},
		{-10, 79}, {-12, 86}, {-13, 90}, {-14, 97},
		/* 276, void */
		{0, 0},
		/* 277 to 307 */
		{-6, 93}, {-6, 84}, {-8, 79}, {0, 66}, {-1, 71},
		{0, 62}, {-2, 60}, {-2, 59}, {-5, 75}, {-3, 62},
		{-4, 58}, {-9, 66}, {-1, 79}, {0, 71}, {3, 68},
		{10, 44}, {-7, 62}, {15, 36}, {14, 40}, {16, 27},
		{12, 29}, {1, 44}, {20, 36}, {18, 32}, {5, 42},
		{1, 48}, {10, 62}, {17, 46}, {9, 64}, {-12, 104},
		{-11, 97},
		/* 308 to 337 */
		{-16, 96}, {-7, 88}, {-8, 85}, {-7, 85}, {-9, 85},
		{-13, 88}, {4, 66}, {-3, 77}, {-3, 76}, {-6, 76},
		{10, 58}, {-1, 76}, {-1, 83}, {-7, 99}, {-14, 95},
		{2, 95}, {0, 76}, {-5, 74}, {0, 70}, {-11, 75},
		{1, 68}, {0, 65}, {-14, 73}, {3, 62}, {4, 62},
		{-1, 68}, {-13, 75}, {11, 55}, {5, 64}, {12, 70},
		/* 338 to 368 */
		{15, 6}, {6, 19}, {7, 16}, {12, 14}, {18, 13},
		{13, 11}, {13, 15}, {15, 16}, {12, 23}, {13, 23},
		{15, 20}, {14, 26}, {14, 44}, {17, 40}, {17, 47},
		{24, 17}, {21, 21}, {25, 22}, {31, 27}, {22, 29},
		{19, 35}, {14, 50}, {10, 57}, {7, 63}, {-2, 77},
		{-4, 82}, {-3, 94}, {9, 69}, {-12, 109}, {36, -35},
		{36, -34},
		/* 369 to 398 */
		{32, -26}, {37, -30}, {44, -32}, {34, -18}, {34, -15},
		{40, -15}, {33, -7}, {35, -5}, {33, 0}, {38, 2},
		{33, 13}, {23, 35}, {13, 58}, {29, -3}, {26, 0},
		{22, 30}, {31, -7}, {35, -15}, {34, -3}, {34, 3},
		{36, -1}, {34, 5}, {32, 11}, {35, 5}, {34, 12},
		{39, 11}, {30, 29}, {34, 26}, {29, 39}, {19, 66},
		/* 399 to 401*/
		{31, 21}, {31, 31}, {25, 50},
		/* 402 to 435 */
		{-17, 120}, {-20, 112}, {-18, 114}, {-11, 85}, {-15, 92},
		{-14, 89}, {-26, 71}, {-15, 81}, {-14, 80}, {0, 68},
		{-14, 70}, {-24, 56}, {-23, 68}, {-24, 50}, {-11, 74},
		{23, -13}, {26, -13}, {40, -15}, {49, -14}, {44, 3},
		{45, 6}, {44, 34}, {33, 54}, {19, 82}, {-3, 75},
		{-1, 23}, {1, 34}, {1, 43}, {0, 54}, {-2, 55},
		{0, 61}, {1, 64}, {0, 68}, {-9, 92},
		/* 436 to 459 */
		{-14, 106}, {-13, 97}, {-15, 90}, {-12, 90}, {-18, 88},
		{-10, 73}, {-9, 79}, {-14, 86}, {-10, 73}, {-10, 70},
		{-10, 69}, {-5, 66}, {-9, 64}, {-5, 58}, {2, 59},
		{21, -10}, {24, -11}, {28, -8}, {28, -1}, {29, 3},
		{29, 9}, {35, 20}, {29, 36}, {14, 67},
	},
	{
		/* 0 to 10 */
		{20, -15}, {2, 54}, {3, 74}, {20, -15},	{2, 54},
		{3, 74}, {-28, 127}, {-23, 104}, {-6, 53}, {-1, 54},
		{7, 51},
		/* 11 to 59 */
		{23, 33}, {23, 2}, {21, 0}, {1, 9}, {0, 49},
		{-37, 118}, {5, 57}, {-13, 78}, {-11, 65}, {1, 62},
		{12, 49}, {-4, 73}, {17, 50}, {18, 64}, {9,43},
		{29, 0}, {26, 67}, {16, 90}, {9, 104}, {-46, 127},
		{-20, 104}, {1, 67}, {-13, 78}, {-11, 65}, {1, 62},
		{-6, 86}, {-17, 95}, {-6, 61}, {9, 45}, {-3, 69},
		{-6, 81}, {-11, 96}, {6, 55}, {7, 67}, {-5, 86},
		{2, 88}, {0, 58}, {-3, 76}, {-10, 94}, {5, 54},
		{4, 69}, {-3, 81}, {0, 88}, {-7, 67}, {-5, 74},
		{-4, 74}, {-5, 80}, {-7, 72}, {1, 58},
		/* 60 to 69 */
		{0, 41}, {0, 63}, {0, 63}, {0, 63}, {-9, 83},
		{4, 86}, {0, 97}, {-7, 72}, {13, 41}, {3, 62},
		/* 70 to 104 */
		{0, 45}, {-4, 78}, {-3, 96}, {-27, 126}, {-28, 98},
		{-25, 101}, {-23, 67}, {-28, 82}, {-20, 94}, {-16, 83},
		{-22, 110}, {-21, 91}, {-18, 102}, {-13, 93}, {-29, 127},
		{-7, 92}, {-5, 89}, {-7, 96}, {-13, 108}, {-3, 46},
		{-1, 65}, {-1, 57}, {-9, 93}, {-3, 74}, {-9, 92},
		{-8, 87}, {-23, 126}, {5, 54}, {6, 60}, {6, 59},
		{6, 69}, {-1, 48}, {0, 68}, {-4, 69}, {-8, 88},
		/* 105 to 226 */
		{-2, 85}, {-6, 78}, {-1, 75}, {-7, 77}, {2, 54},
		{5, 50}, {-3, 68}, {1, 50}, {6, 42}, {-4, 81},
		{1, 63}, {-4, 70}, {0, 67}, {2, 57}, {-2, 76},
		{11, 35}, {4, 64}, {1, 61}, {11, 35}, {18, 25},
		{12, 24}, {13, 29}, {13, 36}, {-10, 93}, {-7, 73},
		{-2, 73}, {13, 46}, {9, 49}, {-7, 100}, {9, 53},
		{2, 53}, {5, 53}, {-2, 61}, {0, 56}, {0, 56},
		{-13, 63}, {-5, 60}, {-1, 62}, {4, 57}, {-6, 69},
		{4, 57}, {14, 39}, {4, 51}, {13, 68}, {3, 64},
		{1, 61}, {9, 63}, {7, 50}, {16, 39}, {5, 44},
		{4, 52}, {11, 48}, {-5, 60}, {-1, 59}, {0, 59},
		{22, 33}, {5, 44}, {14, 43}, {-1, 78}, {0, 60},
		{9, 69}, {11, 28}, {2, 40}, {3, 44}, {0, 49},
		{0, 46}, {2, 44}, {2, 51}, {0, 47}, {4, 39},
		{2, 62}, {6, 46}, {0, 54}, {3, 54}, {2, 58},
		{4, 63}, {6, 51}, {6, 57}, {7, 53}, {6, 52},
		{6, 55}, {11, 45}, {14, 36}, {8, 53}, {-1, 82},
		{7, 55}, {-3, 78}, {15, 46}, {22, 31}, {-1, 84},
		{25, 7}, {30, -7}, {28, 3}, {28, 4}, {32, 0},
		{34, -1}, {30, 6}, {30, 6}, {32, 9}, {31, 19},
		{26, 27}, {26, 30}, {37, 20}, {28, 34}, {17, 70},
		{1, 67}, {5, 59}, {9, 67}, {16, 30}, {18, 32},
		{18, 35}, {22, 29}, {24, 31}, {23, 38}, {18, 43},
		{20, 41}, {11, 63}, {9, 59}, {9, 64}, {-1, 94},
		{-2, 89}, {-9, 108},
		/* 227 to 275 */
		{-6, 76}, {-2, 44}, {0, 45}, {0, 52}, {-3, 64},
		{-2, 59}, {-4, 70}, {-4, 75}, {-8, 82}, {-17, 102},
		{-9, 77}, {3, 24}, {0, 42}, {0, 48}, {0, 55},
		{-6, 59}, {-7, 71}, {-12, 83}, {-11, 87}, {-30, 119},
		{1, 58}, {-3, 29}, {-1, 36}, {1, 38}, {2, 43},
		{-6, 55}, {0, 58}, {0, 64}, {-3, 74}, {-10, 90},
		{0, 70}, {-4, 29}, {5, 31}, {7, 42}, {1, 59},
		{-2, 58}, {-3, 72}, {-3, 81}, {-11, 97}, {0, 58},
		{8, 5}, {10, 14}, {14, 18}, {13, 27}, {2, 40},
		{0, 58}, {-3, 70}, {-6, 79}, {-8, 85},
		/* 276, void */
		{0, 0},
		/* 277 to 337 */
		{-13, 106}, {-16, 106}, {-10, 87}, {-21, 114}, {-18, 110},
		{-14, 98}, {-22, 110}, {-21, 106}, {-18, 103}, {-21, 107},
		{-23, 108}, {-26, 112}, {-10, 96}, {-12, 95}, {-5, 91},
		{-9, 93}, {-22, 94}, {-5, 86}, {9, 67}, {-4, 80},
		{-10, 85}, {-1, 70}, {7, 60}, {9, 58}, {5, 61},
		{12, 50}, {15, 50}, {18, 49}, {17, 54}, {10, 41},
		{7, 46}, {-1, 51}, {7, 49}, {8, 52}, {9, 41},
		{6, 47}, {2, 55}, {13, 41}, {10, 44}, {6, 50},
		{5, 53}, {13, 49}, {4, 63}, {6, 64}, {-2, 69},
		{-2, 59}, {6, 70}, {10, 44}, {9, 31}, {12, 43},
		{3, 53}, {14, 34}, {10, 38}, {-3, 52}, {13, 40},
		{17, 32}, {7, 44}, {7, 38}, {13, 50}, {10, 57},
		{26, 43},
		/* 338 to 435 */
		{14, 11}, {11, 14}, {9, 11}, {18, 11}, {21, 9},
		{23, -2}, {32, -15}, {32, -15}, {34, -21}, {39, -23},
		{42, -33}, {41, -31}, {46, -28}, {38, -12}, {21, 29},
		{45, -24}, {53, -45}, {48, -26}, {65, -43}, {43, -19},
		{39, -10}, {30, 9}, {18, 26}, {20, 27}, {0, 57},
		{-14, 82}, {-5, 75}, {-19, 97}, {-35, 125}, {27, 0},
		{28, 0}, {31, -4}, {27, 6}, {34, 8}, {30, 10},
		{24, 22}, {33, 19}, {22, 32}, {26, 31}, {21, 41},
		{26, 44}, {23, 47}, {16, 65}, {14, 71}, {8, 60},
		{6, 63}, {17, 65}, {21, 24}, {23, 20}, {26, 23},
		{27, 32}, {28, 23}, {28, 24}, {23, 40}, {24, 32},
		{28, 29}, {23, 42}, {19, 57}, {22, 53}, {22, 61},
		{11, 86}, {12, 40}, {11, 51}, {14, 59}, {-4, 79},
		{-7, 71}, {-5, 69}, {-9, 70}, {-8, 66}, {-10, 68},
		{-19, 73}, {-12, 69}, {-16, 70}, {-15, 67}, {-20, 62},
		{-19, 70}, {-16, 66}, {-22, 65}, {-20, 63}, {9, -2},
		{26, -9}, {33, -9}, {39, -7}, {41, -2}, {45, 3},
		{49, 9}, {45, 27}, {36, 59}, {-6, 66}, {-7, 35},
		{-7, 42}, {-8, 45}, {-5, 48}, {-12, 56}, {-6, 60},
		{-5, 62}, {-8, 66}, {-8, 76},
		/* 436 to 459 */
		{-5, 85}, {-6, 81}, {-10, 77}, {-7, 81}, {-17, 80},
		{-18, 73}, {-4, 74}, {-10, 83}, {-9, 71}, {-9, 67},
		{-1, 61}, {-8, 66}, {-14, 66}, {0, 59}, {2, 59},
		{21, -13}, {33, -14}, {39, -7}, {46, -2}, {51, 2},
		{60, 6}, {61, 17}, {55, 34}, {42, 62},
	},
	{
		/* 0 to 10 */
		{20, -15}, {2, 54}, {3, 74}, {20, -15}, {2, 54},
		{3, 74}, {-28, 127}, {-23, 104}, {-6, 53}, {-1, 54},
		{7, 51},
		/* 11 to 104 */
		{22, 25}, {34, 0}, {16, 0}, {-2, 9}, {4, 41},
		{-29, 118}, {2, 65}, {-6, 71}, {-13, 79}, {5, 52},
		{9, 50}, {-3, 70}, {10, 54}, {26, 34}, {19, 22},
		{40, 0}, {57, 2}, {41, 36}, {26, 69}, {-45, 127},
		{-15, 101}, {-4, 76}, {-6, 71}, {-13, 79}, {5, 52},
		{6, 69}, {-13, 90}, {0, 52}, {8, 43}, {-2, 69},
		{-5, 82}, {-10, 96}, {2, 59}, {2, 75}, {-3, 87},
		{-3, 100}, {1, 56}, {-3, 74}, {-6, 85}, {0, 59},
		{-3, 81}, {-7, 86}, {-5, 95}, {-1, 66}, {-1, 77},
		{1, 70}, {-2, 86}, {-5, 72}, {0, 61}, {0, 41},
		{0, 63}, {0, 63},  {0, 63}, {-9, 83}, {4, 86},
		{0, 97},  {-7, 72}, {13, 41}, {3, 62}, {13, 15},
		{7, 51}, {2, 80}, {-39, 127}, {-18, 91}, {-17, 96},
		{-26, 81}, {-35, 98}, {-24, 102}, {-23, 97}, {-27, 119},
		{-24, 99}, {-21, 110}, {-18, 102}, {-36, 127}, {0, 80},
		{-5, 89}, {-7, 94}, {-4, 92}, {0, 39}, {0, 65},
		{-15, 84}, {-35, 127}, {-2, 73}, {-12, 104}, {-9, 91},
		{-31, 127}, {3, 55}, {7, 56}, {7, 55},
		{8, 61}, {-3, 53}, {0, 68}, {-7, 74}, {-9, 88},
		/* 105 to 165 */
		{-13, 103}, {-13, 91}, {-9, 89}, {-14, 92}, {-8, 76},
		{-12, 87}, {-23, 110}, {-24, 105}, {-10, 78}, {-20, 112},
		{-17, 99}, {-78, 127}, {-70, 127}, {-50, 127}, {-46, 127},
		{-4, 66}, {-5, 78}, {-4, 71}, {-8, 72}, {2, 59},
		{-1, 55}, {-7, 70}, {-6, 75}, {-8, 89}, {-34, 119},
		{-3, 75}, {32, 20}, {30, 22}, {-44, 127}, {0, 54},
		{-5, 61}, {0, 58}, {-1, 60}, {-3, 61}, {-8, 67},
		{-25, 84}, {-14, 74}, {-5, 65}, {5, 52}, {2, 57},
		{0, 61}, {-9, 69}, {-11, 70}, {18, 55}, {-4, 71},
		{0, 58}, {7, 61}, {9, 41}, {18, 25}, {9, 32},
		{5, 43}, {9, 47}, {0, 44}, {0, 51}, {2, 46},
		{19, 38}, {-4, 66}, {15, 38}, {12, 42}, {9, 34},
		{0, 89},
		/* 166 to 226 */
		{4, 45}, {10, 28}, {10, 31}, {33, -11}, {52, -43},
		{18, 15}, {28, 0}, {35, -22}, {38, -25}, {34, 0},
		{39, -18}, {32, -12}, {102, -94}, {0, 0}, {56, -15},
		{33, -4}, {29, 10}, {37, -5}, {51, -29}, {39, -9},
		{52, -34}, {69, -58}, {67, -63}, {44, -5}, {32, 7},
		{55, -29}, {32, 1}, {0, 0}, {27, 36}, {33, -25},
		{34, -30}, {36, -28}, {38, -28}, {38, -27}, {34, -18},
		{35, -16}, {34, -14}, {32, -8}, {37, -6}, {35, 0},
		{30, 10}, {28, 18}, {26, 25}, {29, 41}, {0, 75},
		{2, 72}, {8, 77}, {14, 35}, {18, 31}, {17, 35},
		{21, 30}, {17, 45}, {20, 42}, {18, 45}, {27, 26},
		{16, 54}, {7, 66}, {16, 56}, {11, 73}, {10, 67},
		{-10, 116},
		/* 227 - 275 */
		{-23, 112}, {-15, 71}, {-7, 61}, {0, 53}, {-5, 66},
		{-11, 77}, {-9, 80}, {-9, 84}, {-10, 87}, {-34, 127},
		{-21, 101}, {-3, 39}, {-5, 53}, {-7, 61}, {-11, 75},
		{-15, 77}, {-17, 91}, {-25, 107}, {-25, 111}, {-28, 122},
		{-11, 76}, {-10, 44}, {-10, 52}, {-10, 57}, {-9, 58},
		{-16, 72}, {-7, 69}, {-4, 69}, {-5, 74}, {-9, 86},
		{2, 66}, {-9, 34}, {1, 32}, {11, 31}, {5, 52},
		{-2, 55}, {-2, 67}, {0, 73}, {-8, 89}, {3, 52},
		{7, 4}, {10, 8}, {17, 8}, {16, 19}, {3, 37},
		{-1, 61}, {-5, 73}, {-1, 70}, {-4, 78},
		/* 276, void */
		{0, 0},
		/* 277 to 337 */
		{-21, 126}, {-23, 124}, {-20, 110}, {-26, 126}, {-25, 124},
		{-17, 105}, {-27, 121}, {-27, 117}, {-17, 102}, {-26, 117},
		{-27, 116}, {-33, 122}, {-10, 95}, {-14, 100}, {-8, 95},
		{-17, 111}, {-28, 114}, {-6, 89}, {-2, 80}, {-4, 82},
		{-9, 85}, {-8, 81}, {-1, 72}, {5, 64}, {1, 67},
		{9, 56}, {0, 69}, {1, 69}, {7, 69}, {-7, 69},
		{-6, 67}, {-16, 77}, {-2, 64}, {2, 61}, {-6, 67},
		{-3, 64}, {2, 57}, {-3, 65}, {-3, 66}, {0, 62},
		{9, 51}, {-1, 66}, {-2, 71}, {-2, 75}, {-1, 70},
		{-9, 72}, {14, 60}, {16, 37}, {0, 47}, {18, 35},
		{11, 37}, {12, 41}, {10, 41}, {2, 48}, {12, 41},
		{13, 41}, {0, 59}, {3, 50}, {19, 40}, {3, 66},
		{18, 50},
		/* 338 to 398 */
		{19, -6}, {18, -6}, {14, 0}, {26, -12}, {31, -16},
		{33, -25}, {33, -22}, {37, -28}, {39, -30}, {42, -30},
		{47, -42}, {45, -36}, {49, -34}, {41, -17}, {32, 9},
		{69, -71}, {63, -63}, {66, -64}, {77, -74}, {54, -39},
		{52, -35}, {41, -10}, {36, 0}, {40, -1}, {30, 14},
		{28, 26}, {23, 37}, {12, 55}, {11, 65}, {37, -33},
		{39, -36}, {40, -37}, {38, -30}, {46, -33}, {42, -30},
		{40, -24}, {49, -29}, {38, -12}, {40, -10}, {38, -3},
		{46, -5}, {31, 20}, {29, 30}, {25, 44}, {12, 48},
		{11, 49}, {26, 45}, {22, 22}, {23, 22}, {27, 21},
		{33, 20}, {26, 28}, {30, 24}, {27, 34}, {18, 42},
		{25, 39}, {18, 50}, {12, 70}, {21, 54}, {14, 71},
		{11, 83},
		/* 399 to 435 */
		{25, 32}, {21, 49}, {21, 54}, {-5, 85}, {-6, 81},
		{-10, 77}, {-7, 81},{-17, 80}, {-18, 73}, {-4, 74},
		{-10, 83}, {-9, 71}, {-9, 67}, {-1, 61}, {-8, 66},
		{-14, 66}, {0, 59}, {2, 59}, {17, -10}, {32, -13},
		{42, -9}, {49, -5}, {53, 0}, {64, 3}, {68, 10},
		{66, 27}, {47, 57}, {-5, 71}, {0, 24}, {-1, 36},
		{-2, 42}, {-2, 52}, {-9, 57}, {-6, 63}, {-4, 65},
		{-4, 67}, {-7, 82},
		/* 436 to 459 */
		{-3, 81}, {-3, 76}, {-7, 72}, {-6, 78}, {-12, 72},
		{-14, 68}, {-3, 70}, {-6, 76}, {-5, 66}, {-5, 62},
		{0, 57}, {-4, 61}, {-9, 60}, {1, 54}, {2, 58},
		{17, -10}, {32, -13}, {42, -9}, {49, -5}, {53, 0},
		{64, 3}, {68, 10}, {66, 27}, {47, 57},
	},
	{
		/* 0 to 10 */
		{20, -15}, {2, 54}, {3, 74}, {20, -15}, {2, 54},
		{3, 74}, {-28, 127}, {-23, 104}, {-6, 53}, {-1, 54},
		{7, 51},
		/* 11 to 104 */
		{29, 16}, {25, 0}, {14, 0}, {-10, 51}, {-3, 62},
		{-27, 99}, {26, 16}, {-4, 85}, {-24, 102}, {5, 57},
		{6, 57}, {-17, 73}, {14, 57}, {20, 40}, {20, 10},
		{29, 0}, {54, 0}, {37, 42}, {12, 97}, {-32, 127},
		{-22, 117}, {-2, 74}, {-4, 85}, {-24, 102}, {5, 57},
		{-6, 93}, {-14, 88}, {-6, 44}, {4, 55}, {-11, 89},
		{-15, 103}, {-21, 116}, {19, 57}, {20, 58}, {4, 84},
		{6, 96}, {1, 63}, {-5, 85}, {-13, 106}, {5, 63},
		{6, 75}, {-3, 90}, {-1, 101}, {3, 55}, {-4, 79},
		{-2, 75}, {-12, 97}, {-7, 50}, {1, 60}, {0, 41},
		{0, 63}, {0, 63},  {0, 63}, {-9, 83}, {4, 86},
		{0, 97},  {-7, 72}, {13, 41}, {3, 62}, {7, 34},
		{-9, 88}, {-20, 127}, {-36, 127}, {-17, 91}, {-14, 95},
		{-25, 84}, {-25, 86}, {-12, 89}, {-17, 91}, {-31, 127},
		{-14, 76}, {-18, 103}, {-13, 90}, {-37, 127}, {11, 80},
		{5, 76}, {2, 84}, {5, 78}, {-6, 55}, {4, 61},
		{-14, 83}, {-37, 127}, {-5, 79}, {-11, 104}, {-11, 91},
		{-30, 127}, {0, 65}, {-2, 79}, {0, 72}, {-4, 92},
		{-6, 56}, {3, 68}, {-8, 71}, {-13, 98},
		/* 105 to 165 */
		{-4, 86}, {-12, 88}, {-5, 82}, {-3, 72}, {-4, 67},
		{-8, 72}, {-16, 89}, {-9, 69}, {-1, 59}, {5, 66},
		{4, 57}, {-4, 71}, {-2, 71}, {2, 58}, {-1, 74},
		{-4, 44}, {-1, 69}, {0, 62}, {-7, 51}, {-4, 47},
		{-6, 42}, {-3, 41}, {-6, 53}, {8, 76}, {-9, 78},
		{-11, 83}, {9, 52}, {0, 67}, {-5, 90}, {1, 67},
		{-15, 72}, {-5, 75}, {-8, 80}, {-21, 83}, {-21, 64},
		{-13, 31}, {-25, 64}, {-29, 94}, {9, 75}, {17, 63},
		{-8, 74}, {-5, 35}, {-2, 27}, {13, 91}, {3, 65},
		{-7, 69}, {8, 77}, {-10, 66}, {3, 62}, {-3, 68},
		{-20, 81}, {0, 30}, {1, 7}, {-3, 23}, {-21, 74},
		{16, 66}, {-23, 124}, {17, 37}, {44, -18}, {50, -34},
		{-22, 127},
		/* 166 to 275 */
		{4, 39}, {0, 42}, {7, 34}, {11, 29}, {8, 31},
		{6, 37}, {7, 42}, {3, 40}, {8, 33}, {13, 43},
		{13, 36}, {4, 47}, {3, 55}, {2, 58}, {6, 60},
		{8, 44}, {11, 44}, {14, 42}, {7, 48}, {4, 56},
		{4, 52}, {13, 37}, {9, 49}, {19, 58}, {10, 48},
		{12, 45}, {0, 69}, {20, 33}, {8, 63}, {35, -18},
		{33, -25}, {28, -3}, {24, 10}, {27, 0}, {34, -14},
		{52, -44}, {39, -24}, {19, 17}, {31, 25}, {36, 29},
		{24, 33}, {34, 15}, {30, 20}, {22, 73}, {20, 34},
		{19, 31}, {27, 44}, {19, 16}, {15, 36}, {15, 36},
		{21, 28}, {25, 21}, {30, 20}, {31, 12}, {27, 16},
		{24, 42}, {0, 93}, {14, 56}, {15, 57}, {26, 38},
		{-24, 127}, {-24, 115}, {-22, 82}, {-9, 62}, {0, 53},
		{0, 59}, {-14, 85}, {-13, 89}, {-13, 94}, {-11, 92},
		{-29, 127}, {-21, 100}, {-14, 57}, {-12, 67}, {-11, 71},
		{-10, 77}, {-21, 85}, {-16, 88}, {-23, 104}, {-15, 98},
		{-37, 127}, {-10, 82}, {-8, 48}, {-8, 61}, {-8, 66},
		{-7, 70}, {-14, 75}, {-10, 79}, {-9, 83}, {-12, 92},
		{-18, 108}, {-4, 79}, {-22, 69}, {-16, 75}, {-2, 58},
		{1, 58}, {-13, 78}, {-9, 83}, {-4, 81}, {-13, 99},
		{-13, 81}, {-6, 38}, {-13, 62}, {-6, 58}, {-2, 59},
		{-16, 73}, {-10, 76}, {-13, 86}, {-9, 83}, {-10, 87},
		/* 276, void */
		{0, 0},
		/* 277 to 337 */
		{-22, 127}, {-25, 127}, {-25, 120}, {-27, 127},	{-19, 114},
		{-23, 117}, {-25, 118}, {-26, 117}, {-24, 113}, {-28, 118},
		{-31, 120}, {-37, 124}, {-10, 94}, {-15, 102}, {-10, 99},
		{-13, 106}, {-50, 127}, {-5, 92}, {17, 57}, {-5, 86},
		{-13, 94}, {-12, 91}, {-2, 77}, {0, 71}, {-1, 73},
		{4, 64}, {-7, 81}, {5, 64}, {15, 57}, {1, 67},
		{0, 68}, {-10, 67}, {1, 68}, {0, 77}, {2, 64},
		{0, 68}, {-5, 78}, {7, 55}, {5, 59}, {2, 65},
		{14, 54}, {15, 44}, {5, 60}, {2, 70}, {-2, 76},
		{-18, 86}, {12, 70}, {5, 64}, {-12, 70}, {11, 55},
		{5, 56}, {0, 69}, {2, 65}, {-6, 74}, {5, 54},
		{7, 54}, {-6, 76}, {-11, 82}, {-2, 77}, {-2, 77},
		{25, 42},
		/* 338 to 398 */
		{17, -13}, {16, -9}, {17, -12}, {27, -21}, {37, -30},
		{41, -40}, {42, -41}, {48, -47}, {39, -32}, {46, -40},
		{52, -51}, {46, -41}, {52, -39}, {43, -19}, {32, 11},
		{61, -55}, {56, -46}, {62, -50}, {81, -67}, {45, -20},
		{35, -2}, {28, 15}, {34, 1}, {39, 1}, {30, 17},
		{20, 38}, {18, 45}, {15, 54}, {0, 79}, {36, -16},
		{37, -14}, {37, -17}, {32, 1}, {34, 15}, {29, 15},
		{24, 25}, {34, 22}, {31, 16}, {35, 18}, {31, 28},
		{33, 41}, {36, 28}, {27, 47}, {21, 62}, {18, 31},
		{19, 26}, {36, 24}, {24, 23}, {27, 16}, {24, 30},
		{31, 29}, {22, 41}, {22, 42}, {16, 60}, {15, 52},
		{14, 60}, {3, 78}, {-16, 123}, {21, 53}, {22, 56},
		{25, 61},
		/* 399 to 435 */
		{21, 33}, {19, 50}, {17, 61}, {-3, 78}, {-8, 74},
		{-9, 72}, {-10, 72}, {-18, 75}, {-12, 71}, {-11, 63},
		{-5, 70}, {-17, 75}, {-14, 72}, {-16, 67}, {-8, 53},
		{-14, 59}, {-9, 52}, {-11, 68}, {9, -2}, {30, -10},
		{31, -4}, {33, -1}, {33, 7}, {31, 12}, {37, 23},
		{31, 38}, {20, 64}, {-9, 71}, {-7, 37}, {-8, 44},
		{-11, 49}, {-10, 56}, {-12, 59}, {-8, 63}, {-9, 67},
		{-6, 68}, {-10, 79},
		/* 436 to 459 */
		{-3, 78}, {-8, 74}, {-9, 72}, {-10, 72}, {-18, 75},
		{-12, 71}, {-11, 63}, {-5, 70}, {-17, 75}, {-14, 72},
		{-16, 67}, {-8, 53}, {-14, 59}, {-9, 52}, {-11, 68},
		{9, -2}, {30, -10}, {31, -4}, {33, -1}, {33, 7},
		{31, 12}, {37, 23}, {31, 38}, {20, 64}
	}
};

static void init_cabac_context(h264d_cabac_t *cabac, int slice_qp, int idc)
{
	const ctx_idx_mn_t *lut = ctx_idx_mn_IPB[idc];
	int8_t *ctx = cabac->context;
	int i = 460;
	do {
		int pre_state = ((lut->m * slice_qp) >> 4) + lut->n;
		if (pre_state < 64) {
			/* valMPS == 0 */
			pre_state = (pre_state <= 0) ? 1 : pre_state;
			*ctx = (63 - pre_state) * 2; /* -1..-63 ==> 124, 122,..0 */
		} else {
			/* valMPS == 1 */
			pre_state = (126 < pre_state) ? 126 : pre_state;
			*ctx = (pre_state - 64) * 2 + 1; /* 0..62 ==> 1, 3,..125 */
		}
		ctx++;
		lut++;
	} while (--i);
}


static inline void cabac_renorm(h264d_cabac_t *cb, dec_bits *st, int range, int offset)
{
	if (256 <= range) {
		cb->range = range;
		return;
	}
	do {
		range = (uint16_t)(range * 2);
		offset = (uint16_t)(offset * 2 + get_onebit_inline(st));
	} while (range < 256);
	cb->range = range;
	cb->offset = offset;
}

static inline int cabac_decode_decision(h264d_cabac_t *cb, dec_bits *st, int ctxIdx)
{
	static const uint8_t rangeTabLPS[64][4] = {
		{128, 176, 208, 240}, {128, 167, 197, 227},
		{128, 158, 187, 216}, {123, 150, 178, 205},
		{116, 142, 169, 195}, {111, 135, 160, 185},
		{105, 128, 152, 175}, {100, 122, 144, 166},
		{95, 116, 137, 158}, {90, 110, 130, 150},
		{85, 104, 123, 142}, {81, 99, 117, 135},
		{77, 94, 111, 128}, {73, 89, 105, 122},
		{69, 85, 100, 116}, {66, 80, 95, 110},
		{62, 76, 90, 104}, {59, 72, 86, 99},
		{56, 69, 81, 94}, {53, 65, 77, 89},
		{51, 62, 73, 85}, {48, 59, 69, 80},
		{46, 56, 66, 76}, {43, 53, 63, 72},
		{41, 50, 59, 69}, {39, 48, 56, 65},
		{37, 45, 54, 62}, {35, 43, 51, 59},
		{33, 41, 48, 56}, {32, 39, 46, 53},
		{30, 37, 43, 50}, {29, 35, 41, 48},
		{27, 33, 39, 45}, {26, 31, 37, 43},
		{24, 30, 35, 41}, {23, 28, 33, 39},
		{22, 27, 32, 37}, {21, 26, 30, 35},
		{20, 24, 29, 33}, {19, 23, 27, 31},
		{18, 22, 26, 30}, {17, 21, 25, 28},
		{16, 20, 23, 27}, {15, 19, 22, 25},
		{14, 18, 21, 24}, {14, 17, 20, 23},
		{13, 16, 19, 22}, {12, 15, 18, 21},
		{12, 14, 17, 20}, {11, 14, 16, 19},
		{11, 13, 15, 18}, {10, 12, 15, 17},
		{10, 12, 14, 16}, {9 ,11, 13, 15},
		{9, 11, 12, 14}, {8, 10, 12, 14},
		{8, 9, 11, 13}, {7, 9, 11, 12},
		{7, 9, 10, 12}, {7, 8, 10, 11},
		{6, 8, 9, 11}, {6, 7, 9, 10},
		{6, 7, 8, 9}, {2, 2, 2, 2}
	};
	static const int8_t state_trans[2][64] = {
		{
			1, 2, 3, 4, 5, 6, 7, 8,
			9, 10, 11, 12, 13, 14, 15, 16,
			17, 18, 19, 20, 21, 22, 23, 24,
			25, 26, 27, 28, 29, 30, 31, 32,
			33, 34, 35, 36, 37, 38, 39, 40,
			41, 42, 43, 44, 45, 46, 47, 48,
			49, 50, 51, 52, 53, 54, 55, 56,
			57, 58, 59, 60, 61, 62, 62, 63
		},
		{
			0, 0, 1, 2, 2, 4, 4, 5,
			6, 7, 8, 9, 9, 11, 11, 12,
			13, 13, 15, 15, 16, 16, 18, 18,
			19, 19, 21, 21, 22, 22, 23, 24,
			24, 25, 26, 26, 27, 27, 28, 29,
			29, 30, 30, 30, 31, 32, 32, 33,
			33, 33, 34, 34, 35, 35, 35, 36,
			36, 36, 37, 37, 37, 38, 38, 63
		}
	};
	int pStateIdx;
	uint32_t valMPS, binVal;
	uint32_t range, offset, lps;
	int is_lps;

	pStateIdx = cb->context[ctxIdx];
	valMPS = pStateIdx & 1;
	pStateIdx >>= 1;
	range = cb->range;
	offset = cb->offset;
	lps = rangeTabLPS[pStateIdx][(range >> 6) & 3];
	range = range - lps;
	if (offset < range) {
		binVal = valMPS;
		is_lps = 0;
	} else {
		binVal = valMPS ^ 1;
		cb->offset = offset = offset - range;
		range = lps;
		valMPS = pStateIdx ? valMPS : valMPS ^ 1;
		is_lps = 1;
	}
	pStateIdx = state_trans[is_lps][pStateIdx] * 2 + valMPS;
	cb->context[ctxIdx] = pStateIdx;
	cabac_renorm(cb, st, range, offset);
	return binVal;
}

static inline int cabac_decode_bypass(h264d_cabac_t *cb, dec_bits *st)
{
	int range, offset;
	range = cb->range;
	offset = (cb->offset << 1) | get_onebit_inline(st);
	if (offset < range) {
		cb->offset = offset;
		return 0;
	} else {
		cb->offset = offset - range;
		return 1;
	}
}

static inline int cabac_decode_terminate(h264d_cabac_t *cb, dec_bits *st)
{
	int range = cb->range - 2;
	int offset = cb->offset;
	if (range <= offset) {
		cb->range = range;
		return 1;
	} else {
		cabac_renorm(cb, st, range, offset);
		return 0;
	}
}

static int mb_type_cabac_I(h264d_mb_current *mb, dec_bits *st, int avail, int ctx_idx, int slice_type)
{
	h264d_cabac_t *cb = mb->cabac;
	int is_i_slice = (slice_type == I_SLICE);
	int mb_type;

	if (is_i_slice) {
		int add = ((avail & 2) && (mb->top4x4inter->type != MB_INxN)) + ((avail & 1) && (mb->left4x4inter->type != MB_INxN));
		if (!cabac_decode_decision(cb, st, ctx_idx + add)) {
			return MB_INxN;
		}
		ctx_idx = 5;
	} else if (!cabac_decode_decision(cb, st, ctx_idx)) {
		return MB_INxN;
	}
	if (cabac_decode_terminate(cb, st)) {
		return MB_IPCM;
	}
	mb_type = cabac_decode_decision(cb, st, ctx_idx + 1) * 12 + 1;
	if (cabac_decode_decision(cb, st, ctx_idx + 2)) {
		mb_type = mb_type + cabac_decode_decision(cb, st, ctx_idx + 2 + is_i_slice) * 4 + 4;
	}
	mb_type += cabac_decode_decision(cb, st, ctx_idx + 3 + is_i_slice) * 2;
	mb_type += cabac_decode_decision(cb, st, ctx_idx + 3 + is_i_slice * 2);
	return mb_type;
}

static int mb_type_cabac_P(h264d_mb_current *mb, dec_bits *st, int avail, int ctx_idx, int slice_type)
{
	h264d_cabac_t *cb = mb->cabac;
	if (cabac_decode_decision(cb, st, 14)) {
		return 5 + mb_type_cabac_I(mb, st, avail, 17, 0);
	}
	if (cabac_decode_decision(cb, st, 15)) {
		return cabac_decode_decision(cb, st, 17) ? 1 : 2;
	} else {
		return cabac_decode_decision(cb, st, 16) ? 3 : 0;
	}
}

static int mb_type_cabac_B(h264d_mb_current *mb, dec_bits *st, int avail, int ctx_idx, int slice_type)
{
	h264d_cabac_t *cb = mb->cabac;
	int add = ((avail & 1) && (mb->left4x4inter->type != MB_BDIRECT16x16)) + ((avail & 2) && (mb->top4x4inter->type != MB_BDIRECT16x16));
	int mode;

	if (!cabac_decode_decision(cb, st, 27 + add)) {
		return 0;
	}
	if (!cabac_decode_decision(cb, st, 27 + 3)) {
		return 1 + cabac_decode_decision(cb, st, 27 + 5);
	}
	mode = cabac_decode_decision(cb, st, 27 + 4) * 8;
	mode += cabac_decode_decision(cb, st, 27 + 5) * 4;
	mode += cabac_decode_decision(cb, st, 27 + 5) * 2;
	mode += cabac_decode_decision(cb, st, 27 + 5);
	if (mode < 8) {
		return mode + 3;
	} else if (mode < 13) {
		return mode * 2 + cabac_decode_decision(cb, st, 27 + 5) - 4;
	} else if (mode == 13) {
		return 23 + mb_type_cabac_I(mb, st, avail, 32, 0);
	} else if (mode == 14) {
		return 11;
	} else {
		/* mode == 15 */
		return 22;
	}
}

static int mb_skip_cabac(h264d_mb_current *mb, dec_bits *st, int slice_type)
{
	int avail = get_availability(mb);
	int offset = (slice_type == P_SLICE) ? 11 : 24;
	if ((avail & 1) && (mb->left4x4inter->mb_skip == 0)) {
		offset += 1;
	}
	if ((avail & 2) && (mb->top4x4inter->mb_skip == 0)) {
		offset += 1;
	}
	return cabac_decode_decision(mb->cabac, st, offset);
}

struct intra4x4pred_mode_cabac {
	int operator()(int a, int b, dec_bits *st, h264d_cabac_t *cb) const {
		int pred = MIN(a, b);
		if (!cabac_decode_decision(cb, st, 68)) {
			int rem;
			rem = cabac_decode_decision(cb, st, 69);
			rem += cabac_decode_decision(cb, st, 69) * 2;
			rem += cabac_decode_decision(cb, st, 69) * 4;
			pred = (rem < pred) ? rem : rem + 1;
		}
		return pred;
	}
};

struct intra_chroma_pred_mode_cabac {
	uint32_t operator()(h264d_mb_current *mb, dec_bits *st, int avail) const {
		h264d_cabac_t *cb = mb->cabac;
		int ctx_idx = 64 + ((avail & 2) && (mb->top4x4inter->type < MB_IPCM) && mb->top4x4inter->chroma_pred_mode) + ((avail & 1) && (mb->left4x4inter->type < MB_IPCM) && mb->left4x4inter->chroma_pred_mode);
		int pred_mode = cabac_decode_decision(cb, st, ctx_idx);
		if (pred_mode) {
			while ((pred_mode < 3) && cabac_decode_decision(cb, st, 64 + 3)) {
				pred_mode++;
			}
		}
		mb->chroma_pred_mode = pred_mode;
		return pred_mode;
	}
};

struct cbp_cabac {
	uint32_t operator()(h264d_mb_current *mb, dec_bits *st, int avail) const {
		int cbp;
		int inc;
		h264d_cabac_t *cb = mb->cabac;
		int cbp_a = (avail & 1) ? mb->left4x4inter->cbp : 0x0f;
		int cbp_b = (avail & 2) ? mb->top4x4inter->cbp : 0x0f;
		/* luma */
		inc = (!(cbp_a & 2)) + (!(cbp_b & 4)) * 2;
		cbp = cabac_decode_decision(cb, st, 73 + inc);
		inc = !(cbp & 1) + (!(cbp_b & 8)) * 2;
		cbp += cabac_decode_decision(cb, st, 73 + inc) * 2;
		inc = (!(cbp_a & 8)) + !(cbp & 1) * 2;
		cbp += cabac_decode_decision(cb, st, 73 + inc) * 4;
		inc = !(cbp & 4) + !(cbp & 2) * 2;
		cbp += cabac_decode_decision(cb, st, 73 + inc) * 8;
		/* chroma */
		cbp_a >>= 4;
		cbp_b >>= 4;
		inc = (cbp_a != 0) + (cbp_b != 0) * 2;
		if (cabac_decode_decision(cb, st, 77 + inc)) {
			inc = (cbp_a >> 1) + (cbp_b & 2);
			cbp = cbp + cabac_decode_decision(cb, st, 77 + 4 + inc) * 16 + 16;
		}
		return cbp;
	}
};

static inline uint32_t unary_cabac(h264d_cabac_t *cb, dec_bits *st, int limit)
{
	int x = 0;
	int idx = 62;
	do {
		if (cabac_decode_decision(cb, st, idx)) {
			x = x + 1;
			idx = 63;
		} else {
			break;
		}
	} while (--limit);
	return x;
}

struct qp_delta_cabac {
	int operator()(h264d_mb_current *mb, dec_bits *st, int avail) const {
		int ctx_idx = 60 + (mb->prev_qp_delta != 0);
		h264d_cabac_t *cb = mb->cabac;
		int qp_delta = cabac_decode_decision(cb, st, ctx_idx);
		if (qp_delta) {
			qp_delta = unary_cabac(cb, st, 52) + 1;
			qp_delta = (((qp_delta & 1) ? qp_delta : -qp_delta) + 1) >> 1;
		}
		mb->prev_qp_delta = qp_delta;
		return qp_delta;
	}
};

static int ctxidxinc_cbf0(h264d_mb_current *mb, uint32_t cbf, int avail)
{
	int ab;
	if (avail & 1) {
		ab = mb->left4x4inter->cbf & 1;
	} else {
		ab = (mb->type < MB_IPCM);
	}
	if (avail & 2) {
		ab += (mb->top4x4inter->cbf & 1) * 2;
	} else {
		ab += (mb->type < MB_IPCM) * 2;
	}
	return ab;
}

static int ctxidxinc_cbf1(h264d_mb_current *mb, uint32_t cbf, int avail)
{
	int ab = cbf & 1;
	if (avail & 2) {
		ab += mb->top4x4inter->cbf & 2;
	} else {
		ab += (mb->type < MB_IPCM) * 2;
	}
	return ab;
}

static int ctxidxinc_cbf2(h264d_mb_current *mb, uint32_t cbf, int avail)
{
	int ab;
	if (avail & 1) {
		ab = (mb->left4x4inter->cbf >> 1) & 1;
	} else {
		ab = (mb->type < MB_IPCM);
	}
	ab += (cbf * 2) & 2;
	return ab;
}

template <int N>
static int ctxidxinc_cbf_inner3(h264d_mb_current *mb, uint32_t cbf, int avail)
{
	return ((cbf >> (N + 2)) & 1) | ((cbf >> N) & 2);
}

static int ctxidxinc_cbf4(h264d_mb_current *mb, uint32_t cbf, int avail)
{
	int ab = (cbf >> 1) & 1;
	if (avail & 2) {
		ab += (mb->top4x4inter->cbf >> 1) & 2;
	} else {
		ab += (mb->type < MB_IPCM) * 2;
	}
	return ab;
}

static int ctxidxinc_cbf5(h264d_mb_current *mb, uint32_t cbf, int avail)
{
	int ab = (cbf >> 4) & 1;
	if (avail & 2) {
		ab += (mb->top4x4inter->cbf >> 2) & 2;
	} else {
		ab += (mb->type < MB_IPCM) * 2;
	}
	return ab;
}

static int ctxidxinc_cbf6(h264d_mb_current *mb, uint32_t cbf, int avail)
{
	return (cbf >> 3) & 3;
}

static int ctxidxinc_cbf8(h264d_mb_current *mb, uint32_t cbf, int avail)
{
	int ab;
	if (avail & 1) {
		ab = (mb->left4x4inter->cbf >> 2) & 1;
	} else {
		ab = (mb->type < MB_IPCM);
	}
	ab += (cbf >> 1) & 2;
	return ab;
}

static int ctxidxinc_cbf9(h264d_mb_current *mb, uint32_t cbf, int avail)
{
	return ((cbf >> 8) & 1) | ((cbf >> 2) & 2);
}

static int ctxidxinc_cbf10(h264d_mb_current *mb, uint32_t cbf, int avail)
{
	int ab;
	if (avail & 1) {
		ab = (mb->left4x4inter->cbf >> 3) & 1;
	} else {
		ab = (mb->type < MB_IPCM);
	}
	ab += (cbf >> 7) & 2;
	return ab;
}

static int ctxidxinc_cbf12(h264d_mb_current *mb, uint32_t cbf, int avail)
{
	return ((cbf >> 9) & 1) | ((cbf >> 5) & 2);
}

static int ctxidxinc_cbf13(h264d_mb_current *mb, uint32_t cbf, int avail)
{
	return ((cbf >> 12) & 1) | ((cbf >> 6) & 2);
}

static int ctxidxinc_cbf14(h264d_mb_current *mb, uint32_t cbf, int avail)
{
	return (cbf >> 11) & 3;
}

template <int N>
static int ctxidxinc_cbf_chroma_dc(h264d_mb_current *mb, uint32_t cbf, int avail)
{
	int ab;
	if (avail & 1) {
		if (mb->left4x4inter->type == MB_IPCM) {
			ab = 1;
		} else {
			ab = (mb->left4x4inter->cbf >> (4 + N)) & 1;
		}
	} else {
		ab = (mb->type < MB_IPCM);
	}
	if (avail & 2) {
		if (mb->top4x4inter->type == MB_IPCM) {
			ab += 2;
		} else {
			ab += (mb->top4x4inter->cbf >> (3 + N)) & 2;
		}
	} else {
		ab += (mb->type < MB_IPCM) * 2;
	}
	return ab;
}

template <int N>
static int ctxidxinc_cbf_chroma_ac0(h264d_mb_current *mb, uint32_t cbf, int avail)
{
	int ab;
	if (avail & 1) {
		ab = (mb->left4x4inter->cbf >> (6 + N * 2)) & 1;
	} else {
		ab = (mb->type < MB_IPCM);
	}
	if (avail & 2) {
		ab += (mb->top4x4inter->cbf >> (5 + N * 2)) & 2;
	} else {
		ab += (mb->type < MB_IPCM) * 2;
	}
	return ab;
}

template <int N>
static int ctxidxinc_cbf_chroma_ac1(h264d_mb_current *mb, uint32_t cbf, int avail)
{
	int ab = (cbf >> (18 + N * 4)) & 1;
	if (avail & 2) {
		ab += (mb->top4x4inter->cbf >> (6 + N * 2)) & 2;
	} else {
		ab += (mb->type < MB_IPCM) * 2;
	}
	return ab;
}

template <int N>
static int ctxidxinc_cbf_chroma_ac2(h264d_mb_current *mb, uint32_t cbf, int avail)
{
	int ab = (cbf >> (17 + N * 4)) & 2;
	if (avail & 1) {
		ab += (mb->left4x4inter->cbf >> (7 + N * 2)) & 1;
	} else {
		ab += (mb->type < MB_IPCM);
	}
	return ab;
}

static int ctxidxinc_cbf_intra16x16dc(h264d_mb_current *mb, uint32_t cbf, int avail)
{
	int inc;
	if (avail & 1) {
		inc = (mb->left4x4inter->cbf >> 10) & 1;
	} else {
		inc = 1;
	}
	if (avail & 2) {
		inc += (mb->top4x4inter->cbf >> 9) & 2;
	} else {
		inc += 2;
	}
	return inc;
}

static int (* const ctxidxinc_cbf[16 + 2 + 8 + 1])(h264d_mb_current *mb, uint32_t cbf, int avail) = {
	ctxidxinc_cbf0, ctxidxinc_cbf1,
	ctxidxinc_cbf2, ctxidxinc_cbf_inner3<0>,
	ctxidxinc_cbf4, ctxidxinc_cbf5,
	ctxidxinc_cbf6, ctxidxinc_cbf_inner3<4>,

	ctxidxinc_cbf8, ctxidxinc_cbf9,
	ctxidxinc_cbf10, ctxidxinc_cbf_inner3<8>,
	ctxidxinc_cbf12, ctxidxinc_cbf13,
	ctxidxinc_cbf14, ctxidxinc_cbf_inner3<12>,

	ctxidxinc_cbf_chroma_dc<0>, ctxidxinc_cbf_chroma_dc<1>,

	ctxidxinc_cbf_chroma_ac0<0>, ctxidxinc_cbf_chroma_ac1<0>,
	ctxidxinc_cbf_chroma_ac2<0>, ctxidxinc_cbf_inner3<18>,
	ctxidxinc_cbf_chroma_ac0<1>, ctxidxinc_cbf_chroma_ac1<1>,
	ctxidxinc_cbf_chroma_ac2<1>, ctxidxinc_cbf_inner3<22>,

	ctxidxinc_cbf_intra16x16dc
};

static inline int get_coeff_map_cabac(h264d_cabac_t *cb, dec_bits *st, int cat, int num_coeff, int *coeff_map)
{
	static const int16_t significant_coeff_flag_offset[6][2] = {
		{105, 166}, {105 + 15, 166 + 15}, {105 + 29, 166 + 29},
		{105 + 44, 166 + 44}, {105 + 47, 166 + 47}, {402, 417}
	};
	int sigc_offset = significant_coeff_flag_offset[cat][0];
	int last_offset = significant_coeff_flag_offset[cat][1];
	int map_cnt = 0;
	int i;

	for (i = 0; i < num_coeff - 1; ++i) {
		if (cabac_decode_decision(cb, st, sigc_offset + i)) {
			coeff_map[map_cnt++] = i;
			if (cabac_decode_decision(cb, st, last_offset + i)) {
				i = 0;
				break;
			}
		}
	}
	if (i == num_coeff - 1) {
		coeff_map[map_cnt++] = i;
	}
	return map_cnt;
}

static inline void get_coeff_from_map_cabac(h264d_cabac_t *cb, dec_bits *st, int cat, int *coeff_map, int map_cnt, int *coeff, const int16_t *qmat, uint32_t dc_mask)
{
	static const int16_t coeff_abs_level_offset[6] = {
		227, 227 + 10, 227 + 20, 227 + 30, 227 + 39, 426
	};
	static const int8_t coeff_abs_level_ctx[2][8] = {
		{1, 2, 3, 4, 0, 0, 0, 0},
		{5, 5, 5, 5, 6, 7, 8, 9}
	};
	static const int8_t coeff_abs_level_transition[2][8] = {
		{1, 2, 3, 3, 4, 5, 6, 7},
		{4, 4, 4, 4, 5, 6, 7, 7}
	};
	int abs_offset = coeff_abs_level_offset[cat];
	int node_ctx = 0;
	int mp = map_cnt;
	do {
		int ctx = abs_offset + coeff_abs_level_ctx[0][node_ctx];
		int abs_level;
		int idx;
		if (!cabac_decode_decision(cb, st, ctx)) {
			abs_level = 1;
			node_ctx = coeff_abs_level_transition[0][node_ctx];
		} else {
			abs_level = 2;
			ctx = abs_offset + coeff_abs_level_ctx[1][node_ctx];
			node_ctx = coeff_abs_level_transition[1][node_ctx];
			while (abs_level < 15 && cabac_decode_decision(cb, st, ctx)) {
				abs_level++;
			}
			if (abs_level == 15) {
				int left_cnt = 0;
				while (cabac_decode_bypass(cb, st)) {
					left_cnt++;
				}
				abs_level = 1;
				while (left_cnt--) {
					abs_level += abs_level + cabac_decode_bypass(cb, st);
				}
				abs_level += 14;
			}
		}
		idx = coeff_map[--mp];
		coeff[idx] = (cabac_decode_bypass(cb, st) ? -abs_level : abs_level) * qmat[idx & dc_mask];
	} while (mp);
}

struct residual_block_cabac {
	int operator()(h264d_mb_current *mb, int na, int nb, dec_bits *st, int *coeff, int num_coeff, const int16_t *qmat, int avail, int pos4x4, int cat, uint32_t dc_mask) const {
		int coeff_map[8 * 8];
		h264d_cabac_t *cb = mb->cabac;
		if (cat != 5) {
			int coded_block_flag;
			int coded_flag_inc;
			coded_flag_inc = ctxidxinc_cbf[pos4x4](mb, mb->cbf, avail);
			coded_block_flag = cabac_decode_decision(cb, st, 85 + coded_flag_inc + cat * 4);
			if (!coded_block_flag) {
				return 0;
			}
			mb->cbf |= coded_block_flag << pos4x4;
		}
		int coeff_offset = (dc_mask >> 4);
		memset(coeff + coeff_offset, 0, sizeof(*coeff) * num_coeff);
		int map_cnt = get_coeff_map_cabac(cb, st, cat, num_coeff, coeff_map);
		get_coeff_from_map_cabac(cb, st, cat, coeff_map, map_cnt, coeff + coeff_offset, qmat + coeff_offset, dc_mask);
		return map_cnt <= 15 ? map_cnt : 15;
	}
};

static int mb_intra4x4_cabac(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_intra4x4(mb, mbc, st, avail, intra4x4pred_mode_cabac(), intra_chroma_pred_mode_cabac(), cbp_cabac(), qp_delta_cabac(), residual_block_cabac());
} 

static int mb_intra16x16_dconly_cabac(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_intra16x16_dconly(mb, mbc, st, avail, intra_chroma_pred_mode_cabac(), qp_delta_cabac(), residual_block_cabac());
} 

static int mb_intra16x16_acdc_cabac(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_intra16x16_acdc(mb, mbc, st, avail, intra_chroma_pred_mode_cabac(), qp_delta_cabac(), residual_block_cabac());
} 

struct sub_mb_type_p_cabac {
	int operator()(h264d_mb_current *mb, dec_bits *st, int8_t *sub_mb_type, prev8x8_t *curr_blk, int avail) const {
		h264d_cabac_t *cb = mb->cabac;
		for (int i = 0; i < 4; ++i) {
			int t;
			if (cabac_decode_decision(cb, st, 21)) {
				t = 0;
			} else if (!cabac_decode_decision(cb, st, 22)) {
				t = 1;
			} else if (cabac_decode_decision(cb, st, 23)) {
				t = 2;
			} else {
				t = 3;
			}
			sub_mb_type[i] = t;
		}
		return 0;
	}
};


static inline int sub_mb_type_b_one_cabac(h264d_cabac_t *cb, dec_bits *st)
{
	int t;
	if (!cabac_decode_decision(cb, st, 36)) {
		return 0;
	} else if (!cabac_decode_decision(cb, st, 37)) {
		return 1 + cabac_decode_decision(cb, st, 39);
	} else 	if (cabac_decode_decision(cb, st, 38)) {
		if (cabac_decode_decision(cb, st, 39)) {
			return 11 + cabac_decode_decision(cb, st, 39);
		} else {
			t = 7;
		}
	} else {
		t = 3;
	}
	t += cabac_decode_decision(cb, st, 39) * 2;
	return t + cabac_decode_decision(cb, st, 39);
}

struct sub_mb_type_b_cabac {
	int operator()(h264d_mb_current *mb, dec_bits *st) const {
		return sub_mb_type_b_one_cabac(mb->cabac, st);
	}
};

struct sub_mb_types_b_cabac {
	int operator()(h264d_mb_current *mb, dec_bits *st, int8_t *sub_mb_type, prev8x8_t *curr_blk, int avail) const {
		return sub_mb_type_b_base(mb, st, sub_mb_type, curr_blk, avail, sub_mb_type_b_cabac());
	}
};

static inline int mvd_cabac(h264d_mb_current *mb, dec_bits *st, h264d_cabac_t *cb, int ctx, int mva, int mvb)
{
	int mvd;
	int sum;
	int inc;

	sum = ABS(mva) + ABS(mvb);
	if (sum < 3) {
		inc = 0;
	} else if (sum <= 32) {
		inc = 1;
	} else {
		inc = 2;
	}
	if (!cabac_decode_decision(cb, st, ctx + inc)) {
		return 0;
	}
	mvd = 1;
	ctx += 3;
	while (cabac_decode_decision(cb, st, ctx)) {
		ctx += (mvd < 4) ? 1 : 0;
		mvd += 1;
		if (9 <= mvd) {
			unsigned exp = 3;
			while (cabac_decode_bypass(cb, st) && (exp < sizeof(mvd) * 4)) {
				mvd += (1 << exp);
				exp += 1;
			}
			while(exp--) {
				mvd += cabac_decode_bypass(cb, st) << exp;
			}
			break;
		}
	}
	return cabac_decode_bypass(cb, st) ? -mvd : mvd;
}

struct mvd_xy_cabac {
	void operator()(h264d_mb_current *mb, dec_bits *st, int16_t mv[], const int16_t mva[], const int16_t mvb[]) const {
		h264d_cabac_t *cb = mb->cabac;
		mv[0] = mvd_cabac(mb, st, cb, 40, mva[0], mvb[0]);
		mv[1] = mvd_cabac(mb, st, cb, 47, mva[1], mvb[1]);
	}
};

static void sub_mb8x8_mv_cabac(h264d_mb_current *mb, dec_bits *st, int avail, int blk_idx, prev8x8_t *pblk, int lx)
{
	sub_mb8x8_mv(mb, st, avail, blk_idx, pblk, lx, mvd_xy_cabac());
}

static void sub_mb8x4_mv_cabac(h264d_mb_current *mb, dec_bits *st, int avail, int blk_idx, prev8x8_t *pblk, int lx)
{
	sub_mb8x4_mv(mb, st, avail, blk_idx, pblk, lx, mvd_xy_cabac());
}

static void sub_mb4x8_mv_cabac(h264d_mb_current *mb, dec_bits *st, int avail, int blk_idx, prev8x8_t *pblk, int lx)
{
	sub_mb4x8_mv(mb, st, avail, blk_idx, pblk, lx, mvd_xy_cabac());
}

static void sub_mb4x4_mv_cabac(h264d_mb_current *mb, dec_bits *st, int avail, int blk_idx, prev8x8_t *pblk, int lx)
{
	sub_mb4x4_mv(mb, st, avail, blk_idx, pblk, lx, mvd_xy_cabac());
}

static void (* const sub_mb_p_cabac[4])(h264d_mb_current *mb, dec_bits *st, int avail, int blk_idx, prev8x8_t *pblk, int lx) = {
	sub_mb8x8_mv_cabac,
	sub_mb8x4_mv_cabac,
	sub_mb4x8_mv_cabac,
	sub_mb4x4_mv_cabac
};

static void (* const sub_mb_b_cabac[13])(h264d_mb_current *mb, dec_bits *st, int avail, int blk_idx, prev8x8_t *pblk, int lx) = {
	sub_mb8x8_direct,
	sub_mb8x8_mv_cabac,
	sub_mb8x8_mv_cabac,
	sub_mb8x8_mv_cabac,
	sub_mb8x4_mv_cabac,
	sub_mb4x8_mv_cabac,
	sub_mb8x4_mv_cabac,
	sub_mb4x8_mv_cabac,
	sub_mb8x4_mv_cabac,
	sub_mb4x8_mv_cabac,
	sub_mb4x4_mv_cabac,
	sub_mb4x4_mv_cabac,
	sub_mb4x4_mv_cabac
};

struct sub_mbs_p_cabac {
	void operator()(h264d_mb_current *mb, dec_bits *st, int avail, int8_t sub_mb_type[], prev8x8_t curr_blk[], int lx) const {
		if (lx == 0) {
			for (int i = 0; i < 4; ++i) {
				sub_mb_p_cabac[sub_mb_type[i]](mb, st, avail, i, curr_blk, lx);
			}
		}
	}
};

struct sub_mbs_b_cabac {
	void operator()(h264d_mb_current *mb, dec_bits *st, int avail, int8_t sub_mb_type[], prev8x8_t curr_blk[], int lx) const {
		for (int i = 0; i < 4; ++i) {
			sub_mb_b_cabac[sub_mb_type[i]](mb, st, avail, i, curr_blk, lx);
		}
	}
};

static inline int ref_idx_cabac_sub(dec_bits *st, h264d_cabac_t *cb, int inc)
{
	int idx = 0;
	while (cabac_decode_decision(cb, st, 54 + inc)) {
		inc = ((unsigned)inc >> 2) + 4;
		idx += 1;
	}
	return idx;
}

struct ref_idx16x16_cabac {
	int operator()(h264d_mb_current *mb, dec_bits *st, int lx, int avail) const {
		h264d_cabac_t *cb = mb->cabac;
		if (*mb->num_ref_idx_lx_active_minus1[lx]) {
			int inc = ((avail & 1) && !(mb->left4x4inter->direct8x8 & 1) && (0 < mb->left4x4inter->ref[0][lx])) + ((avail & 2) && !(mb->top4x4inter->direct8x8 & 1) && (0 < mb->top4x4inter->ref[0][lx])) * 2;
			return ref_idx_cabac_sub(st, cb, inc);
		} else {
			return 0;
		}
	}
};

struct ref_idx16x8_cabac {
	void operator()(h264d_mb_current *mb, dec_bits *st, int8_t *ref_idx, uint32_t blk_map, int avail) const {
		h264d_cabac_t *cb = mb->cabac;
		int8_t * const *num = mb->num_ref_idx_lx_active_minus1;
		for (int lx = 0; lx < 2; ++lx) {
			int t = *(num[0]);
			ref_idx[0] = (blk_map & 1) ?
				(t ? ref_idx_cabac_sub(st, cb,
							((avail & 1) && !(mb->left4x4inter->direct8x8 & 1) && (0 < mb->left4x4inter->ref[0][lx]))
							+ ((avail & 2) && !(mb->top4x4inter->direct8x8 & 1) && (0 < mb->top4x4inter->ref[0][lx])) * 2) : 0)
				: -1;
			ref_idx[2] = (blk_map & 2) ?
				(t ? ref_idx_cabac_sub(st, cb,
							((avail & 1) && !(mb->left4x4inter->direct8x8 & 2) && (0 < mb->left4x4inter->ref[1][lx]))
							+ (0 < ref_idx[0]) * 2) : 0)
				: -1;
			blk_map >>= 2;
			ref_idx++;
			num++;
		}
	}
};

struct ref_idx8x16_cabac {
	void operator()(h264d_mb_current *mb, dec_bits *st, int8_t *ref_idx, uint32_t blk_map, int avail) const {
		h264d_cabac_t *cb = mb->cabac;
		int8_t * const *num = mb->num_ref_idx_lx_active_minus1;
		for (int lx = 0; lx < 2; ++lx) {
			int t = *(num[0]);
			ref_idx[0] = (blk_map & 1) ?
				(t ? ref_idx_cabac_sub(st, cb,
							((avail & 1) && !(mb->left4x4inter->direct8x8 & 1) && (0 < mb->left4x4inter->ref[0][lx]))
							+ ((avail & 2) && !(mb->top4x4inter->direct8x8 & 1) && (0 < mb->top4x4inter->ref[0][lx])) * 2) : 0)
				: -1;
			ref_idx[2] = (blk_map & 2) ?
				(t ? ref_idx_cabac_sub(st, cb,
							(0 < ref_idx[0])
							+ ((avail & 2) && !(mb->top4x4inter->direct8x8 & 2) && (0 < mb->top4x4inter->ref[1][lx])) * 2) : 0)
				: -1;
			blk_map >>= 2;
			ref_idx++;
			num++;
		}
	}
};

static inline int valid_block(const prev8x8_t *pblk, int blk_idx, int lx, int non_direct) {
	return (0 <= non_direct) && (0 < pblk[blk_idx].ref[lx]);
}

struct ref_idx8x8_cabac {
	void operator()(h264d_mb_current *mb, dec_bits *st, const int8_t *sub_mb_type, prev8x8_t *pblk, int avail, int lx) const {
		int t = (mb->type != MB_P8x8REF0) ? *(mb->num_ref_idx_lx_active_minus1[lx]) : 0;
		int dir = 1 << lx;
		h264d_cabac_t *cb = mb->cabac;
		const int8_t *sub_mb_ref_map = mb->sub_mb_ref_map;
		int sub_dir0, sub_dir1, sub_dir2, sub_dir3;

		sub_dir0 = sub_mb_ref_map[*sub_mb_type++];
		if ((0 <= sub_dir0) && (dir & sub_dir0)) {
			pblk[0].ref[lx] = t ? ref_idx_cabac_sub(st, cb, ((avail & 1) && !(mb->left4x4inter->direct8x8 & 1) && (0 < mb->left4x4inter->ref[0][lx])) + ((avail & 2) && !(mb->top4x4inter->direct8x8 & 1) && (0 < mb->top4x4inter->ref[0][lx])) * 2) : 0;
		}
		sub_dir1 = sub_mb_ref_map[*sub_mb_type++];
		if ((0 <= sub_dir1) && (dir & sub_dir1)) {
			pblk[1].ref[lx] = t ? ref_idx_cabac_sub(st, cb, valid_block(pblk, 0, lx, sub_dir0) + ((avail & 2) && !(mb->top4x4inter->direct8x8 & 2) && (0 < mb->top4x4inter->ref[1][lx])) * 2) : 0;
		}
		sub_dir2 = sub_mb_ref_map[*sub_mb_type++];
		if ((0 <= sub_dir2) && (dir & sub_dir2)) {
			pblk[2].ref[lx] = t ? ref_idx_cabac_sub(st, cb, ((avail & 1) && !(mb->left4x4inter->direct8x8 & 2) && (0 < mb->left4x4inter->ref[1][lx])) + valid_block(pblk, 0, lx, sub_dir0) * 2) : 0;
		}
		sub_dir3 = sub_mb_ref_map[*sub_mb_type];
		if ((0 <= sub_dir3) && (dir & sub_dir3)) {
			pblk[3].ref[lx] = t ? ref_idx_cabac_sub(st, cb, valid_block(pblk, 2, lx, sub_dir2) + valid_block(pblk, 1, lx, sub_dir1) * 2) : 0;
		}
	}
};

static int mb_inter16x16_cabac(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_inter16x16(mb, mbc, st, avail, ref_idx16x16_cabac(), mvd_xy_cabac(), cbp_cabac(), qp_delta_cabac(), residual_block_cabac());
}

static int mb_inter16x8_cabac(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_inter16x8(mb, mbc, st, avail, ref_idx16x8_cabac(), mvd_xy_cabac(), cbp_cabac(), qp_delta_cabac(), residual_block_cabac());
}

static int mb_inter8x16_cabac(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_inter8x16(mb, mbc, st, avail, ref_idx8x16_cabac(), mvd_xy_cabac(), cbp_cabac(), qp_delta_cabac(), residual_block_cabac());
}

static int mb_inter8x8p_cabac(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_inter8x8(mb, mbc, st, avail, sub_mb_type_p_cabac(), ref_idx8x8_cabac(), sub_mbs_p_cabac(), sub_mbs_dec_p(), cbp_cabac(), qp_delta_cabac(), store_direct8x8_info_p(), residual_block_cabac());
}

static int mb_inter8x8b_cabac(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_inter8x8(mb, mbc, st, avail, sub_mb_types_b_cabac(), ref_idx8x8_cabac(), sub_mbs_b_cabac(), sub_mbs_dec_b(), cbp_cabac(), qp_delta_cabac(), store_direct8x8_info_b(), residual_block_cabac());
}

static int mb_bdirect16x16_cabac(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_bdirect16x16(mb, mbc, st, avail, cbp_cabac(), qp_delta_cabac(), residual_block_cabac());
}


static const mb_code mb_decode_cabac[] = {
	{mb_intra4x4_cabac, 0, 0},
	{mb_intra16x16_dconly_cabac, mb_intra16xpred_vert<16>, 0},
	{mb_intra16x16_dconly_cabac, mb_intra16x16pred_horiz, 0},
	{mb_intra16x16_dconly_cabac, mb_intra16x16pred_dc, 0},
	{mb_intra16x16_dconly_cabac, mb_intra16x16pred_planer, 0},
	{mb_intra16x16_dconly_cabac, mb_intra16xpred_vert<16>, 0x10},
	{mb_intra16x16_dconly_cabac, mb_intra16x16pred_horiz, 0x10},
	{mb_intra16x16_dconly_cabac, mb_intra16x16pred_dc, 0x10},
	{mb_intra16x16_dconly_cabac, mb_intra16x16pred_planer, 0x10},
	{mb_intra16x16_dconly_cabac, mb_intra16xpred_vert<16>, 0x20},
	{mb_intra16x16_dconly_cabac, mb_intra16x16pred_horiz, 0x20},
	{mb_intra16x16_dconly_cabac, mb_intra16x16pred_dc, 0x20},
	{mb_intra16x16_dconly_cabac, mb_intra16x16pred_planer, 0x20},
	{mb_intra16x16_acdc_cabac, mb_intra16xpred_vert<16>, 0x0f},
	{mb_intra16x16_acdc_cabac, mb_intra16x16pred_horiz, 0x0f},
	{mb_intra16x16_acdc_cabac, mb_intra16x16pred_dc, 0x0f},
	{mb_intra16x16_acdc_cabac, mb_intra16x16pred_planer, 0x0f},
	{mb_intra16x16_acdc_cabac, mb_intra16xpred_vert<16>, 0x1f},
	{mb_intra16x16_acdc_cabac, mb_intra16x16pred_horiz, 0x1f},
	{mb_intra16x16_acdc_cabac, mb_intra16x16pred_dc, 0x1f},
	{mb_intra16x16_acdc_cabac, mb_intra16x16pred_planer, 0x1f},
	{mb_intra16x16_acdc_cabac, mb_intra16xpred_vert<16>, 0x2f},
	{mb_intra16x16_acdc_cabac, mb_intra16x16pred_horiz, 0x2f},
	{mb_intra16x16_acdc_cabac, mb_intra16x16pred_dc, 0x2f},
	{mb_intra16x16_acdc_cabac, mb_intra16x16pred_planer, 0x2f},
	{mb_intrapcm, 0, 0},
	{mb_inter16x16_cabac, 0, 1},
	{mb_inter16x8_cabac, 0, 3},
	{mb_inter8x16_cabac, 0, 3},
	{mb_inter8x8p_cabac, 0, 0xf},
	{mb_inter8x8p_cabac, 0, 0xf},
	{mb_bdirect16x16_cabac, 0, 0},
	{mb_inter16x16_cabac, 0, 1}, {mb_inter16x16_cabac, 0, 2}, {mb_inter16x16_cabac, 0, 3},
	{mb_inter16x8_cabac, 0, 0x3}, {mb_inter8x16_cabac, 0, 0x3},
	{mb_inter16x8_cabac, 0, 0xc}, {mb_inter8x16_cabac, 0, 0xc},
	{mb_inter16x8_cabac, 0, 0x9}, {mb_inter8x16_cabac, 0, 0x9},
	{mb_inter16x8_cabac, 0, 0x6}, {mb_inter8x16_cabac, 0, 0x6},
	{mb_inter16x8_cabac, 0, 0xb}, {mb_inter8x16_cabac, 0, 0xb},
	{mb_inter16x8_cabac, 0, 0xe}, {mb_inter8x16_cabac, 0, 0xe},
	{mb_inter16x8_cabac, 0, 0x7}, {mb_inter8x16_cabac, 0, 0x7},
	{mb_inter16x8_cabac, 0, 0xd}, {mb_inter8x16_cabac, 0, 0xd},
	{mb_inter16x8_cabac, 0, 0xf}, {mb_inter8x16_cabac, 0, 0xf},
	{mb_inter8x8b_cabac, 0, 0}
};

static inline int macroblock_layer_cabac(h264d_mb_current *mb, int slice_type, dec_bits *st)
{
	static int (* const mb_type_cabac[3])(h264d_mb_current *mb, dec_bits *st, int avail, int ctx_idx, int slice_type) = {
		mb_type_cabac_P,
		mb_type_cabac_B,
		mb_type_cabac_I,
	};
	const mb_code *mbc;
	int mbtype;
	int avail;
	avail = get_availability(mb);
	mb->type = mbtype = adjust_mb_type(mb_type_cabac[slice_type](mb, st, avail, 3, slice_type), slice_type);
	mbc = &mb_decode_cabac[mbtype];
	mbc->mb_dec(mb, mbc, st, avail);
	return 0;
}

static const m2d_func_table_t h264d_func_ = {
	sizeof(h264d_context),
	(int (*)(void *, int, int (*)(void *, int), void *))h264d_init,
	(dec_bits *(*)(void *))h264d_stream_pos,
	(int (*)(void *, m2d_info_t *))h264d_get_info,
	(int (*)(void *, int, m2d_frame_t *, uint8_t *, int))h264d_set_frames,
	(int (*)(void *))h264d_decode_picture,
	(int (*)(void *, m2d_frame_t *, int))h264d_get_decoded_frame
};

const m2d_func_table_t * const h264d_func = &h264d_func_;
