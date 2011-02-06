#include <assert.h>
#include <string.h>
#include <setjmp.h>
#include <stdlib.h>
#include <limits.h>
#include "bitio.h"
#include "mpeg2.h"
#include "h264.h"
#include "h264vld.h"

#define MIN(a, b) ((a) <= (b) ? (a) : (b))
#define ABS(a) ((0 <= (a)) ? (a) : -(a))

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
	return (((ue & 1) ? ue : -ue) + 1) >> 1;
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

static int32_t me_golomb(dec_bits *stream, int is_inter)
{
	uint32_t ue = ue_golomb(stream);
	assert((unsigned)is_inter <= 1);
	return me_golomb_lut[is_inter][(ue < 48) ? ue : 0];
}

static int32_t te_golomb(dec_bits *stream, int range)
{
	return (range == 1) ? (get_onebit_inline(stream) ^ 1) : ue_golomb(stream);
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

static void set_mb_size(h264d_mb_current *mb, int width, int height);

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
		sps->num_ref_frames_in_pic_order_cnt_cycle = ue_golomb(stream);
		for (unsigned i = 0; i < sps->num_ref_frames_in_pic_order_cnt_cycle; ++i) {
			sps->offset_for_ref_frame[i] = se_golomb(stream);
		}
	}
	READ_UE_RANGE(sps->num_ref_frames, stream, 16);
	sps->gaps_in_frame_num_value_allowed_flag = get_onebit(stream);
	sps->pic_width = (ue_golomb(stream) + 1) * 16;
	sps->pic_height = (ue_golomb(stream) + 1) * 16;
	if ((sps->frame_mbs_only_flag = get_onebit(stream)) == 0) {
		sps->mb_adaptive_frame_field_flag = get_onebit(stream);
	}
	sps->direct_8x8_inference_flag = get_onebit(stream);
	if ((sps->frame_cropping_flag = get_onebit(stream)) != 0) {
		for (int i = 0; i < 4; ++i) {
			sps->frame_crop[i] = ue_golomb(stream) * 2;
		}
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

int h264d_init(h264d_context *h2d)
{
	if (!h2d) {
		return -1;
	}
	memset(h2d, 0, sizeof(*h2d));
	h2d->stream = &h2d->stream_i;
	h2d->slice_header = &h2d->slice_header_i;
	h2d->mb_current.frame = &h2d->mb_current.frame_i;
	h2d->mb_current.cabac = &h2d->mb_current.cabac_i;
	h2d->mb_current.num_ref_idx_l0_active_minus1 = &h2d->slice_header->num_ref_idx_l0_active_minus1;
	h2d->mb_current.num_ref_idx_l1_active_minus1 = &h2d->slice_header->num_ref_idx_l1_active_minus1;
	dec_bits_open(h2d->stream, h264d_load_bytes_skip03);
	return 0;
}

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

int h264d_get_info(h264d_context *h2d, h264d_info_t *info)
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
		+ sizeof(deblock_info_t) * ((src_width * info->src_height) >> 8);
	return 0;
}

static void init_mb_buffer(h264d_mb_current *mb, uint8_t *second)
{
	mb->mb_base = (prev_mb_t *)second;
	second += sizeof(*mb->mb_base) * (mb->max_x + 1);
	mb->top4x4pred_base = (int32_t *)second;
	second += sizeof(*mb->top4x4pred_base) * mb->max_x;
	mb->top4x4coef_base = (int32_t *)second;
	second += sizeof(*mb->top4x4coef_base) * mb->max_x;

	mb->deblock_base = (deblock_info_t *)second;
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
	mb->deblock_curr++;
	mb->deblock_curr->idc = 0;
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
	if (0 <= mb->firstline) {
		mb->firstline--;
	}
	return 0;
}

static inline void frames_init(h264d_mb_current *mb, int num_frame, h264d_frame *frame)
{
	h264d_frame_info_t *frm = mb->frame;
	frm->num = num_frame;
	for (int i = 0; i < num_frame; ++i) {
		frm->frames[i] = frame[i];
	}
	memset(frm->lru, 0, sizeof(frm->lru));
}

int h264d_set_frames(h264d_context *h2d, int num_frame, h264d_frame *frame, uint8_t *second_frame)
{
	h264d_mb_current *mb;

	if (!h2d || (num_frame < 3) || !frame || !second_frame) {
		return -1;
	}
	mb = &h2d->mb_current;
	frames_init(mb, num_frame, frame);
	mb->frame->refs = h2d->slice_header->reorder->ref_frames;
	init_mb_buffer(mb, second_frame);
	return 0;
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
	if ((stream->jmp) && setjmp(*stream->jmp) != 0) {
		return -1;
	}
	h2d->slice_header->first_mb_in_slice = UINT_MAX;
	err = 0;
	code_type = 0;
	do {
		if (0 <= (err = m2d_find_mpeg_data(stream))) {
			code_type = get_bits(stream, 8);
			err = h2d_dispatch_one_nal(h2d, code_type);
		} else {
			break;
		}
		VC_CHECK;
	} while (err == 0 || (code_type == SPS_NAL && 0 < err));
#ifdef DUMP_COEF
	print_coefs();
#endif
	return err;
}

int h246d_get_decoded_frame(h264d_context *h2d, uint8_t **luma, uint8_t **chroma)
{
	if (!h2d || !luma || !chroma) {
		return -1;
	}
	*luma = h2d->mb_current.frame->curr_luma;
	*chroma = h2d->mb_current.frame->curr_chroma;
	return 0;
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
static int pred_weight_table(h264d_slice_header *hdr, dec_bits *st);
static int dec_ref_pic_marking(int nal_id, h264d_marking_t *mrk, dec_bits *st);

static int slice_type_adjust(int slice_type)
{
	return (SI_SLICE < slice_type) ? slice_type - SI_SLICE - 1 : slice_type;
}

static inline uint32_t frames_in_use(h264d_ref_frame_t *refs)
{
	int i = 16;
	uint32_t map = 0;
	do {
		if (refs->in_use) {
			map |= 1 << refs->frame_idx;
		}
		refs++;
	} while (--i);
	return map;
}

static inline void find_empty_frame(h264d_mb_current *mb)
{
	h264d_frame_info_t *frm = mb->frame;
	int max_idx = 0;
	int max_val = -1;
	uint32_t in_use = frames_in_use(frm->refs);
	uint32_t mask = 1;
	int *lru = frm->lru;
	for (int i = 0; i < frm->num; ++i) {
		if (!(in_use & mask)) {
			int val = lru[i];
			if (max_val < val) {
				max_val = val;
				max_idx = i;
			}
			lru[i] = val + 1;
		}
		mask <<= 1;
	}
	lru[max_idx] = 0;
	frm->index = max_idx;
	frm->curr_luma = frm->frames[max_idx].luma;
	frm->curr_chroma = frm->frames[max_idx].chroma;
}

static const int8_t normAdjust[6][3] = {
	{10, 16, 13},
	{11, 18, 14},
	{13, 20, 16},
	{14, 23, 18},
	{16, 25, 20},
	{18, 29, 23}
};

static const int8_t qpc_adjust[22] = {
	29, 30, 31, 32, 32, 33, 34, 34,
	35, 35, 36, 36, 37, 37, 37, 38,
	38, 38, 39, 39, 39, 39
};

static void qp_matrix(int16_t *matrix, int scale, int shift)
{
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



static int slice_header(h264d_context *h2d, dec_bits *st)
{
	h264d_slice_header *hdr = h2d->slice_header;
	h264d_sps *sps;
	h264d_pps *pps;
	uint32_t prev_first_mb = hdr->first_mb_in_slice;
	int slice_type;

	if ((hdr->first_mb_in_slice = ue_golomb(st)) <= prev_first_mb) {
		if (prev_first_mb != UINT_MAX) {
			return -2;
		}
		find_empty_frame(&h2d->mb_current);
		memset(h2d->mb_current.deblock_base, 0, sizeof(deblock_info_t) * h2d->mb_current.max_x * h2d->mb_current.max_y);
	}
	READ_UE_RANGE(slice_type, st, 9);
	hdr->slice_type = slice_type_adjust(slice_type);
	READ_UE_RANGE(hdr->pic_parameter_set_id, st, 255);
	pps = &h2d->pps_i[hdr->pic_parameter_set_id];
	sps = &h2d->sps_i[pps->seq_parameter_set_id];
	h2d->mb_current.pps = pps;
	h2d->mb_current.is_constrained_intra = pps->constrained_intra_pred_flag;
	hdr->frame_num = get_bits(st, sps->log2_max_frame_num);
	if (!sps->frame_mbs_only_flag) {
		if ((hdr->field_pic_flag = get_onebit(st)) != 0) {
			hdr->bottom_field_flag = get_onebit(st);
		}
	}
	if ((h2d->id & 31) == 5) {
		READ_UE_RANGE(hdr->idr_pic_id, st, 65535);
	}
	set_mb_size(&h2d->mb_current, sps->pic_width, sps->pic_height);
	set_mb_pos(&h2d->mb_current, hdr->first_mb_in_slice);
	if (sps->poc_type == 0) {
		hdr->poc0.lsb = get_bits(st, sps->log2_max_poc_lsb);
		if (!hdr->field_pic_flag && pps->pic_order_present_flag) {
			hdr->poc0.delta_pic_order_cnt_bottom = se_golomb(st);
		}
	} else if ((sps->poc_type == 1) && !sps->delta_pic_order_always_zero_flag) {
		hdr->poc1.delta_pic_order_cnt[0] = se_golomb(st);
		if (!hdr->field_pic_flag && pps->pic_order_present_flag) {
			hdr->poc1.delta_pic_order_cnt[1] = se_golomb(st);
		}
	}
	if (pps->redundant_pic_cnt_present_flag) {
		hdr->redundant_pic_cnt = ue_golomb(st);
	}
	switch (hdr->slice_type)
	{
	case B_SLICE:
		hdr->direct_spatial_mv_pred_flag = get_onebit(st);
		/* FALLTHROUGH */
	case P_SLICE:
	case SP_SLICE:
		if ((hdr->num_ref_idx_active_override_flag = get_onebit(st)) != 0) {
			READ_UE_RANGE(hdr->num_ref_idx_l0_active_minus1, st, 31);
			if (hdr->slice_type == B_SLICE) {
				READ_UE_RANGE(hdr->num_ref_idx_l1_active_minus1, st, 31);
			}
		} else {
			hdr->num_ref_idx_l0_active_minus1 = pps->num_ref_idx_l0_active_minus1;
			hdr->num_ref_idx_l1_active_minus1 = pps->num_ref_idx_l1_active_minus1;
		}
		if (ref_pic_list_reordering(&hdr->reorder[0], st, sps->num_ref_frames, hdr->frame_num, 1 << sps->log2_max_frame_num)) {
			return -1;
		}
		if (hdr->slice_type == B_SLICE) {
			if (ref_pic_list_reordering(&hdr->reorder[1], st, sps->num_ref_frames, hdr->frame_num, 1 << sps->log2_max_frame_num)) {
				return -1;
			}
			if (pps->weighted_bipred_idc == 1) {
				pred_weight_table(hdr, st);
			}
		} else {
			if (pps->weighted_pred_flag) {
				pred_weight_table(hdr, st);
			}
		}
	}
	if (h2d->id & 0x60) {
		if (dec_ref_pic_marking(h2d->id & 31, &hdr->marking, st) != 0) {
			hdr->frame_num = 0;
		}
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

static void dump_ref_list(h264d_ref_frame_t *refs, int num_ref_frames)
{
	printf("length: %d\n", num_ref_frames);
	if (!num_ref_frames) {
		return;
	}
	do {
		printf("\t%d %d\n", refs->num, refs->in_use);
		refs++;
	} while (--num_ref_frames);
}

static int ref_pic_list_reordering(h264d_reorder_t *rdr, dec_bits *st, int num_ref_frames, int frame_num, int max_frame_num)
{
	assert((unsigned)num_ref_frames <= 16);
	if ((rdr->ref_pic_list_reordering_flag = get_onebit(st)) != 0) {
		h264d_ref_frame_t *refs = rdr->ref_frames;
		int refIdxLx = -1;
		while (++refIdxLx < 16) {
			h264d_ref_frame_t tmp_ref;
			int idc;
			uint32_t num;
			int curr_idx;
			int mode;

			READ_UE_RANGE(idc, st, 3);
			if (idc == 3) {
				break;
			} else if (3 < idc) {
				return -1;
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
			if (!refs[num_ref_frames].in_use) {
				while (!refs[num_ref_frames].in_use && 0 <= num_ref_frames) {
					--num_ref_frames;
				}
				++num_ref_frames;
			}
			curr_idx = refIdxLx;
			while (!(refs[curr_idx].num == num && refs[curr_idx].in_use == mode) && curr_idx < num_ref_frames) {
				curr_idx++;
			}
			if (num_ref_frames <= curr_idx) {
				continue;
			}
			tmp_ref = refs[curr_idx];
			while (refIdxLx < curr_idx) {
				refs[curr_idx] = refs[curr_idx - 1];
				--curr_idx;
			}
			refs[curr_idx] = tmp_ref;
		}
//		dump_ref_list(refs, num_ref_frames);
	}
	return 0;
}

static int pred_weight_table(h264d_slice_header *hdr, dec_bits *st)
{
	return 0;
}

static int dec_ref_pic_marking(int nal_unit_type, h264d_marking_t *mrk, dec_bits *st)
{
	uint32_t tmp = get_onebit(st);
	int op5_detect = 0;

	if (nal_unit_type == 5) {
		mrk->idr.no_output_of_prior_pic_flag = tmp;
		mrk->idr.long_term_reference_flag = get_onebit(st);
	} else {
		mrk->non_idr.adaptive_ref_pic_marking_mode_flag = tmp;
		if (tmp) {
			h264d_mmco *mmco = mrk->non_idr.mmco;
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
	return op5_detect;
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

int m2d_dec_vld_unary(dec_bits *stream, const vlc_t *vld_tab, int bitlen);

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
	int operator()(h264d_mb_current *mb, int na, int nb, dec_bits *st, int *coeff, int num_coeff, const int16_t *qmat, int avail, int pos4x4, int cat, uint32_t dc_mask)
	{
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
	int operator()(int a, int b, dec_bits *st, h264d_cabac_t *cb) {
		int pred = MIN(a, b);
		if (!get_onebit_inline(st)) {
			int rem = get_bits(st, 3);
			pred = (rem < pred) ? rem : rem + 1;
		}
		return pred;
	}
};

int intra4x4pred_dc(uint8_t *dst, int stride, int avail)
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
static int mb_pred_intra4x4(h264d_mb_current *mb, dec_bits *st, int avail, int8_t *pred4x4, h264d_cabac_t *cb, F I4x4PredMode) {
	uint32_t left = mb->left4x4pred;
	uint32_t top = *mb->top4x4pred;
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
	mb->lefttop_ref = mb->top4x4inter->ref[1];
	mb->lefttop_mv[0] = mb->top4x4inter->mv[0][3][0];
	mb->lefttop_mv[1] = mb->top4x4inter->mv[0][3][1];
	memset(mb->left4x4inter->mv, 0, sizeof(mb->left4x4inter->mv));
	memset(mb->top4x4inter->mv, 0, sizeof(mb->top4x4inter->mv));
	mb->left4x4inter->ref[0] = -1;
	mb->left4x4inter->ref[1] = -1;
	mb->top4x4inter->ref[0] = -1;
	mb->top4x4inter->ref[1] = -1;
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

template <typename F0, typename F1, typename F2, typename F3, typename F4>
static inline int mb_intra4x4(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail,
				F0 Intra4x4PredMode,
				F1 IntraChromaPredMode,
				F2 CodedBlockPattern,
				F3 QpDelta,
				F4 ResidualBlock)
{
	int8_t pred4x4[16];
	const int8_t *pr;
	uint32_t intra_chroma_pred_mode;
	int stride;
	uint32_t cbp;
	uint8_t *luma;
	const int *offset;
	int avail_intra;

	avail_intra = avail;
	if (mb->is_constrained_intra) {
		avail_intra &= ~((MB_IPCM < mb->top4x4inter[1].type) * 4 | ((MB_IPCM < mb->top4x4inter->type) * 2) | (MB_IPCM < mb->left4x4inter->type));
	}
	fill_dc_if_unavailable(mb, avail_intra);
	mb_pred_intra4x4(mb, st, avail_intra, pred4x4, mb->cabac, Intra4x4PredMode);
	VC_CHECK;
	intra_chroma_pred_mode = IntraChromaPredMode(mb, st, mb->cabac, avail_intra);
	stride = mb->max_x * 16;
	intra_chroma_pred[intra_chroma_pred_mode](mb->chroma, stride, avail_intra);
	cbp = CodedBlockPattern(mb, st, avail);
	luma = mb->luma;
	offset = mb->offset4x4;
	pr = pred4x4;
	if (cbp) {
		int32_t qp_delta = QpDelta(mb, st, mb->cabac, avail);
		if (qp_delta) {
			set_qp(mb, mb->qp + qp_delta);
		}
	}
	if (cbp & 15) {
		int coeff[16];
		uint32_t top, left;
		int c0, c1, c2, c3, c4, c5;
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
	} else {
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
		mb->prev_qp_delta = 0;
		mb->left4x4coef &= 0xffff0000;
		*mb->top4x4coef &= 0xffff0000;
	}
	store_strength_intra(mb);
	mb_intra_save_info(mb);
	mb->cbp = cbp;
	VC_CHECK;
	return residual_chroma(mb, cbp, st, avail, ResidualBlock);
}

struct intra_chroma_pred_mode_cavlc {
	uint32_t operator()(h264d_mb_current *mb, dec_bits *st, h264d_cabac_t *cb, int avail) {
		uint32_t pred_mode = ue_golomb(st);
		pred_mode = pred_mode <= 3 ? pred_mode : 0;
		mb->chroma_pred_mode = pred_mode;
		return pred_mode;
	}
};

struct cbp_intra_cavlc {
	uint32_t operator()(h264d_mb_current *mb, dec_bits *st, int avail) {
		return me_golomb(st, 0);
	}
};

struct cbp_inter_cavlc {
	uint32_t operator()(h264d_mb_current *mb, dec_bits *st, int avail) {
		return me_golomb(st, 1);
	}
};

struct qp_delta_cavlc {
	int operator()(h264d_mb_current *mb, dec_bits *st, h264d_cabac_t *cb, int avail) {
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
	uint32_t operator()(uint32_t x, uint32_t y) {
		uint32_t msk;
		msk = ((x & y) + (((x ^ y) >> 1) & 0x7f7f7f7f)) & ~0x7f7f7f7f;
		msk = (msk << 1) - (msk >> 7);
		return ((x + y) - msk) | msk;
	}
};

struct SubSaturate {
	uint32_t operator()(uint32_t x, uint32_t y) {
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
	intra_chroma_pred_mode = IntraChromaPredMode(mb, st, mb->cabac, avail_intra);
	intra_chroma_pred[intra_chroma_pred_mode](mb->chroma, stride, avail_intra);
	qp_delta = QpDelta(mb, st, mb->cabac, avail);
	if (qp_delta) {
		set_qp(mb, mb->qp + qp_delta);
	}
	mb->left4x4coef &= 0xffff0000;
	*mb->top4x4coef &= 0xffff0000;
	mb->left4x4pred = 0x22222222;
	*mb->top4x4pred = 0x22222222;
	if (ResidualBlock(mb, avail & 1 ? UNPACK(mb->left4x4coef, 0) : -1, avail & 2 ? UNPACK(*mb->top4x4coef, 0) : -1, st, coeff, 16, mb->qmaty, avail_intra, 26, 0, 0)) {
		intra16x16_dc_transform(coeff, dc);
		offset = mb->offset4x4;
		for (int i = 0; i < 16; ++i) {
			ac4x4transform_dconly(luma + *offset++, dc[i], stride);
		}
	}
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
	intra_chroma_pred_mode = IntraChromaPredMode(mb, st, mb->cabac, avail_intra);
	intra_chroma_pred[intra_chroma_pred_mode](mb->chroma, stride, avail_intra);
	qp_delta = QpDelta(mb, st, mb->cabac, avail);
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


static inline void intrapcm_luma(uint8_t *dst, int stride, dec_bits *st)
{
	int y = 16;
	do {
#ifdef WORDS_BIGENDIAN
		*(uint32_t *)dst = get_bits(st, 32);
		*(uint32_t *)(dst + 4) = get_bits(st, 32);
		*(uint32_t *)(dst + 8) = get_bits(st, 32);
		*(uint32_t *)(dst + 12) = get_bits(st, 32);
#else
		dst[0] = get_bits(st, 8);
		dst[1] = get_bits(st, 8);
		dst[2] = get_bits(st, 8);
		dst[3] = get_bits(st, 8);
		dst[4] = get_bits(st, 8);
		dst[5] = get_bits(st, 8);
		dst[6] = get_bits(st, 8);
		dst[7] = get_bits(st, 8);
		dst[8] = get_bits(st, 8);
		dst[9] = get_bits(st, 8);
		dst[10] = get_bits(st, 8);
		dst[11] = get_bits(st, 8);
		dst[12] = get_bits(st, 8);
		dst[13] = get_bits(st, 8);
		dst[14] = get_bits(st, 8);
		dst[15] = get_bits(st, 8);
#endif
		dst += stride;
	} while (--y);
}

static inline void intrapcm_chroma(uint8_t *dst, int stride, dec_bits *st)
{
	int y = 8;
	do {
		dst[0] = get_bits(st, 8);
		dst[2] = get_bits(st, 8);
		dst[4] = get_bits(st, 8);
		dst[6] = get_bits(st, 8);
		dst[8] = get_bits(st, 8);
		dst[10] = get_bits(st, 8);
		dst[12] = get_bits(st, 8);
		dst[14] = get_bits(st, 8);
		dst += stride;
	} while (--y);
}

static int mb_intrapcm(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	int stride;

	stride = mb->max_x * 16;
	byte_align(st);
	intrapcm_luma(mb->luma, stride, st);
	intrapcm_chroma(mb->chroma, stride, st);
	intrapcm_chroma(mb->chroma + 1, stride, st);
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


static inline void copy_inter(const uint8_t *src, uint8_t *dst, int width, int height, int src_stride, int stride)
{
	width = (unsigned)width >> 2;
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

static inline int inter_pred_mvoffset_luma(int mvint_x, int mvint_y, int stride)
{
	return mvint_y * stride + mvint_x;
}

static inline int inter_pred_mvoffset_chroma(int mvint_x, int mvint_y, int stride)
{
	return (mvint_y >> 1) * stride + (mvint_x >> 1) * 2;
}

typedef struct {
	const uint8_t *src_luma;
	const uint8_t *src_chroma;
	uint8_t *dst_luma;
	uint8_t *dst_chroma;
	int16_t pos_x, pos_y;
} mb_pred_t;

static inline void filter_chroma_horiz(const uint8_t *src, uint8_t *dst, int width, int height, int frac, int stride)
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
		src += stride;
		dst += stride;
	} while (--height);
}

static inline void filter_chroma_vert(const uint8_t *src, uint8_t *dst, int width, int height, int frac, int stride)
{
	int c0 = (8 - frac) * 8;
	int c1 = frac * 8;

	do {
		const uint8_t *s = src++;
		uint8_t *d = dst++;
		int t0 = *s;
		s += stride;
		int y = height;
		do {
			int t1 = *s;
			*d = (t0 * c0 + t1 * c1 + 32) >> 6;
			t0 = t1;
			s += stride;
			d += stride;
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

static inline void chroma_inter_umv(const uint8_t *src, uint8_t *dst, int posx, int posy, int width, int height, int stride, int vert_size, chroma_filter_info_t *filter)
{
	uint8_t buf[18 * 9];

	if (posx < -width) {
		src += -width - posx;
		posx = -width;
	} else if (stride - 2 < posx) {
		src -= posx - stride + 2;
		posx = stride - 2;
	}
	if (posy < -height) {
		src -= (height + posy) * stride;
		posy = -height;
	} else if (vert_size - 1 < posy) {
		src -= (posy - vert_size + 1) * stride;
		posy = vert_size - 1;
	}

	if (filter) {
		width += 2;
		height += 1;
	}
	fill_rect_umv_chroma(src, buf, width, height, stride, vert_size, posx, posy);
	if (filter) {
		width -= 2;
		height -= 1;
		filter_chroma_vert_horiz(buf, dst, width, height, filter->mvx, filter->mvy, width + 2, stride);
	} else {
		copy_inter(buf, dst, width, height, width, stride);
	}
}


static void inter_pred_chroma(const mb_pred_t *pred, int mvx, int mvy, int width, int height, int stride, int vert_stride)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	mvx &= 7;
	mvy &= 7;
	const uint8_t *src_chroma = pred->src_chroma + inter_pred_mvoffset_luma(posx, posy, stride);
	uint8_t *dst = pred->dst_chroma;
	if (mvx || mvy) {
		if ((unsigned)posx <= (unsigned)(stride - width - 2) && (unsigned)posy <= (unsigned)(vert_stride - height - 1)) {
			if (mvy) {
				if (mvx) {
					filter_chroma_vert_horiz(src_chroma, dst, width, height, mvx, mvy, stride, stride);
				} else {
					filter_chroma_vert(src_chroma, dst, width, height, mvy, stride);
				}
			} else {
				filter_chroma_horiz(src_chroma, dst, width, height, mvx, stride);
			}
		} else {
			/* UMV */
			chroma_filter_info_t f = {
				mvx, mvy
			};
			chroma_inter_umv(src_chroma, pred->dst_chroma, posx, posy, width, height, stride, vert_stride, &f);
		}

	} else {
		if ((unsigned)posx <= (unsigned)(stride - width) && (unsigned)posy <= (unsigned)(vert_stride - height)) {
			copy_inter(src_chroma, pred->dst_chroma, width, height, stride, stride);
		} else {
			chroma_inter_umv(src_chroma, pred->dst_chroma, posx, posy, width, height, stride, vert_stride, 0);
		}
	}
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
	int operator()(int c0, int c1, int c2, int c3, int c4, int c5) {
		int t = (((c2 + c3) * 4 - c1 - c4) * 5 + c0 + c5 + 512) >> 10;
		int c = (c2 + 16) >> 5;
		return (CLIP255I(t) + CLIP255I(c) + 1) >> 1;
	}
};

struct PPred22 {
	int operator()(int c0, int c1, int c2, int c3, int c4, int c5) {
		int t = (((c2 + c3) * 4 - c1 - c4) * 5 + c0 + c5 + 512) >> 10;
		return CLIP255C(t);
	}
};

struct PPred32 {
	int operator()(int c0, int c1, int c2, int c3, int c4, int c5) {
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

static void inter_pred_luma_umv(const mb_pred_t *pred, int width, int height, int stride, int vert_size, int fracx, int fracy)
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
	} else if (stride - 1 < posx - 2) {
		src -= posx - stride - 1;
		posx = stride + 1;
	}
	if (posy < 3 - height) {
		src += (3 - height - posy) * stride;
		posy = 3 - height;
	} else if (vert_size - 1 < posy - 2) {
		src -= (posy - vert_size - 1) * stride;
		posy = vert_size + 1;
	}

	posx = posx - 2;
	posy = posy - 2;
	fill_rect_umv_luma(src, buf, width, height, stride, vert_size, posx, posy);
	width -= 6;
	height -= 6;
	inter_pred_luma_filter[fracy][fracx](buf, pred->dst_luma, width, height, width + 6, stride);
}

static void inter_pred_luma_frac00(const mb_pred_t *pred, int width, int height, int stride, int vert_stride)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if ((unsigned)posx < (unsigned)(stride - width) && (unsigned)posy < (unsigned)(vert_stride - height)) {
		const uint8_t *src_luma = pred->src_luma + inter_pred_mvoffset_luma(2, 2, stride);
		copy_inter(src_luma, pred->dst_luma, width, height, stride, stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, 0, 0);
	}
}

static void inter_pred_luma_frac01(const mb_pred_t *pred, int width, int height, int stride, int vert_stride)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && (unsigned)posy < (unsigned)(vert_stride - height)) {
		inter_pred_luma_filter01(pred->src_luma, pred->dst_luma, width, height, stride, stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, 1, 0);
	}
}

static void inter_pred_luma_frac02(const mb_pred_t *pred, int width, int height, int stride, int vert_stride)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && (unsigned)posy < (unsigned)(vert_stride - height)) {
		inter_pred_luma_filter02(pred->src_luma, pred->dst_luma, width, height, stride, stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, 2, 0);
	}
}

static void inter_pred_luma_frac03(const mb_pred_t *pred, int width, int height, int stride, int vert_stride)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && (unsigned)posy < (unsigned)(vert_stride - height)) {
		inter_pred_luma_filter03(pred->src_luma, pred->dst_luma, width, height, stride, stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, 3, 0);
	}
}

static void inter_pred_luma_frac10(const mb_pred_t *pred, int width, int height, int stride, int vert_stride)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if ((unsigned)posx < (unsigned)(stride - width) && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter10(pred->src_luma, pred->dst_luma, width, height, stride, stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, 0, 1);
	}
}

static void inter_pred_luma_frac11(const mb_pred_t *pred, int width, int height, int stride, int vert_stride)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter11(pred->src_luma, pred->dst_luma, width, height, stride, stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, 1, 1);
	}
}

static void inter_pred_luma_frac12(const mb_pred_t *pred, int width, int height, int stride, int vert_stride)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter12(pred->src_luma, pred->dst_luma, width, height, stride, stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, 2, 1);
	}
}

static void inter_pred_luma_frac13(const mb_pred_t *pred, int width, int height, int stride, int vert_stride)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter13(pred->src_luma, pred->dst_luma, width, height, stride, stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, 3, 1);
	}
}

static void inter_pred_luma_frac20(const mb_pred_t *pred, int width, int height, int stride, int vert_stride)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if ((unsigned)posx < (unsigned)(stride - width) && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter20(pred->src_luma, pred->dst_luma, width, height, stride, stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, 0, 2);
	}
}

static void inter_pred_luma_frac21(const mb_pred_t *pred, int width, int height, int stride, int vert_stride)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter21(pred->src_luma, pred->dst_luma, width, height, stride, stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, 1, 2);
	}
}

static void inter_pred_luma_frac22(const mb_pred_t *pred, int width, int height, int stride, int vert_stride)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter22(pred->src_luma, pred->dst_luma, width, height, stride, stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, 2, 2);
	}
}

static void inter_pred_luma_frac23(const mb_pred_t *pred, int width, int height, int stride, int vert_stride)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter23(pred->src_luma, pred->dst_luma, width, height, stride, stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, 3, 2);
	}
}

static void inter_pred_luma_frac30(const mb_pred_t *pred, int width, int height, int stride, int vert_stride)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if ((unsigned)posx < (unsigned)(stride - width) && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter30(pred->src_luma, pred->dst_luma, width, height, stride, stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, 0, 3);
	}
}

static void inter_pred_luma_frac31(const mb_pred_t *pred, int width, int height, int stride, int vert_stride)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter31(pred->src_luma, pred->dst_luma, width, height, stride, stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, 1, 3);
	}
}

static void inter_pred_luma_frac32(const mb_pred_t *pred, int width, int height, int stride, int vert_stride)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter32(pred->src_luma, pred->dst_luma, width, height, stride, stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, 2, 3);
	}
}

static void inter_pred_luma_frac33(const mb_pred_t *pred, int width, int height, int stride, int vert_stride)
{
	int posx = pred->pos_x;
	int posy = pred->pos_y;
	if (2 <= posx && posx < stride - width - 2 && 2 <= posy && posy < vert_stride - height - 2) {
		inter_pred_luma_filter33(pred->src_luma, pred->dst_luma, width, height, stride, stride);
	} else {
		inter_pred_luma_umv(pred, width, height, stride, vert_stride, 3, 3);
	}
}

static void (* const inter_pred_luma[4][4])(const mb_pred_t *pred, int width, int height, int stride, int vert_stride) = {
	{
		inter_pred_luma_frac00,
		inter_pred_luma_frac01,
		inter_pred_luma_frac02,
		inter_pred_luma_frac03,
	},
	{
		inter_pred_luma_frac10,
		inter_pred_luma_frac11,
		inter_pred_luma_frac12,
		inter_pred_luma_frac13,
	},
	{
		inter_pred_luma_frac20,
		inter_pred_luma_frac21,
		inter_pred_luma_frac22,
		inter_pred_luma_frac23,
	},
	{
		inter_pred_luma_frac30,
		inter_pred_luma_frac31,
		inter_pred_luma_frac32,
		inter_pred_luma_frac33,
	},
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

	qp_delta = QpDelta(mb, st, mb->cabac, avail);
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

static const int8_t not_one_hot[8] = {
	1, 0, 0, 1, 0, 1, 1, 1
};

static const int16_t zero_mv[2] = {
	0, 0
};

static inline void determine_pmv(const int16_t *mva, const int16_t *mvb, const int16_t *mvc, int16_t pmv[], int avail, int idx_map)
{
	int pmvx, pmvy;
	if (((avail & 7) == 1) || (idx_map == 1)) {
		pmvx = mva[0];
		pmvy = mva[1];
	} else if (not_one_hot[idx_map]) {
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


static inline void calc_mv16x16(h264d_mb_current *mb, int16_t *pmv, const int16_t *&mvd_a, const int16_t *&mvd_b, int ref_idx, int avail)
{
	prev_mb_t *pmb;
	const int16_t *mva, *mvb, *mvc;
	int idx_map;

	if (avail & 1) {
		pmb = mb->left4x4inter;
		idx_map = (ref_idx == pmb->ref[0]);
		mva = pmb->mv[0][0];
		mvd_a = pmb->mv[1][0];
	} else {
		idx_map = 0;
		mvd_a = mva = zero_mv;
	}
	if (avail & 2) {
		pmb = mb->top4x4inter;
		idx_map |= (ref_idx == pmb->ref[0]) * 2;
		mvb = pmb->mv[0][0];
		mvd_b = pmb->mv[1][0];
	} else {
		mvd_b = mvb = zero_mv;
	}
	if (avail & 4) {
		pmb = mb->top4x4inter + 1;
		idx_map |= (ref_idx == pmb->ref[0]) * 4;
		mvc = pmb->mv[0][0];
	} else if (avail & 8) {
		idx_map |= (ref_idx == mb->lefttop_ref) * 4;
		mvc = mb->lefttop_mv;
	} else {
		mvc = zero_mv;
	}
	determine_pmv(mva, mvb, mvc, pmv, avail, idx_map);
}

static inline void inter_pred8x8(h264d_mb_current *mb, const int16_t mv[], int width, int height, int ref_idx, int offsetx, int offsety);
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
	for (int i = 0; i < 4; ++i) {
		if ((prev4x4 & 0xf) != 0) {
			map |= 2 << (i * 2);
		}
		prev4x4 >>= 4;
	}
	return map;
}

static inline int DIF_SQUARE(int a, int b) {
	int t = a - b;
	return t * t;
}

static inline uint32_t str_mv_calc16x16(uint32_t str, const int16_t mvxy[], int ref_idx, const prev_mb_t *top)
{
	int mvx = mvxy[0];
	int mvy = mvxy[1];
	for (int i = 0; i < 4; ++i) {
		if (!(str & (3 << (i * 2)))) {
			if ((ref_idx != top->ref[i >> 1]) || (16 <= DIF_SQUARE(mvx, top->mv[0][i][0])) || (16 <= DIF_SQUARE(mvy, top->mv[0][i][1]))) {
				str |= 1 << (i * 2);
			}
		}

	}
	return str;
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
	int ref_idx;
	int16_t mv[2];
	int16_t mvd[2];
	int t;
	uint32_t cbp;
	uint32_t str_vert, str_horiz;

	t = *mb->num_ref_idx_l0_active_minus1;
	ref_idx = (0 < t) ? RefIdx16x16(mb, st, t, avail) : 0;
	calc_mv16x16(mb, mv, mvd_a, mvd_b, ref_idx, avail);
	MvdXY(mb, st, mvd, mvd_a, mvd_b);
	mv[0] += mvd[0];
	mv[1] += mvd[1];
	inter_pred8x8(mb, mv, 16, 16, ref_idx, 0, 0);
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

	deblock_info_t *deb = mb->deblock_curr;
	deb->qpy = mb->qp;
	deb->qpc = mb->qp_chroma;
	if (mb->y != 0) {
		if (mb->top4x4inter->type <= MB_IPCM) {
			deb->str4_vert = 1;
			str_vert |= 1;
		}
		str_vert = str_previous_coef(str_vert, top4x4);
		str_vert = str_mv_calc16x16(str_vert, mv, ref_idx, mb->top4x4inter);
	}
	deb->str_vert = str_vert;
	if (mb->x != 0) {
		if (mb->left4x4inter->type <= MB_IPCM) {
			deb->str4_horiz = 1;
			str_horiz |= 1;
		}
		str_horiz = str_previous_coef(str_horiz, left4x4);
		str_horiz = str_mv_calc16x16(str_horiz, mv, ref_idx, mb->left4x4inter);
	}
	deb->str_horiz = str_horiz;

	mb->left4x4pred = 0x22222222;
	*mb->top4x4pred = 0x22222222;

	mb->lefttop_ref = mb->top4x4inter->ref[1];
	mb->lefttop_mv[0] = mb->top4x4inter->mv[0][3][0];

	mb->lefttop_mv[1] = mb->top4x4inter->mv[0][3][1];
	mb->top4x4inter->ref[0] = ref_idx;
	mb->left4x4inter->ref[0] = ref_idx;
	mb->top4x4inter->ref[1] = ref_idx;
	mb->left4x4inter->ref[1] = ref_idx;
	for (int i = 0; i < 4; ++i) {
		mb->top4x4inter->mv[0][i][0] = mv[0];
		mb->top4x4inter->mv[0][i][1] = mv[1];
		mb->left4x4inter->mv[0][i][0] = mv[0];
		mb->left4x4inter->mv[0][i][1] = mv[1];
		mb->top4x4inter->mv[1][i][0] = mvd[0];
		mb->top4x4inter->mv[1][i][1] = mvd[1];
		mb->left4x4inter->mv[1][i][0] = mvd[0];
		mb->left4x4inter->mv[1][i][1] = mvd[1];
	}

	return residual_chroma(mb, cbp, st, avail, ResidualBlock);
}

static inline void calc_mv16x8top(h264d_mb_current *mb, int16_t pmv[], const int16_t *&mvd_a, const int16_t *&mvd_b, int ref_idx, int avail)
{
	prev_mb_t *pmb;
	const int16_t *mva, *mvb, *mvc;
	int idx_map;

	if (avail & 2) {
		pmb = mb->top4x4inter;
		mvd_b = pmb->mv[1][0];
		if (ref_idx == pmb->ref[0]) {
			pmv[0] = pmb->mv[0][0][0];
			pmv[1] = pmb->mv[0][0][1];
			mvd_a = (avail & 1) ? mb->left4x4inter->mv[1][0] : zero_mv;
			return;
		}
		mvb = pmb->mv[0][0];
	} else {
		mvd_b = mvb = zero_mv;
	}
	if (avail & 1) {
		pmb = mb->left4x4inter;
		idx_map = (ref_idx == pmb->ref[0]);
		mva = pmb->mv[0][0];
		mvd_a = pmb->mv[1][0];
	} else {
		mvd_a = mva = zero_mv;
		idx_map = 0;
	}
	if (avail & 4) {
		pmb = mb->top4x4inter + 1;
		idx_map |= (ref_idx == pmb->ref[0]) * 4;
		mvc = pmb->mv[0][0];
	} else if (avail & 8) {
		idx_map |= (ref_idx == mb->lefttop_ref) * 4;
		mvc = mb->lefttop_mv;
	} else {
		mvc = zero_mv;
	}
	determine_pmv(mva, mvb, mvc, pmv, avail, idx_map);
}

static inline void calc_mv16x8bottom(h264d_mb_current *mb, int16_t mv[], const int16_t *&mvd_a, const int16_t *&mvd_b, int ref_idx, int avail, int prev_ref, const int16_t *prev_mv)
{
	prev_mb_t *pmb;
	const int16_t *mva, *mvb, *mvc;
	int idx_map;

	if (avail & 1) {
		pmb = mb->left4x4inter;
		mvd_a = pmb->mv[1][2];
		if (ref_idx == pmb->ref[1]) {
			mv[0] = pmb->mv[0][2][0];
			mv[1] = pmb->mv[0][2][1];
			mvd_b = prev_mv + 2;
			return;
		}
		idx_map = (ref_idx == pmb->ref[0]) * 4;
		mva = pmb->mv[0][2];
		mvc = pmb->mv[0][1];
	} else {
		idx_map = 0;
		mvd_a = mva = zero_mv;
		mvc = zero_mv;
	}
	/* upper block: always exists */
	mvb = prev_mv;
	mvd_b = prev_mv + 2;
	idx_map |= (ref_idx == prev_ref) * 2;
	avail |= 2;
	determine_pmv(mva, mvb, mvc, mv, avail, idx_map);
}

static inline uint32_t str_mv_calc16x8_left(uint32_t str, int ref_idx0, int ref_idx1, const int16_t mvxy0[], const int16_t mvxy1[], const prev_mb_t *top)
{
	int mvx, mvy;
	int i;
	mvx = mvxy0[0];
	mvy = mvxy0[1];
	for (i = 0; i < 2; ++i) {
		if (!(str & (3 << (i * 2))) && ((ref_idx0 != top->ref[0]) || (16 <= DIF_SQUARE(mvx, top->mv[0][i][0])) || (16 <= DIF_SQUARE(mvy, top->mv[0][i][1])))) {
			str |= 1 << (i * 2);
		}
	}
	mvx = mvxy1[0];
	mvy = mvxy1[1];
	for (;i < 4; ++i) {
		int shift = i * 2;
		if (!(str & (3 << shift)) && ((ref_idx1 != top->ref[1]) || (16 <= DIF_SQUARE(mvx, top->mv[0][i][0])) || (16 <= DIF_SQUARE(mvy, top->mv[0][i][1])))) {
			str |= 1 << shift;
		}
	}
	return str;
}

static inline uint32_t str_mv_calc16x8_vert(uint32_t str, int ref_idx0, int ref_idx1, const int16_t mvxy0[], const int16_t mvxy1[])
{
	if ((ref_idx0 != ref_idx1) || (16 <= DIF_SQUARE(mvxy0[0], mvxy1[0])) || (16 <= DIF_SQUARE(mvxy0[1], mvxy1[1]))) {
		for (int i = 0; i < 4; ++i) {
			int shift = i * 2 + 16;
			if (!(str & (3 << shift))) {
				str |= 1 << shift;
			}
		}
	}
	return str;
}

template <typename F0 ,typename F1, typename F2, typename F3, typename F4>
static int mb_inter16x8(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail,
			 F0 RefIdx16x8,
			 F1 MvdXY,
			 F2 CodedBlockPattern,
			 F3 QpDelta,
			 F4 ResidualBlock)
{
	int t;
	uint32_t left4x4, top4x4;
	uint32_t cbp;
	uint32_t str_vert, str_horiz;
	int16_t mv[2][2][2];
	const int16_t *mvd_a, *mvd_b;
	int8_t ref_idx[2];

	t = *mb->num_ref_idx_l0_active_minus1;
	if (0 < t) {
		RefIdx16x8(mb, st, t, ref_idx, avail);
	} else {
		ref_idx[0] = 0;
		ref_idx[1] = 0;
	}
	calc_mv16x8top(mb, mv[0][0], mvd_a, mvd_b, ref_idx[0], avail);
	MvdXY(mb, st, mv[0][1], mvd_a, mvd_b);
	mv[0][0][0] += mv[0][1][0];
	mv[0][0][1] += mv[0][1][1];
	inter_pred8x8(mb, mv[0][0], 16, 8, ref_idx[0], 0, 0);

	calc_mv16x8bottom(mb, mv[1][0], mvd_a, mvd_b, ref_idx[1], avail, ref_idx[0], mv[0][0]);
	MvdXY(mb, st, mv[1][1], mvd_a, mvd_b);
	mv[1][0][0] += mv[1][1][0];
	mv[1][0][1] += mv[1][1][1];
	inter_pred8x8(mb, mv[1][0], 16, 8, ref_idx[1], 0, 8);

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

	deblock_info_t *deb = mb->deblock_curr;
	deb->qpy = mb->qp;
	deb->qpc = mb->qp_chroma;
	if (mb->y != 0) {
		if (mb->top4x4inter->type <= MB_IPCM) {
			deb->str4_vert = 1;
			str_vert |= 1;
		}
		str_vert = str_previous_coef(str_vert, top4x4);
		str_vert = str_mv_calc16x16(str_vert, mv[0][0], ref_idx[0], mb->top4x4inter);
	}
	deb->str_vert = str_mv_calc16x8_vert(str_vert, ref_idx[0], ref_idx[1], mv[0][0], mv[1][0]);
	if (mb->x != 0) {
		if (mb->left4x4inter->type <= MB_IPCM) {
			deb->str4_horiz = 1;
			str_horiz |= 1;
		}
		str_horiz = str_previous_coef(str_horiz, left4x4);
		str_horiz = str_mv_calc16x8_left(str_horiz, ref_idx[0], ref_idx[1], mv[0][0], mv[1][0], mb->left4x4inter);
	}
	deb->str_horiz = str_horiz;

	mb->left4x4pred = 0x22222222;
	*mb->top4x4pred = 0x22222222;
	mb->lefttop_ref = mb->top4x4inter->ref[1];
	mb->lefttop_mv[0] = mb->top4x4inter->mv[0][3][0];
	mb->lefttop_mv[1] = mb->top4x4inter->mv[0][3][1];
	for (int i = 0; i < 4; ++i) {
		mb->top4x4inter->mv[0][i][0] = mv[1][0][0];
		mb->top4x4inter->mv[0][i][1] = mv[1][0][1];
		mb->top4x4inter->mv[1][i][0] = mv[1][1][0];
		mb->top4x4inter->mv[1][i][1] = mv[1][1][1];
	}
	mb->top4x4inter->ref[0] = ref_idx[1];
	mb->top4x4inter->ref[1] = ref_idx[1];
	for (int i = 0; i < 2; ++i) {
		mb->left4x4inter->mv[0][i][0] = mv[0][0][0];
		mb->left4x4inter->mv[0][i][1] = mv[0][0][1];
		mb->left4x4inter->mv[0][2 + i][0] = mv[1][0][0];
		mb->left4x4inter->mv[0][2 + i][1] = mv[1][0][1];
		mb->left4x4inter->mv[1][i][0] = mv[0][1][0];
		mb->left4x4inter->mv[1][i][1] = mv[0][1][1];
		mb->left4x4inter->mv[1][2 + i][0] = mv[1][1][0];
		mb->left4x4inter->mv[1][2 + i][1] = mv[1][1][1];
	}
	mb->left4x4inter->ref[0] = ref_idx[0];
	mb->left4x4inter->ref[1] = ref_idx[1];
	return residual_chroma(mb, cbp, st, avail, ResidualBlock);
}

static inline void calc_mv8x16left(h264d_mb_current *mb, int16_t pmv[], const int16_t *&mvd_a, const int16_t *&mvd_b, int ref_idx, int avail)
{
	prev_mb_t *pmb;
	const int16_t *mva, *mvb, *mvc;
	int idx_map;

	if (avail & 1) {
		pmb = mb->left4x4inter;
		mvd_a = pmb->mv[1][0];
		if (ref_idx == pmb->ref[0]) {
			pmv[0] = pmb->mv[0][0][0];
			pmv[1] = pmb->mv[0][0][1];
			mvd_b = (avail & 2) ? mb->top4x4inter->mv[1][0] : zero_mv;
			return;
		}
		mva = pmb->mv[0][0];
	} else {
		mvd_a = mva = zero_mv;
	}
	idx_map = 0;
	if (avail & 2) {
		pmb = mb->top4x4inter;
		idx_map |= (ref_idx == pmb->ref[0]) * 2;
		idx_map |= (ref_idx == pmb->ref[1]) * 4;
		avail |= 4;
		mvb = pmb->mv[0][0];
		mvd_b = pmb->mv[1][0];
		mvc = pmb->mv[0][2];
	} else {
		mvd_b = mvb = zero_mv;
		avail &= ~4;
		if (avail & 8) {
			idx_map |= (ref_idx == mb->lefttop_ref) * 4;
			mvc = mb->lefttop_mv;
		} else {
			mvc = zero_mv;
		}
	}
	determine_pmv(mva, mvb, mvc, pmv, avail, idx_map);
}

static inline void calc_mv8x16right(h264d_mb_current *mb, int16_t pmv[], const int16_t *&mvd_a, const int16_t *&mvd_b, int ref_idx, int avail, int prev_ref, const int16_t *prev_mv)
{
	prev_mb_t *pmb;
	const int16_t *mva, *mvb, *mvc;
	int idx_map;

	idx_map = 0;
	if (avail & 4) {
		pmb = mb->top4x4inter + 1;
		if (ref_idx == pmb->ref[0]) {
			pmv[0] = pmb->mv[0][0][0];
			pmv[1] = pmb->mv[0][0][1];
			mvd_a = prev_mv + 2;
			mvd_b = (avail & 2) ? mb->top4x4inter->mv[1][2] : zero_mv;
			return;
		}
		mvc = pmb->mv[0][0];
	} else if (avail & 2) {
		pmb = mb->top4x4inter;
		idx_map = (ref_idx == pmb->ref[0]) * 4;
		mvd_b = pmb->mv[1][2];
		if (idx_map) {
			pmv[0] = pmb->mv[0][1][0];
			pmv[1] = pmb->mv[0][1][1];
			mvd_a = prev_mv + 2;
			return;
		} else {
			mvc = pmb->mv[0][1];
		}
	} else {
		mvc = zero_mv;
	}
	/* left block are always available */
	idx_map |= (ref_idx == prev_ref);
	mva = prev_mv;
	mvd_a = prev_mv + 2;
	avail |= 1;
	if (avail & 2) {
		pmb = mb->top4x4inter;
		idx_map |= (ref_idx == pmb->ref[1]) * 2;
		mvb = pmb->mv[0][2];
		mvd_b = pmb->mv[1][2];
	} else {
		mvd_b = mvb = zero_mv;
	}
	determine_pmv(mva, mvb, mvc, pmv, avail, idx_map);
}

template <typename F0 ,typename F1, typename F2, typename F3, typename F4>
static int mb_inter8x16(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail,
			 F0 RefIdx8x16,
			 F1 MvdXY,
			 F2 CodedBlockPattern,
			 F3 QpDelta,
			 F4 ResidualBlock)
{
	int8_t ref_idx[2];
	int16_t mv[2][2][2];
	const int16_t *mvd_a, *mvd_b;
	int t;
	uint32_t cbp;
	uint32_t top4x4, left4x4;
	uint32_t str_vert, str_horiz;

	t = *mb->num_ref_idx_l0_active_minus1;
	if (0 < t) {
		RefIdx8x16(mb, st, t, ref_idx, avail);
	} else {
		ref_idx[0] = 0;
		ref_idx[1] = 0;
	}
	calc_mv8x16left(mb, mv[0][0], mvd_a, mvd_b, ref_idx[0], avail);
	MvdXY(mb, st, mv[0][1], mvd_a, mvd_b);
	mv[0][0][0] += mv[0][1][0];
	mv[0][0][1] += mv[0][1][1];
	inter_pred8x8(mb, mv[0][0], 8, 16, ref_idx[0], 0, 0);

	calc_mv8x16right(mb, mv[1][0], mvd_a, mvd_b, ref_idx[1], avail, ref_idx[0], mv[0][0]);
	MvdXY(mb, st, mv[1][1], mvd_a, mvd_b);
	mv[1][0][0] += mv[1][1][0];
	mv[1][0][1] += mv[1][1][1];
	inter_pred8x8(mb, mv[1][0], 8, 16, ref_idx[1], 8, 0);
	VC_CHECK;

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

	deblock_info_t *deb = mb->deblock_curr;
	deb->qpy = mb->qp;
	deb->qpc = mb->qp_chroma;
	if (mb->y != 0) {
		if (mb->top4x4inter->type <= MB_IPCM) {
			deb->str4_vert = 1;
			str_vert |= 1;
		}
		str_vert = str_previous_coef(str_vert, top4x4);
		str_vert = str_mv_calc16x8_left(str_vert, ref_idx[0], ref_idx[1], mv[0][0], mv[1][0], mb->top4x4inter); /* same as 16x8 left */
	}
	deb->str_vert = str_vert;
	if (mb->x != 0) {
		if (mb->left4x4inter->type <= MB_IPCM) {
			deb->str4_horiz = 1;
			str_horiz |= 1;
		}
		str_horiz = str_previous_coef(str_horiz, left4x4);
		str_horiz = str_mv_calc16x16(str_horiz, mv[0][0], ref_idx[0], mb->left4x4inter);
	}
	deb->str_horiz = str_mv_calc16x8_vert(str_horiz, ref_idx[0], ref_idx[1], mv[0][0], mv[1][0]); /* same as 16x8 vert */

	mb->left4x4pred = 0x22222222;
	*mb->top4x4pred = 0x22222222;
	mb->lefttop_ref = mb->top4x4inter->ref[1];
	mb->lefttop_mv[0] = mb->top4x4inter->mv[0][3][0];
	mb->lefttop_mv[1] = mb->top4x4inter->mv[0][3][1];
	for (int i = 0; i < 2; ++i) {
		mb->top4x4inter->mv[0][i][0] = mv[0][0][0];
		mb->top4x4inter->mv[0][i][1] = mv[0][0][1];
		mb->top4x4inter->mv[0][2 + i][0] = mv[1][0][0];
		mb->top4x4inter->mv[0][2 + i][1] = mv[1][0][1];
		mb->top4x4inter->mv[1][i][0] = mv[0][1][0];
		mb->top4x4inter->mv[1][i][1] = mv[0][1][1];
		mb->top4x4inter->mv[1][2 + i][0] = mv[1][1][0];
		mb->top4x4inter->mv[1][2 + i][1] = mv[1][1][1];
	}
	mb->top4x4inter->ref[0] = ref_idx[0];
	mb->top4x4inter->ref[1] = ref_idx[1];
	for (int i = 0; i < 4; ++i) {
		mb->left4x4inter->mv[0][i][0] = mv[1][0][0];
		mb->left4x4inter->mv[0][i][1] = mv[1][0][1];
		mb->left4x4inter->mv[1][i][0] = mv[1][1][0];
		mb->left4x4inter->mv[1][i][1] = mv[1][1][1];
	}
	mb->left4x4inter->ref[0] = ref_idx[1];
	mb->left4x4inter->ref[1] = ref_idx[1];
	return residual_chroma(mb, cbp, st, avail, ResidualBlock);
}

typedef struct {
	int16_t ref;
	int16_t mv[2][4][2];
} prev8x8_t;

static inline void inter_pred8x8(h264d_mb_current *mb, const int16_t mv[], int width, int height, int ref_idx, int offsetx, int offsety)
{
	h264d_frame_info_t *frmi = mb->frame;
	h264d_frame *frms = &(frmi->frames[frmi->refs[ref_idx].frame_idx]);
	int stride, vert_size;
	int posx, posy;
	int mvx, mvy;
	mb_pred_t pred;

	stride = mb->max_x * 16;
	vert_size = mb->max_y * 16;
	pred.src_luma = frms->luma;
	pred.src_chroma = frms->chroma;
	pred.dst_luma = mb->luma + offsety * stride + offsetx;
	pred.dst_chroma = mb->chroma + (offsety >> 1) * stride + offsetx;
	mvx = mv[0];
	mvy = mv[1];
	pred.pos_x = posx = mb->x * 16 + (mvx >> 2) + offsetx;
	pred.pos_y = posy = mb->y * 16 + (mvy >> 2) + offsety;
	pred.src_luma = pred.src_luma + inter_pred_mvoffset_luma(posx - 2, posy - 2, stride);
	inter_pred_luma[mvy & 3][mvx & 3](&pred, width, height, stride, vert_size);
	pred.pos_x = mb->x * 16 + (mvx >> 3) * 2 + offsetx;
	pred.pos_y = mb->y * 8 + (mvy >> 3) + (offsety >> 1);
	inter_pred_chroma(&pred, mvx, mvy, width, height >> 1, stride, vert_size >> 1);
}

static inline void calc_mv8x8_sub8x8(h264d_mb_current *mb, int16_t *pmv, const int16_t *&mvd_a, const int16_t *&mvd_b, int avail, int ref_idx, int blk_idx, prev8x8_t *pblk)
{
	const int16_t *mva;
	const int16_t *mvb;
	const int16_t *mvc;
	prev_mb_t *pmb;
	int idx_map;

	if (blk_idx & 1) {
		idx_map = (ref_idx == pblk[blk_idx - 1].ref);
		mva = pblk[blk_idx - 1].mv[0][1];
		mvd_a = pblk[blk_idx - 1].mv[1][1];
		avail |= 1;
	} else if (avail & 1) {
		pmb = mb->left4x4inter;
		idx_map = (ref_idx == pmb->ref[blk_idx >> 1]);
		mva = pmb->mv[0][blk_idx];
		mvd_a = pmb->mv[1][blk_idx];
	} else {
		idx_map = 0;
		mvd_a = mva = zero_mv;
	}

	if (blk_idx & 2) {
		idx_map |= (ref_idx == pblk[blk_idx - 2].ref) * 2;
		mvb = pblk[blk_idx - 2].mv[0][2];
		mvd_b = pblk[blk_idx - 2].mv[1][2];
		avail |= 2;
	} else if (avail & 2) {
		pmb = mb->top4x4inter;
		idx_map |= (ref_idx == pmb->ref[blk_idx]) * 2;
		mvb = pmb->mv[0][blk_idx * 2];
		mvd_b = pmb->mv[1][blk_idx * 2];
	} else {
		mvd_b = mvb = zero_mv;
	}

	switch (blk_idx) {
	case 0:
		if (avail & 2) {
			pmb = mb->top4x4inter;
			idx_map |= (ref_idx == pmb->ref[1]) * 4;
			mvc = pmb->mv[0][2];
			avail |= 4;
		} else if (avail & 8) {
			idx_map |= (ref_idx == mb->lefttop_ref) * 4;
			mvc = mb->lefttop_mv;
			avail |= 4;
		} else {
			avail &= ~4;
			mvc = zero_mv;
		}
		break;
	case 1:
		if (avail & 4) {
			pmb = mb->top4x4inter + 1;
			idx_map |= (ref_idx == pmb->ref[0]) * 4;
			mvc = pmb->mv[0][0];
		} else if (avail & 2) {
			pmb = mb->top4x4inter;
			idx_map |= (ref_idx == pmb->ref[0]) * 4;
			mvc = pmb->mv[0][1];
		} else {
			mvc = zero_mv;
		}
		break;
	case 2:
		idx_map |= (ref_idx == pblk[1].ref) * 4;
		mvc = pblk[1].mv[0][2];
		avail |= 4;
		break;
	case 3:
		idx_map |= (ref_idx == pblk[0].ref) * 4;
		mvc = pblk[0].mv[0][3];
		avail |= 4;
		break;
	}
	determine_pmv(mva, mvb, mvc, pmv, avail, idx_map);
}

static inline void calc_mv8x8_sub8x4(h264d_mb_current *mb, int16_t *pmv, const int16_t *&mvd_a, const int16_t *&mvd_b, int avail, int ref_idx, int blk_idx, prev8x8_t *pblk, int y)
{
	const int16_t *mva;
	const int16_t *mvb;
	const int16_t *mvc;
	prev_mb_t *pmb;
	int idx_map;

	if (blk_idx & 1) {
		idx_map = (ref_idx == pblk[blk_idx - 1].ref);
		mva = pblk[blk_idx - 1].mv[0][y * 2 + 1];
		mvd_a = pblk[blk_idx - 1].mv[1][y * 2 + 1];
		avail |= 1;
	} else if (avail & 1) {
		pmb = mb->left4x4inter;
		idx_map = (ref_idx == pmb->ref[blk_idx >> 1]);
		mva = pmb->mv[0][(blk_idx & 2) + y];
		mvd_a = pmb->mv[1][(blk_idx & 2) + y];
	} else {
		idx_map = 0;
		mvd_a = mva = zero_mv;
	}

	if (y != 0) {
		idx_map |= 2;
		mvb = pblk[blk_idx].mv[0][0];
		mvd_b = pblk[blk_idx].mv[1][0];
		avail |= 2;
	} else if (blk_idx & 2) {
		idx_map |= (ref_idx == pblk[blk_idx - 2].ref) * 2;
		mvb = pblk[blk_idx - 2].mv[0][2];
		mvd_b = pblk[blk_idx - 2].mv[1][2];
		avail |= 2;
	} else if (avail & 2) {
		pmb = mb->top4x4inter;
		idx_map |= (ref_idx == pmb->ref[blk_idx & 1]) * 2;
		mvb = pmb->mv[0][blk_idx * 2];
		mvd_b = pmb->mv[1][blk_idx * 2];
	} else {
		mvd_b = mvb = zero_mv;
	}

	switch (blk_idx) {
	case 0:
		if (y == 0) {
			if (avail & 2) {
				pmb = mb->top4x4inter;
				idx_map |= (ref_idx == pmb->ref[1]) * 4;
				avail |= 4;
				mvc = pmb->mv[0][2];
			} else if (avail & 8) {
				idx_map |= (ref_idx == mb->lefttop_ref) * 4;
				avail |= 4;
				mvc = mb->lefttop_mv;
			} else {
				avail &= ~4;
				mvc = zero_mv;
			}
		} else if (avail & 1) {
			pmb = mb->left4x4inter;
			idx_map |= (ref_idx == pmb->ref[0]) * 4;
			mvc = pmb->mv[0][0];
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
				idx_map |= (ref_idx == pmb->ref[0]) * 4;
				mvc = pmb->mv[0][0];
				avail |= 4;
			} else if (avail & 2) {
				pmb = mb->top4x4inter;
				idx_map |= (ref_idx == pmb->ref[0]) * 4;
				mvc = pmb->mv[0][1];
				avail |= 4;
			} else {
				mvc = zero_mv;
			}
		} else {
			idx_map |= (ref_idx == pblk[0].ref) * 4;
			mvc = pblk[0].mv[0][1];
			avail |= 4;
		}
		break;
	case 2:
		if (y == 0) {
			idx_map |= (ref_idx == pblk[1].ref) * 4;
			mvc = pblk[1].mv[0][2];
			avail |= 4;
		} else if (avail & 1) {
			pmb = mb->left4x4inter;
			idx_map |= (ref_idx == pmb->ref[1]) * 4;
			mvc = pmb->mv[0][2];
			avail |= 4;
		} else {
			avail &= ~4;
			mvc = zero_mv;
		}
		break;
	case 3:
		idx_map |= (ref_idx == pblk[y * 2].ref) * 4;
		mvc = pblk[y * 2].mv[0][3 - y * 2];
		avail |= 4;
		break;
	}
	determine_pmv(mva, mvb, mvc, pmv, avail, idx_map);
}

static inline void calc_mv8x8_sub4x8(h264d_mb_current *mb, int16_t *pmv, const int16_t *&mvd_a, const int16_t *&mvd_b, int avail, int ref_idx, int blk_idx, prev8x8_t *pblk, int x)
{
	const int16_t *mva;
	const int16_t *mvb;
	const int16_t *mvc;
	prev_mb_t *pmb;
	int idx_map;

	if (x == 1) {
		idx_map = 1;
		mva = pblk[blk_idx].mv[0][0];
		mvd_a = pblk[blk_idx].mv[1][0];
		avail |= 1;
	} else if (blk_idx & 1) {
		idx_map = (ref_idx == pblk[blk_idx - 1].ref);
		mva = pblk[blk_idx - 1].mv[0][1];
		mvd_a = pblk[blk_idx - 1].mv[1][1];
		avail |= 1;
	} else if (avail & 1) {
		pmb = mb->left4x4inter;
		idx_map = (ref_idx == pmb->ref[blk_idx >> 1]);
		mva = pmb->mv[0][blk_idx]; /* blk_idx shall be 0 or 2 */
		mvd_a = pmb->mv[1][blk_idx];
	} else {
		idx_map = 0;
		mvd_a = mva = zero_mv;
	}

	if (blk_idx & 2) {
		idx_map |= (ref_idx == pblk[blk_idx - 2].ref) * 2;
		mvb = pblk[blk_idx - 2].mv[0][2 + x];
		mvd_b = pblk[blk_idx - 2].mv[1][2 + x];
		avail |= 2;
	} else if (avail & 2) {
		pmb = mb->top4x4inter;
		idx_map |= (ref_idx == pmb->ref[blk_idx & 1]) * 2;
		mvb = pmb->mv[0][blk_idx * 2 + x];
		mvd_b = pmb->mv[1][blk_idx * 2 + x];
	} else {
		mvd_b = mvb = zero_mv;
	}

	switch (blk_idx) {
	case 0:
		if (avail & 2) {
			pmb = mb->top4x4inter;
			idx_map |= (ref_idx == pmb->ref[x]) * 4;
			mvc = pmb->mv[0][x + 1];
			avail |= 4;
		} else {
			avail &= ~4;
			if (x == 0 && (avail & 8)) {
				idx_map |= (ref_idx == mb->lefttop_ref) * 4;
				mvc = mb->lefttop_mv;
			} else {
				mvc = zero_mv;
			}
		}
		break;
	case 1:
		if (x == 0) {
			if (avail & 2) {
				pmb = mb->top4x4inter;
				idx_map |= (ref_idx == pmb->ref[1]) * 4;
				mvc = pmb->mv[0][3];
				avail |= 4;
			} else {
				avail &= ~4;
				mvc = zero_mv;
			}
		} else {
			if (avail & 4) {
				pmb = mb->top4x4inter + 1;
				idx_map |= (ref_idx == pmb->ref[0]) * 4;
				mvc = pmb->mv[0][0];
			} else if (avail & 2) {
				pmb = mb->top4x4inter;
				idx_map |= (ref_idx == pmb->ref[1]) * 4;
				if (0 <= pmb->ref[1]) {
					mvc = pmb->mv[0][2];
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
		if (x == 0) {
			idx_map |= (ref_idx == pblk[0].ref) * 4;
			mvc = pblk[0].mv[0][3];
		} else {
			idx_map |= (ref_idx == pblk[1].ref) * 4;
			mvc = pblk[1].mv[0][2];
		}
		break;
	case 3:
		idx_map |= (ref_idx == pblk[1].ref) * 4;
		avail |= 4;
		if (x == 0) {
			mvc = pblk[1].mv[0][3];
		} else {
			mvc = pblk[1].mv[0][2];
		}
		break;
	}
	determine_pmv(mva, mvb, mvc, pmv, avail, idx_map);
}

static inline void calc_mv8x8_sub4x4(h264d_mb_current *mb, int16_t *pmv, const int16_t *&mvd_a, const int16_t *&mvd_b, int avail, int ref_idx, int blk_idx, prev8x8_t *pblk, int xy)
{
	const int16_t *mva;
	const int16_t *mvb;
	const int16_t *mvc;
	prev_mb_t *pmb;
	int idx_map;

	if (xy & 1) {
		idx_map = 1;
		mva = pblk[blk_idx].mv[0][xy - 1];
		mvd_a = pblk[blk_idx].mv[1][xy - 1];
		avail |= 1;
	} else if (blk_idx & 1) {
		idx_map = (ref_idx == pblk[blk_idx - 1].ref);
		mva = pblk[blk_idx - 1].mv[0][xy + 1]; /* xy shall be 0 or 2 */
		mvd_a = pblk[blk_idx - 1].mv[1][xy + 1];
		avail |= 1;
	} else if (avail & 1) {
		pmb = mb->left4x4inter;
		idx_map = (ref_idx == pmb->ref[blk_idx >> 1]);
		mva = pmb->mv[0][blk_idx + (xy >> 1)]; /* blk_idx shall be 0 or 2 */
		mvd_a = pmb->mv[1][blk_idx + (xy >> 1)];
	} else {
		idx_map = 0;
		mvd_a = mva = zero_mv;
	}

	if (xy & 2) {
		idx_map |= 2;
		mvb = pblk[blk_idx].mv[0][xy - 2];
		mvd_b = pblk[blk_idx].mv[1][xy - 2];
		avail |= 2;
	} else if (blk_idx & 2) {
		idx_map |= (ref_idx == pblk[blk_idx - 2].ref) * 2;
		mvb = pblk[blk_idx - 2].mv[0][2 + (xy & 1)];
		mvd_b = pblk[blk_idx - 2].mv[1][2 + (xy & 1)];
		avail |= 2;
	} else if (avail & 2) {
		pmb = mb->top4x4inter;
		idx_map |= (ref_idx == pmb->ref[blk_idx & 1]) * 2;
		mvb = pmb->mv[0][blk_idx * 2 + (xy & 1)];
		mvd_b = pmb->mv[1][blk_idx * 2 + (xy & 1)];
	} else {
		mvd_b = mvb = zero_mv;
	}

	switch (blk_idx) {
	case 0:
		switch (xy) {
		case 0:
			if (avail & 2) {
				pmb = mb->top4x4inter;
				idx_map |= (ref_idx == pmb->ref[0]) * 4;
				avail |= 4;
				mvc = pmb->mv[0][1];
			} else if (avail & 8) {
				avail &= ~4;
				idx_map |= (ref_idx == mb->lefttop_ref) * 4;
				mvc = mb->lefttop_mv;
			} else {
				avail &= ~4;
				mvc = zero_mv;
			}
			break;
		case 1:
			if (avail & 2) {
				pmb = mb->top4x4inter;
				idx_map |= (ref_idx == pmb->ref[1]) * 4;
				avail |= 4;
				mvc = pmb->mv[0][2];
			} else {
				avail &= ~4;
				mvc = zero_mv;
			}
			break;
		case 2:
			idx_map |= 4;
			avail |= 4;
			mvc = pblk[blk_idx].mv[0][1];
			break;
		case 3:
			idx_map |= 4;
			avail |= 4;
			mvc = pblk[blk_idx].mv[0][0];
			break;
		}
		break;
	case 1:
		switch (xy) {
		case 0:
			if (avail & 2) {
				pmb = mb->top4x4inter;
				idx_map |= (ref_idx == pmb->ref[1]) * 4;
				mvc = pmb->mv[0][3];
				avail |= 4;
			} else {
				avail &= ~4;
				mvc = zero_mv;
			}
			break;
		case 1:
			if (avail & 4) {
				pmb = mb->top4x4inter + 1;
				idx_map |= (ref_idx == pmb->ref[0]) * 4;
				mvc = pmb->mv[0][0];
			} else if (avail & 2) {
				pmb = mb->top4x4inter;
				idx_map |= (ref_idx == pmb->ref[1]) * 4;
				mvc = pmb->mv[0][2];
				avail |= 4;
			} else {
				mvc = zero_mv;
			}
			break;
		case 2:
		case 3:
			idx_map |= 4;
			avail |= 4;
			mvc = pblk[blk_idx].mv[0][3 - xy];
			break;
		}
		break;
	case 2:
		avail |= 4;
		switch (xy) {
		case 0:
		case 1:
			idx_map |= (ref_idx == pblk[xy].ref) * 4;
			mvc = pblk[xy].mv[0][3 - xy];
			break;
		case 2:
		case 3:
			idx_map |= 4;
			mvc = pblk[2].mv[0][3 - xy];
			break;
		}
		break;
	case 3:
		avail |= 4;
		switch (xy) {
		case 0:
		case 1:
			idx_map |= (ref_idx == pblk[1].ref) * 4;
			mvc = pblk[1].mv[0][3 - xy];
			break;
		case 2:
		case 3:
			idx_map |= 4;
			mvc = pblk[3].mv[0][3 - xy];
			break;
		}
		break;
	}
	determine_pmv(mva, mvb, mvc, pmv, avail, idx_map);
}

template<typename F0>
static void sub_mb8x8(h264d_mb_current *mb, dec_bits *st, int avail, int ref_idx, int blk_idx, prev8x8_t *pblk,
		      F0 MvdXY)
{
	const int16_t *mvd_a;
	const int16_t *mvd_b;
	int16_t mv[2];
	int16_t mvd[2];
	int mvx, mvy, mvdx, mvdy;

	calc_mv8x8_sub8x8(mb, mv, mvd_a, mvd_b, avail, ref_idx, blk_idx, pblk);
	MvdXY(mb, st, mvd, mvd_a, mvd_b);
	pblk += blk_idx;
	pblk->ref = ref_idx;
	mvdx = mvd[0];
	mvdy = mvd[1];
	mvx = mv[0] + mvdx;
	mvy = mv[1] + mvdy;
	mv[0] = mvx;
	mv[1] = mvy;
	for (int i = 0; i < 4; ++i) {
		pblk->mv[0][i][0] = mvx;
		pblk->mv[0][i][1] = mvy;
		pblk->mv[1][i][0] = mvdx;
		pblk->mv[1][i][1] = mvdy;
	}
	inter_pred8x8(mb, mv, 8, 8, ref_idx, (blk_idx & 1) * 8, (blk_idx & 2) * 4);
}

template<typename F0>
static void sub_mb8x4(h264d_mb_current *mb, dec_bits *st, int avail, int ref_idx, int blk_idx, prev8x8_t *pblk,
		      F0 MvdXY)

{
	const int16_t *mvd_a;
	const int16_t *mvd_b;
	int16_t mv[2];
	int16_t mvd[2];
	prev8x8_t *p = pblk + blk_idx;

	p->ref = ref_idx;
	for (int y = 0; y < 2; ++y) {
		int mvx, mvy, mvdx, mvdy;
		calc_mv8x8_sub8x4(mb, mv, mvd_a, mvd_b, avail, ref_idx, blk_idx, pblk, y);
		MvdXY(mb, st, mvd, mvd_a, mvd_b);
		mvdx = mvd[0];
		mvdy = mvd[1];
		mvx = mv[0] + mvdx;
		mvy = mv[1] + mvdy;
		mv[0] = mvx;
		mv[1] = mvy;
		p->mv[0][y * 2][0] = mvx;
		p->mv[1][y * 2][0] = mvdx;
		p->mv[0][y * 2][1] = mvy;
		p->mv[1][y * 2][1] = mvdy;
		p->mv[0][y * 2 + 1][0] = mvx;
		p->mv[1][y * 2 + 1][0] = mvdx;
		p->mv[0][y * 2 + 1][1] = mvy;
		p->mv[1][y * 2 + 1][1] = mvdy;
		inter_pred8x8(mb, mv, 8, 4, ref_idx, (blk_idx & 1) * 8, ((blk_idx & 2) + y) * 4);
	}
}

template <typename F0>
static void sub_mb4x8(h264d_mb_current *mb, dec_bits *st, int avail, int ref_idx, int blk_idx, prev8x8_t *pblk,
		      F0 MvdXY)
{
	const int16_t *mvd_a;
	const int16_t *mvd_b;
	int16_t mv[2];
	int16_t mvd[2];
	prev8x8_t *p = pblk + blk_idx;

	p->ref = ref_idx;
	for (int x = 0; x < 2; ++x) {
		int mvx, mvy, mvdx, mvdy;
		calc_mv8x8_sub4x8(mb, mv, mvd_a, mvd_b, avail, ref_idx, blk_idx, pblk, x);
		MvdXY(mb, st, mvd, mvd_a, mvd_b);
		mvdx = mvd[0];
		mvdy = mvd[1];
		mvx = mv[0] + mvdx;
		mvy = mv[1] + mvdy;
		mv[0] = mvx;
		mv[1] = mvy;
		p->mv[0][x][0] = mvx;
		p->mv[1][x][0] = mvdx;
		p->mv[0][x][1] = mvy;
		p->mv[1][x][1] = mvdy;
		p->mv[0][x + 2][0] = mvx;
		p->mv[1][x + 2][0] = mvdx;
		p->mv[0][x + 2][1] = mvy;
		p->mv[1][x + 2][1] = mvdy;
		inter_pred8x8(mb, mv, 4, 8, ref_idx, (blk_idx & 1) * 8 + x * 4, (blk_idx & 2) * 4);
	}
}

template <typename F0>
static void sub_mb4x4(h264d_mb_current *mb, dec_bits *st, int avail, int ref_idx, int blk_idx, prev8x8_t *pblk,
		      F0 MvdXY)
{
	const int16_t *mvd_a;
	const int16_t *mvd_b;
	int16_t mv[2];
	int16_t mvd[2];
	prev8x8_t *p = pblk + blk_idx;

	p->ref = ref_idx;
	for (int xy = 0; xy < 4; ++xy) {
		int mvx, mvy, mvdx, mvdy;
		calc_mv8x8_sub4x4(mb, mv, mvd_a, mvd_b, avail, ref_idx, blk_idx, pblk, xy);
		MvdXY(mb, st, mvd, mvd_a, mvd_b);
		mvdx = mvd[0];
		mvdy = mvd[1];
		mvx = mv[0] + mvdx;
		mvy = mv[1] + mvdy;
		mv[0] = mvx;
		mv[1] = mvy;
		p->mv[0][xy][0] = mvx;
		p->mv[0][xy][1] = mvy;
		p->mv[1][xy][0] = mvdx;
		p->mv[1][xy][1] = mvdy;
		inter_pred8x8(mb, mv, 4, 4, ref_idx, (blk_idx & 1) * 8 + (xy & 1) * 4, (blk_idx & 2) * 4 + (xy & 2) * 2);
	}
}

struct mvd_xy_cavlc {
	void operator()(h264d_mb_current *mb, dec_bits *st, int16_t mv[], const int16_t mva[], const int16_t mvb[]) {
	       mv[0] = se_golomb(st);
	       mv[1] = se_golomb(st);
	}
};

static void sub_mb8x8_cavlc(h264d_mb_current *mb, dec_bits *st, int avail, int ref_idx, int blk_idx, prev8x8_t *pblk)
{
	sub_mb8x8(mb, st, avail, ref_idx, blk_idx, pblk, mvd_xy_cavlc());
}

static void sub_mb8x4_cavlc(h264d_mb_current *mb, dec_bits *st, int avail, int ref_idx, int blk_idx, prev8x8_t *pblk)
{
	sub_mb8x4(mb, st, avail, ref_idx, blk_idx, pblk, mvd_xy_cavlc());
}

static void sub_mb4x8_cavlc(h264d_mb_current *mb, dec_bits *st, int avail, int ref_idx, int blk_idx, prev8x8_t *pblk)
{
	sub_mb4x8(mb, st, avail, ref_idx, blk_idx, pblk, mvd_xy_cavlc());
}

static void sub_mb4x4_cavlc(h264d_mb_current *mb, dec_bits *st, int avail, int ref_idx, int blk_idx, prev8x8_t *pblk)
{
	sub_mb4x4(mb, st, avail, ref_idx, blk_idx, pblk, mvd_xy_cavlc());
}

static void (* const sub_mb_cavlc[4])(h264d_mb_current *mb, dec_bits *st, int avail, int ref_idx, int blk_idx, prev8x8_t *pblk) = {
	sub_mb8x8_cavlc,
	sub_mb8x4_cavlc,
	sub_mb4x8_cavlc,
	sub_mb4x4_cavlc
};

struct sub_mbs_cavlc {
	void operator()(h264d_mb_current *mb, dec_bits *st, int avail, int8_t sub_mb_type[], int8_t ref_idx[], prev8x8_t curr_blk[]) {
		for (int i = 0; i < 4; ++i) {
			sub_mb_cavlc[sub_mb_type[i]](mb, st, avail, ref_idx[i], i, curr_blk);
		}
	}
};

static uint32_t str_mv_calc8x8_top(uint32_t str, const prev8x8_t *p, const prev_mb_t *top)
{
	for (int i = 0; i < 4; ++i) {
		int i8x8 = i >> 1;
		if (!(str & (3 << (i * 2))) && ((p[i8x8].ref != top->ref[i8x8]) || (16 <= DIF_SQUARE(p[i8x8].mv[0][i & 1][0], top->mv[0][i][0])) || (16 <= DIF_SQUARE(p[i8x8].mv[0][i & 1][1], top->mv[0][i][1])))) {
			str |= 1 << (i * 2);
		}
	}
	return str;
}

static uint32_t str_mv_calc8x8_vert(uint32_t str, const prev8x8_t *p)
{
	for (int i = 0; i < 4; ++i) {
		int i8x8 = i >> 1;
		int shift = i * 2 + 8;
		if (!(str & (3 << shift)) && ((16 <= DIF_SQUARE(p[i8x8].mv[0][i & 1][0], p[i8x8].mv[0][(i & 1) + 2][0])) || (16 <= DIF_SQUARE(p[i8x8].mv[0][i & 1][1], p[i8x8].mv[0][(i & 1) + 2][1])))) {
			str |= (1 << shift);
		}
	}
	for (int i = 0; i < 4; ++i) {
		int i8x8 = i >> 1;
		int shift = i * 2 + 16;
		if (!(str & (3 << shift)) && ((p[i8x8].ref != p[i8x8 + 2].ref) || (16 <= DIF_SQUARE(p[i8x8].mv[0][(i & 1) + 2][0], p[i8x8 + 2].mv[0][i & 1][0])) || (16 <= DIF_SQUARE(p[i8x8].mv[0][(i & 1) + 2][1], p[i8x8 + 2].mv[0][i & 1][1])))) {
			str |= (1 << shift);
		}
	}
	for (int i = 0; i < 4; ++i) {
		int i8x8 = (i >> 1) + 2;
		int shift = i * 2 + 24;
		if (!(str & (3 << shift)) && ((16 <= DIF_SQUARE(p[i8x8].mv[0][i & 1][0], p[i8x8].mv[0][(i & 1) + 2][0])) || (16 <= DIF_SQUARE(p[i8x8].mv[0][i & 1][1], p[i8x8].mv[0][(i & 1) + 2][1])))) {
			str |= (1 << shift);
		}
	}
	return str;
}

static uint32_t str_mv_calc8x8_left(uint32_t str, const prev8x8_t *p, const prev_mb_t *left)
{
	for (int i = 0; i < 4; ++i) {
		int i8x8 = i >> 1;
		if (!(str & (3 << (i * 2))) && ((p[i & 2].ref != left->ref[i8x8]) || (16 <= DIF_SQUARE(p[i & 2].mv[0][(i & 1) * 2][0], left->mv[0][i][0])) || (16 <= DIF_SQUARE(p[i & 2].mv[0][(i & 1) * 2][1], left->mv[0][i][1])))) {
			str |= 1 << (i * 2);
		}
	}
	return str;
}


static uint32_t str_mv_calc8x8_horiz(uint32_t str, const prev8x8_t *p)
{
	for (int i = 0; i < 4; ++i) {
		int shift = i * 2 + 8;
		if (!(str & (3 << shift)) && ((16 <= DIF_SQUARE(p[i & 2].mv[0][(i & 1) * 2][0], p[i & 2].mv[0][(i & 1) * 2 + 1][0])) || (16 <= DIF_SQUARE(p[i & 2].mv[0][(i & 1) * 2][1], p[i & 2].mv[0][(i & 1) * 2 + 1][1])))) {
			str |= (1 << shift);
		}
	}
	for (int i = 0; i < 4; ++i) {
		int shift = i * 2 + 16;
		if (!(str & (3 << shift)) && ((p[i & 2].ref != p[(i & 2) + 1].ref) || (16 <= DIF_SQUARE(p[i & 2].mv[0][(i & 1) * 2 + 1][0], p[(i & 2) + 1].mv[0][(i & 1) * 2][0])) || (16 <= DIF_SQUARE(p[i & 2].mv[0][(i & 1) * 2 + 1][1], p[(i & 2) + 1].mv[0][(i & 1) * 2][1])))) {
			str |= (1 << shift);
		}
	}
	for (int i = 0; i < 4; ++i) {
		int shift = i * 2 + 24;
		if (!(str & (3 << shift)) && ((16 <= DIF_SQUARE(p[(i & 2) + 1].mv[0][(i & 1) * 2][0], p[(i & 2) + 1].mv[0][(i & 1) * 2 + 1][0])) || (16 <= DIF_SQUARE(p[(i & 2) + 1].mv[0][(i & 1) * 2][1], p[(i & 2) + 1].mv[0][(i & 1) * 2 + 1][1])))) {
			str |= (1 << shift);
		}
	}
	return str;
}

template <typename F0, typename F1, typename F2, typename F3, typename F4, typename F5>
static int mb_inter8x8(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail,
		       F0 SubMbType,
		       F1 RefIdx8x8,
		       F2 SubMb,
		       F3 CodedBlockPattern,
		       F4 QpDelta,
		       F5 ResidualBlock)
{
	prev8x8_t curr_blk[4];
	int8_t sub_mb_type[4];
	int8_t ref_idx[4];
	int t;
	uint32_t cbp;
	uint32_t str_vert, str_horiz;
	uint32_t left4x4, top4x4;

	if (SubMbType(mb, st, sub_mb_type) < 0) {
		return -1;
	}
	t = *mb->num_ref_idx_l0_active_minus1;
	if (t && (mb->type != MB_P8x8REF0)) {
		RefIdx8x8(mb, st, t, ref_idx, sub_mb_type, avail);
	} else {
		memset(ref_idx, 0, sizeof(ref_idx));
	}
	SubMb(mb, st, avail, sub_mb_type, ref_idx, curr_blk);
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
			str_vert |= 1;
		}
		str_vert = str_previous_coef(str_vert, top4x4);
		str_vert = str_mv_calc8x8_top(str_vert, curr_blk, mb->top4x4inter);
	}
	deb->str_vert = str_mv_calc8x8_vert(str_vert, curr_blk);
	if (mb->x != 0) {
		if (mb->left4x4inter->type <= MB_IPCM) {
			deb->str4_horiz = 1;
			str_horiz |= 1;
		}
		str_horiz = str_previous_coef(str_horiz, left4x4);
		str_horiz = str_mv_calc8x8_left(str_horiz, curr_blk, mb->left4x4inter);
	}
	deb->str_horiz = str_mv_calc8x8_horiz(str_horiz, curr_blk);


	mb->left4x4pred = 0x22222222;
	*mb->top4x4pred = 0x22222222;

	mb->lefttop_ref = mb->top4x4inter->ref[1];
	mb->lefttop_mv[0] = mb->top4x4inter->mv[0][3][0];
	mb->lefttop_mv[1] = mb->top4x4inter->mv[0][3][1];

	mb->top4x4inter->ref[0] = curr_blk[2].ref;
	memcpy(mb->top4x4inter->mv[0][0], curr_blk[2].mv[0][2], sizeof(curr_blk[2].mv[0]) / 2);
	memcpy(mb->top4x4inter->mv[1][0], curr_blk[2].mv[1][2], sizeof(curr_blk[2].mv[1]) / 2);
	mb->top4x4inter->ref[1] = curr_blk[3].ref;
	memcpy(mb->top4x4inter->mv[0][2], curr_blk[3].mv[0][2], sizeof(curr_blk[3].mv[0]) / 2);
	memcpy(mb->top4x4inter->mv[1][2], curr_blk[3].mv[1][2], sizeof(curr_blk[3].mv[1]) / 2);

	mb->left4x4inter->ref[0] = curr_blk[1].ref;
	mb->left4x4inter->mv[0][0][0] = curr_blk[1].mv[0][1][0];
	mb->left4x4inter->mv[0][0][1] = curr_blk[1].mv[0][1][1];
	mb->left4x4inter->mv[0][1][0] = curr_blk[1].mv[0][3][0];
	mb->left4x4inter->mv[0][1][1] = curr_blk[1].mv[0][3][1];
	mb->left4x4inter->mv[1][0][0] = curr_blk[1].mv[1][1][0];
	mb->left4x4inter->mv[1][0][1] = curr_blk[1].mv[1][1][1];
	mb->left4x4inter->mv[1][1][0] = curr_blk[1].mv[1][3][0];
	mb->left4x4inter->mv[1][1][1] = curr_blk[1].mv[1][3][1];
	mb->left4x4inter->ref[1] = curr_blk[3].ref;
	mb->left4x4inter->mv[0][2][0] = curr_blk[3].mv[0][1][0];
	mb->left4x4inter->mv[0][2][1] = curr_blk[3].mv[0][1][1];
	mb->left4x4inter->mv[0][3][0] = curr_blk[3].mv[0][3][0];
	mb->left4x4inter->mv[0][3][1] = curr_blk[3].mv[0][3][1];
	mb->left4x4inter->mv[1][2][0] = curr_blk[3].mv[1][1][0];
	mb->left4x4inter->mv[1][2][1] = curr_blk[3].mv[1][1][1];
	mb->left4x4inter->mv[1][3][0] = curr_blk[3].mv[1][3][0];
	mb->left4x4inter->mv[1][3][1] = curr_blk[3].mv[1][3][1];

	return residual_chroma(mb, cbp, st, avail, ResidualBlock);
}

struct sub_mb_type_cavlc {
	int operator()(h264d_mb_current *mb, dec_bits *st, int8_t *sub_mb_type) {
		for (int i = 0; i < 4; ++i) {
			READ_UE_RANGE(sub_mb_type[i], st, 3);
		}
		return 0;
	}
};

struct ref_idx16x16_cavlc {
	int operator()(h264d_mb_current *mb, dec_bits *st, int limit, int avail) {
		return te_golomb(st, limit);
	}
};

struct ref_idx16x8_cavlc {
	void operator()(h264d_mb_current *mb, dec_bits *st, int limit, int8_t *ref_idx, int avail) {
		ref_idx[0] = te_golomb(st, limit);
		ref_idx[1] = te_golomb(st, limit);
	}
};

struct ref_idx8x8_cavlc {
	void operator()(h264d_mb_current *mb, dec_bits *st, int limit, int8_t *ref_idx, const int8_t *sub_mb_type, int avail) {
		for (int i = 0; i < 4; ++i) {
			ref_idx[i] = te_golomb(st, limit);
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

static int mb_inter8x8_cavlc(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_inter8x8(mb, mbc, st, avail, sub_mb_type_cavlc(), ref_idx8x8_cavlc(), sub_mbs_cavlc(), cbp_inter_cavlc(), qp_delta_cavlc(), residual_block_cavlc());
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
	{mb_intrapcm, 0, 0}, /* FIXME: IPCM */
	{mb_inter16x16_cavlc, 0, 0},
	{mb_inter16x8_cavlc, 0, 0},
	{mb_inter8x16_cavlc, 0, 0},
	{mb_inter8x8_cavlc, 0, 0},
	{mb_inter8x8_cavlc, 0, 0},
};

/** Convert MB type number into unified order:
 * Intra < Inter < Bidirectional
 */
static int adjust_mb_type(int mb_type, int slice_type)
{
	if (slice_type == P_SLICE) {
		if (mb_type < 30) {
			mb_type -= 5;
			return mb_type < 0 ? mb_type + 31 : mb_type;
		} else {
			return -1;
		}
	} else if (slice_type == B_SLICE) {
		mb_type -= 23;
		return mb_type < 0 ? mb_type + 55 : mb_type;
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

	if ((avail & 3) != 3) {
		return;
	}
	pmb = mb->left4x4inter;
	if (pmb->ref[0] == 0 && pmb->mv[0][0][0] == 0 && pmb->mv[0][0][1] == 0) {
		return;
	}
	pmb = mb->top4x4inter;
	if (pmb->ref[0] == 0 && pmb->mv[0][0][0] == 0 && pmb->mv[0][0][1] == 0) {
		return;
	}
	calc_mv16x16(mb, pmv, mvd_a, mvd_b, 0, avail);
	mv[0] += pmv[0];
	mv[1] += pmv[1];
}

static int skip_mbs(h264d_mb_current *mb, uint32_t skip_mb_num, int slice_type)
{
	uint32_t max_mb_run = mb->max_x * mb->max_y - (mb->y * mb->max_x + mb->x);
	uint32_t left4x4, top4x4;

	mb->type = MB_PSKIP;
	skip_mb_num = skip_mb_num < max_mb_run ? skip_mb_num : max_mb_run;
	mb->left4x4pred = 0x22222222;
	left4x4 = mb->left4x4coef;
	mb->left4x4coef = 0;
	mb->cbp = 0;
	mb->cbf = 0;
	do {
		deblock_info_t *deb;
		uint32_t str_vert, str_horiz;
		int16_t mv[2] = {
			0, 0
		};
		calc_mv_pskip(mb, mv, get_availability(mb));
		inter_pred8x8(mb, mv, 16, 16, 0, 0, 0);
		*mb->top4x4pred = 0x22222222;
		top4x4 = *mb->top4x4coef;
		*mb->top4x4coef = 0;

		deb = mb->deblock_curr;
		deb->qpy = mb->qp;
		deb->qpc = mb->qp_chroma;
		str_vert = 0;
		str_horiz = 0;
		if (mb->y != 0) {
			if (mb->top4x4inter->type <= MB_IPCM) {
				deb->str4_vert = 1;
				str_vert = 1;
			}
			str_vert = str_previous_coef(str_vert, top4x4);
			str_vert = str_mv_calc16x16(str_vert, mv, 0, mb->top4x4inter);
		}
		deb->str_vert = str_vert;
		if (mb->x != 0) {
			if (mb->left4x4inter->type <= MB_IPCM) {
				deb->str4_horiz = 1;
				str_horiz = 1;
			}
			str_horiz = str_previous_coef(str_horiz, left4x4);
			str_horiz = str_mv_calc16x16(str_horiz, mv, 0, mb->left4x4inter);
		}
		deb->str_horiz = str_horiz;
		left4x4 = 0;

		mb->prev_qp_delta = 0;
		mb->lefttop_ref = mb->top4x4inter->ref[1];
		mb->lefttop_mv[0] = mb->top4x4inter->mv[0][3][0];
		mb->lefttop_mv[1] = mb->top4x4inter->mv[0][3][1];
		mb->top4x4inter->ref[0] = 0;
		mb->top4x4inter->ref[1] = 0;
		mb->left4x4inter->ref[0] = 0;
		mb->left4x4inter->ref[1] = 0;
		for (int i = 0; i < 4; ++i) {
			mb->top4x4inter->mv[0][i][0] = mv[0];
			mb->top4x4inter->mv[0][i][1] = mv[1];
			mb->left4x4inter->mv[0][i][0] = mv[0];
			mb->left4x4inter->mv[0][i][1] = mv[1];
		}
		memset(mb->top4x4inter->mv[1][0], 0, sizeof(mb->top4x4inter->mv[1]));
		memset(mb->left4x4inter->mv[1][0], 0, sizeof(mb->left4x4inter->mv[1]));

		if (increment_mb_pos(mb) < 0) {
			return -1;
		}
	} while (--skip_mb_num);
	return 0;
}

static int more_rbsp_data(dec_bits *st)
{
	uint32_t dt;
	int bits;

	bits = not_aligned_bits(st);
	if (bits == 0) {
		bits = 8;
	}
	dt = show_bits(st, bits);
	if (dt == (1U << (bits - 1))) {
		/* FIXME */
		dt = show_bits(st, bits + 24);
		if (1 < (dt & 0xffffff)) {
			return 1;
		} else {
			return 0;
		}
	} else {
		return 1;
	}
}

static int post_process(h264d_context *h2d, h264d_mb_current *mb);
void init_cabac_context(h264d_cabac_t *cabac, int slice_qp, int idc);
static void init_cabac_engine(h264d_cabac_t *cb, dec_bits *st);
static inline int cabac_decode_terminate(h264d_cabac_t *cb, dec_bits *st);
static int mb_skip_cabac(h264d_mb_current *mb, dec_bits *st, int slice_type);

static int slice_data(h264d_context *h2d, dec_bits *st)
{
	h264d_slice_header *hdr = h2d->slice_header;
	h264d_pps *pps = &h2d->pps_i[hdr->pic_parameter_set_id];
	h264d_mb_current *mb = &h2d->mb_current;
	int mb_addr;
	int is_ae = pps->entropy_coding_mode_flag;

	if (is_ae) {
		init_cabac_context(mb->cabac, mb->qp, (hdr->slice_type == I_SLICE) ? 0 : hdr->cabac_init_idc + 1);
		byte_align(st);
		init_cabac_engine(mb->cabac, st);
	}
	mb_addr = hdr->first_mb_in_slice;
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
	void operator()(uint8_t *dst, int q0, int q1, int p0, int p1, int a, int b2, int tc0) {
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
	void operator()(uint8_t *dst, int q0, int q1, int p0, int p1, int a, int b2, int tc0) {
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
		str >>= 2;
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
		str >>= 2;
	}
}


template <int N>
struct Strength4v {
	void operator()(uint8_t *dst, int q0, int q1, int p0, int p1, int a, int b2, int tc0, int stride) {
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
	void operator()(uint8_t *dst, int q0, int q1, int p0, int p1, int a, int b2, int tc0, int stride) {
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
		str >>= 2;
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
			if ((y != 0) && (!idc || !mb->firstline) && (str & 255)) {
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
			VC_CHECK;
			curr++;
			luma += 16;
			chroma += 16;
		}
		luma += stride * 15;
		chroma += stride * 7;
	}
}

static inline void marking_sliding_window(h264d_ref_frame_t *refs, int frame_ptr, int frame_num, int max_frame_num, int num_ref_frames)
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
		refs +=	empty_idx;
	} else {
		refs += min_idx;
	}
	refs->in_use = SHORT_TERM;
	refs->frame_idx = frame_ptr;
	refs->num = frame_num;
}

static h264d_ref_frame_t *mmco_search(h264d_ref_frame_t *refs, int in_use, uint32_t target_num)
{
	int i = 16;
	do {
		if (refs->num == target_num && refs->in_use == in_use) {
			break;
		}
		refs++;
	} while (--i);
	return i == 0 ? 0 : refs;
}

static void mmco_discard(h264d_ref_frame_t *refs, int in_use, int target_num)
{
	h264d_ref_frame_t *ref = mmco_search(refs, in_use, target_num);
	if (ref) {
		ref->in_use = NOT_IN_USE;
	}
}

static void mmco_op1(const h264d_mmco *mmco, h264d_ref_frame_t *refs, int frame_ptr, int frame_num, int max_frame_num)
{
	int num = frame_num - mmco->arg1 - 1;
	while (num < 0) {
		num += max_frame_num;
	}
	mmco_discard(refs, SHORT_TERM, num);
}

static void mmco_op2(const h264d_mmco *mmco, h264d_ref_frame_t *refs, int frame_ptr, int frame_num, int max_frame_num)
{
	mmco_discard(refs, LONG_TERM, mmco->arg1);
}

static void mmco_op3(const h264d_mmco *mmco, h264d_ref_frame_t *refs, int frame_ptr, int frame_num, int max_frame_num)
{
	uint32_t long_num = mmco->arg2;
	uint32_t target_num = frame_num - mmco->arg1 - 1;
	int i = 16;

	while ((int)target_num < 0) {
		target_num += max_frame_num;
	}
	do {
		if (refs->in_use == LONG_TERM) {
			if (refs->num == long_num) {
				refs->in_use = NOT_IN_USE;
			}
		} else if (refs->in_use == SHORT_TERM) {
			if (refs->num == target_num) {
				refs->in_use = LONG_TERM;
				refs->num = long_num;
			}
		}
		refs++;
	} while (--i);
}

static void mmco_op4(const h264d_mmco *mmco, h264d_ref_frame_t *refs, int frame_ptr, int frame_num, int max_frame_num)
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

static void mmco_op5(const h264d_mmco *mmco, h264d_ref_frame_t *refs, int frame_ptr, int frame_num, int max_frame_num)
{
	int i = 16;
	do {
		refs->in_use = NOT_IN_USE;
		refs++;
	} while (--i);
}

static h264d_ref_frame_t *find_empty_ref(h264d_ref_frame_t *refs)
{
	int i = 16;
	h264d_ref_frame_t *empty_ref;
	do {
		if (refs->in_use == 0) {
			empty_ref = refs;
			break;
		}
		refs++;
	} while (--i);
	return i ? empty_ref : refs - 1;
}

static void mmco_op6(const h264d_mmco *mmco, h264d_ref_frame_t *refs, int frame_ptr, int frame_nuim, int max_frame_num)
{
	h264d_ref_frame_t *empty_ref = find_empty_ref(refs);
	empty_ref->in_use = LONG_TERM;
	empty_ref->frame_idx = frame_ptr;
	empty_ref->num = mmco->arg1;
}

static void insert_short_ref(h264d_ref_frame_t *refs, int frame_ptr, int frame_num)
{
	h264d_ref_frame_t *empty_ref = find_empty_ref(refs);
	while (empty_ref != refs) {
		*empty_ref = *(empty_ref - 1);
		empty_ref--;
	}
	empty_ref->in_use = SHORT_TERM;
	empty_ref->frame_idx = frame_ptr;
	empty_ref->num = frame_num;
}

static void (* const mmco_ops[6])(const h264d_mmco *mmco, h264d_ref_frame_t *refs, int frame_ptr, int frame_num, int max_frame_num) = {
	mmco_op1, mmco_op2, mmco_op3,
	mmco_op4, mmco_op5, mmco_op6,
};

static inline void marking_mmco(h264d_marking_t *mrk, h264d_ref_frame_t *refs, int frame_ptr, int frame_num, int max_frame_num)
{
	const h264d_mmco *mmco = mrk->non_idr.mmco;
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
		mmco_ops[op - 1](mmco, refs, frame_ptr, frame_num, max_frame_num);
		mmco++;
	} while (--i);
	if (!op6_detect) {
		insert_short_ref(refs, frame_ptr, op5_detect ? 0 : frame_num);
	}
}

static inline void post_ref_pic_marking(h264d_slice_header *hdr, int nal_unit_type, int max_frame_num, int num_ref_frames, h264d_mb_current *mb, h264d_ref_frame_t *refs)
{
	h264d_marking_t *mrk = &hdr->marking;
	int frame_num = hdr->frame_num;
	if (nal_unit_type == SLICE_IDR_NAL) {
		refs[0].in_use = mrk->idr.long_term_reference_flag ? LONG_TERM : SHORT_TERM;
		refs[0].frame_idx = mb->frame->index;
		refs[0].num = frame_num;
		for (int i = 1; i < 16; ++i) {
			refs[i].in_use = 0;
		}
	} else {
		int prev_frame_num = hdr->prev_frame_num;
		int gap = frame_num - prev_frame_num;
		while (gap < 0) {
			gap += max_frame_num;
		}
		if (0 < --gap) {
			if (16 < gap) {
				gap = 16;
				prev_frame_num = frame_num - 17;
			}
			do {
				if (max_frame_num <= ++prev_frame_num) {
					prev_frame_num -= max_frame_num;
				}
				marking_sliding_window(refs, mb->frame->index, prev_frame_num, max_frame_num, num_ref_frames);
			} while (--gap);
		}
		if (mrk->non_idr.adaptive_ref_pic_marking_mode_flag) {
			marking_mmco(mrk, refs, mb->frame->index, frame_num, max_frame_num);
		} else {
			marking_sliding_window(refs, mb->frame->index, frame_num, max_frame_num, num_ref_frames);
		}
	}
}

static inline int ref_pic_eval(const h264d_ref_frame_t& ref, int frame_num, int max_frame_num)
{
	int in_use = ref.in_use;
	int num = ref.num;
	if (in_use == SHORT_TERM) {
		return (frame_num < num) ? num - max_frame_num : num;
	} else if (in_use == LONG_TERM) {
		return -(max_frame_num + num);
	} else {
		return INT_MIN;
	}
}

static inline void swap_ref(h264d_ref_frame_t *left, h264d_ref_frame_t *right) {
	h264d_ref_frame_t tmp = *left;
	*left = *right;
	*right = tmp;
}

static inline void sort_ref(h264d_ref_frame_t *refs, int num_elem, int frame_num, int max_frame_num)
{
	for (int end = num_elem; 0 < end; --end) {
		int ev0 = ref_pic_eval(refs[0], frame_num, max_frame_num);
		int swapped = 0;
		for (int i = 1; i < end; ++i) {
			int ev1 = ref_pic_eval(refs[i], frame_num, max_frame_num);
			if (ev0 < ev1) {
				swap_ref(refs + i - 1, refs + i);
				swapped = 1;
			} else {
				ev0 = ev1;
			}
		}
		if (!swapped) {
			break;
		}
	}
}

static inline void ref_pic_init_ip(h264d_ref_frame_t *ref, int frame_num, int max_frame_num, int num_ref_frames)
{
	sort_ref(ref, num_ref_frames, frame_num, max_frame_num);
	for (int i = num_ref_frames; i < 16; ++i) {
		ref[i].in_use = 0;
	}
}

static int post_process(h264d_context *h2d, h264d_mb_current *mb)
{
	h264d_slice_header *hdr;
	int nal_id;
	int is_filled = (mb->y >= mb->max_y);

	hdr = h2d->slice_header;
	if (is_filled) {
		deblock_pb(mb);
		h264d_sps *sps = &h2d->sps_i[h2d->pps_i[hdr->pic_parameter_set_id].seq_parameter_set_id];
		int max_frame_num = 1 << sps->log2_max_frame_num;
		int num_ref_frames = sps->num_ref_frames;
		nal_id = h2d->id;
		if (nal_id & 0x60) {
			post_ref_pic_marking(hdr, nal_id & 31, max_frame_num, num_ref_frames, &h2d->mb_current, hdr->reorder[0].ref_frames);
			if (hdr->slice_type == B_SLICE) {
				post_ref_pic_marking(hdr, nal_id & 31, max_frame_num, num_ref_frames, &h2d->mb_current, hdr->reorder[1].ref_frames);
			}
		}
		ref_pic_init_ip(hdr->reorder[0].ref_frames, hdr->frame_num, max_frame_num, num_ref_frames);
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

void init_cabac_context(h264d_cabac_t *cabac, int slice_qp, int idc)
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
	int add = !((avail & 1) && (mb->left4x4inter->type == MB_BDIRECT16x16)) + !((avail & 2) && (mb->top4x4inter->type == MB_BDIRECT16x16));
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
	if ((avail & 1) && (mb->left4x4inter->type != MB_PSKIP)) {
		offset += 1;
	}
	if ((avail & 2) && (mb->top4x4inter->type != MB_PSKIP)) {
		offset += 1;
	}
	return cabac_decode_decision(mb->cabac, st, offset);
}

struct intra4x4pred_mode_cabac {
	int operator()(int a, int b, dec_bits *st, h264d_cabac_t *cb)
	{
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
	uint32_t operator()(h264d_mb_current *mb, dec_bits *st, h264d_cabac_t *cb, int avail) {
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
	uint32_t operator()(h264d_mb_current *mb, dec_bits *st, int avail) {
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

static inline uint32_t unary_cabac(h264d_cabac_t *cb, dec_bits *st, int idx, int limit)
{
	int x = 0;
	do {
		if (cabac_decode_decision(cb, st, idx)) {
			x = x + 1;
		} else {
			break;
		}
	} while (--limit);
	return x;
}

struct qp_delta_cabac {
	int operator()(h264d_mb_current *mb, dec_bits *st, h264d_cabac_t *cb, int avail) {
		int qp_delta;
		int ctx_idx = 60 + (mb->prev_qp_delta != 0);
		qp_delta = cabac_decode_decision(cb, st, ctx_idx);
		if (qp_delta) {
			qp_delta = unary_cabac(cb, st, ctx_idx, 52) + 1;
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

static inline int get_coeff_from_map_cabac(h264d_cabac_t *cb, dec_bits *st, int cat, int *coeff_map, int map_cnt, int *coeff, const int16_t *qmat, uint32_t dc_mask)
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
	return map_cnt;
}

struct residual_block_cabac {
	int operator()(h264d_mb_current *mb, int na, int nb, dec_bits *st, int *coeff, int num_coeff, const int16_t *qmat, int avail, int pos4x4, int cat, uint32_t dc_mask)
	{
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
		return get_coeff_from_map_cabac(cb, st, cat, coeff_map, map_cnt, coeff + coeff_offset, qmat + coeff_offset, dc_mask);
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

struct sub_mb_type_cabac {
	int operator()(h264d_mb_current *mb, dec_bits *st, int8_t *sub_mb_type) {
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
	void operator()(h264d_mb_current *mb, dec_bits *st, int16_t mv[], const int16_t mva[], const int16_t mvb[]) {
		h264d_cabac_t *cb = mb->cabac;
		mv[0] = mvd_cabac(mb, st, cb, 40, mva[0], mvb[0]);
		mv[1] = mvd_cabac(mb, st, cb, 47, mva[1], mvb[1]);
	}
};

static void sub_mb8x8_cabac(h264d_mb_current *mb, dec_bits *st, int avail, int ref_idx, int blk_idx, prev8x8_t *pblk)
{
	sub_mb8x8(mb, st, avail, ref_idx, blk_idx, pblk, mvd_xy_cabac());
}

static void sub_mb8x4_cabac(h264d_mb_current *mb, dec_bits *st, int avail, int ref_idx, int blk_idx, prev8x8_t *pblk)
{
	sub_mb8x4(mb, st, avail, ref_idx, blk_idx, pblk, mvd_xy_cabac());
}

static void sub_mb4x8_cabac(h264d_mb_current *mb, dec_bits *st, int avail, int ref_idx, int blk_idx, prev8x8_t *pblk)
{
	sub_mb4x8(mb, st, avail, ref_idx, blk_idx, pblk, mvd_xy_cabac());
}

static void sub_mb4x4_cabac(h264d_mb_current *mb, dec_bits *st, int avail, int ref_idx, int blk_idx, prev8x8_t *pblk)
{
	sub_mb4x4(mb, st, avail, ref_idx, blk_idx, pblk, mvd_xy_cabac());
}

static void (* const sub_mb_cabac[4])(h264d_mb_current *mb, dec_bits *st, int avail, int ref_idx, int blk_idx, prev8x8_t *pblk) = {
	sub_mb8x8_cabac,
	sub_mb8x4_cabac,
	sub_mb4x8_cabac,
	sub_mb4x4_cabac
};

struct sub_mbs_cabac {
	void operator()(h264d_mb_current *mb, dec_bits *st, int avail, int8_t sub_mb_type[], int8_t ref_idx[], prev8x8_t curr_blk[]) {
		for (int i = 0; i < 4; ++i) {
			sub_mb_cabac[sub_mb_type[i]](mb, st, avail, ref_idx[i], i, curr_blk);
		}
	}
};

static inline int ref_idx_cabac_sub(dec_bits *st, h264d_cabac_t *cb, int inc)
{
	int idx = 0;
	while (cabac_decode_decision(cb, st, 54 + inc)) {
		inc = (inc >> 2) + 4;
		idx += 1;
	}
	return idx;
}

struct ref_idx16x16_cabac {
	int operator()(h264d_mb_current *mb, dec_bits *st, int limit, int avail) {
		h264d_cabac_t *cb = mb->cabac;
		int inc;
		inc = ((avail & 1) && (0 < mb->left4x4inter->ref[0])) + ((avail & 2) && (0 < mb->top4x4inter->ref[0])) * 2;
		return ref_idx_cabac_sub(st, cb, inc);
	}
};

struct ref_idx16x8_cabac {
	void operator()(h264d_mb_current *mb, dec_bits *st, int limit, int8_t *ref_idx, int avail) {
		h264d_cabac_t *cb = mb->cabac;
		int inc;
		inc = ((avail & 1) && (0 < mb->left4x4inter->ref[0])) + ((avail & 2) && (0 < mb->top4x4inter->ref[0])) * 2;
		ref_idx[0] = ref_idx_cabac_sub(st, cb, inc);
		inc = ((avail & 1) && (0 < mb->left4x4inter->ref[1])) + (ref_idx[0] != 0) * 2;
		ref_idx[1] = ref_idx_cabac_sub(st, cb, inc);
	}
};

struct ref_idx8x16_cabac {
	void operator()(h264d_mb_current *mb, dec_bits *st, int limit, int8_t *ref_idx, int avail) {
		h264d_cabac_t *cb = mb->cabac;
		int inc;
		inc = ((avail & 1) && (0 < mb->left4x4inter->ref[0])) + ((avail & 2) && (0 < mb->top4x4inter->ref[0])) * 2;
		ref_idx[0] = ref_idx_cabac_sub(st, cb, inc);
		inc = (ref_idx[0] != 0) + ((avail & 2) && (0 < mb->top4x4inter->ref[1])) * 2;
		ref_idx[1] = ref_idx_cabac_sub(st, cb, inc);
	}
};

struct ref_idx8x8_cabac {
	void operator()(h264d_mb_current *mb, dec_bits *st, int limit, int8_t *ref_idx, const int8_t *sub_mb_type, int avail) {
		h264d_cabac_t *cb = mb->cabac;
		int inc;
		inc = ((avail & 1) && (0 < mb->left4x4inter->ref[0])) + ((avail & 2) && (0 < mb->top4x4inter->ref[0])) * 2;
		ref_idx[0] = ref_idx_cabac_sub(st, cb, inc);
		inc = (ref_idx[0] != 0) + ((avail & 2) && (0 < mb->top4x4inter->ref[1])) * 2;
		ref_idx[1] = ref_idx_cabac_sub(st, cb, inc);
		inc = ((avail & 1) && (0 < mb->left4x4inter->ref[1])) + (ref_idx[0] != 0) * 2;
		ref_idx[2] = ref_idx_cabac_sub(st, cb, inc);
		inc = (ref_idx[2] != 0) + (ref_idx[1] != 0) * 2;
		ref_idx[3] = ref_idx_cabac_sub(st, cb, inc);
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

static int mb_inter8x8_cabac(h264d_mb_current *mb, const mb_code *mbc, dec_bits *st, int avail)
{
	return mb_inter8x8(mb, mbc, st, avail, sub_mb_type_cabac(), ref_idx8x8_cabac(), sub_mbs_cabac(), cbp_cabac(), qp_delta_cabac(), residual_block_cabac());
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
	{mb_inter16x16_cabac, 0, 0},
	{mb_inter16x8_cabac, 0, 0},
	{mb_inter8x16_cabac, 0, 0},
	{mb_inter8x8_cabac, 0, 0},
	{mb_inter8x8_cabac, 0, 0},
	{0, 0, 0},
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

