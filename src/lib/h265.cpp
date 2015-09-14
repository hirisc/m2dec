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
#include <limits.h>
#include <cstdlib>
#include <algorithm>
#include "h265.h"
#include "m2d_macro.h"
#if (defined(__GNUC__) && defined(__SSE2__)) || defined(_M_IX86) || defined(_M_AMD64)
#define X86ASM
#endif

typedef enum {
	WRONG_PARAM = -1,
	OUT_OF_RANGE = -2
} h265d_error_t;

typedef enum {
	DIR_LEFT = 1,
	DIR_TOP = 2
} h265d_direction_t;

#define IS_ACTIVE(flag, bit) (((flag) & (1 << (bit))) != 0)
#define NUM_ELEM(arry) (sizeof(arry) / sizeof(arry[0]))
#define READ_CHECK_RANGE(val, dst, max, st) {uint32_t t = (val); (dst) = t; if ((max) < t) error_report(st);}
#define READ_CHECK_RANGE2(val, dst, min, max, st) {int32_t t = (val); (dst) = t; if (t < (min) || (max) < t) error_report(st);}
#define CHECK_RANGE(val, max, st) {if ((max) < (val)) error_report(st);}
#define CHECK_RANGE2(val, min, max, st) {if (((val) < (min)) || ((max) < (val))) error_report(st);}
#define NEGATE(val, sign) (((val) ^ (unsigned)-(signed)(sign)) + (sign))


int h265d_init(h265d_context *h2, int dpb_max, int (*header_callback)(void *, void *), void *arg) {
	if (!h2) {
		return -1;
	}
	h265d_data_t& h2d = *reinterpret_cast<h265d_data_t*>(h2);
	memset(&h2d, 0, sizeof(h2d));
	h2d.coding_tree_unit.cabac.context = reinterpret_cast<int8_t*>(&h2d.coding_tree_unit.context);
	h2d.header_callback = header_callback ? header_callback : header_dummyfunc;
	h2d.header_callback_arg = arg;
	dec_bits_open(&h2d.stream_i, m2d_load_bytes_skip03);
	return 0;
}

dec_bits *h265d_stream_pos(h265d_context *h2) {
	return &reinterpret_cast<h265d_data_t*>(h2)->stream_i;
}

static int set_second_frame(const h265d_sps_t& sps, h265d_ctu_t* ctu, uint8_t* pool) {
	uint8_t* next = pool;
	int col = sps.ctb_info.columns;
	if (ctu) {
		ctu->neighbour_top = reinterpret_cast<h265d_neighbour_t*>(next);
	}
	next += sizeof(ctu->neighbour_left) * col;
	if (ctu) {
		ctu->deblock_topedge = ctu->deblock_topedgebase = reinterpret_cast<h265d_deblocking_strength_t*>(next);
	}
	next += sizeof(ctu->deblock_topedge[0]) * (col << (sps.ctb_info.size_log2 - 3));
	size_t sao_right_len = sizeof(ctu->sao_vlines.right[0][0][0]) << sps.ctb_info.size_log2;
	size_t sao_line_len = (sizeof(ctu->sao_hlines[0][0].bottom[0]) * col) << sps.ctb_info.size_log2;
	size_t sao_map_len = (sizeof(ctu->sao_hlines[0][0].reserved_flag[0]) * col + 7) >> 3;
#ifdef X86ASM
#define ALIGN16(x) (((x) + 15) & ~15)
	next = (uint8_t*)ALIGN16((uintptr_t)next);
	sao_right_len = ALIGN16(sao_right_len);
	sao_line_len = ALIGN16(sao_line_len);
	sao_map_len = ALIGN16(sao_map_len);
#endif
	for (int i = 0; i < 2; ++i) {
		for (int col = 0; col < 2; ++col) {
			h265d_sao_hlines_t& hlines = ctu->sao_hlines[i][col];
			h265d_sao_vlines_t& vlines = ctu->sao_vlines;
			if (ctu) {
				vlines.right[i][col] = reinterpret_cast<uint8_t*>(next);
			}
			next += sao_right_len;
			if (ctu) {
				hlines.bottom = reinterpret_cast<uint8_t*>(next);
			}
			next += sao_line_len;
			if (ctu) {
				hlines.reserved_flag = reinterpret_cast<uint8_t*>(next);
			}
			next += sao_map_len;
		}
	}
	if (ctu) {
		ctu->sao_signbuf = reinterpret_cast<uint8_t*>(next);
	}
	next += (sizeof(ctu->sao_signbuf[0]) * col) << (sps.ctb_info.size_log2 - 2);
	if (ctu) {
		ctu->sao_map = reinterpret_cast<h265d_sao_map_t*>(next);
	}
	next += sizeof(ctu->sao_map[0]) * col * sps.ctb_info.rows;
	return next - pool;
}

int h265d_get_info(h265d_context *h2, m2d_info_t *info) {
	if (!h2 || !info) {
		return -1;
	}
	h265d_data_t& h2d = *reinterpret_cast<h265d_data_t*>(h2);
	h265d_sps_t& sps = h2d.sps[h2d.pps[h2d.slice_header.pps_id].sps_id];
	int width = sps.ctb_info.columns << sps.ctb_info.size_log2;
	int height = sps.ctb_info.rows << sps.ctb_info.size_log2;
	info->src_width = width;
	info->src_height = height;
	info->disp_width = width;
	info->disp_height = height;
	info->frame_num = sps.num_long_term_ref_pics_sps + sps.num_short_term_ref_pic_sets;
	info->crop[0] = sps.cropping[0];
	info->crop[1] = width - sps.pic_width_in_luma_samples + sps.cropping[1];
	info->crop[2] = sps.cropping[2];
	info->crop[3] = height - sps.pic_height_in_luma_samples + sps.cropping[3];
	info->additional_size = set_second_frame(sps, 0, 0);
	return 0;
}
static void init_dpb(h265d_dpb_t& dpb) {
	dpb.max = 16;
	dpb.size = 0;
	dpb.output = -1;
}

static void init_frame_info(h265d_frame_info_t& frame_info, int num_frame, m2d_frame_t *frame) {
	frame_info.num = num_frame;
	std::copy(frame, frame + num_frame, frame_info.frames);
	memset(frame_info.lru, 0, sizeof(frame_info.lru));
	init_dpb(frame_info.dpb);
}

class is_same_frame {
public:
	is_same_frame(int target) : target_(target) {}
	bool operator()(const h265d_dpb_elem_t& frm) const {
		return frm.frame_idx == target_;
	}
private:
	int target_;
};

static inline bool dpb_exist(const h265d_dpb_t& dpb, int frame_idx) {
	const h265d_dpb_elem_t* target = std::find_if(dpb.data, dpb.data + dpb.size, is_same_frame(frame_idx));
	return (target != dpb.data + dpb.size);
}

static inline void find_empty_frame(h265d_frame_info_t& frm) {
	int8_t *lru = frm.lru;
	int lru_len = NUM_ELEM(frm.lru);
	int max_idx = 0;
	int max_val = -1;
	int frm_num = frm.num;
	h265d_dpb_t dpb = frm.dpb;
	for (int i = 0; i < frm_num; ++i) {
		if (dpb_exist(dpb, i)) {
			lru[i] = 0;
		} else {
			lru[i] += 1;
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
	frm.index = max_idx;
	frm.curr_luma = frm.frames[max_idx].luma;
	frm.curr_chroma = frm.frames[max_idx].chroma;
}

int h265d_set_frames(h265d_context *h2, int num_frame, m2d_frame_t *frame, uint8_t *second_frame, int second_frame_size) {
	h265d_data_t& h2d = *reinterpret_cast<h265d_data_t*>(h2);
	if (!h2 || (num_frame < 1) || (NUM_ELEM(h2d.coding_tree_unit.frame_info.frames) < (size_t)num_frame) || !frame || !second_frame) {
		return -1;
	}
	h265d_ctu_t& ctu = h2d.coding_tree_unit;
	init_frame_info(ctu.frame_info, num_frame, frame);
	set_second_frame(h2d.sps[h2d.pps[h2d.slice_header.pps_id].sps_id], &ctu, second_frame);
/*
	h2d->slice_header->reorder[0].ref_frames = mb->frame->refs[0];
	h2d->slice_header->reorder[1].ref_frames = mb->frame->refs[1];
	return init_mb_buffer(mb, second_frame, second_frame_size); */
	return 0;
}

static void sub_layer_info_read(h265d_sub_layer_info_t& dst, uint32_t present, dec_bits& st) {
	if (IS_ACTIVE(present, 15)) {
		dst.sub_layer_profile_first8bit = get_bits(&st, 8);
		dst.sub_layer_profile_compatibility_flag = get_bits32(&st, 32);
		int max = NUM_ELEM(dst.sub_layer_second48bit);
		for (int i = 0; i < max; ++i) {
			dst.sub_layer_second48bit[i] = get_bits(&st, 8);
		}
	}
	if (IS_ACTIVE(present, 14)) {
		dst.sub_layer_level_idc = get_bits(&st, 8);
	}
}

static void profile_tier_level(uint8_t max_sub_layers_minus1, h265d_profile_tier_level_t& dst, dec_bits& st) {
	dst.max_sub_layer = max_sub_layers_minus1 + 1;
	dst.general_profile_first8bit = get_bits(&st, 8);
	dst.general_profile_compatibility_flag = get_bits32(&st, 32);
	for (int i = 0; i < NUM_ELEM(dst.general_second48bit); ++i) {
		dst.general_second48bit[i] = get_bits(&st, 8);
	}
	dst.general_level_idc = get_bits(&st, 8);
	if (max_sub_layers_minus1 != 0) {
		uint32_t present = get_bits(&st, 16);
		dst.sub_layer_profile_level_present_flag = present;
		for (int i = 0; i < max_sub_layers_minus1; ++i) {
			sub_layer_info_read(dst.sub_layer_info[i], present, st);
			present <<= 2;
		}
	}
}

static void error_report(dec_bits& st) {
	longjmp(st.jmp, -1);
}

static void vps_timing_info_read(h265d_vps_timing_info_t& dst, dec_bits& st) {
	dst.num_units_in_tick = get_bits32(&st, 32);
	dst.time_scale = get_bits32(&st, 32);
	if ((dst.poc_proportional_to_timing_flag = get_onebit(&st)) != 0) {
		dst.num_ticks_poc_diff_one_minus1 = ue_golomb(&st);
	}
	READ_CHECK_RANGE(ue_golomb(&st), dst.vps_num_hrd_parameters, 1024, st);
	/* current version omits rest of data */
}

static void sub_layer_reordering_info(h265d_sub_layer_reordering_info_t dst[], uint8_t info_present, uint32_t max_sub_layers_minus1, dec_bits& st) {
	for (uint32_t i = (info_present ? 0 : max_sub_layers_minus1); i <= max_sub_layers_minus1; ++i) {
		dst[i].max_dec_pic_buffering_minus1 = ue_golomb(&st);
		dst[i].max_num_reorder_pic = ue_golomb(&st);
		dst[i].max_latency_increase_plus1 = ue_golomb(&st);
	}
}

static void video_parameter_set(h265d_data_t& h2d, dec_bits& st) {
	h265d_vps_t& vps = h2d.vps;
	vps.id = get_bits(&st, 4);
	skip_bits(&st, 2);
	vps.max_layer = get_bits(&st, 6);
	uint8_t max_sub_layers_minus1 = get_bits(&st, 3);
	CHECK_RANGE(max_sub_layers_minus1, 6, st);
	vps.temporal_id_nesting_flag = get_onebit(&st);
	skip_bits(&st, 16);
	profile_tier_level(max_sub_layers_minus1, vps.profile_tier_level, st);
	vps.sub_layer_ordering_info_present_flag = get_onebit(&st);
	sub_layer_reordering_info(vps.max_buffering, vps.sub_layer_ordering_info_present_flag, max_sub_layers_minus1, st);
	vps.max_layer_id = get_bits(&st, 6);
	READ_CHECK_RANGE(ue_golomb(&st), vps.num_layer_sets_minus1, 1023, st);
	for (int i = 0; i < vps.num_layer_sets_minus1; ++i) {
		skip_bits(&st, vps.max_layer_id + 1);
	}
	vps.timing_info_present_flag = get_onebit(&st);
	if (vps.timing_info_present_flag) {
		vps_timing_info_read(vps.timing_info, st);
	}
}

static void conformance_window_read(h265d_sps_t& dst, dec_bits& st) {
	if ((dst.conformance_window_flag = get_onebit(&st)) != 0) {
		for (int i = 0; i < 4; ++i) {
			dst.cropping[i] = ue_golomb(&st);
		}
	} else {
		memset(dst.cropping, 0, sizeof(dst.cropping));
	}
}

static void scaling_list_read(uint8_t list[], int len, bool upper, dec_bits& st) {
	uint8_t coef = upper ? (se_golomb(&st) + 8) : 8;
	for (int i = 0; i < len; ++i) {
		coef = se_golomb(&st) + coef;
		list[i] = coef;
	}
}

static void scaling_list_read_all(uint8_t list[], int list_num, int list_len, bool upper, dec_bits& st) {
	for (int i = 0; i < list_num; ++i) {
		if (get_onebit(&st) == 0) {
			ue_golomb(&st);
		} else {
			scaling_list_read(list, list_len, upper, st);
		}
		list += list_len;
	}
}

static void scaling_list_data(h265d_scaling_list_data_t& dst, dec_bits& st) {
	assert(0);
	scaling_list_read_all(dst.scale0[0], 6, 16, false, st);
	scaling_list_read_all(dst.scale1[0], 6, 64, false, st);
	scaling_list_read_all(dst.scale2[0], 6, 64, true, st);
	scaling_list_read_all(dst.scale3[0], 2, 64, true, st);
}

static inline void scaling_list_default(h265d_ctu_t& dst, uint32_t size_log2, uint32_t qp) {
	static const uint8_t scaling_default4x4[16] = {
		16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16
	};
	static const uint8_t scaling_default1[64] = {
		16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 16, 17, 16, 17, 18,
		17, 18, 18, 17, 18, 21, 19, 20, 21, 20, 19, 21, 24, 22, 22, 24,
		24, 22, 22, 24, 25, 25, 27, 30, 27, 25, 25, 29, 31, 35, 35, 31,
		29, 36, 41, 44, 41, 36, 47, 54, 54, 47, 65, 70, 65, 88, 88, 115
	};
	static const uint8_t scaling_default2[64] = {
		16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 18,
		18, 18, 18, 18, 18, 20, 20, 20, 20, 20, 20, 20, 24, 24, 24, 24,
		24, 24, 24, 24, 25, 25, 25, 25, 25, 25, 25, 28, 28, 28, 28, 28,
		28, 33, 33, 33, 33, 33, 41, 41, 41, 41, 54, 54, 54, 71, 71, 91
	};
	static const uint8_t* const scaling_default[4][6] = {
		{
			scaling_default4x4,
			scaling_default4x4,
			scaling_default4x4,
			scaling_default4x4,
			scaling_default4x4,
			scaling_default4x4
		},
		{
			scaling_default1,
			scaling_default1,
			scaling_default1,
			scaling_default2,
			scaling_default2,
			scaling_default2
		},
		{
			scaling_default1,
			scaling_default1,
			scaling_default1,
			scaling_default2,
			scaling_default2,
			scaling_default2
		},
		{
			scaling_default1,
			scaling_default2,
			0, 0, 0, 0
		}
	};
	if (dst.pps->pps_scaling_list_data_present_flag) {
	} else {
	}
}

static int short_term_ref_pic_set_nopred_calc(h265d_short_term_ref_pic_elem_t& dst, uint32_t num_pics, bool pos, dec_bits& st) {
	int16_t* delta_poc = dst.delta_poc;
	uint16_t used_flag = 0;
	int val = 0;
	int used_cnt = 0;
	for (uint32_t i = 0; i < num_pics; ++i) {
		uint32_t delta;
		READ_CHECK_RANGE(ue_golomb(&st), delta, 32767, st);
		delta += 1;
		val += pos ? delta : -delta;
		delta_poc[i] = val;
		int used_bit = get_onebit(&st);
		used_flag |= used_bit << i;
		used_cnt += used_bit;
	}
	dst.used_by_curr_pic_flag = used_flag;
	return used_cnt;
}

static void short_term_ref_pic_set_nopred(h265d_short_term_ref_pic_set_t& dst, dec_bits& st) {
	uint32_t neg_pics, pos_pics;
	READ_CHECK_RANGE(ue_golomb(&st), neg_pics, 16, st);
	READ_CHECK_RANGE(ue_golomb(&st), pos_pics, 16 - neg_pics, st);
	dst.ref[0].num_pics = neg_pics;
	dst.ref[1].num_pics = pos_pics;
	int neg_num = short_term_ref_pic_set_nopred_calc(dst.ref[0], neg_pics, false, st);
	int pos_num = short_term_ref_pic_set_nopred_calc(dst.ref[1], pos_pics, true, st);
	dst.total_curr = neg_num + pos_num;
}

static inline int short_term_ref_pic_set_pred_core(int16_t dst_delta_poc[], const int16_t ref_delta_poc[], int delta_rps, uint16_t used_flag, uint16_t use_delta_flag,  uint16_t& used, int idx, bool neg, int j) {
	int dpoc = ref_delta_poc[j] + delta_rps;
	if (((neg ? -dpoc : dpoc) < 0) && (use_delta_flag & (1 << j))) {
		dst_delta_poc[idx] = dpoc;
		if (used_flag & (1 << j)) {
			used |= 1 << idx;
		}
		idx++;
	}
	return idx;
}

static inline int short_term_ref_pic_set_pred_neg(h265d_short_term_ref_pic_elem_t& dst, const h265d_short_term_ref_pic_elem_t& ref, int delta_rps, uint16_t used_flag, uint16_t use_delta_flag, uint16_t& used, int idx, bool neg) {
	for (int j = ref.num_pics - 1; 0 <= j; --j) {
		idx = short_term_ref_pic_set_pred_core(dst.delta_poc, ref.delta_poc, delta_rps, used_flag, use_delta_flag, used, idx, neg, j);
	}
	return idx;
}

static inline int short_term_ref_pic_set_pred_pos(h265d_short_term_ref_pic_elem_t& dst, const h265d_short_term_ref_pic_elem_t& ref, int delta_rps, uint16_t used_flag, uint16_t use_delta_flag, uint16_t& used, int idx, bool neg) {
	for (int j = 0; j < ref.num_pics; ++j) {
		idx = short_term_ref_pic_set_pred_core(dst.delta_poc, ref.delta_poc, delta_rps, used_flag, use_delta_flag, used, idx, neg, j);
	}
	return idx;
}

static void short_term_ref_pic_set_pred_part(h265d_short_term_ref_pic_set_t& dst, const h265d_short_term_ref_pic_set_t& ref, int delta_rps, uint16_t used_flag, uint16_t use_delta_flag, int s0) {
	uint16_t used0 = 0;
	int i = short_term_ref_pic_set_pred_neg(dst.ref[s0], ref.ref[s0 ^ 1], delta_rps, used_flag >> ((s0 != 0) ? 0 : ref.ref[0].num_pics), use_delta_flag >> ((s0 != 0) ? 0 : ref.ref[0].num_pics), used0, 0, s0 != 0);
	int mask = 1 << (ref.ref[0].num_pics + ref.ref[1].num_pics);
	if ((((s0 != 0) ? -delta_rps : delta_rps) < 0) && (use_delta_flag & mask)) {
		dst.ref[s0].delta_poc[i] = delta_rps;
		if (used_flag & mask) {
			used0 |= 1 << i;
		}
		i++;
	}
	dst.ref[s0].num_pics = short_term_ref_pic_set_pred_pos(dst.ref[s0], ref.ref[s0], delta_rps, used_flag >> ((s0 == 0) ? 0 : ref.ref[0].num_pics), use_delta_flag >> ((s0 == 0) ? 0 : ref.ref[0].num_pics), used0, i, s0 != 0);
	dst.ref[s0].used_by_curr_pic_flag = used0;
}

static void short_term_ref_pic_set_pred(h265d_short_term_ref_pic_set_t& dst, const h265d_short_term_ref_pic_set_t& ref, dec_bits& st) {
	int sign = get_onebit(&st);
	int abs_delta_rps = ue_golomb(&st) + 1;
	int delta_rps = sign ? -abs_delta_rps : abs_delta_rps;
	int num = ref.ref[0].num_pics + ref.ref[1].num_pics;
	uint16_t used_flag = 0;
	uint16_t use_delta_flag = 0;
	int used_cnt = 0;
	for (int j = 0; j <= num; ++j) {
		uint32_t used_by = get_onebit(&st);
		uint16_t bit = 1 << j;
		used_cnt += used_by;
		if (used_by) {
			used_flag |= bit;
			use_delta_flag |= bit;
		} else if (get_onebit(&st)) {
			use_delta_flag |= bit;
		}
	}
	short_term_ref_pic_set_pred_part(dst, ref, delta_rps, used_flag, use_delta_flag, 0);
	short_term_ref_pic_set_pred_part(dst, ref, delta_rps, used_flag, use_delta_flag, 1);
	dst.total_curr = used_cnt;
}

static void sps_short_term_ref_pic_set(h265d_short_term_ref_pic_set_t short_term_ref_pic_set[], uint32_t num, dec_bits& st) {
	short_term_ref_pic_set_nopred(short_term_ref_pic_set[0], st);
	for (uint32_t i = 1; i < num; ++i) {
		if (get_onebit(&st)) {
			short_term_ref_pic_set_pred(short_term_ref_pic_set[i], short_term_ref_pic_set[i - 1], st);
		} else {
			short_term_ref_pic_set_nopred(short_term_ref_pic_set[i], st);
		}
	}
}

static void vui_parameters(h265d_vui_parameters_t& dst, dec_bits& st) {
	dst.aspect_ratio_idc = (get_onebit(&st) != 0) ? get_bits(&st, 8) : 0;
	if (dst.aspect_ratio_idc == 255) {
		dst.sar_width = get_bits(&st, 16);
		dst.sar_height = get_bits(&st, 16);
	}
	if ((dst.overscan_info_present_flag = get_onebit(&st)) != 0) {
		dst.overscan_appropriate_flag = get_onebit(&st);
	}
	dst.video_format = 5;
	dst.video_full_range_flag = 0;
	dst.colour_primaries = 2;
	dst.transfer_characteristics = 2;
	dst.matrix_coeffs = 2;
	if (get_onebit(&st) != 0) {
		dst.video_format = get_bits(&st, 3);
		dst.video_full_range_flag = get_onebit(&st);
		if ((dst.colour_description_present_flag = get_onebit(&st)) != 0) {
			dst.colour_primaries = get_bits(&st, 8);
			dst.transfer_characteristics = get_bits(&st, 8);
			dst.matrix_coeffs = get_bits(&st, 8);
		}
	}
}

static inline uint32_t log2ceil(uint32_t num) {
	static const int8_t MultiplyDeBruijnBitPosition[32] = {
		0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
		8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
	};
	num = num | (num >> 1);
	num = num | (num >> 2);
	num = num | (num >> 4);
	num = num | (num >> 8);
	num = num | (num >> 16);
	return 1 + MultiplyDeBruijnBitPosition[(uint32_t)(num * 0x07C4ACDDU) >> 27];
}

static void set_ctb_info(h265d_sps_ctb_info_t& ctb_info, const h265d_sps_t& sps) {
	ctb_info.size_log2_min = sps.log2_min_luma_coding_block_size_minus3 + 3;
	uint32_t ctb_log2 = sps.log2_min_luma_coding_block_size_minus3 + 3 + sps.log2_diff_max_min_luma_coding_block_size;
	ctb_info.size_log2 = ctb_log2;
	ctb_info.pcm_log2_min = sps.pcm_enabled_flag ? sps.log2_min_pcm_luma_coding_block_size_minus3 + 3 : 8;
	ctb_info.pcm_log2 = ctb_info.pcm_log2_min + sps.log2_diff_max_min_pcm_luma_coding_block_size;
	ctb_info.transform_log2_min = sps.log2_min_transform_block_size_minus2 + 2;
	ctb_info.transform_log2 = ctb_info.transform_log2_min + sps.log2_diff_max_min_transform_block_size;
	uint32_t columns = (sps.pic_width_in_luma_samples + (1 << ctb_log2) - 1) >> ctb_log2;
	uint32_t rows = (sps.pic_height_in_luma_samples + (1 << ctb_log2) - 1) >> ctb_log2;
	ctb_info.stride = columns << ctb_log2;
	ctb_info.columns = columns;
	ctb_info.rows = rows;
	ctb_info.num_ctb_log2 = log2ceil(columns * rows);
}

static uint32_t sps_prefix(h265d_sps_prefix_t& dst, dec_bits& st) {
	dst.vps_id = get_bits(&st, 4);
	uint8_t max_sub_layers_minus1 = get_bits(&st, 3);
	dst.max_sub_layers_minus1 = max_sub_layers_minus1;
	dst.temporal_id_nesting_flag = get_onebit(&st);
	profile_tier_level(max_sub_layers_minus1, dst.profile_tier_level, st);
	uint32_t sps_id;
	READ_CHECK_RANGE(ue_golomb(&st), sps_id, 15, st);
	return sps_id;
}

static void sps_residual(h265d_sps_t& dst, const h265d_sps_prefix_t& prefix, dec_bits& st) {
	dst.prefix = prefix;
	READ_CHECK_RANGE(ue_golomb(&st), dst.chroma_format_idc, 3, st);
	if (dst.chroma_format_idc == 3) {
		dst.separate_colour_plane_flag = get_onebit(&st);
	}
	dst.pic_width_in_luma_samples = ue_golomb(&st);
	dst.pic_height_in_luma_samples = ue_golomb(&st);
	conformance_window_read(dst, st);
	READ_CHECK_RANGE(ue_golomb(&st), dst.bit_depth_luma_minus8, 6, st);
	READ_CHECK_RANGE(ue_golomb(&st), dst.bit_depth_chroma_minus8, 6, st);
	READ_CHECK_RANGE(ue_golomb(&st), dst.log2_max_pic_order_cnt_lsb_minus4, 12, st);
	dst.sub_layer_ordering_info_present_flag = get_onebit(&st);
	sub_layer_reordering_info(dst.max_buffering, dst.sub_layer_ordering_info_present_flag, prefix.max_sub_layers_minus1, st);
	READ_CHECK_RANGE(ue_golomb(&st), dst.log2_min_luma_coding_block_size_minus3, 2, st);
	READ_CHECK_RANGE(ue_golomb(&st), dst.log2_diff_max_min_luma_coding_block_size, 3, st);
	READ_CHECK_RANGE(ue_golomb(&st), dst.log2_min_transform_block_size_minus2, 2, st);
	READ_CHECK_RANGE(ue_golomb(&st), dst.log2_diff_max_min_transform_block_size, 3, st);
	READ_CHECK_RANGE(ue_golomb(&st), dst.max_transform_hierarchy_depth_inter, 5, st);
	READ_CHECK_RANGE(ue_golomb(&st), dst.max_transform_hierarchy_depth_intra, 5, st);
	if ((dst.scaling_list_enabled_flag = get_onebit(&st)) != 0) {
		if ((dst.scaling_list_data_present_flag = get_onebit(&st)) != 0) {
			scaling_list_data(dst.scaling_list_data, st);
		}
	}
	dst.amp_enabled_flag = get_onebit(&st);
	dst.sample_adaptive_offset_enabled_flag = get_onebit(&st);
	if ((dst.pcm_enabled_flag = get_onebit(&st)) != 0) {
		dst.pcm_sample_bit_depth_luma_minus1 = get_bits(&st, 4);
		dst.pcm_sample_bit_depth_chroma_minus1 = get_bits(&st, 4);
		READ_CHECK_RANGE(ue_golomb(&st), dst.log2_min_pcm_luma_coding_block_size_minus3, 2, st);
		READ_CHECK_RANGE(ue_golomb(&st), dst.log2_diff_max_min_pcm_luma_coding_block_size, 3, st);
		dst.pcm_loop_filter_disabled_flag = get_onebit(&st);
	} else {
		dst.log2_min_pcm_luma_coding_block_size_minus3 = 8;
	}
	READ_CHECK_RANGE(ue_golomb(&st), dst.num_short_term_ref_pic_sets, 64, st);
	sps_short_term_ref_pic_set(dst.short_term_ref_pic_set, dst.num_short_term_ref_pic_sets, st);
	if ((dst.long_term_ref_pics_present_flag = get_onebit(&st)) != 0) {
		uint32_t num_lt;
		READ_CHECK_RANGE(ue_golomb(&st), num_lt, 32, st);
		dst.num_long_term_ref_pics_sps = num_lt;
		uint32_t max_lsb = (1 << (dst.log2_max_pic_order_cnt_lsb_minus4 + 4)) - 1;
		dst.used_by_curr_pic_lt_sps_flag = 0;
		for (uint32_t i = 0; i < num_lt; ++i) {
			READ_CHECK_RANGE(ue_golomb(&st), dst.lt_ref_pic_poc_lsb_sps[i], max_lsb, st);
			if (get_onebit(&st)) {
				dst.used_by_curr_pic_lt_sps_flag |= (1 << i);
			}
		}
	}
	dst.sps_temporal_mvp_enabled_flag = get_onebit(&st);
	dst.strong_intra_smoothing_enabled_flag = get_onebit(&st);
	if ((dst.vui_parameters_present_flag = get_onebit(&st)) != 0) {
		vui_parameters(dst.vui_parameters, st);
	}
	set_ctb_info(dst.ctb_info, dst);
}

static void seq_parameter_set(h265d_data_t& h2d, dec_bits& st) {
	h265d_sps_prefix_t prefix;
	sps_residual(h2d.sps[sps_prefix(prefix, st)], prefix, st);
}

static void pps_tiles(h265d_tiles_t& dst, dec_bits& st, const h265d_sps_t& sps) {
	uint32_t columns_minus1;
	READ_CHECK_RANGE(ue_golomb(&st), columns_minus1, (uint32_t)(sps.ctb_info.columns - 1), st);
	dst.num_tile_columns_minus1 = columns_minus1;
	uint32_t rows_minus1;
	READ_CHECK_RANGE(ue_golomb(&st), rows_minus1, (uint32_t)(sps.ctb_info.rows - 1), st);
	dst.num_tile_rows_minus1 = rows_minus1;
	if ((dst.uniform_spacing_flag = get_onebit(&st)) == 0) {
		for (uint32_t i = 0; i < columns_minus1; ++i) {
			READ_CHECK_RANGE(ue_golomb(&st), dst.column_width_minus1[i], columns_minus1, st);
		}
		for (uint32_t i = 0; i < rows_minus1; ++i) {
			READ_CHECK_RANGE(ue_golomb(&st), dst.row_height_minus1[i], rows_minus1, st);
		}
	}
	dst.loop_filter_across_tiles_enabled_flag = get_onebit(&st);
}

static void pic_parameter_set(h265d_data_t& h2d, dec_bits& st) {
	uint32_t pps_id;
	READ_CHECK_RANGE(ue_golomb(&st), pps_id, 63, st);
	h265d_pps_t& dst = h2d.pps[pps_id];
	READ_CHECK_RANGE(ue_golomb(&st), dst.sps_id, 15, st);
	const h265d_sps_t& sps = h2d.sps[dst.sps_id];
	dst.dependent_slice_segments_enabled_flag = get_onebit(&st);
	dst.output_flag_present_flag = get_onebit(&st);
	dst.num_extra_slice_header_bits = get_bits(&st, 3);
	dst.sign_data_hiding_enabled_flag = get_onebit(&st);
	dst.cabac_init_present_flag = get_onebit(&st);
	READ_CHECK_RANGE(ue_golomb(&st), dst.num_ref_idx_lx_default_active_minus1[0], 14, st);
	READ_CHECK_RANGE(ue_golomb(&st), dst.num_ref_idx_lx_default_active_minus1[1], 14, st);
	READ_CHECK_RANGE(ue_golomb(&st), dst.init_qp_minus26, 52, st);
	dst.constrained_intra_pred_flag = get_onebit(&st);
	dst.transform_skip_enabled_flag = get_onebit(&st);
	if ((dst.cu_qp_delta_enabled_flag = get_onebit(&st)) != 0) {
		READ_CHECK_RANGE(ue_golomb(&st), dst.diff_cu_qp_delta_depth, 52, st);
	}
	READ_CHECK_RANGE2(se_golomb(&st), dst.pps_cb_qp_offset, -12, 12, st);
	READ_CHECK_RANGE2(se_golomb(&st), dst.pps_cr_qp_offset, -12, 12, st);
	dst.pps_slice_chroma_qp_offsets_present_flag = get_onebit(&st);
	dst.weighted_pred_flag = get_onebit(&st);
	dst.weighted_bipred_flag = get_onebit(&st);
	dst.transquant_bypass_enabled_flag = get_onebit(&st);
	dst.tiles_enabled_flag = get_onebit(&st);
	dst.entropy_coding_sync_enabled_flag = get_onebit(&st);
	if (dst.tiles_enabled_flag) {
		pps_tiles(dst.tiles, st, sps);
	}
	dst.pps_loop_filter_across_slices_enabled_flag = get_onebit(&st);
	if ((dst.deblocking_filter_control_present_flag = get_onebit(&st)) != 0) {
		dst.deblocking_filter_override_enabled_flag = get_onebit(&st);
		if ((dst.pps_deblocking_filter_disabled_flag = get_onebit(&st)) == 0) {
			READ_CHECK_RANGE2(se_golomb(&st), dst.pps_beta_offset_div2, -12, 12, st);
			READ_CHECK_RANGE2(se_golomb(&st), dst.pps_tc_offset_div2, -12, 12, st);
		}
	}
	if ((dst.pps_scaling_list_data_present_flag = get_onebit(&st)) != 0) {
		scaling_list_data(dst.scaling_list_data, st);
	}
	dst.lists_modification_present_flag = get_onebit(&st);
	dst.log2_parallel_merge_level_minus2 = ue_golomb(&st);
	CHECK_RANGE(dst.log2_parallel_merge_level_minus2 + 2, sps.ctb_info.size_log2, st);
	dst.slice_segment_header_extension_present_flag = get_onebit(&st);
	dst.pps_extension_flag = get_onebit(&st);
}

static void au_delimiter(h265d_data_t& h2d, dec_bits& st) {
	h265d_access_unit_delimite_t aud;
	aud.pic_type = get_bits(&st, 3);
}

static void skip_nal(h265d_data_t& h2d, dec_bits& st) {
}

static void entry_points(h265d_entry_point_t& dst, const h265d_pps_t& pps, const h265d_sps_t& sps, dec_bits& st) {
	uint32_t max_num;
	if (!pps.tiles_enabled_flag) {
		max_num = sps.ctb_info.rows - 1;
	} else if (!pps.entropy_coding_sync_enabled_flag) {
		max_num = (pps.tiles.num_tile_columns_minus1 + 1) * (pps.tiles.num_tile_rows_minus1 + 1) - 1;
	} else {
		max_num = (pps.tiles.num_tile_columns_minus1 + 1) * sps.ctb_info.rows - 1;
	}
	READ_CHECK_RANGE(ue_golomb(&st), dst.num_entry_point_offsets, max_num, st);
	uint32_t num_points = dst.num_entry_point_offsets;
	if (0 < num_points) {
		READ_CHECK_RANGE(ue_golomb(&st), dst.offset_len_minus1, 31, st);
		uint32_t offset_bits = dst.offset_len_minus1 + 1;
		uint32_t* points = dst.entry_point_offset_minus1;
		do {
			*points++ = ue_golomb(&st);
		} while (--num_points);
	}
}

static void slice_header_short_term_ref_pic_set(h265d_short_term_ref_pic_set_t& dst, const h265d_short_term_ref_pic_set_t sps_short_term_ref_pic_set[], uint32_t ref_num, dec_bits& st) {
	if (get_onebit(&st)) {
		int delta_idx_minus1;
		READ_CHECK_RANGE(ue_golomb(&st), delta_idx_minus1, ref_num, st);
		short_term_ref_pic_set_pred(dst, sps_short_term_ref_pic_set[ref_num - delta_idx_minus1 - 1], st);
	} else {
		short_term_ref_pic_set_nopred(dst, st);
	}
}

static void init_pic_order_cnt(h265d_poc_t& poc) {
	poc.lsb = 0;
	poc.msb = 0;
}

static void update_pic_order_cnt(h265d_slice_header_body_t& hdr, unsigned curr_lsb, const h265d_sps_t& sps) {
	h265d_poc_t& poc = hdr.slice_pic_order_cnt;
	unsigned prev_lsb = poc.lsb;
	unsigned prev_msb = poc.msb;
	unsigned max_lsb_div2 = 8 << sps.log2_max_pic_order_cnt_lsb_minus4;
	poc.lsb = curr_lsb;
	if ((BLA_W_LP <= hdr.nal_type) && (hdr.nal_type <= BLA_N_LP)) {
		poc.msb = 0;
	} else if ((curr_lsb < prev_lsb) && (max_lsb_div2 <= prev_lsb - curr_lsb)) {
		poc.msb++;
	} else if ((prev_lsb < curr_lsb) && (max_lsb_div2 < curr_lsb - prev_lsb)) {
		poc.msb--;
	}
}

static uint32_t pic_order_cnt(const h265d_poc_t& poc, const h265d_sps_t& sps) {
	return ((16 * poc.msb) << sps.log2_max_pic_order_cnt_lsb_minus4) + poc.lsb;
}

static void slice_header_nonidr(h265d_slice_header_body_t& dst, const h265d_pps_t& pps, const h265d_sps_t& sps, dec_bits& st) {
	unsigned pic_order_cnt_lsb = get_bits(&st, sps.log2_max_pic_order_cnt_lsb_minus4 + 4);
	update_pic_order_cnt(dst, pic_order_cnt_lsb, sps);
	if (get_onebit(&st)) {
		uint32_t idx = 0;
		if (1 < sps.num_short_term_ref_pic_sets) {
			idx = get_bits(&st, log2ceil(sps.num_short_term_ref_pic_sets));
		}
		dst.short_term_ref_pic_set = sps.short_term_ref_pic_set[idx];
	} else {
		slice_header_short_term_ref_pic_set(dst.short_term_ref_pic_set, sps.short_term_ref_pic_set, sps.num_short_term_ref_pic_sets, st);
	}
	if (sps.long_term_ref_pics_present_flag) {
		uint32_t num_long_term_sps = (0 < sps.num_long_term_ref_pics_sps) ? ue_golomb(&st) : 0;
		uint32_t num_long_term_pics = ue_golomb(&st);
		uint32_t num = num_long_term_sps + num_long_term_pics;
		assert(0);
		for (uint32_t i = 0; i < num; ++i) {
			if (i < num_long_term_sps) {
				if (1 < sps.num_long_term_ref_pics_sps) {
				}
			} else {
			}
		}
	}
	dst.slice_temporal_mvp_enabled_flag = sps.sps_temporal_mvp_enabled_flag ? get_onebit(&st) : 0;
}

static void pred_weight_table(h265d_slice_header_body_t& dst, const h265d_pps_t& pps, const h265d_sps_t& sps, dec_bits& st) {
	uint32_t luma_log2_weight_denom;
	assert(0);
}

static int init_ref_pic_list_short_lx(h265d_ref_pic_list_elem_t list[], const h265d_short_term_ref_pic_elem_t& ref, int curr_poc, int rest) {
	unsigned used = ref.used_by_curr_pic_flag;
	int i;
	for (i = 0; (i < ref.num_pics) && (i < rest); ++i) {
		if (used & 1) {
			int poc = curr_poc + ref.delta_poc[i];
			list[i].poc = poc;
			list[i].in_use = true;
			list[i].is_longterm = false;
		}
		used >>= 1;
	}
	return i;
}

static void init_ref_pic_list(h265d_ref_pic_list_elem_t list[][16], const h265d_short_term_ref_pic_set_t& srps, const h265d_slice_header_body_t& hdr, int curr_poc) {
	for (int lx = 0; lx < 2; ++lx) {
		int num_tmp = std::max(hdr.num_ref_idx_lx_active_minus1[lx] + 1, static_cast<int>(srps.total_curr));
		int idx = 0;
		while (idx < num_tmp) {
			idx += init_ref_pic_list_short_lx(list[lx], srps.ref[lx], curr_poc, num_tmp - idx);
			idx += init_ref_pic_list_short_lx(list[lx], srps.ref[lx ^ 1], curr_poc, num_tmp - idx);
		}
	}
}

static void slice_header_nonintra(h265d_slice_header_body_t& dst, const h265d_pps_t& pps, const h265d_sps_t& sps, dec_bits& st) {
	if (get_onebit(&st)) {
		READ_CHECK_RANGE(ue_golomb(&st), dst.num_ref_idx_lx_active_minus1[0], 14, st);
		if (dst.slice_type == 0) {
			READ_CHECK_RANGE(ue_golomb(&st), dst.num_ref_idx_lx_active_minus1[1], 14, st);
		}
	} else {
		memcpy(dst.num_ref_idx_lx_active_minus1, pps.num_ref_idx_lx_default_active_minus1, sizeof(dst.num_ref_idx_lx_active_minus1));
	}
	if (pps.lists_modification_present_flag && (1 < dst.short_term_ref_pic_set.total_curr)) {
		assert(0);
	}
	init_ref_pic_list(dst.ref_list, dst.short_term_ref_pic_set, dst, pic_order_cnt(dst.slice_pic_order_cnt, sps));
	if (dst.slice_type == 0) {
		dst.mvd_l1_zero_flag = get_onebit(&st);
	}
	if (pps.cabac_init_present_flag) {
		dst.cabac_init_flag = get_onebit(&st);
	}
	if (dst.slice_temporal_mvp_enabled_flag) {
		uint32_t col_l0_flag = (dst.slice_type == 0) ? get_onebit(&st) : 1;
		dst.colocated_from_l0_flag = col_l0_flag;
		if (col_l0_flag && (0 < dst.num_ref_idx_lx_active_minus1[0])) {
			READ_CHECK_RANGE(ue_golomb(&st), dst.collocated_ref_idx, dst.num_ref_idx_lx_active_minus1[0], st);
		} else if (!col_l0_flag && (0 < dst.num_ref_idx_lx_active_minus1[1])) {
			READ_CHECK_RANGE(ue_golomb(&st), dst.collocated_ref_idx, dst.num_ref_idx_lx_active_minus1[1], st);
		}
	}
	if (((dst.slice_type == 0) && pps.weighted_bipred_flag) || ((dst.slice_type == 1) && pps.weighted_pred_flag)) {
		pred_weight_table(dst, pps, sps, st);
	}
	uint32_t five_minus_max_num_merge_cand;
	READ_CHECK_RANGE(ue_golomb(&st), five_minus_max_num_merge_cand, 4, st);
	dst.max_num_merge_cand = 5 - five_minus_max_num_merge_cand;
}

static void slice_header_body(h265d_slice_header_body_t& dst, const h265d_pps_t& pps, const h265d_sps_t& sps, dec_bits& st) {
	if (pps.num_extra_slice_header_bits) {
		skip_bits(&st, pps.num_extra_slice_header_bits);
	}
	READ_CHECK_RANGE(ue_golomb(&st), dst.slice_type, 2, st);
	dst.pic_output_flag = (pps.output_flag_present_flag) ? get_onebit(&st) : 1;
	if (sps.separate_colour_plane_flag) {
		dst.colour_plane_id = get_bits(&st, 2);
	}
	if (dst.nal_type != IDR_W_RADL && dst.nal_type != IDR_N_LP) {
		slice_header_nonidr(dst, pps, sps, st);
	} else {
		init_pic_order_cnt(dst.slice_pic_order_cnt);
	}
	if (sps.sample_adaptive_offset_enabled_flag) {
		dst.slice_sao_luma_flag = get_onebit(&st);
		dst.slice_sao_chroma_flag = get_onebit(&st);
	}
	if (dst.slice_type != 2) {
		slice_header_nonintra(dst, pps, sps, st);
	}
	int32_t qp_delta = se_golomb(&st);
	dst.slice_qpy = pps.init_qp_minus26 + qp_delta + 26;
	CHECK_RANGE2(dst.slice_qpy, -sps.bit_depth_luma_minus8 * 6, 51, st);
	int32_t cb_qp_offset = 0;
	int32_t cr_qp_offset = 0;
	if (pps.pps_slice_chroma_qp_offsets_present_flag) {
		READ_CHECK_RANGE2(se_golomb(&st), cb_qp_offset, -12, 12, st);
		READ_CHECK_RANGE2(se_golomb(&st), cr_qp_offset, -12, 12, st);
	}
	cb_qp_offset += pps.pps_cb_qp_offset;
	dst.slice_qpc_delta[0] = cb_qp_offset;
	CHECK_RANGE2(cb_qp_offset, -12, 12, st);
	cr_qp_offset += pps.pps_cr_qp_offset;
	dst.slice_qpc_delta[1] = cr_qp_offset;
	CHECK_RANGE2(cr_qp_offset, -12, 12, st);
	dst.deblocking_filter_disabled_flag = pps.pps_deblocking_filter_disabled_flag;
	dst.deblocking_filter_override_flag = pps.deblocking_filter_override_enabled_flag ? get_onebit(&st) : 0;
	if (dst.deblocking_filter_override_flag) {
		if ((dst.deblocking_filter_disabled_flag = get_onebit(&st)) == 0) {
			READ_CHECK_RANGE2(se_golomb(&st), dst.slice_beta_offset_div2, -6, 6, st);
			READ_CHECK_RANGE2(se_golomb(&st), dst.slice_tc_offset_div2, -6, 6, st);
		}
	}
	if (pps.pps_loop_filter_across_slices_enabled_flag && (dst.slice_sao_luma_flag || dst.slice_sao_chroma_flag || !dst.deblocking_filter_disabled_flag)) {
		dst.slice_loop_filter_across_slices_enabled_flag = get_onebit(&st);
	} else {
		dst.slice_loop_filter_across_slices_enabled_flag = pps.pps_loop_filter_across_slices_enabled_flag;
	}
}

static void slice_header(h265d_slice_header_t& dst, const h265d_pps_t& pps, const h265d_sps_t& sps, dec_bits& st) {
	dst.dependent_slice_segment_flag = 0;
	if (!dst.first_slice_segment_in_pic_flag) {
		if (pps.dependent_slice_segments_enabled_flag) {
			dst.dependent_slice_segment_flag = get_onebit(&st);
		}
		READ_CHECK_RANGE(get_bits(&st, sps.ctb_info.num_ctb_log2), dst.slice_segment_address, (uint32_t)(sps.ctb_info.columns * sps.ctb_info.rows - 1), st);
	} else {
		dst.slice_segment_address = 0;
	}
	if (!dst.dependent_slice_segment_flag) {
		slice_header_body(dst.body, pps, sps, st);
	}
	if (pps.tiles_enabled_flag || pps.entropy_coding_sync_enabled_flag) {
		entry_points(dst.entry_points, pps, sps, st);
	}
	if (pps.slice_segment_header_extension_present_flag) {
		uint32_t ext_len = ue_golomb(&st);
		dst.slice_segment_header_extension_length = ext_len;
		while (ext_len--) {
			get_bits(&st, 8);
		}
	}
	int to_be_skipped = not_aligned_bits(&st);
	skip_bits(&st, to_be_skipped ? to_be_skipped : 8);
}

/** Constants for initialization of CABAC context.
    m, n are to be derived by:
        m = (elem >> 4) * 5 - 45
        n = ((elem & 15) << 3) - 16
 */
static const m2d_cabac_init_mn_t cabac_initial_value[3][157] = {
	{
		{0, 56}, {15, 48}, {-5, 72}, {-5, 88}, {0, 88}, {0, 64}, {0, 0}, {0, 0},
		{0, 0}, {0, 0}, {10, 48}, {0, 0}, {0, 0}, {0, 0}, {10, 48}, {-30, 104},
		{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
		{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 56}, {-5, 64},
		{-5, 64}, {-15, 104}, {-5, 88}, {-20, 96}, {-5, 64}, {10, 32}, {0, 64}, {0, 0},
		{0, 0}, {0, 64}, {0, 64}, {-5, 72}, {-5, 72}, {-15, 96}, {-15, 96}, {-10, 80},
		{-10, 88}, {-5, 80}, {0, 56}, {-10, 88}, {-10, 104}, {-5, 80}, {-15, 88}, {-15, 104},
		{-5, 104}, {-10, 104}, {-15, 104}, {-25, 104}, {-15, 80}, {-10, 72}, {-30, 104}, {-15, 96},
		{-15, 96}, {-10, 80}, {-10, 88}, {-5, 80}, {0, 56}, {-10, 88}, {-10, 104}, {-5, 80},
		{-15, 88}, {-15, 104}, {-5, 104}, {-10, 104}, {-15, 104}, {-25, 104}, {-15, 80}, {-10, 72},
		{-30, 104}, {-20, 72}, {5, 72}, {-5, 32}, {-5, 88}, {-15, 104}, {-15, 104}, {-10, 88},
		{-15, 96}, {-15, 96}, {-20, 96}, {-10, 80}, {-15, 80}, {-10, 80}, {-15, 72}, {-10, 88},
		{-5, 88}, {10, 8}, {0, 56}, {-10, 88}, {-15, 72}, {-10, 88}, {-5, 88}, {10, 8},
		{0, 56}, {-10, 88}, {-15, 72}, {-10, 88}, {-5, 88}, {10, 8}, {0, 56}, {-10, 88},
		{-5, 80}, {-5, 72}, {10, 32}, {10, 32}, {0, 48}, {-5, 48}, {0, 48}, {-5, 48},
		{0, 56}, {-5, 48}, {-5, 72}, {-15, 104}, {-5, 48}, {-5, 72}, {-15, 104}, {-5, 80},
		{-20, 80}, {-5, 56}, {-5, 64}, {-5, 80}, {0, 48}, {-5, 64}, {-5, 72}, {0, 56},
		{-25, 64}, {0, 24}, {-20, 80}, {-5, 72}, {-15, 72}, {-10, 64}, {0, 48}, {-5, 80},
		{10, 8}, {5, 32}, {10, 32}, {-5, 80}, {25, 8}, {-10, 64}, {15, 24}, {-5, 64},
		{0, 56}, {-5, 48}, {5, 40}, {0, 48}, {0, 48}

	},
	{
		{0, 56}, {10, 56}, {-15, 72}, {-5, 72}, {-10, 96}, {0, 64}, {15, 24}, {10, 56},
		{15, 56}, {0, 24}, {0, 64}, {-5, 72}, {0, 64}, {0, 64}, {0, 64}, {0, 48},
		{-25, 104}, {-15, 96}, {-10, 64}, {-20, 104}, {-25, 104}, {-30, 104}, {-40, 104}, {-40, 104},
		{0, 56}, {0, 56}, {0, 56}, {0, 56}, {5, 48}, {5, 48}, {-10, 80}, {-5, 64},
		{-20, 96}, {0, 56}, {-15, 104}, {0, 24}, {-15, 72}, {5, 40}, {0, 64}, {-5, 80},
		{15, 32}, {0, 64}, {0, 64}, {-5, 72}, {-5, 72}, {-10, 88}, {-15, 96}, {-20, 96},
		{-15, 96}, {-20, 104}, {-25, 104}, {-10, 88}, {-15, 104}, {-15, 96}, {-25, 96}, {-15, 96},
		{-15, 104}, {-15, 104}, {-20, 104}, {-20, 96}, {-15, 80}, {-10, 72}, {-15, 80}, {-10, 88},
		{-15, 96}, {-20, 96}, {-15, 96}, {-20, 104}, {-25, 104}, {-10, 88}, {-15, 104}, {-15, 96},
		{-25, 96}, {-15, 96}, {-15, 104}, {-15, 104}, {-20, 104}, {-20, 96}, {-15, 80}, {-10, 72},
		{-15, 80}, {-10, 56}, {-5, 80}, {-30, 88}, {0, 64}, {0, 72}, {0, 64}, {-5, 72},
		{0, 56}, {-5, 72}, {-10, 72}, {-10, 72}, {-30, 104}, {0, 56}, {5, 32}, {10, 40},
		{-5, 80}, {-5, 48}, {0, 56}, {0, 64}, {5, 32}, {10, 40}, {-5, 80}, {-5, 48},
		{0, 56}, {0, 64}, {5, 32}, {10, 40}, {-5, 80}, {-5, 48}, {0, 56}, {0, 64},
		{5, 64}, {0, 56}, {-10, 72}, {-10, 72}, {-15, 72}, {-10, 56}, {-15, 72}, {-10, 56},
		{5, 40}, {0, 40}, {10, 40}, {-5, 80}, {0, 40}, {10, 40}, {-5, 80}, {0, 64},
		{15, 16}, {15, 16}, {5, 40}, {0, 64}, {0, 48}, {5, 40}, {10, 32}, {10, 32},
		{-5, 32}, {0, 24}, {-5, 48}, {0, 56}, {-10, 56}, {-5, 48}, {-5, 56}, {5, 56},
		{15, 0}, {5, 32}, {5, 40}, {0, 64}, {5, 40}, {-5, 56}, {10, 32}, {-15, 72},
		{5, 40}, {-20, 72}, {-10, 64}, {-15, 72}, {5, 40}

	},
	{
		{0, 56}, {5, -16}, {-15, 72}, {-5, 72}, {-10, 96}, {0, 64}, {15, 24}, {10, 56},
		{15, 56}, {-5, 32}, {0, 64}, {-5, 72}, {0, 64}, {0, 64}, {10, 40}, {0, 48},
		{-25, 104}, {0, 64}, {-5, 56}, {-20, 104}, {-25, 104}, {-30, 104}, {-40, 104}, {-40, 104},
		{0, 56}, {0, 56}, {0, 56}, {0, 56}, {5, 48}, {5, 48}, {25, -16}, {5, 40},
		{-10, 64}, {0, 56}, {-15, 104}, {0, 24}, {-20, 80}, {5, 40}, {0, 64}, {5, 56},
		{15, 32}, {0, 64}, {0, 64}, {-5, 72}, {-5, 72}, {-10, 88}, {-15, 96}, {-10, 80},
		{-15, 96}, {-20, 104}, {-20, 96}, {-10, 88}, {-15, 104}, {-15, 104}, {-25, 104}, {-10, 88},
		{-10, 96}, {-15, 104}, {-15, 104}, {-25, 104}, {-15, 80}, {-10, 72}, {-20, 88}, {-10, 88},
		{-15, 96}, {-10, 80}, {-15, 96}, {-20, 104}, {-20, 96}, {-10, 88}, {-15, 104}, {-15, 104},
		{-25, 104}, {-10, 88}, {-10, 96}, {-15, 104}, {-15, 104}, {-25, 104}, {-15, 80}, {-10, 72},
		{-20, 88}, {-10, 56}, {-5, 80}, {-30, 88}, {0, 64}, {5, 64}, {0, 64}, {-5, 72},
		{0, 56}, {-5, 72}, {-10, 72}, {-10, 72}, {-30, 104}, {-10, 80}, {5, 32}, {10, 40},
		{-5, 80}, {-5, 48}, {0, 56}, {0, 64}, {5, 32}, {10, 40}, {-5, 80}, {-5, 48},
		{0, 56}, {0, 64}, {5, 32}, {10, 40}, {-5, 80}, {-5, 48}, {0, 56}, {0, 64},
		{5, 64}, {0, 56}, {-5, 64}, {-5, 64}, {-10, 64}, {-10, 56}, {-10, 64}, {-10, 56},
		{5, 40}, {0, 40}, {10, 40}, {-5, 80}, {0, 40}, {10, 40}, {-5, 80}, {0, 64},
		{15, 16}, {5, 40}, {5, 40}, {0, 64}, {0, 48}, {5, 40}, {10, 32}, {10, 32},
		{-5, 32}, {0, 24}, {-5, 48}, {0, 56}, {-10, 56}, {-5, 48}, {-10, 64}, {5, 56},
		{20, -16}, {5, 32}, {5, 40}, {0, 64}, {0, 48}, {5, 40}, {10, 32}, {-15, 72},
		{5, 40}, {-20, 72}, {-15, 72}, {-15, 72}, {5, 40}

	}
};

static inline uint32_t sao_merge_flag(m2d_cabac_t& cabac, dec_bits& st) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->sao_merge_flag);
}

static inline uint32_t sao_type_idx(m2d_cabac_t& cabac, dec_bits& st) {
	if (!cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->sao_type_idx)) {
		return 0;
	} else {
		return 1 + cabac_decode_bypass(&cabac, &st);
	}
}

static inline uint32_t sao_offset_abs(m2d_cabac_t& cabac, dec_bits& st, uint32_t max_bits) {
	uint32_t bits = max_bits;
	do {
		if (cabac_decode_bypass(&cabac, &st) == 0) {
			break;
		}
	} while (--bits);
	return max_bits - bits;
}

static inline void sao_read_offset(h265d_sao_map_elem_t& dst, uint32_t idx, uint32_t max_bits, m2d_cabac_t& cabac, dec_bits& st) {
	for (int j = 0; j < 4; ++j) {
		dst.offset[j] = sao_offset_abs(cabac, st, max_bits);
	}
}

static inline void sao_read_offset_band(h265d_sao_map_elem_t& map, m2d_cabac_t& cabac, dec_bits& st) {
	int8_t* offset = map.offset;
	for (int j = 0; j < 4; ++j) {
		int ofs = offset[j];
		if (ofs && cabac_decode_bypass(&cabac, &st)) {
			offset[j] = -ofs;
		}
	}
	map.opt.band_pos = cabac_decode_multibypass(&cabac, &st, 5);
}

static inline void sao_eo_fix_offset(h265d_sao_map_elem_t& dst) {
	dst.offset[2] = -dst.offset[2];
	dst.offset[3] = -dst.offset[3];
}

static inline void sao_eo_class(h265d_sao_map_elem_t& dst, m2d_cabac_t& cabac, dec_bits& st) {
	dst.opt.edge = cabac_decode_multibypass(&cabac, &st, 2);
	sao_eo_fix_offset(dst);
}

static void sao_read_block(h265d_sao_map_t& sao_map, h265d_ctu_t& dst, const h265d_slice_header_t& hdr, dec_bits& st) {
	sao_map.luma_idx = 0;
	if (hdr.body.slice_sao_luma_flag) {
		uint32_t idx = sao_type_idx(dst.cabac, st);
		if (idx != 0) {
			uint32_t max_bits = (1 << (dst.sps->bit_depth_luma_minus8 + 3)) - 1;
			sao_map.luma_idx = idx;
			sao_read_offset(sao_map.elem[0], idx, max_bits, dst.cabac, st);
			if (idx == 1) {
				sao_read_offset_band(sao_map.elem[0], dst.cabac, st);
			} else {
				sao_eo_class(sao_map.elem[0], dst.cabac, st);
			}
		}
	}
	sao_map.chroma_idx = 0;
	if (hdr.body.slice_sao_chroma_flag) {
		uint32_t idx = sao_type_idx(dst.cabac, st);
		if (idx != 0) {
			uint32_t max_bits = (1 << (dst.sps->bit_depth_chroma_minus8 + 3)) - 1;
			sao_map.chroma_idx = idx;
			sao_read_offset(sao_map.elem[1], idx, max_bits, dst.cabac, st);
			if (idx == 1) {
				sao_read_offset_band(sao_map.elem[1], dst.cabac, st);
			} else {
				sao_eo_class(sao_map.elem[1], dst.cabac, st);
			}
			sao_read_offset(sao_map.elem[2], idx, max_bits, dst.cabac, st);
			if (idx == 1) {
				sao_read_offset_band(sao_map.elem[2], dst.cabac, st);
			} else {
				sao_map.elem[2].opt.edge = sao_map.elem[1].opt.edge;
				sao_eo_fix_offset(sao_map.elem[2]);
			}
		}
	}
}

static int sao_search_nonmerged_left(const h265d_sao_map_t* sao_map, int max) {
	int res = max;
	do {
		if (!sao_map->merge_left) {
			break;
		}
		sao_map--;
	} while (--res);
	return res - max;
}

static void sao_read(h265d_ctu_t& ctu, const h265d_slice_header_t& hdr, dec_bits& st) {
	h265d_sao_map_t* sao_map = ctu.sao_map + ctu.pos_y * ctu.sps->ctb_info.columns + ctu.pos_x;
	if (ctu.pos_x != 0) {
		if ((sao_map[0].merge_left = sao_merge_flag(ctu.cabac, st)) != 0) {
			return;
		}
	}
	if (ctu.pos_y != 0) {
		if (sao_merge_flag(ctu.cabac, st)) {
			const h265d_sao_map_t* upper = &sao_map[-ctu.sps->ctb_info.columns];
			sao_map[0] = upper[sao_search_nonmerged_left(upper, ctu.pos_x)];
			return;
		}
	}
	sao_read_block(*sao_map, ctu, hdr, st);
}

static void sao_ignore(h265d_ctu_t& dst, const h265d_slice_header_t& hdr, dec_bits& st) {}

static inline uint32_t split_cu_flag(m2d_cabac_t& cabac, dec_bits& st, uint32_t cond) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->split_cu_flag + cond);
}

static inline uint32_t cu_skip_flag(m2d_cabac_t& cabac, dec_bits& st, uint32_t unavail, const h265d_neighbour_t* left, const h265d_neighbour_t* top) {
	int idx = (!(unavail & 1) && (left->pred_mode & 1)) + (!(unavail & 2) && (top->pred_mode & 1));
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->cu_skip_flag + idx);
}

static inline uint32_t merge_idx(m2d_cabac_t& cabac, dec_bits& st) {
	return cabac_decode_multibypass(&cabac, &st, 3);
}

static inline uint32_t merge_flag(m2d_cabac_t& cabac, dec_bits& st) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->merge_flag);
}

static inline int pred_mode_flag(m2d_cabac_t& cabac, dec_bits& st) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->pred_mode_flag);
}

static inline uint32_t part_mode_inter0(m2d_cabac_t& cabac, dec_bits& st, int8_t* ctx) {
	return cabac_decode_decision_raw(&cabac, &st, ctx) ? 0 : 2 - cabac_decode_decision_raw(&cabac, &st, ctx + 1);
}

static inline uint32_t part_mode_inter1(m2d_cabac_t& cabac, dec_bits& st, int8_t* ctx) {
	uint32_t base = part_mode_inter0(cabac, st, ctx);
	if (base == 0) {
		return base;
	} else {
		if (cabac_decode_decision_raw(&cabac, &st, ctx + 3)) {
			return base;
		} else {
			return (base + 1) * 2 + cabac_decode_bypass(&cabac, &st);
		}
	}
}

static inline uint32_t part_mode_inter2(m2d_cabac_t& cabac, dec_bits& st, int8_t* ctx) {
	uint32_t base = part_mode_inter0(cabac, st, ctx);
	if (base < 2) {
		return base;
	} else {
		return base + (cabac_decode_decision_raw(&cabac, &st, ctx + 2) ^ 1);
	}
}

static inline h265d_inter_part_mode_t part_mode_inter(m2d_cabac_t& cabac, dec_bits& st, int size_log2, int min_size_log2, int amp_enabled_flag) {
	int8_t* ctx = reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->part_mode;
	uint32_t mode;
	if (min_size_log2 < size_log2) {
		if (!amp_enabled_flag) {
			mode = part_mode_inter0(cabac, st, ctx);
		} else {
			mode = part_mode_inter1(cabac, st, ctx);
		}
	} else {
		if (size_log2 == 3) {
			mode = part_mode_inter0(cabac, st, ctx);
		} else {
			mode = part_mode_inter2(cabac, st, ctx);
		}
	}
	return static_cast<h265d_inter_part_mode_t>(mode);
}

static inline uint32_t inter_pred_idc(m2d_cabac_t& cabac, dec_bits& st, int width, int height, int depth) {
	if ((width + height != 12) && cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->inter_pred_idc + depth)) {
		return 2;
	} else {
		return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->inter_pred_idc + 4);
	}
}

static inline int ref_idx_lx(m2d_cabac_t& cabac, dec_bits& st, int max, int8_t* ctx) {
	int idx;
	int max_tmp = std::min(max, 2);
	for (idx = 0; idx < max_tmp; ++idx) {
		if (cabac_decode_decision_raw(&cabac, &st, ctx + idx)) {
			break;
		}
	}
	for (;idx < max; ++idx) {
		if (cabac_decode_bypass(&cabac, &st)) {
			break;
		}
	}
	return idx;
}

static inline int ref_idx_lx(m2d_cabac_t& cabac, dec_bits& st, int lx, const uint8_t num_ref_idx_minus1[]) {
	int num = num_ref_idx_minus1[lx];
	return (0 < num) ? ref_idx_lx(cabac, st, num, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->ref_idx_lx[lx]) : 0;
}

static inline uint32_t abs_mvd_greater_flag(m2d_cabac_t& cabac, dec_bits& st, int idx) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->abs_mvd_greater_flag + idx);
}

static inline uint32_t abs_mvd_minus2(m2d_cabac_t& cabac, dec_bits& st) {
	int bits;
	for (bits = 0; cabac_decode_bypass(&cabac, &st); bits++);
	if (bits) {
		return (2 << bits) - 2 + cabac_decode_multibypass(&cabac, &st, bits + 1);
	} else {
		return 0;
	}
}

static inline uint32_t mvd_sign_flag(m2d_cabac_t& cabac, dec_bits& st) {
	return cabac_decode_bypass(&cabac, &st);
}

static inline uint32_t mvp_lx_flag(m2d_cabac_t& cabac, dec_bits& st, int idx) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->mvp_flag + idx);
}

static inline uint32_t rqt_root_cbf(m2d_cabac_t& cabac, dec_bits& st) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->rqt_root_cbf);
}

static inline uint32_t part_mode_intra(m2d_cabac_t& cabac, dec_bits& st) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->part_mode);
}

static inline uint32_t prev_intra_luma_pred_flag(m2d_cabac_t& cabac, dec_bits& st) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->prev_intra_luma_pred_flag);
}

static inline uint32_t mpm_idx(m2d_cabac_t& cabac, dec_bits& st) {
	return cabac_decode_bypass(&cabac, &st) ? 1 + cabac_decode_bypass(&cabac, &st) : 0;
}

static inline uint32_t rem_intra_luma_pred_mode(m2d_cabac_t& cabac, dec_bits& st, uint8_t candidate[]) {
	uint32_t mode = cabac_decode_multibypass(&cabac, &st, 5);
	std::sort(candidate, candidate + 3);
	for (int i = 0; i < 3; ++i) {
		mode += (candidate[i] <= mode);
	}
	return mode;
}

static inline uint32_t split_transform_flag(m2d_cabac_t& cabac, dec_bits& st, uint32_t log2_transform_size) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->split_transform_flag + 5 - log2_transform_size);
}

static inline uint32_t intra_chroma_pred_mode(m2d_cabac_t& cabac, dec_bits& st) {
	if (cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->intra_chroma_pred_mode)) {
		return cabac_decode_multibypass(&cabac, &st, 2);
	} else {
		return 4;
	}
}

static inline uint32_t cbf_chroma(m2d_cabac_t& cabac, dec_bits& st, uint32_t depth) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->cbf_chroma + depth);
}

static inline uint32_t cbf_luma(m2d_cabac_t& cabac, dec_bits& st, uint32_t depth) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->cbf_luma + (depth == 0));
}

static inline int32_t cu_qp_delta(m2d_cabac_t& cabac, dec_bits& st) {
	assert(0);
	return 0;
}

static inline int32_t transform_skip_flag(m2d_cabac_t& cabac, dec_bits& st, int colour) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->transform_skip_flag + ((colour + 1) >> 1));
}

static inline int32_t last_sig_coeff_prefix_luma(m2d_cabac_t& cabac, dec_bits& st, int8_t* ctx, uint32_t shift, uint32_t max) {
	uint32_t idx;
	for (idx = 0; idx < max; ++idx) {
		if (cabac_decode_decision_raw(&cabac, &st, ctx + (idx >> shift)) == 0) {
			break;
		}
	}
	return idx;
}

static inline int32_t last_sig_coeff_suffix_add(m2d_cabac_t& cabac, dec_bits& st, uint32_t prefix) {
	if (prefix < 4) {
		return prefix;
	} else {
		static const int8_t prefix_adj[] = {
			0x04, 0x06, 0x08, 0x0c, 0x10, 0x18
		};
		return prefix_adj[prefix - 4] + cabac_decode_multibypass(&cabac, &st, (prefix >> 1) - 1);
	}
}

static inline uint32_t coded_sub_block_flag(m2d_cabac_t& cabac, dec_bits& st, uint32_t prev_sbf, uint32_t colour) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->coded_sub_block_flag + ((prev_sbf & 1) | (prev_sbf >> 1)) + ((colour + 1) & 2));
}

static inline uint32_t sig_coeff_flag(m2d_cabac_t& cabac, dec_bits& st, uint8_t inc) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->sig_coeff_flag + inc);
}

static inline uint32_t coeff_abs_level_greater1_flag(m2d_cabac_t& cabac, dec_bits& st, uint8_t inc) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->coeff_abs_level_greater1_flag + inc);
}

static inline uint32_t coeff_abs_level_greater2_flag(m2d_cabac_t& cabac, dec_bits& st, uint8_t inc) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->coeff_abs_level_greater2_flag + inc);
}

static inline uint32_t coeff_sign_flags(m2d_cabac_t& cabac, dec_bits& st, uint32_t num_coeff) {
	return cabac_decode_multibypass(&cabac, &st, num_coeff);
}

static inline uint32_t coeff_abs_level_remaining_prefix(m2d_cabac_t& cabac, dec_bits& st, uint8_t rice) {
	uint32_t i = 0;
	do {
		if (cabac_decode_bypass(&cabac, &st) == 0) {
			break;
		}
	} while (++i < 4);
	if (rice != 0) {
		i = (i << rice) + cabac_decode_multibypass(&cabac, &st, rice);
	}
	return i;
}

static inline uint32_t egk_binstring(m2d_cabac_t& cabac, dec_bits& st, uint32_t k) {
	uint32_t quot;
	for (quot = 0; quot < 32; ++quot) {
		if (cabac_decode_bypass(&cabac, &st) == 0) {
			break;
		}
	}
	return (1 << (quot + k)) - (1 << k) + cabac_decode_multibypass(&cabac, &st, quot + k);
}

static inline uint32_t coeff_abs_level_remaining(m2d_cabac_t& cabac, dec_bits& st, uint8_t rice) {
	uint32_t i;
	for (i = 0; i < 20; ++i) {
		if (cabac_decode_bypass(&cabac, &st) == 0) {
			break;
		}
	}
	if (i < 4) {
		return (rice) ? ((i << rice) + cabac_decode_multibypass(&cabac, &st, rice)) : i;
	} else {
		i -= 4;
		return (1 << (i + rice + 1)) + (2 << rice) + cabac_decode_multibypass(&cabac, &st, i + rice + 1);
	}
}

static inline uint32_t end_of_slice_segment_flag(m2d_cabac_t& cabac, dec_bits& st) {
	uint32_t range = cabac.range - 2;
	uint32_t offset = cabac.offset;
	uint32_t bit;
	if (offset < range) {
		if (range < 256) {
			range = range * 2;
			cabac.offset = offset * 2 + get_onebit(&st);
		}
		bit = 0;
	} else {
		bit = 1;
	}
	cabac.range = range;
	return bit;
}

static inline int8_t intra_chroma_pred_dir(int8_t chroma_pred_mode, int8_t mode) {
	switch (chroma_pred_mode) {
	case 0:
		mode = (mode == 0) ? 34 : 0;
		break;
	case 1:
		mode = (mode == 26) ? 34 : 26;
		break;
	case 2:
		mode = (mode == 10) ? 34 : 10;
		break;
	case 3:
		mode = (mode == 1) ? 34 : 1;
		break;
	}
	return mode;
}

static void intra_pred_candidate(uint8_t cand[], uint8_t candA, uint8_t candB) {
	if (candA == candB) {
		if (candA <= INTRA_DC) {
			cand[0] = INTRA_PLANAR;
			cand[1] = INTRA_DC;
			cand[2] = INTRA_ANGULAR26;
		} else {
			cand[0] = candA;
			cand[1] = ((candA - 3) & 31) + 2;
			cand[2] = ((candA - 1) & 31) + 2;
		}
	} else {
		cand[0] = candA;
		cand[1] = candB;
		uint8_t candC;
		if ((candA != INTRA_PLANAR) && (candB != INTRA_PLANAR)) {
			candC = INTRA_PLANAR;
		} else if ((candA != INTRA_DC) && (candB != INTRA_DC)) {
			candC = INTRA_DC;
		} else {
			candC = INTRA_ANGULAR26;
		}
		cand[2] = candC;
	}
}

static const int8_t h265d_scan_order2x2diag_inverse[2 * 2] = {
	0x00, 0x02, 0x01, 0x03
};

static const int8_t h265d_scan_order2x2diag[2 * 2] = {
	0x00, 0x02, 0x01, 0x03
};

static const int8_t h265d_scan_order8x8diagonal_subblock[4 * 4] = {
	0x00, 0x08, 0x01, 0x10, 0x09, 0x02, 0x18, 0x11,
	0x0a, 0x03, 0x19, 0x12, 0x0b, 0x1a, 0x13, 0x1b
};

static const int8_t h265d_scan_order2x2vertical[2 * 2] = {
	0x00, 0x02, 0x01, 0x03
};

static const int8_t h265d_scan_order8x8vertical_subblock[4 * 4] = {
	0x00, 0x08, 0x10, 0x18, 0x01, 0x09, 0x11, 0x19,
	0x02, 0x0a, 0x12, 0x1a, 0x03, 0x0b, 0x13, 0x1b
};

static const int8_t h265d_scan_order4x4diag_inverse[4 * 4] = {
	0x00, 0x02, 0x05, 0x09, 0x01, 0x04, 0x08, 0x0c,
	0x03, 0x07, 0x0b, 0x0e, 0x06, 0x0a, 0x0d, 0x0f
};

static const int8_t h265d_scan_order4x4diag[4 * 4] = {
	0x00, 0x04, 0x01, 0x08, 0x05, 0x02, 0x0c, 0x09,
	0x06, 0x03, 0x0d, 0x0a, 0x07, 0x0e, 0x0b, 0x0f
};

static const int8_t h265d_scan_order16x16diagonal_subblock[4 * 4] = {
	0x00, 0x10, 0x01, 0x20, 0x11, 0x02, 0x30, 0x21,
	0x12, 0x03, 0x31, 0x22, 0x13, 0x32, 0x23, 0x33
};

static const int8_t h265d_scan_order4x4vertical[4 * 4] = {
	0x00, 0x04, 0x08, 0x0c, 0x01, 0x05, 0x09, 0x0d,
	0x02, 0x06, 0x0a, 0x0e, 0x03, 0x07, 0x0b, 0x0f
};

static const int8_t h265d_scan_order16x16vertical_subblock[4 * 4] = {
	0x00, 0x10, 0x20, 0x30, 0x01, 0x11, 0x21, 0x31,
	0x02, 0x12, 0x22, 0x32, 0x03, 0x13, 0x23, 0x33
};

static const int8_t h265d_scan_order8x8diag_inverse[8 * 8] = {
	0x00, 0x02, 0x05, 0x09, 0x0e, 0x14, 0x1b, 0x23,
	0x01, 0x04, 0x08, 0x0d, 0x13, 0x1a, 0x22, 0x2a,
	0x03, 0x07, 0x0c, 0x12, 0x19, 0x21, 0x29, 0x30,
	0x06, 0x0b, 0x11, 0x18, 0x20, 0x28, 0x2f, 0x35,
	0x0a, 0x10, 0x17, 0x1f, 0x27, 0x2e, 0x34, 0x39,
	0x0f, 0x16, 0x1e, 0x26, 0x2d, 0x33, 0x38, 0x3c,
	0x15, 0x1d, 0x25, 0x2c, 0x32, 0x37, 0x3b, 0x3e,
	0x1c, 0x24, 0x2b, 0x31, 0x36, 0x3a, 0x3d, 0x3f
};

static const int8_t h265d_scan_order8x8diag[8 * 8] = {
	0x00, 0x08, 0x01, 0x10, 0x09, 0x02, 0x18, 0x11,
	0x0a, 0x03, 0x20, 0x19, 0x12, 0x0b, 0x04, 0x28,
	0x21, 0x1a, 0x13, 0x0c, 0x05, 0x30, 0x29, 0x22,
	0x1b, 0x14, 0x0d, 0x06, 0x38, 0x31, 0x2a, 0x23,
	0x1c, 0x15, 0x0e, 0x07, 0x39, 0x32, 0x2b, 0x24,
	0x1d, 0x16, 0x0f, 0x3a, 0x33, 0x2c, 0x25, 0x1e,
	0x17, 0x3b, 0x34, 0x2d, 0x26, 0x1f, 0x3c, 0x35,
	0x2e, 0x27, 0x3d, 0x36, 0x2f, 0x3e, 0x37, 0x3f
};

static const int8_t h265d_scan_order32x32diagonal_subblock[4 * 4] = {
	0x00, 0x20, 0x01, 0x40, 0x21, 0x02, 0x60, 0x41,
	0x22, 0x03, 0x61, 0x42, 0x23, 0x62, 0x43, 0x63
};

static const int8_t h265d_scan_order8x8vertical[8 * 8] = {
	0x00, 0x08, 0x10, 0x18, 0x20, 0x28, 0x30, 0x38,
	0x01, 0x09, 0x11, 0x19, 0x21, 0x29, 0x31, 0x39,
	0x02, 0x0a, 0x12, 0x1a, 0x22, 0x2a, 0x32, 0x3a,
	0x03, 0x0b, 0x13, 0x1b, 0x23, 0x2b, 0x33, 0x3b,
	0x04, 0x0c, 0x14, 0x1c, 0x24, 0x2c, 0x34, 0x3c,
	0x05, 0x0d, 0x15, 0x1d, 0x25, 0x2d, 0x35, 0x3d,
	0x06, 0x0e, 0x16, 0x1e, 0x26, 0x2e, 0x36, 0x3e,
	0x07, 0x0f, 0x17, 0x1f, 0x27, 0x2f, 0x37, 0x3f
};

static const int8_t h265d_scan_order32x32vertical_subblock[4 * 4] = {
	0x00, 0x20, 0x40, 0x60, 0x01, 0x21, 0x41, 0x61,
	0x02, 0x22, 0x42, 0x62, 0x03, 0x23, 0x43, 0x63
};

static const int8_t h265d_scan_order8x8horizontal[8 * 8] = {
	0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
	0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
	0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
	0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
	0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27,
	0x28, 0x29, 0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f,
	0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
	0x38, 0x39, 0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f
};

static const int8_t h265d_scan_order8x8horizontal_subblock[4 * 4] = {
	0x00, 0x01, 0x02, 0x03, 0x08, 0x09, 0x0a, 0x0b,
	0x10, 0x11, 0x12, 0x13, 0x18, 0x19, 0x1a, 0x1b
};

static const int8_t h265d_scan_order16x16horizontal_subblock[4 * 4] = {
	0x00, 0x01, 0x02, 0x03, 0x10, 0x11, 0x12, 0x13,
	0x20, 0x21, 0x22, 0x23, 0x30, 0x31, 0x32, 0x33
};

static const int8_t h265d_scan_order32x32horizontal_subblock[4 * 4] = {
	0x00, 0x01, 0x02, 0x03, 0x20, 0x21, 0x22, 0x23,
	0x40, 0x41, 0x42, 0x43, 0x60, 0x61, 0x62, 0x63
};

struct residual_scan_order_t {
	const int8_t* sub_block_num;
	const int8_t* sub_block_pos;
	const int8_t* inner_pos_num;
	const int8_t* inner_xy_pos;
	const int8_t* macro_xy_pos;
};

template <int OFFSET>
static uint8_t sig_coeff_flag_inc4(uint32_t sxy, int pxy, uint32_t prev_sbf) {
	static const int8_t inc_map[] = {
		0, 1, 4, 5, 2, 3, 4, 5, 6, 6, 8, 8, 7, 7, 8
	};
	return inc_map[pxy] + OFFSET;
}

static const int8_t sig_coeff_inc_init[64] = {
	2, 2, 2, 2, 1, 2, 1, 2, 1, 2, 0, 2, 0, 2, 0, 2,
	1, 1, 2, 2, 1, 1, 1, 2, 0, 1, 0, 2, 0, 1, 0, 2,
	1, 0, 2, 2, 0, 0, 1, 2, 0, 0, 0, 2, 0, 0, 0, 2,
	0, 0, 2, 2, 0, 0, 1, 2, 0, 0, 0, 2, 0, 0, 0, 2
};

template <int MOD>
static uint8_t sig_coeff_flag_inc_luma(uint32_t sxy, int pxy, uint32_t prev_sbf) {
	uint8_t inc;
	if (sxy == 0) {
		if (pxy == 0) {
			return 0;
		}
		inc = 0;
	} else {
		inc = 3;
	}
	return sig_coeff_inc_init[pxy * 4 + prev_sbf] + inc + MOD;
}

template <int MOD>
static uint8_t sig_coeff_flag_inc_chroma(uint32_t sxy, int pxy, uint32_t prev_sbf) {
	if ((sxy == 0) && (pxy == 0)) {
		return 27;
	}
	return sig_coeff_inc_init[pxy * 4 + prev_sbf] + MOD + 27;
}

static const residual_scan_order_t residual_scan_order[3][4] = {
	{
		{h265d_scan_order2x2diag_inverse, h265d_scan_order2x2diag, h265d_scan_order4x4diag_inverse, h265d_scan_order4x4diag, h265d_scan_order4x4diag},
		{h265d_scan_order2x2diag_inverse, h265d_scan_order2x2diag, h265d_scan_order4x4diag_inverse, h265d_scan_order4x4diag, h265d_scan_order8x8diagonal_subblock},
		{h265d_scan_order4x4diag_inverse, h265d_scan_order4x4diag, h265d_scan_order4x4diag_inverse, h265d_scan_order4x4diag, h265d_scan_order16x16diagonal_subblock},
		{h265d_scan_order8x8diag_inverse, h265d_scan_order8x8diag, h265d_scan_order4x4diag_inverse, h265d_scan_order4x4diag, h265d_scan_order32x32diagonal_subblock}
	},
	{
		{h265d_scan_order8x8horizontal, h265d_scan_order8x8horizontal, h265d_scan_order8x8horizontal, h265d_scan_order8x8horizontal, h265d_scan_order8x8horizontal},
		{h265d_scan_order8x8horizontal, h265d_scan_order8x8horizontal, h265d_scan_order8x8horizontal, h265d_scan_order8x8horizontal, h265d_scan_order8x8horizontal_subblock},
		{h265d_scan_order8x8horizontal, h265d_scan_order8x8horizontal, h265d_scan_order8x8horizontal, h265d_scan_order8x8horizontal, h265d_scan_order16x16horizontal_subblock},
		{h265d_scan_order8x8horizontal, h265d_scan_order8x8horizontal, h265d_scan_order8x8horizontal, h265d_scan_order8x8horizontal, h265d_scan_order32x32horizontal_subblock}
	},
	{
		{h265d_scan_order2x2vertical, h265d_scan_order2x2vertical, h265d_scan_order4x4vertical, h265d_scan_order4x4vertical, h265d_scan_order4x4vertical},
		{h265d_scan_order2x2vertical, h265d_scan_order2x2vertical, h265d_scan_order4x4vertical, h265d_scan_order4x4vertical, h265d_scan_order8x8vertical_subblock},
		{h265d_scan_order4x4vertical, h265d_scan_order4x4vertical, h265d_scan_order4x4vertical, h265d_scan_order4x4vertical, h265d_scan_order16x16vertical_subblock},
		{h265d_scan_order8x8vertical, h265d_scan_order8x8vertical, h265d_scan_order4x4vertical, h265d_scan_order4x4vertical, h265d_scan_order32x32vertical_subblock}
	}
};

typedef uint8_t (*sig_coeff_flag_inc_func_t)(uint32_t sxy, int pxy, uint32_t prev_sbf);

static const sig_coeff_flag_inc_func_t sig_coeff_flag_inc_func[3][4][2] = {
	{
		{sig_coeff_flag_inc4<0>, sig_coeff_flag_inc4<27>},
		{sig_coeff_flag_inc_luma<9>, sig_coeff_flag_inc_chroma<9>},
		{sig_coeff_flag_inc_luma<21>, sig_coeff_flag_inc_chroma<12>},
		{sig_coeff_flag_inc_luma<21>, sig_coeff_flag_inc_chroma<12>}
	},
	{
		{sig_coeff_flag_inc4<0>, sig_coeff_flag_inc4<27>},
		{sig_coeff_flag_inc_luma<15>, sig_coeff_flag_inc_chroma<9>},
		{sig_coeff_flag_inc_luma<21>, sig_coeff_flag_inc_chroma<12>},
		{sig_coeff_flag_inc_luma<21>, sig_coeff_flag_inc_chroma<12>}
	},
	{
		{sig_coeff_flag_inc4<0>, sig_coeff_flag_inc4<27>},
		{sig_coeff_flag_inc_luma<15>, sig_coeff_flag_inc_chroma<9>},
		{sig_coeff_flag_inc_luma<21>, sig_coeff_flag_inc_chroma<12>},
		{sig_coeff_flag_inc_luma<21>, sig_coeff_flag_inc_chroma<12>}
	}
};

static inline uint32_t scan_order_index(uint32_t x, uint32_t y, uint32_t size_log2) {
	return (y << size_log2) + x;
}

typedef struct {
	int8_t pos;
	int8_t val;
} h265d_sigcoeff_t;

#define FILL_SIGCOEFF(coeffs, i, p, v) do {h265d_sigcoeff_t& dst = coeffs[++i]; dst.pos = p; dst.val = v;} while (0)

static inline uint32_t sig_coeff_flags_read(m2d_cabac_t& cabac, dec_bits& st, sig_coeff_flag_inc_func_t sig_coeff_flag_inc, const int8_t pix_pos[], bool is_last, int32_t pos, h265d_sigcoeff_t coeffs[], uint32_t sxy, uint32_t prev_sbf) {
	int32_t idx = -1;
	if (is_last) {
		FILL_SIGCOEFF(coeffs, idx, pos, 1);
		pos--;
	}
	while (0 < pos) {
		if (sig_coeff_flag(cabac, st, sig_coeff_flag_inc(sxy, pix_pos[pos], prev_sbf))) {
			FILL_SIGCOEFF(coeffs, idx, pos, 1);
		}
		pos--;
	}
	if ((pos == 0) && (((idx < 0) && sxy) || sig_coeff_flag(cabac, st, sig_coeff_flag_inc(sxy, 0, prev_sbf)))) {
		FILL_SIGCOEFF(coeffs, idx, 0, 1);
	}
	return static_cast<uint32_t>(idx + 1);
}

static inline uint32_t sig_coeff_greater(int colour, int subblock_idx, uint8_t& greater1ctx, h265d_sigcoeff_t coeff[], uint32_t num_coeff, m2d_cabac_t& cabac, dec_bits& st) {
	uint32_t ctxset = (((colour == 0) && (subblock_idx != 0)) ? 2 : 0) + (greater1ctx == 0);
	uint32_t greater1offset = ctxset * 4 + ((colour == 0) ? 0 : 16);
	greater1ctx = 1;
	uint32_t max_flags = 0;
	int last_greater1_idx = -1;
	uint32_t max_num = std::min(num_coeff, 8U);
	for (uint32_t j = 0; j < max_num; ++j) {
		if (coeff_abs_level_greater1_flag(cabac, st, greater1offset + greater1ctx)) {
			greater1ctx = 0;
			coeff[j].val = 2;
			if (0 <= last_greater1_idx) {
				max_flags |= 1 << j;
			} else {
				last_greater1_idx = j;
			}
		} else if ((uint32_t)(greater1ctx - 1) < 2) {
			greater1ctx++;
		}
	}
	if (0 <= last_greater1_idx) {
		if (coeff_abs_level_greater2_flag(cabac, st, (colour == 0) ? ctxset : ctxset + 4)) {
			coeff[last_greater1_idx].val = 3;
			max_flags |= 1 << last_greater1_idx;
		}
	}
	if (8 < num_coeff) {
		max_flags |= ((1 << (num_coeff - 8)) - 1) << 8;
	}
	return max_flags;
}

static inline uint32_t sig_coeff_writeback(uint32_t num_coeff, h265d_sigcoeff_t coeff[], uint32_t max_coeffs, uint8_t hidden_sign, uint16_t sign_flags, const int8_t xy_pos[], int16_t* dst, uint32_t write_pos, const h265d_scaling_func_t scaling, const h265d_scaling_info_t& scale, m2d_cabac_t& cabac, dec_bits& st) {
	uint32_t rice_param = 0;
	uint32_t sign_mask = 1 << (num_coeff - 1 - hidden_sign);
	uint32_t last_write_pos;
	int32_t level_sum = 0;
	uint32_t xy_pos_sum = 0;
	for (uint32_t i = 0; i < num_coeff; ++i) {
		uint32_t abs_level = coeff[i].val;
		if (max_coeffs & 1) {
			abs_level += coeff_abs_level_remaining(cabac, st, rice_param);
			rice_param = std::min(rice_param + ((static_cast<uint32_t>(3) << rice_param) < abs_level), static_cast<uint32_t>(4));
		}
		level_sum += abs_level;
		last_write_pos = write_pos + xy_pos[coeff[i].pos];
		xy_pos_sum |= last_write_pos;
		uint32_t sign = (sign_flags & sign_mask) != 0;
		dst[last_write_pos] = scaling(NEGATE(abs_level, sign), scale, last_write_pos);
		sign_mask >>= 1;
		max_coeffs >>= 1;
	}
	if (hidden_sign && (level_sum & 1)) {
		dst[last_write_pos] = -dst[last_write_pos];
	}
	return xy_pos_sum;
}

class sub_block_flags_t {
	uint8_t log2size_;
	uint8_t pos_max_;
	uint8_t sx_, sy_;
	uint8_t flags[9];
public:
	sub_block_flags_t(uint32_t log2size) : log2size_(log2size - 2), pos_max_((1 << (log2size - 2)) - 1) {
		memset(flags, 0, sizeof(flags));
	}
	uint32_t prev_flags(uint32_t block_pos) {
		uint32_t sx = block_pos & pos_max_;
		uint32_t sy = block_pos >> log2size_;
		sx_ = sx;
		sy_ = sy;
		return ((flags[sy] >> (sx + 1)) & 1) + ((flags[sy + 1] >> sx) & 1) * 2;
	}
	void set_flag() {
		flags[sy_] |= 1 << sx_;
	}
	uint32_t offset() const {
		return ((sy_ << (log2size_ + 2)) + sx_) * 4;
	}
};

#ifdef __arm__
#define SATURATE16BIT(x) __sat((x), 16)
#else
#define SATURATE16BIT(x) (((x) <= 32767) ? ((-32768 <= (x)) ? (x) : -32768) : 32767)
#endif

template <int LOG2>
static int16_t scaling_default_base(int32_t val, const h265d_scaling_info_t& scale, int idx) {
	return SATURATE16BIT((val * scale.scale + (1 << (LOG2 - 2))) >> (LOG2 - 1));
}

static const h265d_scaling_func_t scaling_default_func[4] = {
	scaling_default_base<2>,
	scaling_default_base<3>,
	scaling_default_base<4>,
	scaling_default_base<5>
};

template<int N, int SHIFT, int GAP, typename T>
static inline void NxNtransform_dconly(uint8_t *dst, int16_t* coeff, int stride) {
	acNxNtransform_dconly<N, SHIFT, GAP - 1, T>(dst, coeff[0], stride);
}

template <int LOG2>
struct sat16 {
	int16_t operator()(int val) const {
		val = (val + (1 << (LOG2 - 1))) >> LOG2;
		return SATURATE16BIT(val);
	}
};

template <int LOG2, int N>
static inline void add_transformed_coeff_line(uint8_t* dst, const int16_t* coeff) {
	for (int x = 0; x < (1 << LOG2); ++x) {
		int v = dst[x * N] + coeff[x];
		dst[x * N] = CLIP255C(v);
	}
}

template <int LOG2>
static inline void transform_horiz_pretruncate(int16_t *dst) {
	for (int i = 0; i < (1 << LOG2); ++i) {
		dst[i] = (dst[i] + 1) >> 1;
	}
}

template <int LOG2>
static inline void transform_vert_pretruncate(int16_t *dst) {
	for (int i = 0; i < (1 << LOG2); ++i) {
		dst[i << LOG2] = (dst[i << LOG2] + 1) >> 1;
	}
}

template <typename T, typename F>
static inline void transformdst_dc_line(int dc, T& d0, T& d1, T& d2, T& d3, F Saturate) {
	d0 = Saturate(dc * 29);
	d1 = Saturate(dc * 55);
	d2 = Saturate(dc * 74);
	d3 = Saturate(dc * (55 + 29));
}

static inline void transformdst_dconly(uint8_t *dst, int16_t* coeff, int stride) {
	int d0, d1, d2, d3;
	int e0, e1, e2, e3;
	transformdst_dc_line(coeff[0], d0, d1, d2, d3, sat16<7>());
	transformdst_dc_line(d0, e0, e1, e2, e3, sat16<12>());
	dst[0] = CLIP255C(dst[0] + e0);
	dst[1] = CLIP255C(dst[1] + e1);
	dst[2] = CLIP255C(dst[2] + e2);
	dst[3] = CLIP255C(dst[3] + e3);
	dst += stride;
	transformdst_dc_line(d1, e0, e1, e2, e3, sat16<12>());
	dst[0] = CLIP255C(dst[0] + e0);
	dst[1] = CLIP255C(dst[1] + e1);
	dst[2] = CLIP255C(dst[2] + e2);
	dst[3] = CLIP255C(dst[3] + e3);
	dst += stride;
	transformdst_dc_line(d2, e0, e1, e2, e3, sat16<12>());
	dst[0] = CLIP255C(dst[0] + e0);
	dst[1] = CLIP255C(dst[1] + e1);
	dst[2] = CLIP255C(dst[2] + e2);
	dst[3] = CLIP255C(dst[3] + e3);
	dst += stride;
	transformdst_dc_line(d3, e0, e1, e2, e3, sat16<12>());
	dst[0] = CLIP255C(dst[0] + e0);
	dst[1] = CLIP255C(dst[1] + e1);
	dst[2] = CLIP255C(dst[2] + e2);
	dst[3] = CLIP255C(dst[3] + e3);
}

#ifdef X86ASM
void transformdst_ac4x4(uint8_t* dst, int16_t* coeff, int stride);
void transform_ac4x4(uint8_t* dst, int16_t* coeff, int stride);
void transform_ac8x8(uint8_t* dst, int16_t* src, int stride);
void transform_ac16x16(uint8_t* dst, int16_t* src, int stride);
void transform_ac32x32(uint8_t* dst, int16_t* src, int stride);
void transform_ac4x4chroma(uint8_t* dst, int16_t* coeff, int stride);
void transform_ac8x8chroma(uint8_t* dst, int16_t* coeff, int stride);
void transform_ac16x16chroma(uint8_t* dst, int16_t* src, int stride);
void transform_ac32x32chroma(uint8_t* dst, int16_t* src, int stride);
#else
template <typename F>
static inline void transformdst_line4(int16_t *dst, const int16_t* coeff, F Saturate) {
	int c0 = coeff[0];
	int c1 = coeff[4];
	int c2 = coeff[8];
	int c3 = coeff[12];
	int d0 = c0 + c2;
	int d1 = c2 + c3;
	int d2 = c0 - c3;
	int d3 = c1 * 74;
	dst[0] = Saturate(d0 * 29 + d1 * 55 + d3);
	dst[1] = Saturate(d2 * 55 - d1 * 29 + d3);
	dst[2] = Saturate((c0 - c2 + c3) * 74);
	dst[3] = Saturate(d0 * 55 + d2 * 29 - d3);
}

static inline void transformdst_ac4x4(uint8_t* dst, int16_t* coeff, int stride) {
	int16_t* tmp = coeff + 16;
	for (int x = 0; x < 4; ++x) {
		transformdst_line4(tmp + x * 4, coeff + x, sat16<7>());
	}
	for (int y = 0; y < 4; ++y) {
		transformdst_line4(coeff, tmp, sat16<12>());
		add_transformed_coeff_line<2, 1>(dst, coeff);
		dst += stride;
		tmp++;
	}
}
#endif

template <int LOG2, typename T0, typename F>
static inline void transform_line4(T0 *dst, const int16_t* coeff, F Saturate) {
	int c0 = coeff[0];
	int c1 = coeff[1 << LOG2];
	int c2 = coeff[2 << LOG2];
	int c3 = coeff[3 << LOG2];
	int odd0 = c1 * 83 + c3 * 36;
	int even0 = (c0 + c2) * 64;
	int odd1 = c1 * 36 - c3 * 83;
	int even1 = (c0 - c2) * 64;
	dst[0] = Saturate(even0 + odd0);
	dst[1] = Saturate(even1 + odd1);
	dst[2] = Saturate(even1 - odd1);
	dst[3] = Saturate(even0 - odd0);
}

struct nosat16 {
	int operator()(int val) const {
		return val;
	}
};

template <int LOG2, typename T0, typename F>
static inline void transform_line8(T0 *dst, const int16_t* coeff, F Saturate) {
	int32_t tmp[4];
	transform_line4<(LOG2 + 1)>(tmp, coeff, nosat16());
	int c1 = coeff[1 << LOG2];
	int c3 = coeff[3 << LOG2];
	int c5 = coeff[5 << LOG2];
	int c7 = coeff[7 << LOG2];
	int eo0 = 89 * c1 + 75 * c3 + 50 * c5 + 18 * c7;
	int eo1 = 75 * c1 - 18 * c3 - 89 * c5 - 50 * c7;
	int eo2 = 50 * c1 - 89 * c3 + 18 * c5 + 75 * c7;
	int eo3 = 18 * c1 - 50 * c3 + 75 * c5 - 89 * c7;
	dst[0] = Saturate(tmp[0] + eo0);
	dst[7] = Saturate(tmp[0] - eo0);
	dst[1] = Saturate(tmp[1] + eo1);
	dst[6] = Saturate(tmp[1] - eo1);
	dst[2] = Saturate(tmp[2] + eo2);
	dst[5] = Saturate(tmp[2] - eo2);
	dst[3] = Saturate(tmp[3] + eo3);
	dst[4] = Saturate(tmp[3] - eo3);
}

static const int8_t transform_oddc8[] = {
	90, 87, 80, 70, 57, 43, 25, 9,
	87, 57, 9, -43, -80, -90, -70, -25,
	80, 9, -70, -87, -25, 57, 90, 43,
	70, -43, -87, 9, 90, 25, -80, -57,
	57, -80, -25, 90, -9, -87, 43, 70,
	43, -90, 57, 25, -87, 70, 9, -80,
	25, -70, 90, -80, 43, 9, -57, 87,
	9, -25, 43, -57, 70, -80, 87, -90
};

static const int8_t transform_oddc16[] = {
	90, 90, 88, 85, 82, 78, 73, 67, 61, 54, 46, 38, 31, 22, 13, 4,
	90, 82, 67, 46, 22, -4, -31, -54, -73, -85, -90, -88, -78, -61, -38, -13,
	88, 67, 31, -13, -54, -82, -90, -78, -46, -4, 38, 73, 90, 85, 61, 22,
	85, 46, -13, -67, -90, -73, -22, 38, 82, 88, 54, -4, -61, -90, -78, -31,
	82, 22, -54, -90, -61, 13, 78, 85, 31, -46, -90, -67, 4, 73, 88, 38,
	78, -4, -82, -73, 13, 85, 67, -22, -88, -61, 31, 90, 54, -38, -90, -46,
	73, -31, -90, -22, 78, 67, -38, -90, -13, 82, 61, -46, -88, -4, 85, 54,
	67, -54, -78, 38, 85, -22, -90, 4, 90, 13, -88, -31, 82, 46, -73, -61,
	61, -73, -46, 82, 31, -88, -13, 90, -4, -90, 22, 85, -38, -78, 54, 67,
	54, -85, -4, 88, -46, -61, 82, 13, -90, 38, 67, -78, -22, 90, -31, -73,
	46, -90, 38, 54, -90, 31, 61, -88, 22, 67, -85, 13, 73, -82, 4, 78,
	38, -88, 73, -4, -67, 90, -46, -31, 85, -78, 13, 61, -90, 54, 22, -82,
	31, -78, 90, -61, 4, 54, -88, 82, -38, -22, 73, -90, 67, -13, -46, 85,
	22, -61, 85, -90, 73, -38, -4, 46, -78, 90, -82, 54, -13, -31, 67, -88,
	13, -38, 61, -78, 88, -90, 85, -73, 54, -31, 4, 22, -46, 67, -82, 90,
	4, -13, 22, -31, 38, -46, 54, -61, 67, -73, 78, -82, 85, -88, 90, -90
};

template <int LOG2, typename T0, typename F>
static inline void transform_line16(T0 *dst, const int16_t* coeff, F Saturate) {
	int32_t even[8];
	transform_line8<(LOG2 + 1)>(even, coeff, nosat16());
	const int8_t* c = transform_oddc8;
	int c0 = coeff[1 << LOG2];
	int c1 = coeff[3 << LOG2];
	int c2 = coeff[5 << LOG2];
	int c3 = coeff[7 << LOG2];
	int c4 = coeff[9 << LOG2];
	int c5 = coeff[11 << LOG2];
	int c6 = coeff[13 << LOG2];
	int c7 = coeff[15 << LOG2];
	for (int i = 0; i < 8; ++i) {
		int sum = c0 * *c++;
		sum += c1 * *c++;
		sum += c2 * *c++;
		sum += c3 * *c++;
		sum += c4 * *c++;
		sum += c5 * *c++;
		sum += c6 * *c++;
		sum += c7 * *c++;
		int ev = even[i];
		dst[i] = Saturate(ev + sum);
		dst[15 - i] = Saturate(ev - sum);
	}
}

template <int LOG2, typename F>
static inline void transform_line32(int16_t *dst, const int16_t* coeff, F Saturate) {
	int32_t even[16];
	transform_line16<(LOG2 + 1)>(even, coeff, nosat16());
	const int8_t* c = transform_oddc16;
	const int16_t* co = coeff + (1 << LOG2);
	for (int i = 0; i < 16; ++i) {
		int sum = 0;
		for (int j = 0; j < 16; ++j) {
			sum += co[j << (LOG2 + 1)] * *c++;
		}
		int ev = even[i];
		dst[i] = Saturate(ev + sum);
		dst[31 - i] = Saturate(ev - sum);
	}
}

template <int LOG2, int N>
static inline void transform_horiz(uint8_t* dst, int16_t* coeff, int stride) {
	transform_horiz_pretruncate<LOG2>(coeff);
	int16_t* tmp = coeff + (1 << LOG2);
	if (LOG2 == 2) {
		transform_line4<0>(tmp, coeff, sat16<12>());
	} else if (LOG2 == 3) {
		transform_line8<0>(tmp, coeff, sat16<12>());
	} else if (LOG2 == 4) {
		transform_line16<0>(tmp, coeff, sat16<12>());
	} else {
		transform_line32<0>(tmp, coeff, sat16<12>());
	}
	for (int y = 0; y < (1 << LOG2); ++y) {
		add_transformed_coeff_line<LOG2, N>(dst, tmp);
		dst += stride;
	}
}

template <int LOG2, int N>
static inline void transform_vert(uint8_t* dst, int16_t* coeff, int stride) {
	int16_t* tmp = coeff + (1 << (LOG2 * 2));
	if (LOG2 == 2) {
		transform_line4<2>(tmp, coeff, sat16<7>());
	} else if (LOG2 == 3) {
		transform_line8<3>(tmp, coeff, sat16<7>());
	} else if (LOG2 == 4) {
		transform_line16<4>(tmp, coeff, sat16<7>());
	} else {
		transform_line32<5>(tmp, coeff, sat16<7>());
	}
	for (int y = 0; y < (1 << LOG2); ++y) {
		int diff = (tmp[y] + 32) >> 6;
		for (int x = 0; x < (1 << LOG2); ++x) {
			int v = dst[x * N] + diff;
			dst[x * N] = CLIP255C(v);
		}
		dst += stride;
	}
}

template <int LOG2, int N>
static inline void transform_acNxN(uint8_t* dst, int16_t* coeff, int stride) {
	int16_t* tmp = coeff + (1 << (LOG2 * 2));
	for (int x = 0; x < (1 << LOG2); ++x) {
		if (LOG2 == 2) {
			transform_line4<2>(tmp + x * 4, coeff + x, sat16<7>());
		} else if (LOG2 == 3) {
			transform_line8<3>(tmp + x * 8, coeff + x, sat16<7>());
		} else if (LOG2 == 4) {
			transform_line16<4>(tmp + x * 16, coeff + x, sat16<7>());
		} else {
			transform_line32<5>(tmp + x * 32, coeff + x, sat16<7>());
		}
	}
	for (int y = 0; y < (1 << LOG2); ++y) {
		if (LOG2 == 2) {
			transform_line4<2>(coeff, tmp, sat16<12>());
		} else if (LOG2 == 3) {
			transform_line8<3>(coeff, tmp, sat16<12>());
		} else if (LOG2 == 4) {
			transform_line16<4>(coeff, tmp, sat16<12>());
		} else {
			transform_line32<5>(coeff, tmp, sat16<12>());
		}
		add_transformed_coeff_line<LOG2, N>(dst, coeff);
		dst += stride;
		tmp++;
	}
}

static void (* const transform_func[4][2][4])(uint8_t *dst, int16_t* coeff, int stride) = {
	{
		{
			transformdst_dconly,
			NxNtransform_dconly<8, 7, 1, uint64_t>,
			NxNtransform_dconly<16, 7, 1, uint64_t>,
			NxNtransform_dconly<32, 7, 1, uint64_t>
		},
		{
			NxNtransform_dconly<4, 7, 2, uint64_t>,
			NxNtransform_dconly<8, 7, 2, uint64_t>,
			NxNtransform_dconly<16, 7, 2, uint64_t>,
			NxNtransform_dconly<32, 7, 2, uint64_t>
		}
	},
	{
		{
			transformdst_ac4x4,
			transform_horiz<3, 1>,
			transform_horiz<4, 1>,
			transform_horiz<5, 1>
		},
		{
			transform_horiz<2, 2>,
			transform_horiz<3, 2>,
			transform_horiz<4, 2>,
			transform_horiz<5, 2>
		}
	},
	{
		{
			transformdst_ac4x4,
			transform_vert<3, 1>,
			transform_vert<4, 1>,
			transform_vert<5, 1>
		},
		{
			transform_vert<2, 2>,
			transform_vert<3, 2>,
			transform_vert<4, 2>,
			transform_vert<5, 2>
		}
	},
	{
		{
			transformdst_ac4x4,
#ifdef X86ASM
			transform_ac8x8,
			transform_ac16x16,
			transform_ac32x32
#else
			transform_acNxN<3, 1>,
			transform_acNxN<4, 1>,
			transform_acNxN<5, 1>
#endif
		},
		{
#ifdef X86ASM
			transform_ac4x4chroma,
			transform_ac8x8chroma,
			transform_ac16x16chroma,
			transform_ac32x32chroma
#else
			transform_acNxN<2, 2>,
			transform_acNxN<3, 2>,
			transform_acNxN<4, 2>,
			transform_acNxN<5, 2>
#endif
		}
	}
};

static inline void transform(int16_t coeff[], uint32_t size_log2, uint8_t* frame, uint32_t stride, int colour, uint32_t xy_pos_sum) {
	uint32_t size = 1 << size_log2;
	transform_func[(size <= xy_pos_sum) * 2 + ((xy_pos_sum & (size - 1)) != 0)][(colour + 1) >> 1][size_log2 - 2](frame + (colour >> 1), coeff, stride);
}

template <int N>
static void skip_transform_core(const int16_t* coeff, uint8_t* frame, int stride) {
	for (int y = 0; y < 4; ++y) {
		for (int x = 0; x < 4; ++x) {
			int v = frame[y * stride + x * N] + ((coeff[y * 4 + x] + 16) >> 5);
			frame[y * stride + x * N] = CLIP255C(v);
		}
	}
}

static inline void skip_transform(int16_t coeff[], uint8_t* frame, int stride, int colour, uint32_t xy_pos_sum) {
	static void (* const transform_func[3])(const int16_t* coeff, uint8_t* dst, int stride) = {
		skip_transform_core<1>, skip_transform_core<2>, skip_transform_core<2>
	};
	if (!xy_pos_sum) {
		int v = frame[colour >> 1] + ((coeff[0] + 16) >> 5);
		frame[colour >> 1] = CLIP255C(v);
	} else {
		transform_func[colour](coeff, frame + (colour >> 1), stride);
	}
}

static void residual_coding(h265d_ctu_t& dst, dec_bits& st, uint32_t size_log2, int colour, int pred_idx, uint8_t* frame) {
	static const uint8_t scan_order[35] = {
		0, 0, 0, 0, 0, 0,
		2, 2, 2, 2, 2, 2,
		2, 2, 2, 0, 0, 0,
		0, 0, 0, 0, 1, 1,
		1, 1, 1, 1, 1, 1,
		1, 0, 0, 0, 0
	};
	bool transform_skip = ((size_log2 == 2) && dst.pps->transform_skip_enabled_flag && transform_skip_flag(dst.cabac, st, colour));
	uint32_t max = size_log2 * 2 - 1;
	uint32_t offset, shift;
	if (colour == 0) {
		offset = (size_log2 - 2) * 3 + ((size_log2 - 1) >> 2);
		shift = (size_log2 + 1) >> 2;
	} else {
		offset = 15;
		shift = size_log2 - 2;
	}
	uint32_t x = last_sig_coeff_prefix_luma(dst.cabac, st, reinterpret_cast<h265d_cabac_context_t*>(dst.cabac.context)->last_sig_coeff_x_prefix + offset, shift, max);
	uint32_t y = last_sig_coeff_prefix_luma(dst.cabac, st, reinterpret_cast<h265d_cabac_context_t*>(dst.cabac.context)->last_sig_coeff_y_prefix + offset, shift, max);
	uint32_t last_sig_coeff_x = last_sig_coeff_suffix_add(dst.cabac, st, x);
	uint32_t last_sig_coeff_y = last_sig_coeff_suffix_add(dst.cabac, st, y);
	if ((last_sig_coeff_x != 0) || (last_sig_coeff_y != 0)) {
		memset(dst.coeff_buf, 0, sizeof(dst.coeff_buf[0]) << (size_log2 * 2));
	}
	sub_block_flags_t sub_block_info(size_log2);
	uint32_t order_idx = 0;
	if (colour == 0) {
		if (size_log2 <= 3) {
			order_idx = scan_order[dst.order_luma[pred_idx]];
		}
	} else if (size_log2 == 2) {
		order_idx = scan_order[dst.order_chroma];
	}
	if (order_idx == 2) {
		std::swap(last_sig_coeff_x, last_sig_coeff_y);
	}
	const residual_scan_order_t& order = residual_scan_order[order_idx][size_log2 - 2];
	sig_coeff_flag_inc_func_t sig_coeff_flag_inc = sig_coeff_flag_inc_func[order_idx][size_log2 - 2][(colour + 1) >> 1];
	uint32_t last_subblock_pos = order.sub_block_num[scan_order_index(last_sig_coeff_x >> 2, last_sig_coeff_y >> 2, size_log2 - 2)];
	int i = last_subblock_pos;
	uint8_t greater1ctx = 1;
	uint32_t num = order.inner_pos_num[scan_order_index(last_sig_coeff_x & 3, last_sig_coeff_y & 3, 2)];
	const h265d_scaling_func_t scaling = dst.scaling_func[size_log2 - 2];
	const h265d_scaling_info_t& scaling_info = dst.qp_scale[colour];
	uint32_t xy_pos_sum = 0;
	do {
		uint32_t sxy = order.sub_block_pos[i];
		uint32_t prev_sbf = sub_block_info.prev_flags(sxy);
		if (((unsigned)(last_subblock_pos - 1) <= (unsigned)(i - 1)) || coded_sub_block_flag(dst.cabac, st, prev_sbf, colour)) {
			sub_block_info.set_flag();
			h265d_sigcoeff_t sig_coeff_flags[4 * 4];
			uint32_t num_coeff = sig_coeff_flags_read(dst.cabac, st, sig_coeff_flag_inc, order.inner_xy_pos, ((uint32_t)i == last_subblock_pos), num, sig_coeff_flags, sxy, prev_sbf);
			if (num_coeff == 0) {
				break;
			}
			uint32_t max_coeffs = sig_coeff_greater(colour, i, greater1ctx, sig_coeff_flags, num_coeff, dst.cabac, st);
			uint8_t hidden_bit = (dst.pps->sign_data_hiding_enabled_flag && (3 < sig_coeff_flags[0].pos - sig_coeff_flags[num_coeff - 1].pos));
			uint16_t sign_flags = coeff_sign_flags(dst.cabac, st, num_coeff - hidden_bit);
			xy_pos_sum |= sig_coeff_writeback(num_coeff, sig_coeff_flags, max_coeffs, hidden_bit, sign_flags, order.macro_xy_pos, dst.coeff_buf, sub_block_info.offset(), scaling, scaling_info, dst.cabac, st);
		}
		num = 15;
	} while (0 <= --i);
	if (!transform_skip) {
		transform(dst.coeff_buf, size_log2, frame, dst.sps->ctb_info.stride, colour, xy_pos_sum);
	} else {
		skip_transform(dst.coeff_buf, frame, dst.sps->ctb_info.stride, colour, xy_pos_sum);
	}
}

static void transform_unit(h265d_ctu_t& dst, dec_bits& st, uint32_t size_log2, uint8_t cbf, int idx, int pred_idx, uint32_t offset_x, uint32_t offset_y) {
	uint32_t stride = dst.sps->ctb_info.stride;
	if (cbf & 1) {
		residual_coding(dst, st, size_log2, 0, pred_idx, dst.luma + offset_y * stride + offset_x);
	}
	if (cbf & 6) {
		if (2 < size_log2) {
			size_log2 -= 1;
		} else if (idx != 3) {
			return;
		} else {
			offset_x -= 4;
			offset_y -= 4;
		}
		int offset_chroma = (offset_y >> 1) * stride + offset_x;
		if (cbf & 4) {
			residual_coding(dst, st, size_log2, 1, pred_idx, dst.chroma + offset_chroma);
		}
		if (cbf & 2) {
			residual_coding(dst, st, size_log2, 2, pred_idx, dst.chroma + offset_chroma);
		}
	}
}

struct load_1pix {
	uint32_t operator()(const uint8_t* src) const {
		return src[0];
	}
};

struct store_1pix {
	void operator()(uint8_t* dst, uint32_t val) const {
		dst[0] = val;
	}
};

struct load_2pix {
	uint32_t operator()(const uint8_t* src) const {
		return src[1] * 65536 + src[0];
	}
};

struct store_2pix {
	void operator()(uint8_t* dst, uint32_t val) const {
		dst[0] = val;
		dst[1] = val >> 16;
	}
};

template <typename F>
uint32_t sum_strided(const uint8_t* src, int len, int stride, F Load) {
	uint32_t sum = 0;
	for (int i = 0; i < len; ++i) {
		sum += Load(src + i * stride);
	}
	return sum;
}

template <typename F>
uint32_t sum_edge_base(uint8_t* dst, int len, int valid_main, int valid_sub, int stride, int sub_stride, F Load) {
	if (len <= valid_main) {
		return sum_strided(dst - sub_stride, len, stride, Load);
	} else if (0 < valid_main) {
		return sum_strided(dst - sub_stride, valid_main, stride, Load) + Load(dst + (valid_main - 1) * stride - sub_stride) * (len - valid_main);
	} else if (0 < valid_sub) {
		return Load(dst - stride) * len;
	} else {
		static const uint8_t default_val[] = {
			128, 128
		};
		return Load(default_val) * len;
	}
}

template <int N>
uint32_t sum_edge(uint8_t* dst, int len, int valid_main, int valid_sub, int stride, int sub_stride) {
	return sum_edge_base(dst, len, valid_main, valid_sub, stride, sub_stride, load_1pix());
}

template <>
uint32_t sum_edge<2>(uint8_t* dst, int len, int valid_main, int valid_sub, int stride, int sub_stride) {
	return sum_edge_base(dst, len, valid_main, valid_sub, stride, sub_stride, load_2pix());
}

template <int N>
void fill_dc(uint8_t* dst, int size_log2, int stride, uint32_t dc) {
	if (N == 1) {
		dc *= 0x01010101;
	} else {
		dc = ((dc & 255) | ((dc >> 8) & 0xff00)) * 0x00010001;
	}
	uint32_t y = 1 << size_log2;
	uint32_t cnt = 1 << (size_log2 - 3 + N);
	do {
		for (uint32_t x = 0; x < cnt; ++x) {
			reinterpret_cast<uint32_t*>(dst)[x] = dc;
		}
		dst += stride;
	} while (--y);
}

static inline void intra_dc_filter_left(uint8_t* dst, uint32_t cnt, int stride, uint32_t dc) {
	dc = dc * 3 + 2;
	uint8_t* s = dst - 1;
	do {
		s[1] = (s[0] + dc) >> 2;
		s += stride;
	} while (--cnt);
}

static inline void intra_dc_filter_top(uint8_t* dst, uint32_t cnt, int stride, uint32_t dc) {
	dc = dc * 3 + 2;
	const uint8_t* s = dst - stride;
	for (uint32_t i = 0; i < cnt; ++i) {
		dst[i] = (s[i] + dc) >> 2;
	}
}

static inline void intra_dc_filter_leftonly(uint8_t* dst, uint32_t cnt, int stride, uint32_t dc) {
	uint32_t left = dst[-1];
	dst[0] = (left + dc + 1) >> 1;
	uint32_t dc1 = (left + dc * 3 + 2) >> 2;
	memset(dst + 1, dc1, cnt - 1);
	intra_dc_filter_left(dst + stride, cnt - 1, stride, dc);
}

static inline void intra_dc_filter_toponly(uint8_t* dst, uint32_t cnt, int stride, uint32_t dc) {
	intra_dc_filter_top(dst + 1, cnt - 1, stride, dc);
	uint32_t top = *(dst - stride);
	uint32_t dc1 = (top + dc * 3 + 2) >> 2;
	dst[0] = (top + dc + 1) >> 1;
	do {
		dst += stride;
		dst[0] = dc1;
	} while (--cnt);
}

static inline void intra_dc_filter_both(uint8_t* dst, uint32_t cnt, int stride, uint32_t dc) {
	dst[0] = (*(dst - stride) + dst[-1] + dc * 2 + 2) >> 2;
	intra_dc_filter_top(dst + 1, cnt - 1, stride, dc);
	intra_dc_filter_left(dst + stride, cnt - 1, stride, dc);
}

template <int N>
static inline void intra_pred_dc(uint8_t* dst, int size_log2, int stride, int valid_x, int valid_y) {
	uint32_t size = 1 << size_log2;
	uint32_t dc = sum_edge<N>(dst, size, valid_x, valid_y, N, stride) + sum_edge<N>(dst, size, valid_y, valid_x, stride, N) + (((N == 1) ? 1 : 0x00010001) << size_log2);
	dc = dc >> (size_log2 + 1);
	fill_dc<N>(dst, size_log2, stride, dc);
	if ((N == 1) && (size < 32)) {
		if (0 < valid_x) {
			if (0 < valid_y) {
				intra_dc_filter_both(dst, size, stride, dc);
			} else {
				intra_dc_filter_toponly(dst, size, stride, dc);
			}
		} else if (0 < valid_y) {
			intra_dc_filter_leftonly(dst, size, stride, dc);
		}
	}
}

template <int N>
static void intra_pred_planar_core(uint8_t* dst, const uint8_t* src, int size_log2, int stride) {
	uint32_t size = 1 << size_log2;
	const uint8_t* top = src + (size + 1) * N;
	uint32_t left_bottom = src[size * N];
	uint32_t right_top = top[size * N];
	uint32_t topscale = size;
	uint32_t vleft = 0;
	for (uint32_t y = 0; y < size; ++y) {
		uint32_t left = src[y * N];
		topscale--;
		vleft += left_bottom;
		uint32_t xinc = right_top - left;
		uint32_t base = (left << size_log2) + vleft;
		for (uint32_t x = 0; x < size; ++x) {
			base += xinc;
			dst[x * N] = (base + top[x * N] * topscale + size) >> (size_log2 + 1);
		}
		dst += stride;
	}
}

static inline bool abs_smaller(int val, int threshold) {
	return val * val < threshold * threshold;
}

static inline bool intra_pred_detect_strong_onedir(int lt, const uint8_t* src, int valid_len, int stride) {
	if (64 <= valid_len) {
		return abs_smaller(lt + src[64 * stride] - src[32 * stride] * 2, 8);
	} else if (32 <= valid_len) {
		return abs_smaller(lt - src[32 * stride], 8);
	} else {
		return true;
	}
}

static bool intra_pred_detect_strong_filter(bool strong_filter_enabled, const uint8_t* src, int size_log2, int stride, int valid_x, int valid_y) {
	if (strong_filter_enabled && (size_log2 == 5)) {
		if (0 < valid_x) {
			if (0 < valid_y) {
				src = src - stride - 1;
				int lt = src[0];
				return intra_pred_detect_strong_onedir(lt, src, valid_x, 1) && intra_pred_detect_strong_onedir(lt, src, valid_y, stride);
			} else {
				src -= stride;
				return intra_pred_detect_strong_onedir(src[0], src - 1, valid_x, 1);
			}
		} else if (0 < valid_y) {
			src -= 1;
			return intra_pred_detect_strong_onedir(src[0], src - stride, valid_y, stride);
		} else {
			return false;
		}
	}
	return false;
}

#include "intrapos.h"

template <int N>
struct get_pix_raw {
	uint32_t operator()(const uint8_t* src, int offset, int offset_min, int offset_max, int stride) const {
		int ofs = ((offset_min <= offset) ? ((offset < offset_max) ? offset : offset_max - 1) : offset_min) * stride;
		uint32_t pix = src[ofs];
		for (int i = 1; i < N; ++i) {
			pix |= src[ofs + i] << (8 * i);
		}
		return pix;
	}
};

template <int N>
static inline void multipix_fill_gap(uint8_t* dst, int midlen, int pregap, int postgap) {
	if (0 < pregap) {
		if (N == 1) {
			memset(dst, dst[pregap], pregap);
		} else {
			std::fill((uint16_t*)dst, (uint16_t*)dst + pregap, ((const uint16_t*)dst)[pregap]);
		}
	}
	if (0 < postgap) {
		if (N == 1) {
			memset(dst + pregap + midlen, dst[pregap + midlen - 1], postgap);
		} else {
			std::fill((uint16_t*)dst + pregap + midlen, (uint16_t*)dst + pregap + midlen + postgap, ((const uint16_t*)dst)[pregap + midlen - 1]);
		}
	}
}

template <typename T>
static inline void memcpy_from_vertical(T* dst, const void* src, size_t size, int stride) {
	T* dst_end = dst + size;
	const uint8_t* s0 = reinterpret_cast<const uint8_t*>(src);
	while (dst != dst_end) {
		*dst++ = *reinterpret_cast<const T*>(s0);
		s0 += stride;
	}
}

template <int N>
static inline void get_multipix_raw_core(uint8_t* dst, const uint8_t* src, int offset, int offset_min, int offset_max, int size_log2, int len, int stride, int sub_stride) {
	int pregap;
	if (offset_min <= offset) {
		pregap = 0;
		src += stride * offset;
	} else {
		pregap = offset_min - offset;
		src += stride * offset_min;
	}
	int midlen = std::min(offset_max - offset - pregap, len);
	int postgap = len - pregap - midlen;
	uint8_t* dst_ofs = dst + pregap * N;
	if (stride == N) {
		memcpy(dst_ofs, src, midlen * N);
	} else {
		if (N == 1) {
			memcpy_from_vertical(dst_ofs, src, midlen, stride);
		} else {
			memcpy_from_vertical((uint16_t*)dst_ofs, src, midlen, stride);
		}
	}
	multipix_fill_gap<N>(dst, midlen, pregap, postgap);
}

template <int N>
struct get_multipix_raw {
	void operator()(uint8_t* dst, const uint8_t* src, int offset, int offset_min, int offset_max, int size_log2, int len, int stride, int sub_stride) const {
		get_multipix_raw_core<N>(dst, src, offset, offset_min, offset_max, size_log2, len, stride, sub_stride);
	}
};

struct get_pix_filtered_strong {
	uint32_t operator()(const uint8_t* src, int offset, int offset_min, int offset_max, int stride) const {
		return ((63 - offset) * src[(offset_min < 0) ? -stride : 0] + (offset + 1) * src[stride * std::min(63, offset_max - 1)] + 32) >> 6;
	}
};

static inline void get_multipix_filtered_strong_core(uint8_t* dst, const uint8_t* src, int offset, int offset_min, int offset_max, int size_log2, int len, int stride, int sub_stride) {
	uint32_t c0 = src[(offset_min < 0) ? -stride : 0];
	uint32_t c1 = src[stride * std::min(63, offset_max - 1)];
	for (int i = 0; i < len; ++i) {
		dst[i] = ((63 - offset) * c0 + (offset + 1) * c1 + 32) >> 6;
		offset++;
	}
}

struct get_multipix_filtered_strong {
	void operator()(uint8_t* dst, const uint8_t* src, int offset, int offset_min, int offset_max, int size_log2, int len, int stride, int sub_stride) const {
		get_multipix_filtered_strong_core(dst, src, offset, offset_min, offset_max, size_log2, len, stride, sub_stride);
	}
};

struct get_pix_filtered {
	uint8_t operator()(const uint8_t* src, int offset, int offset_min, int offset_max, int stride) const {
		int c1 = src[offset * stride];
		if (offset_min < offset) {
			int c0 = src[(offset - 1) * stride];
			if (offset < offset_max - 1) {
				return (c0 + c1 * 2 + src[(offset + 1) * stride] + 2) >> 2;
			} else {
				return (c0 + c1 * 3 + 2) >> 2;
			}
		} else {
			return (c1 * 3 + src[(offset + 1) * stride] + 2) >> 2;
		}
	}
};

static inline void get_multipix_filtered_core(uint8_t* dst, const uint8_t* src, int offset, int offset_min, int offset_max, int size_log2, int len, int stride, int sub_stride) {
	int c0, c1;
	if (offset_min < offset) {
		c0 = src[(offset - 1) * stride];
		c1 = src[offset * stride];
	} else if (offset_min == offset) {
		c1 = src[offset * stride];
		if (offset_min < 0) {
			c0 = src[sub_stride - stride];
		} else {
			c0 = c1;
		}
	} else {
		c0 = c1 = src[(offset + 1) * stride];
	}
	src += offset * stride;
	int midlen = std::min(offset_max - offset - 1, len);
	for (int i = 0; i < midlen; ++i) {
		src += stride;
		uint32_t c2 = src[0];
		dst[i] = (c0 + c1 * 2 + c2 + 2) >> 2;
		c0 = c1;
		c1 = c2;
	}
	while (midlen < len) {
		dst[midlen++] = (c0 + c1 * 3 + 2) >> 2;
		c0 = c1;
	}
	if (((2 << size_log2) <= offset_max) && (offset + len == (2 << size_log2))) {
		dst[len - 1] = c1;
	}
}

struct get_multipix_filtered {
	void operator()(uint8_t* dst, const uint8_t* src, int offset, int offset_min, int offset_max, int size_log2, int len, int stride, int sub_stride) const {
		get_multipix_filtered_core(dst, src, offset, offset_min, offset_max, size_log2, len, stride, sub_stride);
	}
};

template <int N>
void fillmem(uint8_t* dst, const uint8_t* src, int len) {}

template <>
void fillmem<1>(uint8_t* dst, const uint8_t* src, int len) {
	memset(dst, *src, len);
}

template <>
void fillmem<2>(uint8_t* dst, const uint8_t* src, int len) {
	std::fill(reinterpret_cast<uint16_t*>(dst), reinterpret_cast<uint16_t*>(dst) + len, *reinterpret_cast<const uint16_t*>(src));
}

template <int N>
static void intra_pred_planar(uint8_t* dst, int size_log2, int stride, int valid_x, int valid_y, bool strong_enabled) {
	uint8_t neighbour[2 * 34];
	if ((valid_x <= 0) && (valid_y <= 0)) {
		fill_dc<N>(dst, size_log2, stride, (N == 1) ? 128 : 0x00800080);
	} else {
		void (*multipix)(uint8_t* dst, const uint8_t* src, int offset, int offset_min, int offset_max, int size_log2, int len, int stride, int sub_stride);
		if ((N == 1) && (3 <= size_log2)) {
			if (intra_pred_detect_strong_filter(strong_enabled, dst, size_log2, stride, valid_x, valid_y)) {
				multipix = get_multipix_filtered_strong_core;
			} else {
				multipix = get_multipix_filtered_core;
			}
		} else {
			multipix = get_multipix_raw_core<N>;
		}
		if (0 < valid_y) {
			multipix(&neighbour[0], dst - N, 0, (0 < valid_x) ? -1 : 0, valid_y, size_log2, (1 << size_log2) + 1, stride, N);
		} else {
			fillmem<N>(neighbour, dst - stride, (1 << size_log2) + 1);
		}
		if (0 < valid_x) {
			multipix(&neighbour[((1 << size_log2) + 1) * N], dst - stride, 0, (0 < valid_y) ? -1 : 0, valid_x, size_log2, (1 << size_log2) + 1, N, stride);
		} else {
			fillmem<N>(&neighbour[((1 << size_log2) + 1) * N], dst - N, (1 << size_log2) + 1);
		}
		for (int i = 0; i < N; ++i) {
			intra_pred_planar_core<N>(dst + i, neighbour + i, size_log2, stride);
		}
	}
}

template <int N>
static void memfill_pix(uint8_t* dst, uint32_t pat, size_t len) {
	for (size_t i = 0; i < len; ++i) {
		for (size_t j = 0; j < N; ++j) {
			dst[i * N + j] = pat >> (j * 8);
		}
	}
}

template <>
void memfill_pix<1>(uint8_t* dst, uint32_t pat, size_t len) {
	memset(dst, pat, len);
}

template <int N, typename F0>
static inline void intra_pred_get_ref_extra(uint8_t dst[], const uint8_t src[], uint32_t size_log2, int base_len, int extra_len, const int8_t* pos_tbl, int main_stride, int sub_stride, int valid_main, int valid_sub, F0 GetPix) {
	if (0 < valid_sub) {
		src -= main_stride;
		int ofs_min = (0 < valid_main) ? -1 : 0;
		for (int i = 0; i < extra_len; ++i) {
			memfill_pix<N>(dst + i * N, GetPix(src, pos_tbl[i], ofs_min, valid_sub, sub_stride), 1);
		}
	} else {
		if (0 < valid_main) {
			uint32_t pix = src[-sub_stride];
			pix = (N == 1) ? pix : pix | (src[-sub_stride + 1] << 8);
			memfill_pix<N>(dst, pix, extra_len);
		} else {
			memset(dst, 128, extra_len * N);
		}
	}
}

template <int N, typename F0, typename F1>
static void intra_pred_get_ref(uint8_t dst[], const uint8_t src[], int size_log2, int main_stride, int sub_stride, int valid_main, int valid_sub, const int8_t* pos_tbl, F0 GetPix, F1 GetMultiPix) {
	int extra_len = *pos_tbl++;
	int base_len = pos_tbl[extra_len + 1];
	if (extra_len != 0) {
		intra_pred_get_ref_extra<N>(dst, src, size_log2, base_len, extra_len, pos_tbl, main_stride, sub_stride, valid_main, valid_sub, GetPix);
		dst = dst + extra_len * N;
	}
	int base_pos = pos_tbl[extra_len];
	if (0 < valid_main) {
		GetMultiPix(dst, src - sub_stride, base_pos, (0 < valid_sub) ? -1 : 0, std::min(2 << size_log2, valid_main), size_log2, base_len, main_stride, sub_stride);
	} else {
		if (0 < valid_sub) {
			memfill_pix<N>(dst, get_pix_raw<N>()(src - main_stride, 0, 0, 4, sub_stride), base_len);
		} else {
			memset(dst, 128, base_len * N);
		}
	}
}

static void intra_pred_get_ref_filtered_strong(uint8_t dst[], const uint8_t src[], int size_log2, int stride, int valid_x, int valid_y, uint32_t mode) {
	const int8_t* pos_tbl = intra_pred_pos[mode][size_log2 - 2];
	if (mode < 16) {
		intra_pred_get_ref<1>(dst, src, size_log2, stride, 1, valid_y, valid_x, pos_tbl, get_pix_filtered_strong(), get_multipix_filtered_strong());
	} else {
		intra_pred_get_ref<1>(dst, src, size_log2, 1, stride, valid_x, valid_y, pos_tbl, get_pix_filtered_strong(), get_multipix_filtered_strong());
	}
}

static void intra_pred_get_ref_filtered(uint8_t dst[], const uint8_t src[], int size_log2, int stride, int valid_x, int valid_y, uint32_t mode) {
	const int8_t* pos_tbl = intra_pred_pos[mode][size_log2 - 2];
	if (mode < 16) {
		intra_pred_get_ref<1>(dst, src, size_log2, stride, 1, valid_y, valid_x, pos_tbl, get_pix_filtered(), get_multipix_filtered());
	} else {
		intra_pred_get_ref<1>(dst, src, size_log2, 1, stride, valid_x, valid_y, pos_tbl, get_pix_filtered(), get_multipix_filtered());
	}
}

template <int N>
static void intra_pred_get_ref_raw(uint8_t dst[], const uint8_t src[], int size_log2, int stride, int valid_x, int valid_y, uint32_t mode) {
	const int8_t* pos_tbl = intra_pred_pos[mode][size_log2 - 2];
	if (mode < 16) {
		intra_pred_get_ref<N>(dst, src, size_log2, stride, N, valid_y, valid_x, pos_tbl, get_pix_raw<N>(), get_multipix_raw<N>());
	} else {
		intra_pred_get_ref<N>(dst, src, size_log2, N, stride, valid_x, valid_y, pos_tbl, get_pix_raw<N>(), get_multipix_raw<N>());
	}
}

template <int N, typename F0, typename F1>
static void intra_pred_angular_filter(uint8_t* dst, const uint8_t* ref, uint32_t size_log2, uint32_t inner_inc, uint32_t outer_inc, uint32_t mode, F0 Load, F1 Store) {
	const int8_t* coef = intra_pred_coef[mode][0];
	const int8_t* inc = intra_pred_coef[mode][1];
	uint32_t size = 1 << size_log2;
	const uint8_t* src = ref + (*inc++ >> (5 - size_log2)) * N;
	for (uint32_t y = 0; y < size; ++y) {
		int c1 = coef[y];
		int c0 = 32 - c1;
		uint32_t d0 = Load(src);
		const uint8_t* s = src + N;
		for (uint32_t x = 0; x < size; ++x) {
			uint32_t d1 = Load(s + x * N);
			Store(dst + x * inner_inc, (d0 * c0 + d1 * c1 + 0x00100010) >> 5);
			d0 = d1;
		}
		src += inc[y] * N;
		dst += outer_inc;
	}
}

static inline void intra_pred_post_edge_filter(uint8_t* dst, const uint8_t* ref, uint32_t len, uint32_t stride) {
	uint32_t d0 = dst[0];
	uint32_t c0 = *ref++;
	for (uint32_t x = 0; x < len; ++x) {
		int t0 = d0 + ((int)((ref[x]) - c0) >> 1);
		dst[x * stride] = CLIP255C(t0);
	}
}

template <int N>
static void intra_pred_diagonal(uint8_t* dst, const uint8_t* ref, uint32_t size_log2, uint32_t stride, uint32_t mode) {
	const int8_t* inc = intra_pred_coef[mode][1];
	uint32_t size = N << size_log2;
	uint32_t y = 1 << size_log2;
	const uint8_t* src = ref + (inc[0] >> (5 - size_log2)) * N;
	int src_inc = inc[1] * N;
	do {
		memcpy(dst, src, size);
		src += src_inc;
		dst += stride;
	} while (--y);
}

template <int N, typename F0, typename F1>
static void intra_pred_angular(uint8_t* dst, int size_log2, int stride, int valid_x, int valid_y, uint32_t mode, bool strong_enabled, F0 Load, F1 Store) {
	static const int8_t filter_thr[16] = {
		56, 48, 48, 48, 48, 48, 48, 32, 0, 32, 48, 48, 48, 48, 48, 48
	};
	uint8_t neighbour[64];
	if ((N == 1) && (filter_thr[mode & 15] & (1 << size_log2))) {
		if (intra_pred_detect_strong_filter(strong_enabled, dst, size_log2, stride, valid_x, valid_y)) {
			intra_pred_get_ref_filtered_strong(neighbour, dst, size_log2, stride, valid_x, valid_y, mode);
		} else {
			intra_pred_get_ref_filtered(neighbour, dst, size_log2, stride, valid_x, valid_y, mode);
		}
	} else {
		intra_pred_get_ref_raw<N>(neighbour, dst, size_log2, stride, valid_x, valid_y, mode);
	}
	if (mode & 7) {
		if (mode < 16) {
			intra_pred_angular_filter<N>(dst, neighbour, size_log2, stride, N, mode, Load, Store);
		} else {
			intra_pred_angular_filter<N>(dst, neighbour, size_log2, N, stride, mode, Load, Store);
		}
	} else {
		intra_pred_diagonal<N>(dst, neighbour, size_log2, stride, mode);
	}
}

static inline void intra_pred_postfilter(uint8_t* dst, int len, int stride, int sub_stride, uint32_t c0) {
	uint32_t d0 = dst[0];
	const uint8_t* src = dst - sub_stride;
	for (int x = 0; x < len; ++x) {
		int t0 = d0 + ((int)((src[x * stride]) - c0) >> 1);
		dst[x * stride] = CLIP255C(t0);
	}
}

template <int N>
static inline void intra_pred_horizontal_raw(uint8_t* dst, int size_log2, int stride) {
	int height = 1 << size_log2;
	int xcnt = N << (size_log2 - 2);
	for (int y = 0; y < height; ++y) {
		uint32_t pat;
		if (N == 1) {
			pat = dst[-1] * 0x01010101;
		} else {
			pat = reinterpret_cast<const uint16_t*>(dst)[-1] * 0x00010001;
		}
		for (int x = 0; x < xcnt; ++x) {
			reinterpret_cast<uint32_t*>(dst)[x] = pat;
		}
		dst += stride;
	}
}

template <int N>
static void intra_pred_horizontal_mode(uint8_t* dst, int size_log2, int stride, int valid_x, int valid_y) {
	if (0 < valid_y) {
		intra_pred_horizontal_raw<N>(dst, size_log2, stride);
		if ((N == 1) && (size_log2 < 5) && (0 < valid_x)) {
			intra_pred_postfilter(dst, 1 << size_log2, 1, stride, dst[-1 - stride]);
		}
	} else {
		uint32_t dc = (0 < valid_x) ? ((N == 1) ? *(dst - stride) : load_2pix()(dst - stride)) : ((N == 1) ? 128 : 0x00800080);
		fill_dc<N>(dst, size_log2, stride, dc);
		if ((N == 1) && (size_log2 < 5) && (0 < valid_x)) {
			intra_pred_postfilter(dst, 1 << size_log2, 1, stride, dc);
		}
	}
}

template <int N>
static void intra_pred_vertical_raw(uint8_t* dst, int size_log2, int stride) {
	int height = (1 << size_log2) + 1;
	int xcnt = N << (size_log2 - 2);
	dst -= stride;
	for (int x = 0; x < xcnt; ++x) {
		uint32_t pat = reinterpret_cast<const uint32_t*>(dst)[x];
		for (int y = 1; y < height; ++y) {
			reinterpret_cast<uint32_t*>(dst + y * stride)[x] = pat;
		}
	}
}

template <int N>
static void intra_pred_vertical_mode(uint8_t* dst, int size_log2, int stride, int valid_x, int valid_y) {
	if (0 < valid_x) {
		intra_pred_vertical_raw<N>(dst, size_log2, stride);
		if ((N == 1) && (size_log2 < 5) && (0 < valid_y)) {
			intra_pred_postfilter(dst, 1 << size_log2, stride, 1, dst[-1 - stride]);
		}
	} else {
		uint32_t dc = (0 < valid_y) ? ((N == 1) ? *(dst - N) : load_2pix()(dst - N)) : ((N == 1) ? 128 : 0x00800080);
		fill_dc<N>(dst, size_log2, stride, dc);
		if ((N == 1) && (size_log2 < 5) && (0 < valid_y)) {
			intra_pred_postfilter(dst, 1 << size_log2, stride, 1, dc);
		}
	}
}

template <int N, typename F0, typename F1>
static inline void intra_prediction_dispatch(uint8_t* dst, int size_log2, int stride, int valid_x, int valid_y, uint32_t mode, bool strong_enabled, F0 Load, F1 Store) {
	switch (mode) {
	case 0:
		intra_pred_planar<N>(dst, size_log2, stride, valid_x, valid_y, strong_enabled);
		break;
	case 1:
		intra_pred_dc<N>(dst, size_log2, stride, valid_x, valid_y);
		break;
	case 8 + 2:
		intra_pred_horizontal_mode<N>(dst, size_log2, stride, valid_x, valid_y);
		break;
	case 24 + 2:
		intra_pred_vertical_mode<N>(dst, size_log2, stride, valid_x, valid_y);
		break;
	default:
		intra_pred_angular<N>(dst, size_log2, stride, valid_x, valid_y, mode - 2, strong_enabled, Load, Store);
		break;
	}
}

static inline void intra_prediction(const h265d_ctu_t& dst, int size_log2, int offset_x, int valid_x, int offset_y, int valid_y, int pred_idx) {
	uint32_t stride = dst.sps->ctb_info.stride;
	intra_prediction_dispatch<1>(dst.luma + offset_y * stride + offset_x, size_log2, stride, valid_x, valid_y, dst.order_luma[pred_idx], dst.sps->strong_intra_smoothing_enabled_flag, load_1pix(), store_1pix());
	if (size_log2 == 2) {
		return;
	}
	intra_prediction_dispatch<2>(dst.chroma + (offset_y >> 1) * stride + offset_x, size_log2 - 1, stride, valid_x >> 1, valid_y >> 1, dst.order_chroma, false, load_2pix(), store_2pix());
}

static inline void qpy_fill(uint8_t left[], uint8_t top[], uint32_t qpy, int num) {
	for (int i = 0; i < num; ++i) {
		left[i] = qpy;
		top[i] = qpy;
	}
}

static inline int qpi_to_qpc(int qpi) {
	static const int8_t qpcadj[52] = {
		0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
		16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 29, 30,
		31, 32, 33, 33, 34, 34, 35, 35, 36, 36, 37, 37, 38, 39, 40, 41,
		42, 43, 44, 45
	};
	return qpcadj[qpi % 52];
}

static inline void qp_to_scale(h265d_scaling_info_t dst[], uint32_t qpy, const int8_t* cbcr_delta) {
	static const int16_t qp_scale[] = {
		40, 45, 51, 57, 64, 72,
		80, 90, 102, 114, 128, 144,
		160, 180, 204, 228, 256, 288,
		320, 360, 408, 456, 512, 576,
		640, 720, 816, 912, 1024, 1152,
		1280, 1440, 1632, 1824, 2048, 2304,
		2560, 2880, 3264, 3648, 4096, 4608,
		5120, 5760, 6528, 7296, 8192, 9216,
		10240, 11520, 13056, 14592
	};
	dst[0].scale = qp_scale[qpy];
	dst[1].scale = qp_scale[qpi_to_qpc(qpy + cbcr_delta[0])];
	dst[2].scale = qp_scale[qpi_to_qpc(qpy + cbcr_delta[1])];
}

static inline void qpy_update(h265d_ctu_t& dst, uint8_t qp_left[], uint8_t qp_top[], int32_t qp_delta) {
	uint32_t qpy_org = dst.qpy;
	uint32_t left = qp_left[0];
	uint32_t top = qp_top[0];
	uint32_t qpy = (((left + top + 1) >> 1) + qp_delta + 52) % 52;
	if (qpy_org != qpy) {
		dst.qpy = qpy;
		qp_to_scale(dst.qp_scale, qpy, dst.qpc_delta);
	}
}

static void record_block_boundary_onedir(h265d_ctu_t& ctu, int dir, int offset_x, int offset_y, int xgap, int ygap, int len, int edgemax, int strength) {
	if (offset_x & 7) {
		return;
	}
	int org_x = offset_x >> 3;
	int org_y = offset_y >> 2;
	h265d_deblocking_strength_t* boundary = ctu.deblock_boundary[dir] + org_x * xgap + (org_y + 1) * ygap;
	int qp = ctu.qpy + 1;
	const uint8_t* prev_qp = &ctu.qp_history[dir][org_y];
	for (int n = 0; n < len; ++n) {
		boundary[n * ygap].qp = (qp + prev_qp[n]) >> 1;
		boundary[n * ygap].str = strength;
	}
}

static void record_block_boundary(h265d_ctu_t& ctu, int size_log2, int offset_x, int offset_y, uint32_t unavail, int str) {
	if (ctu.slice_header->body.deblocking_filter_disabled_flag) {
		return;
	}
	int edgemax = 1 << (ctu.sps->ctb_info.size_log2 - 3);
	int len = 1 << (size_log2 - 2);
	if (offset_x || !(unavail & 1)) {
		record_block_boundary_onedir(ctu, 0, offset_x, offset_y, 1, edgemax, len, edgemax, str);
	}
	if (offset_y || !(unavail & 2)) {
		record_block_boundary_onedir(ctu, 1, offset_y, offset_x, edgemax * 2 + 1, 1, len, edgemax, str);
	}
}

static void transform_tree(h265d_ctu_t& dst, dec_bits& st, int size_log2, uint32_t unavail, uint32_t depth, uint8_t upper_cbf_cbcr, int offset_x, int valid_x, int offset_y, int valid_y, int idx, int pred_idx, bool is_intra) {
	uint32_t split;
	if (dst.sps->ctb_info.transform_log2 < size_log2) {
		split = 1;
	} else if ((depth == 0) && dst.intra_split) {
		split = 2;
	} else {
		if ((dst.sps->ctb_info.transform_log2_min < size_log2) && (depth < dst.sps->max_transform_hierarchy_depth_intra)) {
			split = split_transform_flag(dst.cabac, st, size_log2);
		} else {
			split = 0;
		}
	}
	uint8_t cbf;
	if (2 < size_log2) {
		cbf = 0;
		if (upper_cbf_cbcr & 2) {
			cbf |= cbf_chroma(dst.cabac, st, depth) * 2;
		}
		if (upper_cbf_cbcr & 1) {
			cbf |= cbf_chroma(dst.cabac, st, depth);
		}
	} else {
		cbf = upper_cbf_cbcr;
	}
	if (split) {
		int pi0, pi1, pi2, pi3;
		if (split == 2) {
			pi0 = 0;
			pi1 = 1;
			pi2 = 2;
			pi3 = 3;
		} else {
			pi0 = pred_idx;
			pi1 = pred_idx;
			pi2 = pred_idx;
			pi3 = pred_idx;
		}
		size_log2 -= 1;
		if (is_intra && (size_log2 == 2)) {
			uint32_t stride = dst.sps->ctb_info.stride;
			intra_prediction_dispatch<2>(dst.chroma + (offset_y >> 1) * stride + offset_x, size_log2, stride, (unavail & 2) ? -1 : (valid_x >> 1), (unavail & 1) ? -1 : (valid_y >> 1), dst.order_chroma, false, load_2pix(), store_2pix());
		}
		depth += 1;
		uint32_t block_len = 1 << size_log2;
		transform_tree(dst, st, size_log2, unavail, depth, cbf, offset_x, valid_x, offset_y, valid_y, 0, pi0, is_intra);
		transform_tree(dst, st, size_log2, unavail & ~1, depth, cbf, offset_x + block_len, valid_x - block_len, offset_y, std::min(static_cast<uint32_t>(valid_y), block_len), 1, pi1, is_intra);
		transform_tree(dst, st, size_log2, unavail & ~2, depth, cbf, offset_x, std::min(static_cast<uint32_t>(valid_x), block_len * 2), offset_y + block_len, valid_y - block_len, 2, pi2, is_intra);
		transform_tree(dst, st, size_log2, 0, depth, cbf, offset_x + block_len, std::min(static_cast<uint32_t>(valid_x - block_len), block_len), offset_y + block_len, std::min(static_cast<uint32_t>(valid_y - block_len), block_len), 3, pi3, is_intra);
	} else {
		if (is_intra) {
			intra_prediction(dst, size_log2, offset_x, (unavail & 2) ? -1 : valid_x, offset_y, (unavail & 1) ? -1 : valid_y, pred_idx);
		}
		cbf = cbf * 2;
		if (is_intra || depth || cbf) {
			cbf |= cbf_luma(dst.cabac, st, depth);
		}
		int32_t qp_delta = 0;
		if (dst.qp_delta_req) {
			dst.qp_delta_req = 0;
			qp_delta = cu_qp_delta(dst.cabac, st);
		}
		qpy_update(dst, &dst.qp_history[0][offset_y >> 2], &dst.qp_history[1][offset_x >> 2], qp_delta);
		if (cbf) {
			transform_unit(dst, st, size_log2, cbf, idx, pred_idx, offset_x, offset_y);
		}
		record_block_boundary(dst, size_log2, offset_x, offset_y, unavail, 2);
		qpy_fill(&dst.qp_history[0][offset_y >> 2], &dst.qp_history[1][offset_x >> 2], dst.qpy, 1 << (size_log2 - 2));
	}
}

static inline void cu_pred_mode_fill(h265d_neighbour_t dst0[], h265d_neighbour_t dst1[], uint32_t mode, uint32_t num) {
	for (uint32_t i = 0; i < num; ++i) {
		dst0[i].pred_mode = mode;
		dst1[i].pred_mode = mode;
	}
}

static inline void intra_depth_fill(h265d_neighbour_t dst0[], h265d_neighbour_t dst1[], uint32_t size_log2) {
	int num = 1 << (size_log2 - 2);
	uint8_t depth = 6 - size_log2;
	for (int i = 0; i < num; ++i) {
		dst0[i].depth = depth;
		dst1[i].depth = depth;
	}
}

static void prediction_unit_merge(h265d_ctu_t& dst, dec_bits& st, int offset_x, int offset_y, int width, int height) {
	if (1 < dst.slice_header->body.max_num_merge_cand) {
		merge_idx(dst.cabac, st);
	}
}

static inline int mvd_coding_suffix(m2d_cabac_t& cabac, dec_bits& st, int mvd) {
	if (mvd) {
		if (1 < mvd) {
			mvd += abs_mvd_minus2(cabac, st);
		}
		mvd = mvd_sign_flag(cabac, st) ? -mvd : mvd;
	}
	return mvd;
}

static inline void mvd_coding(m2d_cabac_t& cabac, dec_bits& st, int16_t mvd[]) {
	int mvd0 = abs_mvd_greater_flag(cabac, st, 0);
	int mvd1 = abs_mvd_greater_flag(cabac, st, 0);
	mvd0 = (mvd0) ? mvd0 + abs_mvd_greater_flag(cabac, st, 1) : mvd0;
	mvd1 = (mvd1) ? mvd1 + abs_mvd_greater_flag(cabac, st, 1) : mvd1;
	mvd[0] = mvd_coding_suffix(cabac, st, mvd0);
	mvd[1] = mvd_coding_suffix(cabac, st, mvd1);
}

static const int16_t* pred_mv(uint32_t unavail, int offset_y, int height, h265d_neighbour_t* left, int lx, const pred_info_t& pred, int offset) {
	bool is_scaled = true;
	bool available = false;
	unavail >>= offset;
	if (!(unavail & 5)) {
		int pos;
		if (!(unavail & 4)) {
			pos = height >> 2;
		} else {
			pos = 0;
		}
		int a_ref_idx;
		if ((left[pos].pred.idc == lx) || (left[pos].pred.idc == 2)) {
			a_ref_idx = (lx == 0) ? left[pos].pred.ref_idx0 : left[pos].pred.ref_idx1;
		} else {
			a_ref_idx = (lx == 0) ? left[pos].pred.ref_idx1 : left[pos].pred.ref_idx0;
		}
		int ref_idx = (lx == 0) ? pred.ref_idx0 : pred.ref_idx1;
		if (a_ref_idx == ref_idx) {
			return left[pos].pred.mvd[lx];
		}
//		if (ref_list[a_ref_idx].poc == ref_list[a_ref_idx].poc) {
	} else {
		is_scaled = false;
	}
}

static void pred_block(h265d_ctu_t& dst, uint32_t unavail, int offset_x, int offset_y, int width, int height, h265d_neighbour_t* left, h265d_neighbour_t* top, int lx, const pred_info_t& pred) {
	pred_mv(unavail, offset_y, height, left, lx, pred, 0);
	pred_mv(unavail, offset_x, width, top, lx, pred, 1);
}

static bool prediction_unit(h265d_ctu_t& dst, dec_bits& st, int size_log2, uint32_t unavail, int offset_x, int offset_y, int width, int height, h265d_neighbour_t* left, h265d_neighbour_t* top) {
	if (merge_flag(dst.cabac, st)) {
		prediction_unit_merge(dst, st, offset_x, offset_y, width, height);
		return true;
	} else {
		int pred_idc;
		if (dst.slice_header->body.slice_type == 0) {
			int depth = dst.sps->ctb_info.size_log2 - size_log2;
			pred_idc = inter_pred_idc(dst.cabac, st, width, height, depth);
		} else {
			pred_idc = 0;
		}
		pred_info_t pred;
		pred.idc = pred_idc;
		if (pred_idc != 1) {
			pred.ref_idx0 = ref_idx_lx(dst.cabac, st, 0, dst.slice_header->body.num_ref_idx_lx_active_minus1);
			mvd_coding(dst.cabac, st, pred.mvd[0]);
			pred.lx_flag0 = mvp_lx_flag(dst.cabac, st, 0);
			pred_block(dst, unavail, offset_x, offset_y, width, height, left, top, 0, pred);
		}
		if (pred_idc != 0) {
			pred.ref_idx1 = ref_idx_lx(dst.cabac, st, 1, dst.slice_header->body.num_ref_idx_lx_active_minus1);
			if ((pred_idc == 1) || !dst.slice_header->body.mvd_l1_zero_flag) {
				mvd_coding(dst.cabac, st, pred.mvd[1]);
			} else {
				memset(pred.mvd, 0, sizeof(pred.mvd[1]));
			}
			pred.lx_flag1 = mvp_lx_flag(dst.cabac, st, 1);
			pred_block(dst, unavail, offset_x, offset_y, width, height, left, top, 1, pred);
		}
		return false;
	}
}

static bool prediction_unit_cases(h265d_ctu_t& dst, dec_bits& st, int size_log2, uint32_t unavail, int offset_x, int offset_y, int valid_x, int valid_y, h265d_neighbour_t* left, h265d_neighbour_t* top) {
	h265d_inter_part_mode_t mode = part_mode_inter(dst.cabac, st, size_log2, dst.sps->ctb_info.size_log2_min, dst.sps->amp_enabled_flag);
	int len = 1 << size_log2;
	int len_s;
	switch (mode) {
	case PART_2Nx2N:
		prediction_unit(dst, st, size_log2, unavail, offset_x, offset_y, len, len, left, top);
		return false;
		break;
	case PART_2NxN:
		len_s = len >> 1;
		prediction_unit(dst, st, size_log2, unavail, offset_x, offset_y, len, len_s, left, top);
		prediction_unit(dst, st, size_log2, unavail, offset_x, offset_y + len_s, len, len_s, left, top);
		break;
	case PART_Nx2N:
		len_s = len >> 1;
		prediction_unit(dst, st, size_log2, unavail, offset_x, offset_y, len_s, len, left, top);
		prediction_unit(dst, st, size_log2, unavail, offset_x + len_s, offset_y, len_s, len, left, top);
		break;
	case PART_NxN:
		len_s = len >> 1;
		prediction_unit(dst, st, size_log2, unavail, offset_x, offset_y, len_s, len_s, left, top);
		prediction_unit(dst, st, size_log2, unavail, offset_x + len_s, offset_y, len_s, len_s, left, top);
		prediction_unit(dst, st, size_log2, unavail, offset_x, offset_y + len_s, len_s, len_s, left, top);
		prediction_unit(dst, st, size_log2, unavail, offset_x + len_s, offset_y + len_s, len_s, len_s, left, top);
		break;
	case PART_2NxnU:
		len_s = len >> 2;
		prediction_unit(dst, st, size_log2, unavail, offset_x, offset_y, len, len_s, left, top);
		prediction_unit(dst, st, size_log2, unavail, offset_x, offset_y + len_s, len, len - len_s, left, top);
		break;
	case PART_2NxnD:
		len_s = len >> 2;
		prediction_unit(dst, st, size_log2, unavail, offset_x, offset_y, len, len - len_s, left, top);
		prediction_unit(dst, st, size_log2, unavail, offset_x, offset_y + len - len_s, len, len_s, left, top);
		break;
	case PART_nLx2N:
		len_s = len >> 2;
		prediction_unit(dst, st, size_log2, unavail, offset_x, offset_y, len_s, len, left, top);
		prediction_unit(dst, st, size_log2, unavail, offset_x + len_s, offset_y, len - len_s, len, left, top);
		break;
	case PART_nRx2N:
		len_s = len >> 2;
		prediction_unit(dst, st, size_log2, unavail, offset_x, offset_y, len - len_s, len, left, top);
		prediction_unit(dst, st, size_log2, unavail, offset_x + len - len_s, offset_y, len_s, len, left, top);
		break;
	}
	return true;
}

static void cu_header_intra(h265d_ctu_t& dst, dec_bits& st, int size_log2, h265d_neighbour_t* left, h265d_neighbour_t* top) {
	int part_num = 1;
	dst.intra_split = 0;
	if ((dst.sps->ctb_info.size_log2_min == size_log2) && (part_mode_intra(dst.cabac, st) == 0)) {
		dst.intra_split = 1;
		part_num = 4;
	}
	if ((dst.sps->ctb_info.pcm_log2_min <= size_log2) && (size_log2 <= dst.sps->ctb_info.pcm_log2)) {
		assert(0);
	}
	uint32_t pred_flag = 0;
	for (int i = 0; i < part_num; ++i) {
		pred_flag |= prev_intra_luma_pred_flag(dst.cabac, st) << i;
	}
	uint32_t neighbour_num = 1 << (size_log2 - 2 - (part_num == 4));
	for (int i = 0; i < part_num; ++i) {
		uint8_t candidate[3];
		h265d_neighbour_t* left_tmp = left + (i >> 1);
		h265d_neighbour_t* top_tmp = top + (i & 1);
		intra_pred_candidate(candidate, left_tmp->pred_mode, top_tmp->pred_mode);
		uint8_t mode;
		if (pred_flag & 1) {
			mode = candidate[mpm_idx(dst.cabac, st)];
		} else {
			mode = rem_intra_luma_pred_mode(dst.cabac, st, candidate);
		}
		dst.order_luma[i] = mode;
		pred_flag >>= 1;
		cu_pred_mode_fill(left_tmp, top_tmp, mode, neighbour_num);
	}
	if (part_num != 4) {
		memset(dst.order_luma + 1, dst.order_luma[0], sizeof(dst.order_luma[0]) * 3);
	}
	uint32_t chroma_mode_idx = intra_chroma_pred_mode(dst.cabac, st);
	for (int i = 0; i < 1; ++i) {
		dst.order_chroma = intra_chroma_pred_dir(chroma_mode_idx, dst.order_luma[i]);
	}
}

static void pred_inter(h265d_ctu_t& dst, dec_bits& st, int size_log2, uint32_t unavail, int offset_x, int offset_y, int valid_x, int valid_y, h265d_neighbour_t* left, h265d_neighbour_t* top) {
	uint32_t skip = cu_skip_flag(dst.cabac, st, unavail, left, top);
	cu_pred_mode_fill(left, top, skip, 1 << (size_log2 - 2));
	if (skip) {
		int len = 1 << size_log2;
		prediction_unit_merge(dst, st, 0, 0, len, len);
	} else {
		if (pred_mode_flag(dst.cabac, st) != 0) {
			cu_header_intra(dst, st, size_log2, left, top);
			return;
		}
		bool rqt_root_cbf_needed = prediction_unit_cases(dst, st, size_log2, unavail, offset_x, offset_y, valid_x, valid_y, left, top);
		if (!rqt_root_cbf_needed || rqt_root_cbf(dst.cabac, st)) {
			transform_tree(dst, st, size_log2, unavail, 0, 3, offset_x, valid_x, offset_y, valid_y, 0, 0, false);
		}
	}
}

static void coding_unit_header(h265d_ctu_t& dst, dec_bits& st, int size_log2, uint32_t unavail, h265d_neighbour_t* left, h265d_neighbour_t* top) {
	intra_depth_fill(left, top, size_log2);
	if (dst.pps->cu_qp_delta_enabled_flag) {
		dst.qp_delta_req = 1;
	}
	if (dst.pps->transquant_bypass_enabled_flag) {
		assert(0);
	}
}

static inline bool frame_boundary(int size_log2, int valid_x, int valid_y) {
	return (valid_x < (1 << size_log2)) || (valid_y < (1 << size_log2));
}

static void quad_tree(h265d_ctu_t& dst, dec_bits& st, int size_log2, uint32_t unavail, int offset_x, int valid_x, int offset_y, int valid_y, h265d_neighbour_t* left, h265d_neighbour_t* top) {
	if ((valid_x <= 0) || (valid_y <= 0)) {
		return;
	}
	if ((dst.sps->ctb_info.size_log2_min < size_log2) && (frame_boundary(size_log2, valid_x, valid_y) || split_cu_flag(dst.cabac, st, (6 < size_log2 + left->depth) + (6 < size_log2 + top->depth)))) {
		size_log2 -= 1;
		uint32_t block_len = 1 << size_log2;
		uint32_t info_offset = 1 << (size_log2 - 2);
		quad_tree(dst, st, size_log2, unavail, offset_x, valid_x, offset_y, valid_y, left, top);
		quad_tree(dst, st, size_log2, (unavail & ~1) | 4, offset_x + block_len, valid_x - block_len, offset_y, std::min(static_cast<uint32_t>(valid_y), block_len), left, top + info_offset);
		quad_tree(dst, st, size_log2, unavail & ~10, offset_x, std::min(static_cast<uint32_t>(valid_x), block_len * 2), offset_y + block_len, valid_y - block_len, left + info_offset, top);
		quad_tree(dst, st, size_log2, 12, offset_x + block_len, std::min(static_cast<uint32_t>(valid_x - block_len), block_len), offset_y + block_len, std::min(static_cast<uint32_t>(valid_y - block_len), block_len), left + info_offset, top + info_offset);
	} else {
		coding_unit_header(dst, st, size_log2, unavail, left, top);
		if (dst.slice_header->body.slice_type < 2) {
			pred_inter(dst, st, size_log2, unavail, offset_x, offset_y, valid_x, valid_y, left, top);
		} else {
			cu_header_intra(dst, st, size_log2, left, top);
			transform_tree(dst, st, size_log2, unavail, 0, 3, offset_x, valid_x, offset_y, valid_y, 0, 0, true);
		}
	}
}

static inline int clip2(int val, int lim) {
	return (val < 0) ? 0 : ((lim < val) ? lim : val);
}

static inline bool dsam0(int dpq0, int p3, int p0, int q0, int q3, int beta, int tc) {
	int dpq = dpq0 * 2;
	if ((beta >> 2) <= dpq) {
		return false;
	}
	if (((5 * tc + 1) >> 1) <= std::abs(p0 - q0)) {
		return false;
	}
	if ((beta >> 3) <= std::abs(p3 - p0) + std::abs(q0 - q3)) {
		return false;
	}
	return true;
}

static inline int deblock_weakpq1flag(int dp0, int dp3, int dq0, int dq3, int beta) {
	int dp = dp0 + dp3;
	int dq = dq0 + dq3;
	int beta2 = (beta + (beta >> 1)) >> 3;
	return (dp < beta2) * 2 + (dq < beta2);
}

static inline int clip3delta(int dif, int lim) {
	return (-lim <= dif) ? ((dif <= lim) ? dif : lim) : -lim;
}

static inline int clip3tc(int src, int dst, int tc) {
	return src + clip3delta(dst - src, tc);
}

static inline int deblock_weakpq1(int p0, int p1, int p2, int delta, int tc) {
	int d1 = p1 + clip3delta((((p2 + p0 + 1) >> 1) - p1 + delta) >> 1, tc >> 1);
	return CLIP255C(d1);
}

static inline void deblock_filter1_line(uint8_t* dst, int tc, int depq, int xofs) {
	int p1 = dst[xofs * 2];
	int p0 = dst[xofs * 3];
	int q0 = dst[xofs * 4];
	int q1 = dst[xofs * 5];
	int delta = (9 * (q0 - p0) - 3 * (q1 - p1) + 8) >> 4;
	if (std::abs(delta) < tc * 10) {
		delta = clip3delta(delta, tc);
		dst[xofs * 3] = CLIP255C(p0 + delta);
		dst[xofs * 4] = CLIP255C(q0 - delta);
		if (depq & 2) {
			dst[xofs * 2] = deblock_weakpq1(p0, p1, dst[xofs * 1], delta, tc);
		}
		if (depq & 1) {
			dst[xofs * 5] = deblock_weakpq1(q0, q1, dst[xofs * 6], -delta, tc);
		}
	}
}

static inline void deblock_filter1(uint8_t* dst, int tc, int depq, int stride, int xofs) {
	for (int y = 0; y < 4; ++y) {
		deblock_filter1_line(dst, tc, depq, xofs);
		dst += stride;
	}
}

static inline void deblock_filter2_line(uint8_t* dst, int tc2, int xofs) {
	int p2 = dst[xofs * 1];
	int p1 = dst[xofs * 2];
	int p0 = dst[xofs * 3];
	int q0 = dst[xofs * 4];
	dst[xofs * 1] = clip3tc(p2, (2 * dst[0] + 3 * p2 + p1 + p0 + q0 + 4) >> 3, tc2);
	dst[xofs * 2] = clip3tc(p1, (p2 + p1 + p0 + q0 + 2) >> 2, tc2);
	int q1 = dst[xofs * 5];
	dst[xofs * 3] = clip3tc(p0, (p2 + 2 * p1 + 2 * p0 + 2 * q0 + q1 + 4) >> 3, tc2);
	int q2 = dst[xofs * 6];
	dst[xofs * 4] = clip3tc(q0, (p1 + 2 * p0 + 2 * q0 + 2 * q1 + q2 + 4) >> 3, tc2);
	dst[xofs * 5] = clip3tc(q1, (p0 + q0 + q1 + q2 + 2) >> 2, tc2);
	dst[xofs * 6] = clip3tc(q2, (p0 + q0 + q1 + 3 * q2 + 2 * dst[xofs * 7] + 4) >> 3, tc2);
}

static inline void deblock_filter2(uint8_t* dst, int tc, int stride, int xofs) {
	int tc2 = tc * 2;
	for (int y = 0; y < 4; ++y) {
		deblock_filter2_line(dst, tc2, xofs);
		dst += stride;
	}
}

static const int8_t q_thr[54 - 16][2] = {
	{6, 0}, {7, 0}, {8, 1}, {9, 1}, {10, 1}, {11, 1}, {12, 1}, {13, 1},
	{14, 1}, {15, 1}, {16, 1}, {17, 2}, {18, 2}, {20, 2}, {22, 2}, {24, 3},
	{26, 3}, {28, 3}, {30, 3}, {32, 4}, {34, 4}, {36, 4}, {38, 5}, {40, 5},
	{42, 6}, {44, 6}, {46, 7}, {48, 8}, {50, 9}, {52, 10}, {54, 11}, {56, 13},
	{58, 14}, {60, 16}, {62, 18}, {64, 20}, {64, 22}, {64, 24}
};

static inline void deblock_edge_luma(uint8_t* dst, int beta_qpy, int tc_qpy, int xofs, int stride) {
	const uint8_t* src3 = dst + stride * 3;
	int dp0 = dst[xofs * 1] - 2 * dst[xofs * 2] + dst[xofs * 3];
	int dq0 = dst[xofs * 4] - 2 * dst[xofs * 5] + dst[xofs * 6];
	int dp3 = src3[xofs * 1] - 2 * src3[xofs * 2] + src3[xofs * 3];
	int dq3 = src3[xofs * 4] - 2 * src3[xofs * 5] + src3[xofs * 6];
	dp0 = std::abs(dp0);
	dp3 = std::abs(dp3);
	dq0 = std::abs(dq0);
	dq3 = std::abs(dq3);
	int dpq0 = dp0 + dq0;
	int dpq3 = dp3 + dq3;
	int d = dpq0 + dpq3;
	int beta = q_thr[beta_qpy][0];
	if (d < beta) {
		int tc = q_thr[tc_qpy][1];
		if (dsam0(dpq0, dst[0], dst[xofs * 3], dst[xofs * 4], dst[xofs * 7], beta, tc) && dsam0(dpq3, src3[0], src3[xofs * 3], src3[xofs * 4], src3[xofs * 7], beta, tc)) {
			deblock_filter2(dst, tc, stride, xofs);
		} else {
			deblock_filter1(dst, tc, deblock_weakpq1flag(dp0, dp3, dq0, dq3, beta), stride, xofs);
		}
	}
}

static inline void deblocking_edge_luma_block(const h265d_deblocking_strength_t& strength, int beta_offset, int tc_offset, uint8_t* dst, int xstride, int ystride) {
	int str = strength.str;
	if (str != 0) {
		int qp = strength.qp;
		int beta_qp = (beta_offset ? clip2(qp + beta_offset, 51) : qp) - 16;
		if (0 <= beta_qp) {
			int ofs = tc_offset + (str & 2);
			int tc_qp = (ofs ? clip2(qp + ofs, 51) : qp) - 16;
			if (0 <= tc_qp) {
				deblock_edge_luma(dst, beta_qp, tc_qp, xstride, ystride);
			}
		}
	}
}

static void deblocking_vert_edge_luma_lines(const h265d_deblocking_strength_t str[], int beta_offset, int tc_offset, uint8_t* dst, int stride, int blknum, int edgenum) {
	for (int y = 0; y < blknum; ++y) {
		for (int x = 0; x < edgenum; ++x) {
			deblocking_edge_luma_block(str[x], beta_offset, tc_offset, dst + x * 8, 1, stride);
		}
		str += edgenum;
		dst += stride * 4;
	}
}

static void deblocking_horiz_edge_luma_lines(const h265d_deblocking_strength_t str[], int beta_offset, int tc_offset, uint8_t* dst, int stride, int blknum, int edgenum) {
	for (int y = 0; y < edgenum; ++y) {
		for (int x = 0; x < blknum; ++x) {
			deblocking_edge_luma_block(str[x], beta_offset, tc_offset, dst + x * 4, stride, 1);
		}
		str = str + edgenum * 2 + 1;
		dst += stride * 8;
	}
}

static inline int qpi_to_qpc_deb(int qpi) {
	static const int8_t qpcadj_12[] = {
		-12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1,
		0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
		16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 29, 30,
		31, 32, 33, 33, 34, 34, 35, 35, 36, 36, 37, 37, 38, 39, 40, 41,
		42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59
	};
	return qpcadj_12[qpi + 12];
}

static inline int tc_qpc(int qp, int qpc_offset, int tc_offset) {
	qp = qpi_to_qpc_deb(qp + qpc_offset);
	return clip2(qp + 2 + tc_offset, 53);
}

static inline void deblocking_edge_chroma_block(int qp, int pps_offset, int tc_offset, uint8_t* dst, int xofs, int stride) {
	qp = tc_qpc(qp, pps_offset, tc_offset) - 16;
	if (qp < 0) {
		return;
	}
	int tc = q_thr[qp][1];
	for (int x = 0; x < 2; ++x) {
		int p1 = dst[0];
		int p0 = dst[xofs * 1];
		int q0 = dst[xofs * 2];
		int q1 = dst[xofs * 3];
		int delta = clip3delta(((q0 - p0) * 4 + p1 - q1 + 4) >> 3, tc);
		if (delta) {
			dst[xofs * 1] = CLIP255C(p0 + delta);
			dst[xofs * 2] = CLIP255C(q0 - delta);
		}
		dst += stride;
	}
}

static void deblocking_vert_edge_chroma_lines(const h265d_deblocking_strength_t str[], int cb_offset, int cr_offset, int tc_offset, uint8_t* dst, int stride, int blknum, int edgenum) {
	for (int y = 0; y < blknum; ++y) {
		for (int x = 0; x < edgenum; ++x) {
			if (str[x * 2].str == 2) {
				int qp = str[x * 2].qp;
				deblocking_edge_chroma_block(qp, cb_offset, tc_offset, dst + x * 16, 2, stride);
				deblocking_edge_chroma_block(qp, cr_offset, tc_offset, dst + x * 16 + 1, 2, stride);
			}
		}
		str += edgenum * 2;
		dst += stride * 2;
	}
}

static void deblocking_horiz_edge_chroma_lines(const h265d_deblocking_strength_t str[], int cb_offset, int cr_offset, int tc_offset, uint8_t* dst, int stride, int blknum, int edgenum) {
	for (int y = 0; y < edgenum; ++y) {
		for (int x = 0; x < blknum; ++x) {
			if (str[x].str == 2) {
				int qp = str[x].qp;
				deblocking_edge_chroma_block(qp, cb_offset, tc_offset, dst + x * 4, stride, 2);
				deblocking_edge_chroma_block(qp, cr_offset, tc_offset, dst + x * 4 + 1, stride, 2);
			}
		}
		dst += stride * 8;
		str = str + (edgenum * 4 + 1) * 2;
	}
}

static void pre_deblocking_strength(h265d_ctu_t& ctu, int edgenum) {
	memcpy(ctu.deblock_boundary[0], ctu.deblock_topedge, edgenum * sizeof(ctu.deblock_topedge[0]));
}

static inline void clear_deblocking_strength_left(h265d_ctu_t& ctu, int edgenum) {
	h265d_deblocking_strength_t* left = &ctu.deblock_boundary[1][0];
	int len = edgenum * 2;
	for (int n = 0; n < edgenum; ++n) {
		left[0] = left[len];
		memset(left + 1, 0, len * sizeof(left[0]));
		left = left + len + 1;
	}
}

static void post_deblocking_strength(h265d_ctu_t& ctu, int edgenum) {
	if (ctu.pos_x < ctu.sps->ctb_info.columns - 1) {
		clear_deblocking_strength_left(ctu, edgenum);
	} else {
		memset(ctu.deblock_boundary[1], 0, sizeof(ctu.deblock_boundary[1]));
	}
	int store_size = edgenum * sizeof(ctu.deblock_topedge[0]);
	memcpy(ctu.deblock_topedge, ctu.deblock_boundary[0] + edgenum * edgenum * 2, store_size);
	memset(ctu.deblock_boundary[0] + store_size, 0, sizeof(ctu.deblock_boundary[0]) - store_size);
}

static void sanityedge(const h265d_deblocking_strength_t* bound, int inc, int edgenum, int stride) {
	for (int e = 0; e < edgenum; ++e) {
		if ((bound[e * inc].str != 0) || (bound[e * inc + stride].str != 0)) {
			throw -1;
		}
	}
}

static void sanity(const h265d_ctu_t& ctu) {
#if 0
	int edgenum = 1 << (ctu.sps->ctb_info.size_log2 - 3);
	if (ctu.pos_y == 0) {
		sanityedge(ctu.deblock_boundary[0], 1, edgenum, 0);
	}
	if (ctu.pos_x == 0) {
		sanityedge(ctu.deblock_boundary[1], edgenum * 2 + 2, edgenum, 0);
	}
#endif
}

static void deblock_ctu(h265d_ctu_t& ctu) {
	if (ctu.slice_header->body.deblocking_filter_disabled_flag) {
		return;
	}
	int edgenum = 1 << (ctu.sps->ctb_info.size_log2 - 3);
	sanity(ctu);
	pre_deblocking_strength(ctu, edgenum);
	sanity(ctu);
	int stride = ctu.sps->ctb_info.stride;
	int beta_offset = ctu.slice_header->body.slice_beta_offset_div2 * 2;
	int tc_offset = ctu.slice_header->body.slice_tc_offset_div2 * 2;
	uint8_t* luma = ctu.luma - stride * 4 - 4;
	deblocking_vert_edge_luma_lines(ctu.deblock_boundary[0], beta_offset, tc_offset, luma, stride, edgenum * 2 + (ctu.pos_y == ctu.sps->ctb_info.rows - 1), edgenum);
	deblocking_horiz_edge_luma_lines(ctu.deblock_boundary[1], beta_offset, tc_offset, luma, stride, edgenum * 2 + (ctu.pos_x == ctu.sps->ctb_info.columns - 1), edgenum);
	uint8_t* chroma = ctu.chroma - stride * 2 - 4;
	int cb_offset = ctu.pps->pps_cb_qp_offset;
	int cr_offset = ctu.pps->pps_cr_qp_offset;
	deblocking_vert_edge_chroma_lines(ctu.deblock_boundary[0], cb_offset, cr_offset, tc_offset, chroma, stride, edgenum * 2 + (ctu.pos_y == ctu.sps->ctb_info.rows - 1), edgenum >> 1);
	deblocking_horiz_edge_chroma_lines(ctu.deblock_boundary[1], cb_offset, cr_offset, tc_offset, chroma, stride, edgenum * 2 + (ctu.pos_x == ctu.sps->ctb_info.columns - 1), edgenum >> 1);
	post_deblocking_strength(ctu, edgenum);
	VC_CHECK;
}

static inline int signe(int x) {
	return ((x >> (sizeof(int) * 8 - 1)) & 2) | ((x + 65535) >> 16);
}

static inline uint8_t set_sign_elem(uint8_t dst, int pos, int sign) {
	static const uint8_t mask[] = {
		0xfc, 0xf3, 0xcf, 0x3f
	};
	pos &= 3;
	return (dst & mask[pos]) | (sign << (pos * 2));
}

static inline void set_sign(uint8_t* dst, int pos, int sign) {
	dst += pos >> 2;
	*dst = set_sign_elem(*dst, pos, sign);
}

/* Assume pos are to be incremented sequentially. */
static inline int set_sign_buffered(uint8_t* dst, int pos, int sign, int tmp) {
	if ((pos & 3) == 0) {
		dst[(pos - 1) >> 2] = tmp;
	}
	return set_sign_elem(tmp, pos, sign);
}

static inline int get_sign(const uint8_t* dst, int pos) {
	return (dst[pos >> 2] >> ((pos & 3) * 2)) & 3;
}

static inline int sao_eo_index(int sign0, int sign2) {
	static const int8_t idx[16] = {
		-1, 2, 1, -1, 2, 3, -1, 2, 1, -1, 0, 1, -1, 2, 1, -1
	};
	return idx[sign0 * 4 + sign2];
}

static inline void sao_edge0_line(uint8_t* dst, int width, const int8_t* offset, int xgap) {
	int d1 = dst[0];
	int sign0 = signe(d1 - dst[-xgap]);
	for (int x = 0; x < width; ++x) {
		int d2 = dst[(x + 1) * xgap];
		int sign2 = signe(d1 - d2);
		int idx = sao_eo_index(sign2, sign0);
		if (0 <= idx) {
			dst[x * xgap] = clip2(d1 + offset[idx], 255);
		}
		d1 = d2;
		sign0 = sign2 ^ 3;
	}
}

static inline void sao_edge0(uint8_t* dst, const int8_t offset[], int width, int height, int xgap, int stride, uint32_t unavail) {
	if (unavail & 1) {
		dst += xgap;
		width -= 1;
	}
	if (unavail & 4) {
		width -= 1;
	}
	for (int y = 0; y < height; ++y) {
		sao_edge0_line(dst, width, offset, xgap);
		dst += stride;
	}
}

static inline void fill_sign(const uint8_t* s0, const uint8_t* s1, int len, int xgap, uint8_t* dst) {
	for (int x = 0; x < len; ++x) {
		set_sign(dst, x, signe(s1[x * xgap] - s0[x * xgap]));
	}
}

/*
  XDELTA = -1: left-top to right-bottom
  XDELTA = +1: right-top to left-bottom
 */
template <int XDELTA>
static inline void sao_diag_edge(uint8_t* dst, const int8_t offset[], int width, int height, uint8_t* signbuf, int xgap, int stride, uint32_t unavail) {
	if (XDELTA) {
		if (unavail & 1) {
			dst += xgap;
			width -= 1;
		}
		if (unavail & 4) {
			width -= 1;
		}
	}
	if (unavail & 2) {
		dst += stride;
		height -= 1;
	}
	if (unavail & 8) {
		height -= 1;
	}
	fill_sign(dst - stride + xgap * XDELTA, dst, width, xgap, signbuf);
	for (int y = 0; y < height; ++y) {
		const uint8_t* s2 = dst + stride - xgap * XDELTA;
		int tmp = 0;
		for (int x = 0; x < width; ++x) {
			int d0 = dst[x * xgap];
			int sign0 = get_sign(signbuf, x);
			int sign2 = signe(d0 - s2[x * xgap]);
			int idx = sao_eo_index(sign2, sign0);
			if (0 <= idx) {
				dst[x * xgap] = clip2(d0 + offset[idx], 255);
			}
			if (!XDELTA || ((0 < XDELTA) && (x != 0))) {
				set_sign(signbuf, x - XDELTA, sign2 ^ 3);
			} else if (XDELTA < 0) {
				tmp = set_sign_buffered(signbuf, x - XDELTA, sign2 ^ 3, tmp);
			}
		}
		if (XDELTA < 0) {
			if ((width & 3) != 0) {
				signbuf[width >> 2] = tmp;
			}
			set_sign(signbuf, 0, signe(s2[-xgap] - dst[-xgap]));
		} else if (0 < XDELTA) {
			set_sign(signbuf, width - 1, signe(s2[width * xgap] - dst[width * xgap]));
		}
		dst = dst + stride;
	}
}

static void sao_eo_block(const h265d_sao_map_t& sao_map, int cidx, uint8_t* dst, h265d_ctu_t& ctu, int width, int height, int xgap, int stride, uint32_t unavail) {
	int edge = sao_map.elem[(cidx + 1) >> 1].opt.edge;
	const int8_t* offset = sao_map.elem[cidx].offset;
	switch (edge) {
	case 0:
		sao_edge0(dst, offset, width, height, xgap, stride, unavail);
		break;
	case 1:
		sao_diag_edge<0>(dst, offset, width, height, ctu.sao_signbuf, xgap, stride, unavail);
		break;
	case 2:
		sao_diag_edge<-1>(dst, offset, width, height, ctu.sao_signbuf, xgap, stride, unavail);
		break;
	case 3:
		sao_diag_edge<1>(dst, offset, width, height, ctu.sao_signbuf, xgap, stride, unavail);
		break;
	}
}

static void sao_bo_block(const h265d_sao_map_t& sao_map, int cidx, uint8_t* dst, int width, int height, int xgap, int stride) {
	int bandtop = sao_map.elem[cidx].opt.band_pos;
	const int8_t* offset = sao_map.elem[cidx].offset;
	int shift = sizeof(*dst) * 8 - 5;
	unsigned lim = 4U << shift;
	bandtop <<= shift;
	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			int d0 = dst[x * xgap];
			int dif = d0 - bandtop;
			if ((unsigned)dif < lim) {
				dst[x * xgap] = clip2(d0 + offset[dif >> shift], 255);
			}
		}
		dst += stride;
	}
}

struct swap_hline {
	void operator()(uint8_t* a, uint8_t* b, int len) const {
		uint32_t* a32 = reinterpret_cast<uint32_t*>(a);
		uint32_t* b32 = reinterpret_cast<uint32_t*>(b);
		int cnt = len >> 2;
		for (int x = 0; x < cnt; ++x) {
			uint32_t tmp = a32[x];
			a32[x] = b32[x];
			b32[x] = tmp;
		}
	}
};

template <typename F>
static void sao_swap_hline(const uint8_t* record_bits, uint8_t* src, uint8_t* mod, int xpos, int xnum, int width, F Modify) {
	record_bits += xpos >> 3;
	int bits = *record_bits++;
	for (int x = 0; x < xnum; ++x) {
		if (bits & (1 << (xpos++ & 7))) {
			Modify(src + x * width, mod + x * width, width);
		}
		if ((xpos & 7) == 0) {
			bits = *record_bits++;
		}
	}
}

static inline void write_hline(uint8_t* dst, const uint8_t* src, int len) {
	const uint32_t* a32 = reinterpret_cast<const uint32_t*>(src);
	uint32_t* b32 = reinterpret_cast<uint32_t*>(dst);
	int cnt = len >> 2;
	for (int x = 0; x < cnt; ++x) {
		b32[x] = a32[x];
	}
}

static void sao_reserve_hline(h265d_sao_hlines_t& lines, const uint8_t* src, int xpos, int xnum, int width) {
	uint8_t* record_bits = lines.reserved_flag;
	write_hline(lines.bottom + xpos * width, src, xnum * width);
	record_bits += xpos >> 3;
	int xend = xpos + xnum;
	int bits = *record_bits;
	for (int x = xpos; x < xend; ++x) {
		bits |= 1 << (x & 7);
		if ((x & 7) == 7) {
			*record_bits++ = bits;
			bits = 0;
		}
	}
	if (xend & 7) {
		*record_bits = bits;
	}
}

static inline int vline_flag(int phase, int cidx) {
	return 1 << (cidx * 2 + (phase & 1));
}

template <typename T>
static void sao_swap_vline(h265d_sao_vlines_t& lines, int phase, int cidx, T* src0, T* src1, int stride0, int stride1, int height) {
	if (lines.reserved_flag & vline_flag(phase, cidx)) {
		for (int y = 0; y < height; ++y) {
			T tmp = src0[y * stride0];
			src0[y * stride0] = src1[y * stride1];
			src1[y * stride1] = tmp;
		}
	}
}

template <typename T>
static void sao_reserve_vline(h265d_sao_vlines_t& lines, int phase, int cidx, const T* src, int height, int stride) {
	lines.reserved_flag |= vline_flag(phase, cidx);
	memcpy_from_vertical(reinterpret_cast<T*>(lines.right[phase & 1][cidx]), reinterpret_cast<const T*>(src), height, stride);
}

static void sao_clear_vline(h265d_sao_vlines_t& lines, int phase, int cidx) {
	lines.reserved_flag &= ~vline_flag(phase, cidx);
}

static inline int sao_merged_num(const h265d_sao_map_t* sao_map, int max) {
	int i;
	for (i = 1; i < max; ++i) {
		if (!sao_map[i].merge_left) {
			break;
		}
	}
	return i;
}

static int sao_region(h265d_ctu_t& ctu, const h265d_sao_map_t* sao_map, int width, int height, int stride, uint32_t unavail, int max, int xphase, int pos_y, int valid_width) {
	int x_num;
	sao_clear_vline(ctu.sao_vlines, xphase ^ 1, 0);
	sao_clear_vline(ctu.sao_vlines, xphase ^ 1, 1);
	x_num = sao_merged_num(sao_map, max);
	int hlen = std::min(width * x_num, valid_width);
	int idx = sao_map[0].luma_idx;
	if (idx) {
		if (x_num < max) {
			if (sao_map[x_num].luma_idx == 2) {
				sao_reserve_vline(ctu.sao_vlines, xphase ^ 1, 0, ctu.luma + hlen - 1, height, stride);
			}
		} else {
			unavail |= 4;
		}
		sao_reserve_hline(ctu.sao_hlines[(pos_y ^ 1) & 1][0], ctu.luma + stride * (height - 1), ctu.pos_x, x_num, width);
		if (idx == 1) {
			sao_bo_block(sao_map[0], 0, ctu.luma, hlen, height, 1, stride);
		} else {
			sao_swap_vline(ctu.sao_vlines, xphase, 0, ctu.sao_vlines.right[xphase & 1][0], ctu.luma - 1, 1, stride, height);
			sao_eo_block(sao_map[0], 0, ctu.luma, ctu, hlen, height, 1, stride, unavail);
			sao_swap_vline(ctu.sao_vlines, xphase, 0, ctu.luma - 1, ctu.sao_vlines.right[xphase & 1][0], stride, 1, height);
		}
	}
	idx = sao_map[0].chroma_idx;
	if (idx) {
		if (x_num < max) {
			if (sao_map[x_num].chroma_idx == 2) {
				sao_reserve_vline(ctu.sao_vlines, xphase ^ 1, 1, reinterpret_cast<uint16_t*>(ctu.chroma + hlen - 2), height >> 1, stride);
			}
		} else {
			unavail |= 4;
		}
		sao_reserve_hline(ctu.sao_hlines[(pos_y ^ 1) & 1][1], ctu.chroma + stride * ((height >> 1) - 1), ctu.pos_x, x_num, width);
		if (idx == 1) {
			sao_bo_block(sao_map[0], 1, ctu.chroma, hlen >> 1, height >> 1, 2, stride);
			sao_bo_block(sao_map[0], 2, ctu.chroma + 1, hlen >> 1, height >> 1, 2, stride);
		} else if (1 < idx) {
			sao_swap_vline(ctu.sao_vlines, xphase, 1, reinterpret_cast<uint16_t*>(ctu.sao_vlines.right[xphase & 1][1]), reinterpret_cast<uint16_t*>(ctu.chroma - 2), 1, stride >> 1, height >> 1);
			sao_eo_block(sao_map[0], 1, ctu.chroma, ctu, hlen >> 1, height >> 1, 2, stride, unavail);
			sao_eo_block(sao_map[0], 2, ctu.chroma + 1, ctu, hlen >> 1, height >> 1, 2, stride, unavail);
			sao_swap_vline(ctu.sao_vlines, xphase, 1, reinterpret_cast<uint16_t*>(ctu.chroma - 2), reinterpret_cast<uint16_t*>(ctu.sao_vlines.right[xphase & 1][1]), stride >> 1, 1, height >> 1);
		}
	}
	return x_num;
}

static void sao_oneframe(h265d_ctu_t& ctu) {
	if (!ctu.slice_header->body.slice_sao_luma_flag && !ctu.slice_header->body.slice_sao_chroma_flag) {
		return;
	}
	int ymax = ctu.sps->ctb_info.rows;
	int xmax = ctu.sps->ctb_info.columns;
	int stride = ctu.sps->ctb_info.stride;
	int len = 1 << ctu.sps->ctb_info.size_log2;
	int width = ctu.sps->pic_width_in_luma_samples;
	const h265d_sao_map_t* sao_map = ctu.sao_map;
	uint32_t unavail = 3;
	for (int y = 0; y < ymax; ++y) {
		ctu.pos_y = y;
		ctu.luma = ctu.frame_info.curr_luma + stride * len * y;
		ctu.chroma = ctu.frame_info.curr_chroma + stride * (len >> 1) * y;
		if (y != 0) {
			h265d_sao_hlines_t* prev = ctu.sao_hlines[y & 1];
			sao_swap_hline(prev[0].reserved_flag, prev[0].bottom, ctu.luma - stride, 0, xmax, len, swap_hline());
			sao_swap_hline(prev[1].reserved_flag, prev[1].bottom, ctu.chroma - stride, 0, xmax, len, swap_hline());
		}
		h265d_sao_hlines_t* next = ctu.sao_hlines[(y ^ 1) & 1];
		memset(next[0].reserved_flag, 0, (sizeof(next[0].reserved_flag[0]) * xmax + 7) >> 3);
		memset(next[1].reserved_flag, 0, (sizeof(next[0].reserved_flag[0]) * xmax + 7) >> 3);
		int vlen = (y < ymax - 1) ? len : ((ctu.sps->pic_height_in_luma_samples - 1) & (len - 1)) + 1;
		int x = 0;
		int phase = 0;
		int valid_width = width;
		ctu.sao_vlines.reserved_flag = 0;
		while (x < xmax) {
			ctu.pos_x = x;
			int run = sao_region(ctu, sao_map + x, len, vlen, stride, unavail, xmax - x, phase++, y, valid_width);
			x += run;
			valid_width -= len * run;
			ctu.luma += len * run;
			ctu.chroma += len * run;
			unavail &= ~1;
		}
		if (y != 0) {
			h265d_sao_hlines_t* prev = ctu.sao_hlines[y & 1];
			ctu.luma = ctu.frame_info.curr_luma + stride * len * y;
			ctu.chroma = ctu.frame_info.curr_chroma + stride * (len >> 1) * y;
			sao_swap_hline(prev[0].reserved_flag, ctu.luma - stride, prev[0].bottom, 0, xmax, len, memcpy);
			sao_swap_hline(prev[1].reserved_flag, ctu.chroma - stride, prev[1].bottom, 0, xmax, len, memcpy);
		}
		sao_map += xmax;
		unavail = (y < ymax - 2) ? 1 : 9;
	}
}

static inline uint32_t get_avail(int xpos, int ypos, int cols, int rows) {
	return (xpos == 0) + (ypos == 0) * 2 + (xpos == cols - 1) * 4 + (ypos == rows - 1) * 8;
}

static void coding_tree_unit(h265d_ctu_t& dst, dec_bits& st) {
	dst.sao_read(dst, *dst.slice_header, st);
	uint32_t idx_in_slice = dst.idx_in_slice;
	uint32_t unavail = (!dst.pos_y || (idx_in_slice < dst.sps->ctb_info.columns)) * 10 + (!dst.pos_x || !idx_in_slice) * 5;
	quad_tree(dst, st, dst.sps->ctb_info.size_log2, unavail, 0, dst.valid_x, 0, dst.valid_y, dst.neighbour_left, dst.neighbour_top + dst.pos_x * NUM_ELEM(dst.neighbour_left));
	deblock_ctu(dst);
}

static inline void neighbour_init(h265d_neighbour_t neighbour[], size_t num) {
	static const h265d_neighbour_t neighbour_init = {
		INTRA_DC, 0
	};
	std::fill(neighbour, neighbour + num, neighbour_init);
}

static void ctu_init(h265d_ctu_t& dst, h265d_data_t& h2d, const h265d_pps_t& pps, const h265d_sps_t& sps) {
	const h265d_slice_header_t& hdr = h2d.slice_header;
	const h265d_slice_header_body_t& header = hdr.body;
	int slice_type = header.slice_type;
	int idc = (slice_type < 2) ? (2 - (slice_type ^ header.cabac_init_flag)) : 0;
	init_cabac_context(&dst.cabac, header.slice_qpy, cabac_initial_value[idc], NUM_ELEM(cabac_initial_value[idc]));
	dst.sao_read = (hdr.body.slice_sao_luma_flag || hdr.body.slice_sao_chroma_flag) ? sao_read : sao_ignore;
	dst.sps = &sps;
	if (sps.scaling_list_enabled_flag) {
		dst.scaling_func = 0;
		if (sps.scaling_list_data_present_flag) {
		}
	} else {
		dst.scaling_func = scaling_default_func;
	}
	dst.coeff_buf = reinterpret_cast<int16_t*>((reinterpret_cast<uintptr_t>(dst.coeff_buffer) + 15) & ~15);
	int ctu_address = hdr.slice_segment_address;
	dst.idx_in_slice = 0;
	dst.pos_y = ctu_address / sps.ctb_info.columns;
	uint32_t y_offset = sps.ctb_info.columns * dst.pos_y;
	dst.pos_x = ctu_address - y_offset;
	dst.valid_x = sps.pic_width_in_luma_samples - (dst.pos_x << sps.ctb_info.size_log2);
	dst.valid_y = std::min(sps.pic_height_in_luma_samples - (dst.pos_y << sps.ctb_info.size_log2), 1U << sps.ctb_info.size_log2);
	uint32_t luma_offset = ((y_offset << sps.ctb_info.size_log2) + dst.pos_x) << sps.ctb_info.size_log2;
	m2d_frame_t& frm = dst.frame_info.frames[dst.frame_info.index];
	frm.width = sps.ctb_info.stride;
	frm.height = sps.ctb_info.rows << sps.ctb_info.size_log2;
	frm.crop[0] = sps.cropping[0];
	frm.crop[1] = sps.cropping[1] + frm.width - sps.pic_width_in_luma_samples;
	frm.crop[2] = sps.cropping[2];
	frm.crop[3] = sps.cropping[3] + frm.height - sps.pic_height_in_luma_samples;
	dst.luma = frm.luma + luma_offset;
	dst.chroma = frm.chroma + (luma_offset >> 1);
	dst.deblock_topedge = dst.deblock_topedgebase + (dst.pos_x << (sps.ctb_info.size_log2 - 3));
	dst.pps = &pps;
	dst.slice_header = &hdr;
	if (dst.qpy != header.slice_qpy) {
		dst.qpy = header.slice_qpy;
		qp_to_scale(dst.qp_scale, header.slice_qpy, header.slice_qpc_delta);
		dst.qpc_delta[0] = header.slice_qpc_delta[0];
		dst.qpc_delta[1] = header.slice_qpc_delta[1];
	}
	neighbour_init(dst.neighbour_left, NUM_ELEM(dst.neighbour_left));
	size_t neighbour_flags_top_len = sps.ctb_info.columns * NUM_ELEM(dst.neighbour_left);
	neighbour_init(dst.neighbour_top, neighbour_flags_top_len);
	memset(dst.qp_history[0], dst.qpy, sizeof(dst.qp_history));
	memset(dst.deblock_boundary, 0, sizeof(dst.deblock_boundary));
	memset(dst.deblock_topedgebase, 0, sizeof(dst.deblock_topedgebase[0]) * (sps.ctb_info.columns << (sps.ctb_info.size_log2 - 3)));
}

static uint32_t ctu_pos_increment(h265d_ctu_t& dst) {
	uint32_t pos_x = dst.pos_x + 1;
	uint32_t size_log2 = dst.sps->ctb_info.size_log2;
	if (dst.sps->ctb_info.columns <= pos_x) {
		neighbour_init(dst.neighbour_left, NUM_ELEM(dst.neighbour_left));
		dst.pos_y++;
		dst.valid_x = dst.sps->pic_width_in_luma_samples;
		if (dst.pos_y == dst.sps->ctb_info.rows - 1) {
			dst.valid_y = std::min(dst.sps->pic_height_in_luma_samples - (dst.pos_y << size_log2), 1U << size_log2);
		}
		pos_x = 0;
		uint32_t y_offset = dst.sps->ctb_info.columns << (size_log2 * 2);
		uint32_t line_offset = dst.sps->ctb_info.columns << size_log2;
		dst.luma = dst.luma + y_offset - line_offset;
		dst.chroma = dst.chroma + (y_offset >> 1) - line_offset;
	} else {
		dst.valid_x -= 1 << size_log2;
	}
	dst.luma += 1 << size_log2;
	dst.chroma += 1 << size_log2;
	dst.deblock_topedge = dst.deblock_topedgebase + (pos_x << (size_log2 - 3));
	dst.pos_x = pos_x;
	dst.idx_in_slice++;
	uint32_t num = NUM_ELEM(dst.neighbour_left);
	h265d_neighbour_t* top = dst.neighbour_top + pos_x * num;
	for (uint32_t i = 0; i < num; ++i) {
		top[i].pred_mode = INTRA_DC;
	}
	memset(dst.qp_history[0], dst.qpy, sizeof(dst.qp_history));
	return (dst.sps->ctb_info.rows <= dst.pos_y);
}

static void slice_data(h265d_ctu_t& dst, h265d_data_t& h2d, const h265d_pps_t& pps, const h265d_sps_t& sps, dec_bits& st) {
	ctu_init(dst, h2d, pps, sps);
	init_cabac_engine(&dst.cabac, &st);
	do {
		coding_tree_unit(dst, st);
		if (ctu_pos_increment(dst)) {
			break;
		}
	} while (!end_of_slice_segment_flag(dst.cabac, st));
}

static void insert_dpb(h265d_dpb_t& dpb, int frame_idx, uint32_t poc, bool is_idr);

static void slice_layer(h265d_data_t& h2d, dec_bits& st) {
	h265d_slice_header_t& header = h2d.slice_header;
	header.body.nal_type = h2d.current_nal;
	if ((header.first_slice_segment_in_pic_flag = get_onebit(&st)) != 0) {
		find_empty_frame(h2d.coding_tree_unit.frame_info);
	}
	if ((BLA_W_LP <= h2d.current_nal) &&  (h2d.current_nal <= RSV_IRAP_VCL23)) {
		header.no_output_of_prior_pics_flag = get_onebit(&st);
	}
	READ_CHECK_RANGE(ue_golomb(&st), header.pps_id, 63, st);
	const h265d_pps_t& pps = h2d.pps[header.pps_id];
	const h265d_sps_t& sps = h2d.sps[pps.sps_id];
	slice_header(header, pps, sps, st);
	slice_data(h2d.coding_tree_unit, h2d, pps, sps, st);
	sao_oneframe(h2d.coding_tree_unit);
	insert_dpb(h2d.coding_tree_unit.frame_info.dpb, h2d.coding_tree_unit.frame_info.index, pic_order_cnt(h2d.slice_header.body.slice_pic_order_cnt, sps), (header.body.nal_type == IDR_W_RADL) || (header.body.nal_type == IDR_N_LP));
}

static int dispatch_one_nal(h265d_data_t& h2d, uint32_t nalu_header) {
	int err = 0;
	dec_bits& st = h2d.stream_i;
	switch (h2d.current_nal = static_cast<h265d_nal_t>((nalu_header >> 9) & 63)) {
	case TRAIL_N:
	case TRAIL_R:
	case IDR_W_RADL:
		slice_layer(h2d, st);
		err = 1;
		break;
	case VPS_NAL:
		video_parameter_set(h2d, st);
		break;
	case SPS_NAL:
		seq_parameter_set(h2d, st);
		h2d.header_callback(h2d.header_callback_arg, st.id);
		break;
	case PPS_NAL:
		pic_parameter_set(h2d, st);
		break;
	case AUD_NAL:
		au_delimiter(h2d, st);
		break;
	default:
		skip_nal(h2d, st);
		break;
	}
	return err;
}

int h265d_decode_picture(h265d_context *h2) {
	if (!h2) {
		return -1;
	}
	h265d_data_t& h2d = *reinterpret_cast<h265d_data_t*>(h2);
	dec_bits* stream = &h2d.stream_i;
	if (setjmp(stream->jmp) != 0) {
		return -2;
	}
/*	h2d->slice_header->first_mb_in_slice = UINT_MAX;*/
	int err = 0;
	uint32_t nalu_header = 0;
	do {
		if (0 <= (err = m2d_find_mpeg_data(stream))) {
			nalu_header = get_bits(stream, 16);
			err = dispatch_one_nal(h2d, nalu_header);
		} else {
			error_report(*stream);
		}
		VC_CHECK;
	} while (err == 0);// || (code_type == SPS_NAL && 0 < err));
	return err;
}

static void insert_dpb(h265d_dpb_t& dpb, int frame_idx, uint32_t poc, bool is_idr) {
	int size = dpb.size;
	int max = dpb.max;
	if (max <= size) {
		size -= 1;
		dpb.output = dpb.data[size].frame_idx;
	} else {
		dpb.output = -1;
	}
	if (0 < size) {
		memmove(&dpb.data[1], &dpb.data[0], sizeof(dpb.data[0]) * size);
	}
	dpb.data[0].frame_idx = frame_idx;
	dpb.data[0].poc = poc;
	dpb.data[0].is_idr = is_idr;
	dpb.size = size + 1;
}

static int peek_decoded_frame(h265d_dpb_t& dpb) {
	int size = dpb.size;
	if (size <= 0) {
		return -1;
	}
	return (0 <= dpb.output) ? dpb.output : dpb.data[0].frame_idx;
}

static void force_pop_dpb(h265d_dpb_t& dpb) {
	int size = dpb.size;
	if (0 < size) {
		memmove(&dpb.data[0], &dpb.data[1], sizeof(dpb.data[0]) * size);
		dpb.size = size - 1;
		dpb.output = -1;
	}
}

static int peek_decoded_frame(h265d_frame_info_t& frm, m2d_frame_t *frame, int bypass_dpb) {
	int idx = peek_decoded_frame(frm.dpb);
	if (idx < 0) {
		return 0;
	}
	*frame = frm.frames[idx];
	return 1;
}

int h265d_peek_decoded_frame(h265d_context *h2, m2d_frame_t *frame, int bypass_dpb)
{
	if (!h2 || !frame) {
		return -1;
	}
	return peek_decoded_frame(reinterpret_cast<h265d_data_t*>(h2)->coding_tree_unit.frame_info, frame, bypass_dpb);
}

int h265d_get_decoded_frame(h265d_context *h2, m2d_frame_t *frame, int bypass_dpb)
{
	int idx = h265d_peek_decoded_frame(h2, frame, bypass_dpb);
	if (idx < 0) {
		return -1;
	}
	force_pop_dpb(reinterpret_cast<h265d_data_t*>(h2)->coding_tree_unit.frame_info.dpb);
	return idx;
}

static const m2d_func_table_t h265d_func_ = {
	sizeof(h265d_data_t),
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
