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
#include <algorithm>
#include "h265.h"

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
	info->additional_size = (sizeof(h2d.coding_tree_unit.sao_map[0]) * sps.ctb_info.rows) * sps.ctb_info.columns
		+ sizeof(h2d.coding_tree_unit.neighbour_left) * sps.ctb_info.columns;
	return 0;
}

int h265d_set_frames(h265d_context *h2, int num_frame, m2d_frame_t *frame, uint8_t *second_frame, int second_frame_size) {
	h265d_data_t& h2d = *reinterpret_cast<h265d_data_t*>(h2);
	if (!h2 || (num_frame < 1) || (NUM_ELEM(h2d.coding_tree_unit.frames) < (size_t)num_frame) || !frame || !second_frame) {
		return -1;
	}
	h265d_ctu_t& ctu = h2d.coding_tree_unit;
	ctu.num_frames = num_frame;
	std::copy(frame, frame + num_frame, ctu.frames);
	memset(ctu.lru, 0, sizeof(ctu.lru));
	ctu.neighbour_top = reinterpret_cast<h265d_neighbour_t*>(second_frame);
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
		for (int i = 0; i < NUM_ELEM(dst.sub_layer_second48bit); ++i) {
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

static void short_term_ref_pic_set_nopred_calc(h265d_short_term_ref_pic_set_t& dst, uint32_t num_pics, dec_bits& st) {
	int16_t* delta_poc = dst.delta_poc;
	uint16_t used_flag = 0;
	int val = 0;
	for (uint32_t i = 0; i < num_pics; ++i) {
		uint32_t delta;
		READ_CHECK_RANGE(ue_golomb(&st), delta, 32767, st);
		val = val - (delta + 1);
		delta_poc[i] = val;
		used_flag |= get_onebit(&st) << i;
	}
	dst.used_by_curr_pic_flag = used_flag;
}

static void short_term_ref_pic_set_nopred(h265d_short_term_ref_pic_set_t dst[], dec_bits& st) {
	uint32_t neg_pics, pos_pics;
	READ_CHECK_RANGE(ue_golomb(&st), neg_pics, 16, st);
	READ_CHECK_RANGE(ue_golomb(&st), pos_pics, 16 - neg_pics, st);
	dst[0].num_pics = neg_pics;
	dst[1].num_pics = pos_pics;
	short_term_ref_pic_set_nopred_calc(dst[0], neg_pics, st);
	short_term_ref_pic_set_nopred_calc(dst[1], pos_pics, st);
}

static void sps_short_term_ref_pic_set(h265d_sps_t& dst, uint32_t num, dec_bits& st) {
	short_term_ref_pic_set_nopred(dst.short_term_ref_pic_set[0], st);
	for (uint32_t i = 1; i < num; ++i) {
		if (get_onebit(&st)) {
			assert(0);
		} else {
			short_term_ref_pic_set_nopred(dst.short_term_ref_pic_set[i], st);
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
	return MultiplyDeBruijnBitPosition[(uint32_t)(num * 0x07C4ACDDU) >> 27];
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
	sps_short_term_ref_pic_set(dst, dst.num_short_term_ref_pic_sets, st);
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
	READ_CHECK_RANGE(ue_golomb(&st), dst.num_ref_idx_l0_default_active_minus1, 14, st);
	READ_CHECK_RANGE(ue_golomb(&st), dst.num_ref_idx_l1_default_active_minus1, 14, st);
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
		assert(0);
	}
	if (sps.sample_adaptive_offset_enabled_flag) {
		dst.slice_sao_luma_flag = get_onebit(&st);
		dst.slice_sao_chroma_flag = get_onebit(&st);
	}
	if (dst.slice_type != 2) {
		assert(0);
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
	CHECK_RANGE2(cb_qp_offset, -12, 12, st);
	cr_qp_offset += pps.pps_cr_qp_offset;
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
static const m2d_cabac_init_mn_t cabac_initial_value[3][154] = {
	{
		{0, 56}, {15, 48}, {-5, 72}, {-5, 88}, {0, 88}, {0, 64}, {0, 0}, {0, 0},
		{0, 0}, {0, 0}, {10, 48}, {0, 0}, {0, 0}, {0, 0}, {10, 48}, {-30, 104},
		{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
		{0, 0}, {0, 0}, {0, 0}, {0, 56}, {-5, 64}, {-5, 64}, {-15, 104}, {-5, 88},
		{-20, 96}, {-5, 64}, {10, 32}, {0, 64}, {0, 0}, {0, 0}, {0, 64}, {0, 64},
		{-5, 72}, {-5, 72}, {-15, 96}, {-15, 96}, {-10, 80}, {-10, 88}, {-5, 80}, {0, 56},
		{-10, 88}, {-10, 104}, {-5, 80}, {-15, 88}, {-15, 104}, {-5, 104}, {-10, 104}, {-15, 104},
		{-25, 104}, {-15, 80}, {-10, 72}, {-30, 104}, {-15, 96}, {-15, 96}, {-10, 80}, {-10, 88},
		{-5, 80}, {0, 56}, {-10, 88}, {-10, 104}, {-5, 80}, {-15, 88}, {-15, 104}, {-5, 104},
		{-10, 104}, {-15, 104}, {-25, 104}, {-15, 80}, {-10, 72}, {-30, 104}, {-20, 72}, {5, 72},
		{-5, 32}, {-5, 88}, {-15, 104}, {-15, 104}, {-10, 88}, {-15, 96}, {-15, 96}, {-20, 96},
		{-10, 80}, {-15, 80}, {-10, 80}, {-15, 72}, {-10, 88}, {-5, 88}, {10, 8}, {0, 56},
		{-10, 88}, {-15, 72}, {-10, 88}, {-5, 88}, {10, 8}, {0, 56}, {-10, 88}, {-15, 72},
		{-10, 88}, {-5, 88}, {10, 8}, {0, 56}, {-10, 88}, {-5, 80}, {-5, 72}, {10, 32},
		{10, 32}, {0, 48}, {-5, 48}, {0, 48}, {-5, 48}, {0, 56}, {-5, 48}, {-5, 72},
		{-15, 104}, {-5, 48}, {-5, 72}, {-15, 104}, {-5, 80}, {-20, 80}, {-5, 56}, {-5, 64},
		{-5, 80}, {0, 48}, {-5, 64}, {-5, 72}, {0, 56}, {-25, 64}, {0, 24}, {-20, 80},
		{-5, 72}, {-15, 72}, {-10, 64}, {0, 48}, {-5, 80}, {10, 8}, {5, 32}, {10, 32},
		{-5, 80}, {25, 8}, {-10, 64}, {15, 24}, {-5, 64}, {0, 56}, {-5, 48}, {5, 40},
		{0, 48}, {0, 48},
	},
};

static uint32_t sao_type_idx(m2d_cabac_t& cabac, dec_bits& st) {
	if (!cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->sao_type_idx)) {
		return 0;
	} else {
		return 1 + cabac_decode_bypass(&cabac, &st);
	}
}

static uint32_t sao_offset_abs(m2d_cabac_t& cabac, dec_bits& st) {
	assert(0);
	return 0;
}

static uint32_t sao_offset_sign(m2d_cabac_t& cabac, dec_bits& st) {
	assert(0);
	return 0;
}

static void sao_read_elem(h265d_sao_map_elem_t& dst, uint32_t idx, m2d_cabac_t& cabac, dec_bits& st) {
	for (int j = 0; j < 4; ++j) {
		dst.offset[j] = sao_offset_abs(cabac, st);
	}
	if (idx == 1) {
		for (int j = 0; j < 4; ++j) {
			uint32_t sign = sao_offset_sign(cabac, st);
			dst.offset[j] = NEGATE(dst.offset[j], sign);
		}
	}
}

static void sao_read(h265d_ctu_t& dst, const h265d_slice_header_t& hdr, dec_bits& st) {
	uint32_t merge_flag = 0;
	if (dst.pos_x != 0) {
		merge_flag = DIR_LEFT;
	}
	if (dst.pos_y != 0) {
		merge_flag |= DIR_TOP;
	}
	h265d_sao_map_t& sao_map = dst.sao_map[0];
	if (!merge_flag) {
		sao_map.idx = 0;
		if (hdr.body.slice_sao_luma_flag) {
			uint32_t idx = sao_type_idx(dst.cabac, st);
			if (idx != 0) {
				sao_map.idx = idx;
				sao_read_elem(sao_map.elem[0], idx, dst.cabac, st);
			}
		}
		if (hdr.body.slice_sao_chroma_flag) {
			uint32_t idx = sao_type_idx(dst.cabac, st);
			if (idx != 0) {
				sao_map.idx |= idx << 4;
				for (int i = 1; i < 3; ++i) {
					sao_read_elem(sao_map.elem[i], idx, dst.cabac, st);
				}
			}
		}
	}
}

static void sao_ignore(h265d_ctu_t& dst, const h265d_slice_header_t& hdr, dec_bits& st) {}

static inline uint32_t split_cu_flag(m2d_cabac_t& cabac, dec_bits& st, uint32_t cond) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->split_cu_flag + cond);
}

struct pred_mode_flag_ipic {
	uint8_t operator()(m2d_cabac_t& cabac, dec_bits& st) const {
		return 1;
	}
};

static inline uint32_t prev_intra_luma_pred_flag(m2d_cabac_t& cabac, dec_bits& st) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->prev_intra_luma_pred_flag);
}

static inline uint32_t mpm_idx(m2d_cabac_t& cabac, dec_bits& st) {
	return cabac_decode_bypass(&cabac, &st) ? 1 + cabac_decode_bypass(&cabac, &st) : 0;
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
	if (prefix <= 3) {
		return prefix;
	} else {
		uint32_t max = 1 << ((prefix >> 1) - 1);
		return max * (2 + (prefix & 1)) + cabac_decode_multibypass(&cabac, &st, max);
	}
}

static inline uint32_t coded_sub_block_flag(m2d_cabac_t& cabac, dec_bits& st, uint8_t prev_sbf) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->coded_sub_block_flag + ((prev_sbf & 1) | (prev_sbf >> 1)));
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

static inline uint32_t coeff_sign_flags(m2d_cabac_t& cabac, dec_bits& st, uint32_t first_pos, uint32_t hidden_bit) {
	return cabac_decode_multibypass(&cabac, &st, first_pos + (hidden_bit ^ 1)) << hidden_bit;
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

static inline uint32_t intra_chroma_pred_dir(uint32_t chroma_pred_mode, uint32_t mode) {
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

static void intra_pred_candidate(h265d_ctu_t& dst, uint8_t cand[], uint8_t candA, uint8_t candB) {
	if (candA == candB) {
		if (candA <= INTRA_DC) {
			cand[0] = INTRA_PLANAR;
			cand[1] = INTRA_DC;
			cand[2] = INTRA_ANGULAR26;
		} else {
			cand[0] = candA;
			cand[1] = (candA - 3) & 31;
			cand[2] = (candA - 1) & 31;
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

static const uint8_t h265d_scan_order2x2[3][2][4] = {
	{
		{
			0x0, 0x2, 0x1, 0x3,
		},
		{
			0x0, 0x2, 0x1, 0x3,
		},
	},
	{
		{
			0x0, 0x1, 0x2, 0x3,
		},
		{
			0x0, 0x1, 0x2, 0x3,
		},
	},
	{
		{
			0x0, 0x2, 0x1, 0x3,
		},
		{
			0x0, 0x2, 0x1, 0x3,
		},
	},
};

static const uint8_t h265d_scan_order4x4[3][2][16] = {
	{
		{
			0x0, 0x2, 0x5, 0x9, 0x1, 0x4, 0x8, 0xc,
			0x3, 0x7, 0xb, 0xe, 0x6, 0xa, 0xd, 0xf,
		},
		{
			0x0, 0x4, 0x1, 0x8, 0x5, 0x2, 0xc, 0x9,
			0x6, 0x3, 0xd, 0xa, 0x7, 0xe, 0xb, 0xf,
		},
	},
	{
		{
			0x0, 0x1, 0x2, 0x3, 0x4, 0x5, 0x6, 0x7,
			0x8, 0x9, 0xa, 0xb, 0xc, 0xd, 0xe, 0xf,
		},
		{
			0x0, 0x1, 0x2, 0x3, 0x4, 0x5, 0x6, 0x7,
			0x8, 0x9, 0xa, 0xb, 0xc, 0xd, 0xe, 0xf,
		},
	},
	{
		{
			0x0, 0x4, 0x8, 0xc, 0x1, 0x5, 0x9, 0xd,
			0x2, 0x6, 0xa, 0xe, 0x3, 0x7, 0xb, 0xf,
		},
		{
			0x0, 0x4, 0x8, 0xc, 0x1, 0x5, 0x9, 0xd,
			0x2, 0x6, 0xa, 0xe, 0x3, 0x7, 0xb, 0xf,
		},
	},
};

static const uint8_t h265d_scan_order8x8[3][2][64] = {
	{
		{
			0x0, 0x2, 0x5, 0x9, 0xe, 0x14, 0x1b, 0x23,
			0x1, 0x4, 0x8, 0xd, 0x13, 0x1a, 0x22, 0x2a,
			0x3, 0x7, 0xc, 0x12, 0x19, 0x21, 0x29, 0x30,
			0x6, 0xb, 0x11, 0x18, 0x20, 0x28, 0x2f, 0x35,
			0xa, 0x10, 0x17, 0x1f, 0x27, 0x2e, 0x34, 0x39,
			0xf, 0x16, 0x1e, 0x26, 0x2d, 0x33, 0x38, 0x3c,
			0x15, 0x1d, 0x25, 0x2c, 0x32, 0x37, 0x3b, 0x3e,
			0x1c, 0x24, 0x2b, 0x31, 0x36, 0x3a, 0x3d, 0x3f,
		},
		{
			0x0, 0x8, 0x1, 0x10, 0x9, 0x2, 0x18, 0x11,
			0xa, 0x3, 0x20, 0x19, 0x12, 0xb, 0x4, 0x28,
			0x21, 0x1a, 0x13, 0xc, 0x5, 0x30, 0x29, 0x22,
			0x1b, 0x14, 0xd, 0x6, 0x38, 0x31, 0x2a, 0x23,
			0x1c, 0x15, 0xe, 0x7, 0x39, 0x32, 0x2b, 0x24,
			0x1d, 0x16, 0xf, 0x3a, 0x33, 0x2c, 0x25, 0x1e,
			0x17, 0x3b, 0x34, 0x2d, 0x26, 0x1f, 0x3c, 0x35,
			0x2e, 0x27, 0x3d, 0x36, 0x2f, 0x3e, 0x37, 0x3f,
		},
	},
	{
		{
			0x0, 0x1, 0x2, 0x3, 0x4, 0x5, 0x6, 0x7,
			0x8, 0x9, 0xa, 0xb, 0xc, 0xd, 0xe, 0xf,
			0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
			0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
			0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27,
			0x28, 0x29, 0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f,
			0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
			0x38, 0x39, 0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f,
		},
		{
			0x0, 0x1, 0x2, 0x3, 0x4, 0x5, 0x6, 0x7,
			0x8, 0x9, 0xa, 0xb, 0xc, 0xd, 0xe, 0xf,
			0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
			0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
			0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27,
			0x28, 0x29, 0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f,
			0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
			0x38, 0x39, 0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f,
		},
	},
	{
		{
			0x0, 0x8, 0x10, 0x18, 0x20, 0x28, 0x30, 0x38,
			0x1, 0x9, 0x11, 0x19, 0x21, 0x29, 0x31, 0x39,
			0x2, 0xa, 0x12, 0x1a, 0x22, 0x2a, 0x32, 0x3a,
			0x3, 0xb, 0x13, 0x1b, 0x23, 0x2b, 0x33, 0x3b,
			0x4, 0xc, 0x14, 0x1c, 0x24, 0x2c, 0x34, 0x3c,
			0x5, 0xd, 0x15, 0x1d, 0x25, 0x2d, 0x35, 0x3d,
			0x6, 0xe, 0x16, 0x1e, 0x26, 0x2e, 0x36, 0x3e,
			0x7, 0xf, 0x17, 0x1f, 0x27, 0x2f, 0x37, 0x3f,
		},
		{
			0x0, 0x8, 0x10, 0x18, 0x20, 0x28, 0x30, 0x38,
			0x1, 0x9, 0x11, 0x19, 0x21, 0x29, 0x31, 0x39,
			0x2, 0xa, 0x12, 0x1a, 0x22, 0x2a, 0x32, 0x3a,
			0x3, 0xb, 0x13, 0x1b, 0x23, 0x2b, 0x33, 0x3b,
			0x4, 0xc, 0x14, 0x1c, 0x24, 0x2c, 0x34, 0x3c,
			0x5, 0xd, 0x15, 0x1d, 0x25, 0x2d, 0x35, 0x3d,
			0x6, 0xe, 0x16, 0x1e, 0x26, 0x2e, 0x36, 0x3e,
			0x7, 0xf, 0x17, 0x1f, 0x27, 0x2f, 0x37, 0x3f,
		},
	},
};

struct residual_scan_order_t {
	const uint8_t* sub_block_num;
	const uint8_t* xy_pos;
	uint8_t (*sig_coeff_flag_inc)(uint8_t sx, uint8_t sy, uint8_t px, uint8_t py, uint8_t prev_sbf);
};

template <int OFFSET>
static uint8_t sig_coeff_flag_inc4(uint8_t sx, uint8_t sy, uint8_t px, uint8_t py, uint8_t prev_sbf) {
	static const int8_t inc_map[] = {
		0, 1, 4, 5, 2, 3, 4, 5, 6, 6, 8, 8, 7, 7, 8
	};
	return inc_map[py * 4 + px] + OFFSET;
}

static const int8_t sig_coeff_inc_init[64] = {
	2, 2, 2, 2, 1, 2, 1, 2, 1, 2, 0, 2, 0, 2, 0, 2,
	1, 1, 2, 2, 1, 1, 1, 2, 0, 1, 0, 2, 0, 1, 0, 2,
	1, 0, 2, 2, 0, 0, 1, 2, 0, 0, 0, 2, 0, 0, 0, 2,
	0, 0, 2, 2, 0, 0, 1, 2, 0, 0, 0, 2, 0, 0, 0, 2
};

template <int MOD>
static uint8_t sig_coeff_flag_inc_luma(uint8_t sx, uint8_t sy, uint8_t px, uint8_t py, uint8_t prev_sbf) {
	uint8_t inc = 3;
	if ((sx == 0) && (sy == 0)) {
		if ((px == 0) && (py == 0)) {
			return 0;
		}
		inc = 0;
	}
	return sig_coeff_inc_init[(py << 4) | (px << 2) | prev_sbf] + inc + MOD;
}

template <int MOD>
static uint8_t sig_coeff_flag_inc_chroma(uint8_t sx, uint8_t sy, uint8_t px, uint8_t py, uint8_t prev_sbf) {
	return sig_coeff_inc_init[(py << 4) | (px << 2) | prev_sbf] + MOD + 27;
}

static const residual_scan_order_t residual_scan_order[3][3][2] = {
	{
		{
			{h265d_scan_order2x2[0][0], h265d_scan_order2x2[0][1], sig_coeff_flag_inc4<0>},
			{h265d_scan_order2x2[0][0], h265d_scan_order2x2[0][1], sig_coeff_flag_inc4<27>}
		},
		{
			{h265d_scan_order4x4[0][0], h265d_scan_order4x4[0][1], sig_coeff_flag_inc_luma<9>},
			{h265d_scan_order4x4[0][0], h265d_scan_order4x4[0][1], sig_coeff_flag_inc_chroma<9>}
		},
		{
			{h265d_scan_order8x8[0][0], h265d_scan_order8x8[0][1], sig_coeff_flag_inc_luma<21>},
			{h265d_scan_order8x8[0][0], h265d_scan_order8x8[0][1], sig_coeff_flag_inc_chroma<12>}
		}
	},
	{
		{
			{h265d_scan_order2x2[1][0], h265d_scan_order2x2[1][1], sig_coeff_flag_inc4<0>},
			{h265d_scan_order2x2[1][0], h265d_scan_order2x2[1][1], sig_coeff_flag_inc4<27>}
		},
		{
			{h265d_scan_order4x4[1][0], h265d_scan_order4x4[1][1], sig_coeff_flag_inc_luma<15>},
			{h265d_scan_order4x4[1][0], h265d_scan_order4x4[1][1], sig_coeff_flag_inc_chroma<9>}
		},
		{
			{h265d_scan_order8x8[1][0], h265d_scan_order8x8[1][1], sig_coeff_flag_inc_luma<21>},
			{h265d_scan_order8x8[1][0], h265d_scan_order8x8[1][1], sig_coeff_flag_inc_chroma<12>}
		}
	},
	{
		{
			{h265d_scan_order2x2[2][0], h265d_scan_order2x2[2][1], sig_coeff_flag_inc4<0>},
			{h265d_scan_order2x2[2][0], h265d_scan_order2x2[2][1], sig_coeff_flag_inc4<27>}
		},
		{
			{h265d_scan_order4x4[2][0], h265d_scan_order4x4[2][1], sig_coeff_flag_inc_luma<15>},
			{h265d_scan_order4x4[2][0], h265d_scan_order4x4[2][1], sig_coeff_flag_inc_chroma<9>}
		},
		{
			{h265d_scan_order8x8[2][0], h265d_scan_order8x8[2][1], sig_coeff_flag_inc_luma<21>},
			{h265d_scan_order8x8[2][0], h265d_scan_order8x8[2][1], sig_coeff_flag_inc_chroma<12>}
		}
	}
};

static inline uint32_t scan_order_index(uint32_t x, uint32_t y, uint32_t size_log2) {
	return (y << size_log2) + x;
}

static inline uint32_t sig_coeff_flags_read(m2d_cabac_t& cabac, dec_bits& st, uint8_t (*sig_coeff_flag_inc)(uint8_t, uint8_t, uint8_t, uint8_t, uint8_t), const uint8_t pix_pos[], bool is_last, int32_t pos, int8_t flags[], uint8_t sx, uint8_t sy, uint8_t prev_sbf) {
	int32_t idx = -1;
	if (is_last) {
		flags[++idx] = pos;
		pos--;
	}
	while (0 < pos) {
		uint8_t e = pix_pos[pos];
		uint8_t px = e & 3;
		uint8_t py = e >> 2;
		if (sig_coeff_flag(cabac, st, sig_coeff_flag_inc(sx, sy, px, py, prev_sbf))) {
			flags[++idx] = pos;
		}
		pos--;
	}
	if ((pos == 0) && ((sx + sy == 0) || (idx < 0) || sig_coeff_flag(cabac, st, sig_coeff_flag_inc(sx, sy, 0, 0, prev_sbf)))) {
		flags[++idx] = 0;
	}
	uint32_t first_pos = (0 <= idx) ? flags[idx] : 0;
	flags[++idx] = -1;
	return first_pos;
}

static void residual_coding(h265d_ctu_t& dst, dec_bits& st, uint32_t size_log2, int colour) {
	uint8_t sub_block_flags[8];
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
	memset(sub_block_flags, 0, sizeof(sub_block_flags));
	uint32_t order_idx;
	if (colour == 0) {
		order_idx = dst.order_luma;
	} else {
		order_idx = dst.order_chroma;
	}
	const residual_scan_order_t& suborder = residual_scan_order[order_idx][1][0];
	const residual_scan_order_t& order = residual_scan_order[order_idx][size_log2 - 2][(colour + 1) >> 1];
	uint32_t last_subblock_pos = order.sub_block_num[scan_order_index(last_sig_coeff_x >> 2, last_sig_coeff_y >> 2, size_log2 - 2)];
	uint8_t sub_pos_max = (1 << (size_log2 - 2)) - 1;
	int i = last_subblock_pos;
	uint8_t greater1ctx = 1;
	uint32_t num = suborder.sub_block_num[scan_order_index(last_sig_coeff_x & 3, last_sig_coeff_y & 3, 2)];
	bool is_last = true;
	do {
		uint8_t sx = order.xy_pos[i] & ((1 << size_log2) - 1);
		uint8_t sy = order.xy_pos[i] >> size_log2;
		uint8_t prev_sbf = (sx < sub_pos_max) ? sub_block_flags[sy] & (1 << (sx + 1)) : 0;
		prev_sbf |= (sy < sub_pos_max) ? (sub_block_flags[sy + 1] & (1 << sx)) << 1 : 0;
		if (((unsigned)(last_subblock_pos - 1) <= (unsigned)(i - 1)) || coded_sub_block_flag(dst.cabac, st, prev_sbf)) {
			sub_block_flags[sy] |= 1 << sx;
			int8_t sig_coeff_flags[4 * 4 + 1];
			uint32_t first_pos = sig_coeff_flags_read(dst.cabac, st, order.sig_coeff_flag_inc, suborder.xy_pos, is_last, num, sig_coeff_flags, sx, sy, prev_sbf);
			uint32_t num_greater1 = 0;
			uint32_t ctxset = (((colour == 0) || (i == 0)) ? 0 : 2) + (greater1ctx == 0);
			uint32_t greater1offset = ctxset * 4 + ((colour == 0) ? 0 : 16);
			greater1ctx = 1;
			uint8_t greater1_flags = 0;
			int last_greater1_pos = -1;
			const int8_t* sig_coeff = sig_coeff_flags;
			int8_t pos;
			while (0 <= (pos = *sig_coeff++)) {
				greater1_flags = (greater1_flags << 1) + coeff_abs_level_greater1_flag(dst.cabac, st, greater1offset + greater1ctx);
				if (greater1_flags & 1) {
					greater1ctx = 0;
					if (last_greater1_pos < 0) {
						last_greater1_pos = pos;
					}
				} else if ((uint32_t)(greater1ctx - 1) < 2) {
					greater1ctx++;
				}
				if (8 <= ++num_greater1) {
					break;
				}
			}
			uint32_t greater2_flag = (0 <= last_greater1_pos) ? coeff_abs_level_greater2_flag(dst.cabac, st, (colour == 0) ? ctxset : ctxset + 4) : 0;
			uint8_t hidden_bit;
			if (dst.pps->sign_data_hiding_enabled_flag && (3 < sig_coeff_flags[0] - sig_coeff_flags[first_pos])) {
				hidden_bit = 1;
			} else {
				hidden_bit = 0;
			}
			uint16_t sign_flags = coeff_sign_flags(dst.cabac, st, first_pos, hidden_bit);
			sig_coeff = sig_coeff_flags;
			uint32_t abs_level = 0;
			uint32_t rice_param = 0;
			uint32_t shift = first_pos;
			int16_t* write_pos = dst.coeff_buf + (sy * (1 << size_log2) + sx) * 4;
			int32_t level_sum = 0;
			while (0 <= (pos = *sig_coeff++)) {
				uint32_t level = ((pos == last_greater1_pos) ? greater2_flag + 1 : 1) + ((greater1_flags >> shift) & 1);
				uint32_t remain_thr = (first_pos < 7) ? ((pos == last_greater1_pos) ? 3 : 2) : 1;
				if (level == remain_thr) {
					rice_param = std::min(rice_param + ((static_cast<uint32_t>(3) << rice_param) < abs_level), static_cast<uint32_t>(4));
					abs_level = coeff_abs_level_remaining(dst.cabac, st, rice_param);
					level += abs_level;
				}
				level_sum += level;
				uint8_t cx = suborder.xy_pos[i] & 3;
				uint8_t cy = suborder.xy_pos[i] >> 2;
				write_pos[cy * (1 << size_log2) + cx] = NEGATE(level, (sign_flags >> shift) & 1);
				shift--;
			}
		}
		num = 15;
	} while (0 <= --i);
}

static void transform_unit(h265d_ctu_t& dst, dec_bits& st, uint32_t size_log2, uint8_t cbf, int idx) {
	if (dst.qp_delta_req) {
		dst.qp_delta_req = 0;
		dst.qp += cu_qp_delta(dst.cabac, st);
	}
	if (cbf & 1) {
		residual_coding(dst, st, size_log2, 0);
	}
	if (cbf & 6) {
		if (2 < size_log2) {
			size_log2 -= 1;
		} else if (idx != 3) {
			return;
		}
		if (cbf & 4) {
			residual_coding(dst, st, size_log2, 1);
		}
		if (cbf & 2) {
			residual_coding(dst, st, size_log2, 2);
		}
	}
}

static void transform_tree(h265d_ctu_t& dst, dec_bits& st, uint32_t size_log2, uint32_t depth, uint8_t upper_cbf_cbcr, int idx) {
	uint32_t split;
	if (size_log2 <= dst.sps->ctb_info.transform_log2) {
		if ((dst.sps->ctb_info.transform_log2_min < size_log2) && (depth < dst.sps->max_transform_hierarchy_depth_intra)) {
			split = split_transform_flag(dst.cabac, st, size_log2);
		} else {
			split = 0;
		}
	} else {
		split = 1;
	}
	uint8_t cbf = 0;
	if (2 < size_log2) {
		if (upper_cbf_cbcr & 2) {
			cbf |= cbf_chroma(dst.cabac, st, depth) * 2;
		}
		if (upper_cbf_cbcr & 1) {
			cbf |= cbf_chroma(dst.cabac, st, depth);
		}
	}
	if (split) {
		size_log2 -= 1;
		depth += 1;
		transform_tree(dst, st, size_log2, depth, cbf, 0);
		transform_tree(dst, st, size_log2, depth, cbf, 1);
		transform_tree(dst, st, size_log2, depth, cbf, 2);
		transform_tree(dst, st, size_log2, depth, cbf, 3);
	} else {
		cbf = cbf * 2;
		if (dst.is_intra || depth || cbf) {
			cbf |= cbf_luma(dst.cabac, st, depth);
		}
		if (cbf) {
			transform_unit(dst, st, size_log2, cbf, idx);
		}
	}
}

static inline void intra_pred_mode_fill(h265d_neighbour_t dst0[], h265d_neighbour_t dst1[], uint32_t mode, uint32_t size_log2) {
	int num = 1 << (size_log2 - 3);
	for (int i = 0; i < num; ++i) {
		dst0[i].pred_mode = mode;
		dst1[i].pred_mode = mode;
	}
}

static inline void intra_depth_fill(h265d_neighbour_t dst0[], h265d_neighbour_t dst1[], uint32_t size_log2) {
	int num = 1 << (size_log2 - 3);
	uint8_t depth = 6 - size_log2;
	for (int i = 0; i < num; ++i) {
		dst0[i].depth = depth;
		dst1[i].depth = depth;
	}
}

template <typename F>
static void coding_unit(h265d_ctu_t& dst, dec_bits& st, uint32_t size_log2, h265d_neighbour_t* left, h265d_neighbour_t* top, F PredModeFlag) {
	static const uint8_t scan_order[35] = {
		0, 0, 0, 0, 0, 0,
		2, 2, 2, 2, 2, 2,
		2, 2, 2, 0, 0, 0,
		0, 0, 0, 0, 1, 1,
		1, 1, 1, 1, 1, 1,
		1, 0, 0, 0, 0
	};
	if (dst.pps->transquant_bypass_enabled_flag) {
		assert(0);
	}
	if (((dst.is_intra = PredModeFlag(dst.cabac, st)) == 0) || (dst.sps->ctb_info.size_log2_min == size_log2)) {
		assert(0);
	}
	if ((dst.sps->ctb_info.pcm_log2_min <= size_log2) && (size_log2 <= dst.sps->ctb_info.pcm_log2)) {
		assert(0);
	}
	uint32_t pred_flag = 0;
	for (int i = 0; i < 1; ++i) {
		pred_flag |= prev_intra_luma_pred_flag(dst.cabac, st) << i;
	}
	uint8_t candidate[3];
	uint8_t luma_mode;
	intra_pred_candidate(dst, candidate, left->pred_mode, top->pred_mode);
	if (pred_flag & 1) {
		luma_mode = candidate[mpm_idx(dst.cabac, st)];
	} else {
		assert(0);
	}
	dst.order_luma = scan_order[luma_mode];
	pred_flag >>= 1;
	intra_pred_mode_fill(left, top, luma_mode, size_log2);
	uint32_t chroma_mode_idx = intra_chroma_pred_mode(dst.cabac, st);
	for (int i = 0; i < 1; ++i) {
		uint8_t chroma_mode = intra_chroma_pred_dir(chroma_mode_idx, luma_mode);
		dst.order_chroma = scan_order[chroma_mode];
	}
	transform_tree(dst, st, size_log2, 0, 3, 0);
}

static void quad_tree_normal(h265d_ctu_t& dst, dec_bits& st, uint32_t size_log2, uint8_t offset_x, uint8_t offset_y, h265d_neighbour_t* left, h265d_neighbour_t* top) {
	if ((dst.sps->ctb_info.size_log2_min < size_log2) && split_cu_flag(dst.cabac, st, (6 < size_log2 + left->depth) + (6 < size_log2 + top->depth))) {
		size_log2 -= 1;
		uint32_t block_len = 1 << size_log2;
		uint32_t info_offset = 1 << (size_log2 - 2);
		quad_tree_normal(dst, st, size_log2, offset_x, offset_y, left, top);
		quad_tree_normal(dst, st, size_log2, offset_x + block_len, offset_y, left, top + info_offset);
		quad_tree_normal(dst, st, size_log2, offset_x, offset_y + block_len, left + info_offset, top);
		quad_tree_normal(dst, st, size_log2, offset_x + block_len, offset_y + block_len, left + info_offset, top + info_offset);
	} else {
		intra_depth_fill(left, top, size_log2);
		if (dst.pps->cu_qp_delta_enabled_flag) {
			dst.qp_delta_req = 1;
		}
		return coding_unit(dst, st, size_log2, left, top, pred_mode_flag_ipic());
	}
}

static void coding_tree_unit(h265d_ctu_t& dst, dec_bits& st) {
	dst.sao_read(dst, *dst.slice_header, st);
	if ((dst.sps->ctb_info.columns - dst.pos_x < 2) || (dst.sps->ctb_info.rows - dst.pos_y < 2)) {
		assert(0);
	} else {
		quad_tree_normal(dst, st, dst.sps->ctb_info.size_log2, 0, 0, dst.neighbour_left, dst.neighbour_top + dst.pos_x * sizeof(dst.neighbour_left));
	}
}

static void ctu_init(h265d_ctu_t& dst, const h265d_slice_header_t& hdr, const h265d_pps_t& pps, const h265d_sps_t& sps) {
	const h265d_slice_header_body_t& header = hdr.body;
	int slice_type = header.slice_type;
	int idc = (slice_type < 2) ? ((slice_type ^ header.cabac_init_flag) + 1) : 0;
	init_cabac_context(&dst.cabac, header.slice_qpy, cabac_initial_value[idc], NUM_ELEM(cabac_initial_value[idc]));
	dst.sao_read = (hdr.body.slice_sao_luma_flag || hdr.body.slice_sao_chroma_flag) ? sao_read : sao_ignore;
	dst.sps = &sps;
	int ctu_address = hdr.slice_segment_address;
	dst.idx_in_slice = 0;
	dst.pos_y = (ctu_address / sps.ctb_info.columns) << sps.ctb_info.size_log2;
	dst.pos_x = (ctu_address - sps.ctb_info.columns * dst.pos_y) << sps.ctb_info.size_log2;
	dst.pps = &pps;
	dst.slice_header = &hdr;
	dst.qp = pps.init_qp_minus26 + 26;
	memset(dst.neighbour_left, INTRA_DC, sizeof(dst.neighbour_left));
	uint32_t neighbour_flags_top_len = sps.ctb_info.columns * sizeof(dst.neighbour_left);
	memset(dst.neighbour_top, INTRA_DC, neighbour_flags_top_len);
	dst.sao_map = reinterpret_cast<h265d_sao_map_t*>(dst.neighbour_top + neighbour_flags_top_len);
}

static void ctu_pos_increment(h265d_ctu_t& dst) {
	uint32_t pos_x = dst.pos_x + 1;
	if (dst.sps->ctb_info.columns <= pos_x) {
		dst.pos_y++;
		memset(dst.neighbour_left, INTRA_DC, sizeof(dst.neighbour_left));
		pos_x = 0;
	}
	dst.pos_x = pos_x;
	dst.idx_in_slice++;
	h265d_neighbour_t* top = dst.neighbour_top + pos_x * 8;
	for (int i = 0; i < 8; ++i) {
		top[i].pred_mode = INTRA_DC;
	}
}

static void slice_data(h265d_ctu_t& dst, const h265d_slice_header_t& hdr, const h265d_pps_t& pps, const h265d_sps_t& sps, dec_bits& st) {
	ctu_init(dst, hdr, pps, sps);
	init_cabac_engine(&dst.cabac, &st);
	int end_of_slice;
	do {
		coding_tree_unit(dst, st);
		ctu_pos_increment(dst);
		end_of_slice = 1;
	} while (!end_of_slice);
}

static void slice_layer(h265d_data_t& h2d, dec_bits& st) {
	h265d_slice_header_t& header = h2d.slice_header;
	header.body.nal_type = h2d.current_nal;
	header.first_slice_segment_in_pic_flag = get_onebit(&st);
	if ((BLA_W_LP <= h2d.current_nal) &&  (h2d.current_nal <= RSV_IRAP_VCL23)) {
		header.no_output_of_prior_pics_flag = get_onebit(&st);
	}
	READ_CHECK_RANGE(ue_golomb(&st), header.pps_id, 63, st);
	const h265d_pps_t& pps = h2d.pps[header.pps_id];
	const h265d_sps_t& sps = h2d.sps[pps.sps_id];
	slice_header(header, pps, sps, st);
	slice_data(h2d.coding_tree_unit, header, pps, sps, st);
}

static int dispatch_one_nal(h265d_data_t& h2d, uint32_t nalu_header) {
	int err = 0;
	dec_bits& st = h2d.stream_i;
	switch (h2d.current_nal = static_cast<h265d_nal_t>((nalu_header >> 9) & 63)) {
	case IDR_W_RADL:
		slice_layer(h2d, st);
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

int h265d_peek_decoded_frame(h265d_context *h2, m2d_frame_t *frame, int bypass_dpb)
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

int h265d_get_decoded_frame(h265d_context *h2, m2d_frame_t *frame, int bypass_dpb)
{
//	h265d_data_t& h2d = *reinterpret_cast<h265d_data_t*>(h2);
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
