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
		+ sizeof(h2d.coding_tree_unit.neighbour_flags_top[0]) * sps.ctb_info.columns;
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
	ctu.neighbour_flags_top = reinterpret_cast<uint8_t*>(second_frame);
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
	ctb_info.pcm_log2_min = sps.log2_min_pcm_luma_coding_block_size_minus3 + 3;
	ctb_info.pcm_log2 = ctb_info.pcm_log2_min + sps.log2_diff_max_min_pcm_luma_coding_block_size;
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
	set_ctb_info(dst.ctb_info, dst);
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

#define NEGATE(val, sign) (((val) ^ (unsigned)-(signed)(sign)) + (sign))

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

static inline uint32_t prev_intra_luma_pred_flag(m2d_cabac_t& cabac, dec_bits& st) {
	return cabac_decode_decision_raw(&cabac, &st, reinterpret_cast<h265d_cabac_context_t*>(cabac.context)->prev_intra_luma_pred_flag);
}

static inline uint32_t mpm_idx(m2d_cabac_t& cabac, dec_bits& st) {
	return cabac_decode_bypass(&cabac, &st) ? 1 + cabac_decode_bypass(&cabac, &st) : 0;
}

typedef union {
	struct {
		uint8_t right, bottom;
	} elem;
	uint16_t all;
} h265d_qtree_neighbour_t;

typedef union {
	struct {
		uint8_t offset_x, offset_y;
		uint8_t neighbour_left, neighbour_top;
	} elem;
	uint32_t all;
} h265d_qtree_info_t;

static void intra_pred_candidate(h265d_ctu_t& dst, uint8_t cand[], h265d_qtree_info_t info) {
	const uint8_t* cand_pos = &dst.intra_pred_record[info.elem.offset_x & ~7];
	uint8_t candA = *cand_pos;
	uint8_t candB = (info.elem.offset_y == 0) ? INTRA_DC : cand_pos[info.elem.offset_y >> 3];
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

static void coding_unit_ipic(h265d_ctu_t& dst, dec_bits& st, uint32_t size_log2, h265d_qtree_info_t info) {
	if (dst.pps->transquant_bypass_enabled_flag) {
		assert(0);
	}
	if ((dst.sps->ctb_info.pcm_log2_min <= size_log2) && (size_log2 <= dst.sps->ctb_info.pcm_log2)) {
		assert(0);
	}
	uint32_t pred_flag = 0;
	for (int i = 0; i < 1; ++i) {
		pred_flag |= prev_intra_luma_pred_flag(dst.cabac, st) << i;
	}
	for (int i = 0; i < 1; ++i) {
		uint8_t candidate[3];
		intra_pred_candidate(dst, candidate, info);
		uint32_t mode;
		if (pred_flag & 1) {
			mode = candidate[mpm_idx(dst.cabac, st)];
		} else {
			assert(0);
		}
		uint8_t* right_edge = &dst.intra_pred_record[(info.elem.offset_x + (1 << size_log2)) & ~7];
		pred_flag >>= 1;
	}
}

static uint16_t quad_tree_normal(h265d_ctu_t& dst, dec_bits& st, uint32_t size_log2, h265d_qtree_info_t info) {
	if ((dst.sps->ctb_info.size_log2_min < size_log2) && split_cu_flag(dst.cabac, st, (info.elem.neighbour_left & 1) + (info.elem.neighbour_top & 1))) {
		uint32_t shift = (size_log2 & 3) | 1;
		size_log2 -= 1;
		uint32_t len = 1 << size_log2;
		info.elem.neighbour_left >>= 1;
		info.elem.neighbour_top >>= 1;
		h265d_qtree_neighbour_t qtree0, qtree1, qtree2, qtree3;
		uint32_t org = info.all;
		qtree0.all = quad_tree_normal(dst, st, size_log2, info);
		info.elem.offset_x += len;
		info.elem.neighbour_left = qtree0.elem.right;
		info.elem.neighbour_top >>= shift;
		qtree1.all = quad_tree_normal(dst, st, size_log2, info);
		info.all = org;
		info.elem.offset_y += len;
		info.elem.neighbour_top = qtree0.elem.bottom;
		info.elem.neighbour_left >>= shift;
		qtree2.all = quad_tree_normal(dst, st, size_log2, info);
		info.elem.offset_x += len;
		info.elem.neighbour_top = qtree1.elem.bottom;
		info.elem.neighbour_left = qtree2.elem.right;
		qtree3.all = quad_tree_normal(dst, st, size_log2, info);
		qtree1.elem.right = qtree1.elem.right | (qtree3.elem.right << shift);
		qtree1.elem.bottom = qtree2.elem.bottom | (qtree3.elem.bottom << shift);
		return qtree1.all;
	} else {
		coding_unit_ipic(dst, st, size_log2, info);
		return 0;
	}
}

static void coding_tree_unit(h265d_ctu_t& dst, dec_bits& st) {
	dst.sao_read(dst, *dst.slice_header, st);
	h265d_qtree_info_t info;
	info.elem.offset_x = 0;
	info.elem.offset_y = 0;
	info.elem.neighbour_left = dst.neighbour_flags_left;
	info.elem.neighbour_top = dst.neighbour_flags_top[dst.pos_x];
	h265d_qtree_neighbour_t neighbour;
	if ((dst.sps->ctb_info.columns - dst.pos_x < 2) || (dst.sps->ctb_info.rows - dst.pos_y < 2)) {
		assert(0);
	} else {
		neighbour.all = quad_tree_normal(dst, st, dst.sps->ctb_info.size_log2, info);
	}
	dst.neighbour_flags_left = neighbour.elem.right;
	dst.neighbour_flags_top[dst.pos_x] = neighbour.elem.bottom;
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
	dst.neighbour_flags_left = 0;
	memset(dst.intra_pred_record, 0, 8);
	uint32_t neighbour_flags_top_len = (sps.ctb_info.columns * sizeof(dst.neighbour_flags_top[0]) + 1) >> 1;
	memset(dst.neighbour_flags_top, 0, neighbour_flags_top_len);
	dst.sao_map = reinterpret_cast<h265d_sao_map_t*>(dst.neighbour_flags_top + neighbour_flags_top_len);
}

static void ctu_pos_increment(h265d_ctu_t& dst) {
	uint32_t pos_x = dst.pos_x + 1;
	if (dst.sps->ctb_info.columns <= pos_x) {
		dst.pos_y++;
		dst.neighbour_flags_left = 0;
		memset(dst.intra_pred_record, 0, 8);
		pos_x = 0;
	}
	dst.pos_x = pos_x;
	dst.idx_in_slice++;
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
