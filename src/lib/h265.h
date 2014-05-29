#ifndef __H265DEC_H__
#define __H265DEC_H__

#ifdef _MSC_VER
#include <crtdbg.h>
#define VC_CHECK assert(_CrtCheckMemory());
#else
#define VC_CHECK
#endif
#include "m2d.h"

typedef struct {
	uint8_t sub_layer_profile_first8bit;
	uint8_t sub_layer_level_idc;
	uint8_t sub_layer_second48bit[6];
	uint32_t sub_layer_profile_compatibility_flag;
} h265d_sub_layer_info_t;

typedef struct {
	uint8_t max_sub_layer;
	uint8_t general_profile_first8bit;
	uint8_t general_level_idc;
	uint8_t general_second48bit[6];
	uint32_t general_profile_compatibility_flag;
	uint16_t sub_layer_profile_level_present_flag;
	h265d_sub_layer_info_t sub_layer_info[8];
} h265d_profile_tier_level_t;

typedef struct {
	uint32_t num_units_in_tick;
	uint32_t time_scale;
	uint8_t poc_proportional_to_timing_flag;
	uint32_t num_ticks_poc_diff_one_minus1;
	uint16_t vps_num_hrd_parameters;
} h265d_vps_timing_info_t;

typedef struct {
	uint32_t max_dec_pic_buffering_minus1;
	uint32_t max_num_reorder_pic;
	uint32_t max_latency_increase_plus1;
} h265d_sub_layer_reordering_info_t;

typedef struct {
	uint8_t id;
	uint8_t max_layer;
	uint8_t max_layer_id;
	uint8_t temporal_id_nesting_flag;
	uint8_t sub_layer_ordering_info_present_flag;
	uint8_t timing_info_present_flag;
	uint16_t num_layer_sets_minus1;
	uint8_t layer_id_included_flag[2];
	h265d_profile_tier_level_t profile_tier_level;
	h265d_sub_layer_reordering_info_t max_buffering[8];
	h265d_vps_timing_info_t timing_info;
} h265d_vps_t;

typedef struct {
	uint32_t left_offset;
	uint32_t right_offset;
	uint32_t top_offset;
	uint32_t bottom_offset;
} h265d_conformance_window_t;

typedef struct {
	uint8_t scaling_list_pred_mode_flag[4][6];
	uint32_t scaling_list_pred_matrix_id_delta[4][6];
	uint8_t scale0[6][16];
	uint8_t scale1[6][64];
	uint8_t scale2[6][64];
	uint8_t scale3[2][64];
} h265d_scaling_list_data_t;

typedef struct {
	uint8_t num_negative_pics;
	uint8_t num_positive_pics;
	uint8_t delta_poc_s0_minus1;
	uint8_t used_by_curr_pic_s0_flag;
	uint8_t delta_poc_s1_minus1;
	uint8_t used_by_curr_pic_s1_flag;
} h265d_short_term_ref_pic_set_nopred_t;

typedef struct {
	uint8_t aspect_ratio_idc;
	uint8_t video_format;
	uint8_t colour_primaries;
	uint8_t transfer_characteristics;
	uint8_t matrix_coeffs;
	uint8_t overscan_info_present_flag;
	uint8_t overscan_appropriate_flag;
	uint8_t video_full_range_flag;
	uint8_t colour_description_present_flag;
	uint16_t sar_width;
	uint16_t sar_height;
} h265d_vui_parameters_t;

typedef struct {
	uint8_t vps_id;
	uint8_t max_sub_layers_minus1;
	uint8_t id;
	uint8_t chroma_format_idc;
	uint32_t pic_width_in_luma_samples;
	uint32_t pic_height_in_luma_samples;
	uint8_t bit_depth_luma_minus8;
	uint8_t bit_depth_chroma_minus8;
	uint8_t log2_max_pic_order_cnt_lsb_minus4;
	uint8_t log2_min_luma_coding_block_size_minus3;
	uint8_t log2_diff_max_min_luma_coding_block_size;
	uint8_t log2_min_transform_block_size_minus2;
	uint8_t log2_diff_max_min_transform_block_size;
	uint8_t max_transform_hierarchy_depth_inter;
	uint8_t max_transform_hierarchy_depth_intra;
	uint8_t pcm_sample_bit_depth_luma_minus1;
	uint8_t pcm_sample_bit_depth_chroma_minus1;
	uint8_t log2_min_pcm_luma_coding_block_size_minus3;
	uint8_t log2_diff_max_min_pcm_luma_coding_block_size;
	uint8_t num_short_term_ref_pic_sets;
	uint8_t num_long_term_ref_pics_sps;
	uint8_t temporal_id_nesting_flag;
	uint8_t separate_colour_plane_flag;
	uint8_t conformance_window_flag;
	uint8_t sub_layer_ordering_info_present_flag;
	uint8_t scaling_list_enabled_flag;
	uint8_t scaling_list_data_present_flag;
	uint8_t amp_enabled_flag;
	uint8_t sample_adaptive_offset_enabled_flag;
	uint8_t pcm_enabled_flag;
	uint8_t pcm_loop_filter_disabled_flag;
	uint8_t long_term_ref_pics_present_flag;
	uint8_t sps_temporal_mvp_enabled_flag;
	uint8_t strong_intra_smoothing_enabled_flag;
	uint8_t vui_parameters_present_flag;
	uint32_t used_by_curr_pic_lt_sps_flag;
	uint16_t lt_ref_pic_poc_lsb_sps[32];
	h265d_profile_tier_level_t profile_tier_level;
	h265d_conformance_window_t conf_win;
	h265d_sub_layer_reordering_info_t max_buffering[8];
	h265d_scaling_list_data_t scaling_list_data;
	h265d_vui_parameters_t vui_parameters;
} h265d_sps_t;

typedef struct {
	int id;
	dec_bits *stream;
	int (*header_callback)(void *arg, void *seq_id);
	void *header_callback_arg;
	dec_bits stream_i;
	h265d_vps_t vps;
	h265d_sps_t sps;
} h265d_context;

#endif /* __H265DEC_H__ */
