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
	uint16_t sub_layer_profile_level_present_flag;
	uint32_t general_profile_compatibility_flag;
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
	uint8_t layer_id_included_flag[2];
	uint16_t num_layer_sets_minus1;
	uint32_t temporal_id_nesting_flag : 1;
	uint32_t sub_layer_ordering_info_present_flag : 1;
	uint32_t timing_info_present_flag : 1;
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
	uint8_t num_pics;
	int16_t delta_poc[16];
	uint16_t used_by_curr_pic_flag;
} h265d_short_term_ref_pic_set_t;

typedef struct {
	uint8_t aspect_ratio_idc;
	uint8_t video_format;
	uint8_t colour_primaries;
	uint8_t transfer_characteristics;
	uint8_t matrix_coeffs;
	uint16_t sar_width;
	uint16_t sar_height;
	uint32_t overscan_info_present_flag : 1;
	uint32_t overscan_appropriate_flag : 1;
	uint32_t video_full_range_flag : 1;
	uint32_t colour_description_present_flag : 1;
} h265d_vui_parameters_t;

typedef struct {
	uint8_t vps_id;
	uint8_t max_sub_layers_minus1;
	uint8_t temporal_id_nesting_flag;
	h265d_profile_tier_level_t profile_tier_level;
} h265d_sps_prefix_t;

typedef struct {
	uint8_t size_log2;
	uint8_t num_ctb_log2;
	uint16_t columns;
	uint16_t rows;
} h265d_sps_ctb_info_t;

typedef struct {
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
	uint32_t separate_colour_plane_flag : 1;
	uint32_t conformance_window_flag : 1;
	uint32_t sub_layer_ordering_info_present_flag : 1;
	uint32_t scaling_list_enabled_flag : 1;
	uint32_t scaling_list_data_present_flag : 1;
	uint32_t amp_enabled_flag : 1;
	uint32_t sample_adaptive_offset_enabled_flag : 1;
	uint32_t pcm_enabled_flag : 1;
	uint32_t pcm_loop_filter_disabled_flag : 1;
	uint32_t long_term_ref_pics_present_flag : 1;
	uint32_t sps_temporal_mvp_enabled_flag : 1;
	uint32_t strong_intra_smoothing_enabled_flag : 1;
	uint32_t vui_parameters_present_flag : 1;
	uint32_t used_by_curr_pic_lt_sps_flag;
	uint16_t lt_ref_pic_poc_lsb_sps[32];
	h265d_sps_ctb_info_t ctb_info;
	h265d_sps_prefix_t prefix;
	h265d_conformance_window_t conf_win;
	h265d_sub_layer_reordering_info_t max_buffering[8];
	h265d_scaling_list_data_t scaling_list_data;
	h265d_short_term_ref_pic_set_t short_term_ref_pic_set[64][2];
	h265d_vui_parameters_t vui_parameters;
} h265d_sps_t;

typedef struct {
	uint16_t num_tile_columns_minus1;
	uint16_t num_tile_rows_minus1;
	uint16_t* column_width_minus1;
	uint16_t* row_height_minus1;
	uint32_t uniform_spacing_flag : 1;
	uint32_t loop_filter_across_tiles_enabled_flag : 1;
} h265d_tiles_t;

typedef struct {
	uint8_t sps_id;
	uint8_t num_ref_idx_l0_default_active_minus1;
	uint8_t num_ref_idx_l1_default_active_minus1;
	uint8_t init_qp_minus26;
	uint8_t diff_cu_qp_delta_depth;
	int8_t pps_cb_qp_offset;
	int8_t pps_cr_qp_offset;
	int8_t pps_beta_offset_div2;
	int8_t pps_tc_offset_div2;
	uint8_t log2_parallel_merge_level_minus2;
	uint32_t dependent_slice_segments_enabled_flag : 1;
	uint32_t output_flag_present_flag : 1;
	uint32_t num_extra_slice_header_bits : 3;
	uint32_t sign_data_hiding_enabled_flag : 1;
	uint32_t cabac_init_present_flag : 1;
	uint32_t constrained_intra_pred_flag : 1;
	uint32_t transform_skip_enabled_flag : 1;
	uint32_t cu_qp_delta_enabled_flag : 1;
	uint32_t pps_slice_chroma_qp_offsets_present_flag : 1;
	uint32_t weighted_pred_flag : 1;
	uint32_t weighted_bipred_flag : 1;
	uint32_t transquant_bypass_enabled_flag : 1;
	uint32_t tiles_enabled_flag : 1;
	uint32_t entropy_coding_sync_enabled_flag : 1;
	uint32_t pps_loop_filter_across_slices_enabled_flag : 1;
	uint32_t deblocking_filter_control_present_flag : 1;
	uint32_t deblocking_filter_override_enabled_flag : 1;
	uint32_t pps_deblocking_filter_disabled_flag : 1;
	uint32_t pps_scaling_list_data_present_flag : 1;
	uint32_t lists_modification_present_flag : 1;
	uint32_t slice_segment_header_extension_present_flag : 1;
	uint32_t pps_extension_flag : 1;
	h265d_tiles_t tiles;
	h265d_scaling_list_data_t scaling_list_data;
} h265d_pps_t;

typedef struct {
	uint8_t pic_type;
} h265d_access_unit_delimite_t;

typedef enum {
	BLA_W_LP = 16,
	IDR_W_RADL = 19,
	IDR_N_LP = 20,
	RSV_IRAP_VCL23 = 23,
	VPS_NAL = 32,
	SPS_NAL = 33,
	PPS_NAL = 34,
	AUD_NAL = 35,
	EOS_NAL = 36,
	EOB_NAL = 37,
	FILLER_NAL = 38,
	PREFIX_SEI = 39,
	SUFFIX_SEI = 40
} h265d_nal_t;

typedef struct {
	int8_t sao_merge_flag[3];
	int8_t sao_type_idx[3];
	int8_t split_cu_flag[9];
	int8_t cu_transquant_bypass_flag[3];
	int8_t cu_skip_flag[6];
	int8_t pred_mode_flag[2];
	int8_t part_mode[9];
	int8_t prev_intra_luma_pred_flag[3];
	int8_t intra_chroma_pred_mode[3];
	int8_t rqt_root_cbf[2];
	int8_t merge_flag[2];
	int8_t merge_idx[2];
	int8_t inter_pred_idc[10];
	int8_t ref_idx[4];
	int8_t mvp_flag[2];
	int8_t split_transform_flag[9];
	int8_t cbf_luma[6];
	int8_t cbf_chroma[12];
	int8_t abs_mvd_greater_flag[4];
	int8_t cu_qp_delta_abs[6];
	int8_t transform_skip_flag[6];
} context;

typedef struct {
	uint32_t range;
	uint32_t offset;
	h265d_cabac_context_t context;
} h265d_cabac_t;

typedef struct {
	h265d_nal_t nal_type;
	uint8_t slice_type;
	int8_t slice_qpy;
	int8_t slice_qpcb;
	int8_t slice_qpcr;
	int8_t slice_beta_offset_div2;
	int8_t slice_tc_offset_div2;
	uint32_t pic_output_flag : 1;
	uint32_t colour_plane_id : 2;
	uint32_t slice_sao_luma_flag : 1;
	uint32_t slice_sao_chroma_flag : 1;
	uint32_t deblocking_filter_override_flag : 1;
	uint32_t deblocking_filter_disabled_flag : 1;
	uint32_t slice_loop_filter_across_slices_enabled_flag : 1;
	h265d_cabac_t cabac;
} h265d_slice_header_body_t;

typedef struct {
	uint8_t offset_len_minus1;
	uint32_t num_entry_point_offsets;
	uint32_t* entry_point_offset_minus1;
} h265d_entry_point_t;

typedef struct {
	uint8_t pps_id;
	uint32_t slice_segment_address;
	uint32_t slice_segment_header_extension_length;
	uint8_t first_slice_segment_in_pic_flag;
	uint8_t no_output_of_prior_pics_flag;
	uint8_t dependent_slice_segment_flag;
	h265d_slice_header_body_t body;
	h265d_entry_point_t entry_points;
} h265d_slice_header_t;

typedef struct {
	h265d_nal_t current_nal;
	int (*header_callback)(void *arg, void *seq_id);
	void *header_callback_arg;
	dec_bits stream_i;
	h265d_slice_header_t slice_header;
	h265d_vps_t vps;
	h265d_sps_t sps[16];
	h265d_pps_t pps[64];
} h265d_data_t;

struct h265d_context;
typedef struct h265d_context h265d_context;

#endif /* __H265DEC_H__ */
