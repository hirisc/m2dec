/** Yet Another H.265 decoder
 *  Copyright 2016 Takayuki Minegishi
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

#ifndef __H265MODULES_H__
#define __H265MODULES_H__

#ifdef _MSC_VER
#include <crtdbg.h>
#define VC_CHECK assert(_CrtCheckMemory());
#else
#define VC_CHECK
#endif
#include "h265.h"

#define MINV(a, b) (((a) <= (b)) ? (a) : (b))
#define MAXV(a, b) (((a) >= (b)) ? (a) : (b))
#define CLIP3(low, high, val) (((low) <= (val)) ? (((val) <= (high)) ? (val) : (high)) : (low))

static const int H265D_MAX_FRAME_NUM = 8;
static const int H265D_NEIGHBOUR_NUM = 16;

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
	uint8_t scaling_list_pred_mode_flag[4][6];
	uint32_t scaling_list_pred_matrix_id_delta[4][6];
	uint8_t scale0[6][16];
	uint8_t scale1[6][64];
	uint8_t scale2[6][64];
	uint8_t scale3[2][64];
} h265d_scaling_list_data_t;

typedef struct {
	uint8_t num_pics;
	uint16_t used_by_curr_pic_flag;
	int16_t delta_poc[16];
} h265d_short_term_ref_pic_elem_t;

typedef struct {
	uint8_t total_curr;
	h265d_short_term_ref_pic_elem_t ref[2];
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
	uint8_t size_log2_min;
	uint8_t pcm_log2;
	uint8_t pcm_log2_min;
	uint8_t transform_log2;
	uint8_t transform_log2_min;
	uint8_t num_ctb_log2;
	uint16_t stride;
	uint16_t columns;
	uint16_t rows;
} h265d_ctb_info_t;

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
	uint16_t cropping[4];
	h265d_ctb_info_t ctb_info;
	h265d_sps_prefix_t prefix;
	h265d_sub_layer_reordering_info_t max_buffering[8];
	h265d_scaling_list_data_t scaling_list_data;
	h265d_short_term_ref_pic_set_t short_term_ref_pic_set[64];
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
	uint8_t num_ref_idx_lx_default_active_minus1[2];
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
	TRAIL_N = 0,
	TRAIL_R = 1,
	BLA_W_LP = 16,
	BLA_N_LP = 18,
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

typedef enum {
	INTRA_PLANAR = 0,
	INTRA_DC = 1,
	INTRA_ANGULAR2 = 2,
	INTRA_ANGULAR3,
	INTRA_ANGULAR4,
	INTRA_ANGULAR5,
	INTRA_ANGULAR6,
	INTRA_ANGULAR7,
	INTRA_ANGULAR8,
	INTRA_ANGULAR9,
	INTRA_ANGULAR10,
	INTRA_ANGULAR11,
	INTRA_ANGULAR12,
	INTRA_ANGULAR13,
	INTRA_ANGULAR14,
	INTRA_ANGULAR15,
	INTRA_ANGULAR16,
	INTRA_ANGULAR17,
	INTRA_ANGULAR18,
	INTRA_ANGULAR19,
	INTRA_ANGULAR20,
	INTRA_ANGULAR21,
	INTRA_ANGULAR22,
	INTRA_ANGULAR23,
	INTRA_ANGULAR24,
	INTRA_ANGULAR25,
	INTRA_ANGULAR26,
	INTRA_ANGULAR27,
	INTRA_ANGULAR28,
	INTRA_ANGULAR29,
	INTRA_ANGULAR30,
	INTRA_ANGULAR31,
	INTRA_ANGULAR32,
	INTRA_ANGULAR33,
	INTRA_ANGULAR34,
	INTRA_STRONG_FILTER = 64
} h265d_intra_pred_mode_t;

typedef enum {
	PART_2Nx2N = 0,
	PART_2NxN = 1,
	PART_Nx2N = 2,
	PART_NxN = 3,
	PART_2NxnU = 4,
	PART_2NxnD = 5,
	PART_nLx2N = 6,
	PART_nRx2N = 7
} h265d_inter_part_mode_t;

typedef struct {
	int8_t sao_merge_flag[1];
	int8_t sao_type_idx[1];
	int8_t split_cu_flag[3];
	int8_t cu_transquant_bypass_flag[1];
	int8_t cu_skip_flag[3];
	int8_t pred_mode_flag[1];
	int8_t part_mode[4];
	int8_t prev_intra_luma_pred_flag[1];
	int8_t intra_chroma_pred_mode[1];
	int8_t rqt_root_cbf[1];
	int8_t merge_flag[1];
	int8_t merge_idx[1];
	int8_t inter_pred_idc[5];
	int8_t ref_idx_lx[2];
	int8_t mvp_flag[1];
	int8_t split_transform_flag[3];
	int8_t cbf_luma[2];
	int8_t cbf_chroma[4];
	int8_t abs_mvd_greater_flag[2];
	int8_t cu_qp_delta_abs[2];
	int8_t transform_skip_flag[2];
	int8_t last_sig_coeff_x_prefix[18];
	int8_t last_sig_coeff_y_prefix[18];
	int8_t coded_sub_block_flag[4];
	int8_t sig_coeff_flag[42];
	int8_t coeff_abs_level_greater1_flag[24];
	int8_t coeff_abs_level_greater2_flag[6];
} h265d_cabac_context_t;

typedef struct {
	m2d_cabac_t cabac;
	h265d_cabac_context_t* context;
} h265d_cabac_t;

typedef struct {
	int8_t offset[4];
	union {
		uint8_t band_pos;
		uint8_t edge;
	} opt;
} h265d_sao_map_elem_t;

typedef struct h265d_sao_map_t {
	uint8_t merge_left : 1;
	uint8_t luma_idx : 2;
	uint8_t chroma_idx : 2;
	h265d_sao_map_elem_t elem[3];
} h265d_sao_map_t;

typedef struct {
	uint8_t* bottom;
	uint8_t* reserved_flag;
} h265d_sao_hlines_t;

typedef struct {
	uint8_t* right[2][2];
	uint8_t reserved_flag;
} h265d_sao_vlines_t;

typedef struct {
	uint16_t lsb;
	uint16_t msb;
	int32_t poc;
} h265d_poc_t;

typedef struct {
	int32_t poc;
	int32_t num;
	int16_t frame_idx;
	bool in_use;
	bool is_longterm;
} h265d_ref_pic_list_elem_t;

typedef struct {
	h265d_nal_t nal_type;
	uint8_t slice_type;
	int8_t slice_qpy;
	int8_t slice_qpc_delta[2];
	int8_t slice_beta_offset_div2;
	int8_t slice_tc_offset_div2;
	uint8_t num_ref_idx_lx_active_minus1[2];
	uint8_t collocated_ref_idx;
	int8_t max_num_merge_cand;
	uint32_t mvd_l1_zero_flag : 1;
	uint32_t colocated_from_l0_flag : 1;
	uint32_t pic_output_flag : 1;
	uint32_t colour_plane_id : 2;
	uint32_t slice_temporal_mvp_enabled_flag : 1;
	uint32_t slice_sao_luma_flag : 1;
	uint32_t slice_sao_chroma_flag : 1;
	uint32_t cabac_init_flag : 1;
	uint32_t deblocking_filter_override_flag : 1;
	uint32_t deblocking_filter_disabled_flag : 1;
	uint32_t slice_loop_filter_across_slices_enabled_flag : 1;
	h265d_poc_t slice_pic_order_cnt;
	h265d_ref_pic_list_elem_t ref_list[2][16];
	h265d_short_term_ref_pic_set_t short_term_ref_pic_set;
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

struct pred_info_t {
	int16_t mvd[2][2];
	int8_t ref_idx[2];
};

typedef struct {
	uint16_t pu_intra : 1;
	uint16_t pu_nonzero_coef : 1;
	uint16_t tu_intra : 1;
	uint16_t tu_nonzero_coef : 1;
	uint16_t skip : 1;
	uint16_t pred_mode : 6;
	uint16_t depth : 4;
	pred_info_t pred;
} h265d_neighbour_t;

typedef struct {
	int32_t scale;
	uint8_t** matrix;
} h265d_scaling_info_t;

typedef struct {
	uint8_t qp : 6;
	uint8_t str : 2;
} h265d_deblocking_strength_t;

typedef int16_t (* h265d_scaling_func_t)(int32_t val, const h265d_scaling_info_t& scale, int idx);

typedef struct {
	int poc;
	int8_t frame_idx;
	uint8_t in_use : 1;
	uint8_t is_longterm : 1;
	uint8_t is_idr : 1;
	uint8_t is_terminal : 1;
} h265d_dpb_elem_t;

typedef struct {
	int8_t size;
	int8_t max;
	int8_t output;
	int8_t is_ready;
	h265d_dpb_elem_t data[16];
} h265d_dpb_t;

typedef struct {
	uint8_t *curr_luma;
	uint8_t *curr_chroma;
	int8_t num;
	int8_t index;
	int32_t poc[H265D_MAX_FRAME_NUM];
	m2d_frame_t frames[H265D_MAX_FRAME_NUM];
	int8_t lru[H265D_MAX_FRAME_NUM];
	h265d_dpb_t dpb;
} h265d_frame_info_t;

class h265d_deblocking_t {
	int8_t leftgap[2];
	int8_t topgap[2];
	int8_t edgemax;
	bool disabled_;
	const uint8_t (*qp_history_)[16];
	const h265d_ref_pic_list_elem_t (*ref_list_)[16];
	h265d_deblocking_strength_t* topedge_;
	h265d_deblocking_strength_t* topedge_base_;
	h265d_deblocking_strength_t boundary_[2][8 * 17];
	static h265d_deblocking_strength_t* boundary_to_fill(h265d_deblocking_strength_t* deb, int offset_x, int offset_y, int xgap, int ygap, int& org_y) {
		int org_x = offset_x >> 3;
		org_y = offset_y >> 2;
		return deb + org_x * xgap + (org_y + 1) * ygap;
	}
	void record_tu_intra_onedir(int qpy, int dir, int offset_x, int offset_y, int unavail, int xgap, int ygap, int len, int edgemax) {
		if ((offset_x & 7) || ((offset_x == 0) && ((unavail >> dir) & 1))) {
			return;
		}
		int org_y;
		h265d_deblocking_strength_t* boundary = boundary_to_fill(boundary_[dir], offset_x, offset_y, xgap, ygap, org_y);
		int qp = qpy + 1;
		const uint8_t* prev_qp = &qp_history_[dir][org_y];
		for (int n = 0; n < len; ++n) {
			boundary[n * ygap].qp = (qp + prev_qp[n]) >> 1;
			boundary[n * ygap].str = 2;
		}
	}
	static int ref_idx_num(const pred_info_t& pred) {
		return (pred.ref_idx[1] < 0) ? 1 : (1 + (0 <= pred.ref_idx[0]));
	}
	static int strength_tu(const h265d_neighbour_t& neighbour) {
		return neighbour.tu_intra ? 2 : (neighbour.tu_nonzero_coef ? 1 : 0);
	}
	void record_tu_onedir(int qpy, int dir, int offset_x, int offset_y, int unavail, int xgap, int ygap, int len, int edgemax, int strength, const h265d_neighbour_t neighbour[]) {
		if ((offset_x & 7) || ((offset_x == 0) && ((unavail >> dir) & 1))) {
			return;
		}
		int org_y;
		h265d_deblocking_strength_t* boundary = boundary_to_fill(boundary_[dir], offset_x, offset_y, xgap, ygap, org_y);
		int qp = qpy + 1;
		const uint8_t* prev_qp = &qp_history_[dir][org_y];
		int prev_str = -1;
		for (int n = 0; n < len; ++n) {
			boundary[n * ygap].qp = (qp + prev_qp[n]) >> 1;
			int str_tu = MAXV(strength, strength_tu(neighbour[n]));
			boundary[n * ygap].str = MAXV(boundary[n * ygap].str, str_tu);
		}
	}
	static inline int DIF_SQUARE(int a, int b) {
		int t = a - b;
		return t * t;
	}
	static bool mv_diff_is_large(const int16_t mvxy_a[], const int16_t mvxy_b[]) {
		return (16 <= DIF_SQUARE(mvxy_a[0], mvxy_b[0])) || (16 <= DIF_SQUARE(mvxy_a[1], mvxy_b[1]));
	}
	static int inter_strength(int nfrm0, int nfrm1, int cfrm0, int cfrm1, const int16_t n_mvxy[][2], const int16_t c_mvxy[][2], int n_swapped, int c_swapped) {
		if ((nfrm0 != cfrm0) || (nfrm1 != cfrm1)) {
			return 1;
		} else {
			if (nfrm0 == nfrm1) {
				return (mv_diff_is_large(n_mvxy[0], c_mvxy[0]) || mv_diff_is_large(n_mvxy[1], c_mvxy[1])) && (mv_diff_is_large(n_mvxy[0], c_mvxy[1]) || mv_diff_is_large(n_mvxy[0], c_mvxy[1]));
			} else {
				return ((0 <= nfrm0) && mv_diff_is_large(n_mvxy[n_swapped], c_mvxy[c_swapped])) || ((0 <= nfrm1) && mv_diff_is_large(n_mvxy[n_swapped ^ 1], c_mvxy[c_swapped ^ 1]));
			}
		}
	}
	int refidx_to_frameidx(int refidx, int lx) {
		return (0 <= refidx) ? ref_list_[lx][refidx].frame_idx : -1;
	}
	int strength_inter_pu_onedir(int qpy, const h265d_neighbour_t& neighbour, const int16_t mvxy[][2], int frm0, int frm1, int c_swapped) {
		if (neighbour.pu_intra) {
			return 2;
		} else if (neighbour.pu_nonzero_coef) {
			return 1;
		} else {
			int nfrm0 = refidx_to_frameidx(neighbour.pred.ref_idx[0], 0);
			int nfrm1 = refidx_to_frameidx(neighbour.pred.ref_idx[1], 1);
			int n_swapped = 0;
			if (nfrm0 < nfrm1) {
				std::swap(nfrm0, nfrm1);
				n_swapped = 1;
			}
			return inter_strength(nfrm0, nfrm1, frm0, frm1, neighbour.pred.mvd, mvxy, c_swapped, n_swapped);
		}
	}
	void record_pu_onedir(int qpy, int dir, int offset_x, int offset_y, uint32_t unavail, int xgap, int ygap, int len, int edgemax, const h265d_neighbour_t* neighbour, int refidx0, int refidx1, const int16_t mvxy[][2]) {
		if ((offset_x & 7) || ((offset_x == 0) && ((unavail >> dir) & 1))) {
			return;
		}
		int frm0 = refidx_to_frameidx(refidx0, 0);
		int frm1 = refidx_to_frameidx(refidx1, 1);
		int c_swapped = 0;
		if (frm0 < frm1) {
			std::swap(frm0, frm1);
			c_swapped = 1;
		}
		int org_y;
		h265d_deblocking_strength_t* boundary = boundary_to_fill(boundary_[dir], offset_x, offset_y, xgap, ygap, org_y);
		int qp = qpy + 1;
		const uint8_t* prev_qp = &qp_history_[dir][org_y];
		len >>= 2;
		for (int i = 0; i < len; ++i) {
			boundary[i * ygap].qp = (qp + prev_qp[i]) >> 1;
			boundary[i * ygap].str = strength_inter_pu_onedir(qpy, neighbour[i], mvxy, frm0, frm1, c_swapped);
		}
	}
	void clear_left() {
		h265d_deblocking_strength_t* left = &boundary_[1][0];
		int edgenum = edgemax;
		int len = edgenum * 2;
		for (int n = 0; n < edgenum; ++n) {
			left[0] = left[len];
			memset(left + 1, 0, len * sizeof(left[0]));
			left = left + len + 1;
		}
	}
public:
	void init(uint8_t* buf) {
		topedge_ = topedge_base_ = reinterpret_cast<h265d_deblocking_strength_t*>(buf);
	}
	size_t allocate_size(int log2size, int columns) {
		return sizeof(topedge_[0]) * (columns << (log2size - 3));
	}
	void set_ctu(bool disable, const uint8_t qp_history[][16], const h265d_ref_pic_list_elem_t ref_list[][16], int log2size, int columns, int pos_x) {
		disabled_ = disable;
		qp_history_ = qp_history;
		ref_list_ = ref_list;
		edgemax = 1 << (log2size - 3);
		leftgap[0] = 1;
		leftgap[1] = edgemax;
		topgap[0] = edgemax * 2 + 1;
		topgap[1] = 1;
		topedge_ = topedge_base_ + (pos_x << (log2size - 3));
		memset(boundary_, 0, sizeof(boundary_));
		memset(topedge_base_, 0, sizeof(topedge_base_[0]) * (columns << (log2size - 3)));
	}

	void set_pos(int log2size, int pos_x) {
		topedge_ = topedge_base_ + (pos_x << (log2size - 3));
	}

	void record_tu_intra(int qpy, int size_log2, int offset_x, int offset_y, uint32_t unavail) {
		if (disabled_) {
			return;
		}
		int len = 1 << (size_log2 - 2);
		record_tu_intra_onedir(qpy, 0, offset_x, offset_y, unavail, leftgap[0], leftgap[1], len, edgemax);
		record_tu_intra_onedir(qpy, 1, offset_y, offset_x, unavail, topgap[0], topgap[1], len, edgemax);
	}

	void record_tu(int qpy, int size_log2, int offset_x, int offset_y, uint32_t unavail, int str, const h265d_neighbour_t* left, const h265d_neighbour_t* top) {
		if (disabled_) {
			return;
		}
		int len = 1 << (size_log2 - 2);
		record_tu_onedir(qpy, 0, offset_x, offset_y, unavail, leftgap[0], leftgap[1], len, edgemax, str, left);
		record_tu_onedir(qpy, 1, offset_y, offset_x, unavail, topgap[0], topgap[1], len, edgemax, str, top);
	}

	void record_pu(int qpy, int width, int height, int offset_x, int offset_y, uint32_t unavail, const h265d_neighbour_t* left, const h265d_neighbour_t* top, int refidx0, int refidx1, const int16_t mvxy[][2]) {
		if (disabled_) {
			return;
		}
		record_pu_onedir(qpy, 0, offset_x, offset_y, unavail, leftgap[0], leftgap[1], height, edgemax, left, refidx0, refidx1, mvxy);
		record_pu_onedir(qpy, 1, offset_y, offset_x, unavail, topgap[0], topgap[1], width, edgemax, top, refidx0, refidx1, mvxy);
	}

	void pre_deblocking() {
		memcpy(boundary_[0], topedge_, edgemax * sizeof(topedge_[0]));
	}

	void post_deblocking(int pos_x, int columns) {
		if (pos_x < columns - 1) {
			clear_left();
		} else {
			memset(boundary_[1], 0, sizeof(boundary_[1]));
		}
		int store_size = edgemax * sizeof(topedge_[0]);
		memcpy(topedge_, boundary_[0] + edgemax * edgemax * 2, store_size);
		memset(boundary_[0] + store_size, 0, sizeof(boundary_[0]) - store_size);
	}

	const h265d_deblocking_strength_t* boundary(int dir) const {
		return boundary_[dir];
	}
};

class frameidx_record_t {
	uint64_t frameidx_[2];
public:
	int frameidx(int lx, int refidx) const {
		return (frameidx_[lx] >> (refidx * 4)) & 7;
	}

	const frameidx_record_t& operator=(const h265d_ref_pic_list_elem_t reflist[][16]) {
		uint64_t v0 = 0;
		uint64_t v1 = 0;
		for (int i = 0; i < 16; ++i) {
			v0 |= reflist[0][i].frame_idx << (i * 4);
			v1 |= reflist[1][i].frame_idx << (i * 4);
		}
		frameidx_[0] = v0;
		frameidx_[1] = v1;
		return *this;
	}
};

class temporal_mvscale_t {
	static int scale(int poc0, int refpoc0, int poc1, int refpoc1) {
		int diff1 = poc1 - refpoc1;
		int diff0 = poc0 - refpoc0;
		if (diff1 == 0) {
			return 4096;
		}
		int td = CLIP3(-128, 127, diff1);
		int tb = CLIP3(-128, 127, diff0);
		int tx = (16384 + (std::abs(td) >> 1)) / td;
		int scale = (tb * tx + 32) >> 6;
		return CLIP3(-4096, 4095, scale);
	}
public:
	int16_t scales_[8][8];
	void init(const h265d_frame_info_t& frm, int poc0, int poc1) {
		for (int i = 0; i < 8; ++i) {
			int refpoc0 = frm.poc[i];
			for (int j = 0; j < 8; ++j) {
				int refpoc1 = frm.poc[j];
				scales_[i][j] = scale(poc0, refpoc0, poc1, refpoc1);
			}
		}
	}
};

class temporal_mvscale_index_t {
	const frameidx_record_t* ref_curr_;
	const frameidx_record_t* ref_col_;
	temporal_mvscale_t colmv_scale_, tmv_scale_;
public:
	void init(const h265d_frame_info_t& frm, const frameidx_record_t& ref0, const frameidx_record_t& ref1, int poc0, int poc1) {
		ref_curr_ = &ref0;
		ref_col_ = &ref1;
		colmv_scale_.init(frm, poc0, poc1);
		tmv_scale_.init(frm, poc0, poc0);
	}

	int colmv_scale(int lx_a, int refidx_a, int lx_b, int refidx_b) const {
		return colmv_scale_.scales_[ref_curr_->frameidx(lx_a, refidx_a)][ref_col_->frameidx(lx_b, refidx_b)];
	}

	int tmv_scale(int lx_a, int refidx_a, int lx_b, int refidx_b) const {
		return tmv_scale_.scales_[ref_curr_->frameidx(lx_a, refidx_a)][ref_curr_->frameidx(lx_b, refidx_b)];
	}
};

class colpics_t {
	h265d_neighbour_t* curr_top;
	const h265d_neighbour_t* ref_top;
	h265d_neighbour_t* curr_pos;
	int8_t ctu_log2size_, ctu_size_;
	int16_t stride_;
	int width_, height_;
	h265d_neighbour_t* colpics[H265D_MAX_FRAME_NUM];
	frameidx_record_t frame_indices[H265D_MAX_FRAME_NUM];
	temporal_mvscale_index_t tmv_scale_;
public:
	void register_reflist(int frame_idx, const h265d_ref_pic_list_elem_t reflist[][16]) {
		frame_indices[frame_idx] = reflist;
	}
	size_t colpic_size(int width, int height) const {
		return sizeof(*colpics[0]) * ((width + 15) >> 4) * ((height + 15) >> 4);
	}
	void init(const h265d_slice_header_body_t& header, const h265d_frame_info_t& frm, int ctu_log2size, int width, int height, int pos_x, int pos_y, int poc) {
		curr_top = colpics[frm.index];
		const h265d_ref_pic_list_elem_t& col_pic = header.ref_list[header.colocated_from_l0_flag ^ 1][header.collocated_ref_idx];
		int col_frmidx = col_pic.frame_idx;
		ref_top = colpics[col_frmidx];
		ctu_log2size_ = ctu_log2size;
		ctu_size_ = 1 << (ctu_log2size_ - 4);
		stride_ = (width + 15) >> 4;
		width_ = width;
		height_ = height;
		set_curr_pos(pos_x, pos_y);
		register_reflist(frm.index, header.ref_list);
		if (header.slice_type < 2) {
			tmv_scale_.init(frm, frame_indices[frm.index], frame_indices[col_frmidx], poc, col_pic.poc);
		}
	}
	void set_colpic(int idx, uint8_t* buf) {
		colpics[idx] = reinterpret_cast<h265d_neighbour_t*>(buf);
	}
	void set_curr_pos(int pos_x, int pos_y) {
		curr_pos = curr_top + (pos_y << (ctu_log2size_ - 4)) * stride_ + (pos_x << (ctu_log2size_ - 4));
	}
	void inc_curr_pos() {
		curr_pos = curr_pos + ctu_size_;
	}

	int ref_offset(int base_x, int base_y, int offset_x, int offset_y) const {
		return ((base_y + offset_y) >> 4) * stride_ + ((base_x + offset_x) >> 4);
	}

	const h265d_neighbour_t* get_ref(int pos_x, int pos_y, int offset_x, int offset_y, int width, int height) const {
		int bottom_right_x = offset_x + width;
		int bottom_right_y = offset_y + height;
		int base_x = pos_x << ctu_log2size_;
		int base_y = pos_y << ctu_log2size_;
		if (!(bottom_right_y >> ctu_log2size_) && (base_x + bottom_right_x < width_) && (base_y + bottom_right_y < height_)) {
			const h265d_neighbour_t* ref = ref_top + ref_offset(base_x, base_y, bottom_right_x, bottom_right_y);
			if (!ref->pu_intra) {
				return ref;
			}
		}
		bottom_right_x = offset_x + (width >> 1);
		bottom_right_y = offset_y + (height >> 1);
		return ref_top + ref_offset(base_x, base_y, bottom_right_x, bottom_right_y);
	}

	const temporal_mvscale_index_t& tmv_scale() const {
		return tmv_scale_;
	}

	class fill_intra {
	public:
		void operator()(h265d_neighbour_t& dst) {
			dst.pu_intra = 1;
		}
	};

	class fill_inter {
		const int16_t (*mvxy_)[2];
		int ref0_;
		int ref1_;
	public:
		fill_inter(const int16_t (*mvxy)[2], int ref0, int ref1) : mvxy_(mvxy), ref0_(ref0), ref1_(ref1) {}
		void operator()(h265d_neighbour_t& dst) {
			dst.pu_intra = 0;
			dst.pred.ref_idx[0] = ref0_;
			dst.pred.ref_idx[1] = ref1_;
			memcpy(dst.pred.mvd, mvxy_, sizeof(dst.pred.mvd));
		}
	};

	template <typename F>
	void fill(int offset_x, int offset_y, int width, int height, F Fill) {
		h265d_neighbour_t* dst = curr_pos;
		int xmax = (offset_x + width) >> 2;
		int xorg = offset_x >> 2;
		int ymax = (offset_y + height) >> 2;
		for (int y = offset_y >> 2; y < ymax; ++y) {
			if (!(y & 3)) {
				dst = curr_pos + (y >> 2) * stride_;
				for (int x = xorg; x < xmax; ++x) {
					if (!(x & 3)) {
						Fill(dst[x >> 2]);
					}
				}
			}
		}
	}
};

typedef struct h265d_ctu_t {
	h265d_cabac_t cabac;
	uint16_t pos_x, pos_y;
	uint16_t pos_in_slice;
	uint32_t valid_x, valid_y;
	uint32_t idx_in_slice;
	int8_t qp_delta_req, qpy;
	int8_t qpc_delta[2];
	h265d_scaling_info_t qp_scale[3];
	uint8_t intra_split : 1;
	uint8_t not_first_row : 1;
	int8_t order_luma[4], order_chroma;
	const h265d_ctb_info_t *size;
	const h265d_sps_t* sps;
	int16_t *coeff_buf;
	const h265d_slice_header_t* slice_header;
	const h265d_pps_t* pps;
	h265d_neighbour_t neighbour_left[H265D_NEIGHBOUR_NUM + 2];
	h265d_neighbour_t* neighbour_top; // use 16 bytes for each CTU
	uint8_t qp_history[2][16];
	const h265d_scaling_func_t* scaling_func;
	h265d_sao_map_t* sao_map;
	void (*sao_read)(struct h265d_ctu_t& dst, const h265d_slice_header_t& hdr, dec_bits& st);
	uint8_t* luma;
	uint8_t* chroma;
	colpics_t colpics;
	h265d_frame_info_t frame_info;
	h265d_cabac_context_t context;
	h265d_deblocking_t deblocking;
	h265d_sao_vlines_t sao_vlines;
	h265d_sao_hlines_t sao_hlines[2][2];
	uint8_t* sao_signbuf;
	int16_t pred_buffer[2][32 * 32];
	int16_t coeff_buffer[32 * 32 * 2 + 7];
} h265d_ctu_t;

typedef struct {
	h265d_nal_t current_nal;
	int (*header_callback)(void *arg, void *seq_id);
	void *header_callback_arg;
	dec_bits stream_i;
	h265d_ctu_t coding_tree_unit;
	h265d_slice_header_t slice_header;
	h265d_vps_t vps;
	h265d_sps_t sps[16];
	h265d_pps_t pps[64];
} h265d_data_t;


#endif /* __H265MODULES_H__ */
