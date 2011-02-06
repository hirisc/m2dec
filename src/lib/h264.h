#ifndef __H264_H__
#define __H264_H__

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _M_IX86
#include <crtdbg.h>
#define VC_CHECK assert(_CrtCheckMemory());

#ifdef _WINDLL
#include <Windows.h>
#endif
#else
#define VC_CHECK
#endif

#ifndef NDEBUG
#include <stdio.h>
#endif

#include "types.h"
#define __LIBH264DEC_API

enum {
	NOT_IN_USE = 0,
	SHORT_TERM = 1,
	LONG_TERM = 2,
	UNDEFINED_NAL = 0,
	SLICE_NONIDR_NAL = 1,
	SLICE_PARTA_NAL = 2,
	SLICE_PARTB_NAL = 3,
	SLICE_PARTC_NAL = 4,
	SLICE_IDR_NAL = 5,
	SEI_NAL = 6,
	SPS_NAL = 7,
	PPS_NAL = 8,
	AUDELIM_NAL = 9,
	EOS_NAL = 10,
	EOSTR_NAL = 11,
	FILLER_NAL = 12,
	P_SLICE = 0,
	B_SLICE = 1,
	I_SLICE = 2,
	SP_SLICE = 3,
	SI_SLICE = 4,
	MB_INxN = 0,
	MB_I16x16 = 1,
	MB_IPCM = 25,
	MB_P16x16 = 26,
	MB_P8x8 = 29,
	MB_P8x8REF0 = 30,
	MB_PSKIP = 31,
	MB_BDIRECT16x16 = 32,
	EXTENDED_SAR = 255
};

typedef struct {
	int8_t cpb_cnt_minus1; /* [0, 31] */
	uint8_t bit_rate_scale;
	uint8_t cpb_size_scale;
	uint8_t initial_cpb_removal_delay_length_minus1;
	uint8_t cpb_removal_delay_length_minus1;
	uint8_t dpb_output_delay_length_minus1;
	uint8_t time_offset_length;
	uint32_t bit_rate_value_minus1[32];
	uint32_t cpb_size_value_minus1[32];
	uint32_t cbr_flag; /* bitfield mapping, LSB-first */
} hrd_parameters_t;

typedef struct {
	unsigned aspect_ratio_info_present_flag : 1;
	unsigned overscan_info_present_flag : 1;
	unsigned overscan_appropriate_flag : 1;
	unsigned video_signal_type_present_flag : 1;
	unsigned video_full_range_flag : 1;
	unsigned colour_description_present_flag : 1;
	unsigned chroma_loc_info_present_flag : 1;
	unsigned timing_info_present_flag : 1;
	unsigned fixed_frame_rate_flag : 1;
	unsigned nal_hrd_parameters_present_flag : 1;
	unsigned vcl_hrd_parameters_present_flag : 1;
	unsigned low_delay_hrd_flag : 1;
	unsigned pic_struct_present_flag : 1;
	unsigned bitstream_restriction_flag : 1;
	unsigned motion_vectors_over_pic_boundaries_flag : 1;
	uint8_t aspect_ratio_idc;
	uint8_t video_format;
	uint8_t colour_primaries;
	uint8_t transfer_characteristics;
	uint8_t matrix_coefficients;
	uint16_t sar_width;
	uint16_t sar_height;
	uint32_t chroma_sample_loc_type_top_field;
	uint32_t chroma_sample_loc_type_bottom_field;
	uint32_t num_units_in_tick;
	uint32_t time_scale;
	uint32_t max_bytes_per_pic_denom;
	uint32_t max_bits_per_mb_denom;
	uint32_t log2_max_mv_length_horizontal;
	uint32_t log2_max_mv_length_vertical;
	uint32_t num_reorder_frames;
	uint32_t max_dec_frame_buffering;
	hrd_parameters_t nal_hrd_parameters, vcl_hrd_parameters;
} vui_parameters_t;
	
typedef struct {
	uint8_t profile_idc;
	uint8_t level_idc;
	uint8_t poc_type;
	int8_t log2_max_frame_num;
	int8_t log2_max_poc_lsb;
	uint8_t num_ref_frames;
	int16_t pic_width;
	int16_t pic_height;
	int16_t frame_crop[4];
	unsigned constraint_set_flag : 8; /* reserved_zero_bits included */
	unsigned delta_pic_order_always_zero_flag : 1;
	unsigned gaps_in_frame_num_value_allowed_flag : 1;
	unsigned frame_mbs_only_flag : 1;
	unsigned mb_adaptive_frame_field_flag : 1;
	unsigned direct_8x8_inference_flag : 1;
	unsigned frame_cropping_flag : 1;
	unsigned vui_parameters_present_flag : 1;
	int32_t offset_for_non_ref_pic;
	int32_t offset_for_top_to_bottom_field;
	uint32_t num_ref_frames_in_pic_order_cnt_cycle;
	vui_parameters_t vui;
	int32_t offset_for_ref_frame[256];
} h264d_sps;

typedef struct {
	unsigned entropy_coding_mode_flag : 1;
	unsigned pic_order_present_flag : 1;
	unsigned weighted_pred_flag : 1;
	unsigned deblocking_filter_control_present_flag : 1;
	unsigned constrained_intra_pred_flag : 1;
	unsigned redundant_pic_cnt_present_flag : 1;
	unsigned transform_8x8_mode_flag : 1;
	unsigned pic_scaling_matrix_present_flag : 1;
	int8_t num_ref_idx_l0_active_minus1;
	int8_t num_ref_idx_l1_active_minus1;
	int8_t weighted_bipred_idc;
	int8_t pic_init_qp;
	int8_t pic_init_qs;
	int8_t chroma_qp_index[2];
	uint8_t pic_scaling_list_present_flag; /* bitfield, LSB first */
	uint32_t seq_parameter_set_id;
	uint32_t num_slice_groups_minus1;
/*	uint32_t slice_group_map_type; */
} h264d_pps;

typedef struct {
	uint8_t *luma;
	uint8_t *chroma;
	int32_t luma_len;
} h264d_frame;

typedef struct {
	int16_t in_use; /* 0, 1, 2 = unused, short_term, long_term */
	int16_t frame_idx;
	uint32_t num;
} h264d_ref_frame_t;

typedef struct {
	int8_t ref_pic_list_reordering_flag;
	h264d_ref_frame_t ref_frames[16];
} h264d_reorder_t;

typedef struct {
	int8_t op;
	uint32_t arg1, arg2;
} h264d_mmco;

typedef struct {
	union {
		struct {
			unsigned no_output_of_prior_pic_flag : 1;
			unsigned long_term_reference_flag : 1;
		} idr;
		struct {
			int8_t adaptive_ref_pic_marking_mode_flag;
			h264d_mmco mmco[16];
		} non_idr;
	};
} h264d_marking_t;

typedef struct {
	unsigned field_pic_flag : 1;
	unsigned bottom_field_flag : 1;
	unsigned direct_spatial_mv_pred_flag : 1;
	unsigned num_ref_idx_active_override_flag : 1;
	unsigned sp_for_switch_flag : 1;
	uint8_t pic_parameter_set_id;
	int8_t slice_type;
	int8_t num_ref_idx_l0_active_minus1;
	int8_t num_ref_idx_l1_active_minus1;
	int8_t cabac_init_idc;
	int8_t qp_delta;
	int8_t qs_delta;
	int8_t disable_deblocking_filter_idc;
	int8_t slice_alpha_c0_offset_div2;
	int8_t slice_beta_offset_div2;
	uint16_t idr_pic_id;
	uint32_t first_mb_in_slice;
	uint32_t frame_num;
	uint32_t prev_frame_num;
	union {
		struct {
			uint32_t lsb;
			uint32_t msb;
			int32_t delta_pic_order_cnt_bottom;
		} poc0;
		struct {
			int32_t delta_pic_order_cnt[2];
		} poc1;
	};
	uint32_t poc;
	uint32_t redundant_pic_cnt;
	h264d_reorder_t reorder[2];
	h264d_marking_t marking;
} h264d_slice_header;

typedef struct {
	uint8_t *curr_luma;
	uint8_t *curr_chroma;
	int num;
	int index;
	h264d_ref_frame_t *refs;
	h264d_frame frames[32];
	int lru[32];
} h264d_frame_info_t;

typedef struct {
	int8_t type;
	int8_t chroma_pred_mode;
	int8_t cbp;
	int8_t ref[2];
	uint16_t cbf;
	int16_t mv[2][4][2];
} prev_mb_t;

typedef struct {
	int8_t idc, qpy, qpc, slicehdr;
	uint8_t str4_vert, str4_horiz;
	uint32_t str_vert, str_horiz;
} deblock_info_t;

typedef struct {
	uint32_t range;
	uint32_t offset;
	int8_t context[460];
} h264d_cabac_t;

#define ENC_SLICEHDR(hdr, a, b) (hdr = ((((b) + 6) << 4) | ((a) + 6)))
#define DEC_SLICEHDR(hdr, a, b) (a = ((hdr & 15) - 6) * 2), (b = (((uint8_t)hdr >> 4) - 6) * 2)

typedef struct {
	int8_t is_constrained_intra;
	int8_t type;
	int8_t qp, qp_chroma;
	int8_t lefttop_ref;
	int8_t prev_qp_delta;
	int8_t chroma_pred_mode;
	int16_t x, y;
	int16_t max_x, max_y;
	int16_t firstline; /* # of first line of MBs */
	int16_t lefttop_mv[2];
	uint8_t *luma; /* current destination point */
	uint8_t *chroma;
	int32_t left4x4pred;
	int32_t left4x4coef;
	int32_t *top4x4pred;
	int32_t *top4x4coef;
	prev_mb_t *left4x4inter;
	prev_mb_t *top4x4inter;
	uint32_t cbp, cbf;
	deblock_info_t *deblock_curr;
	deblock_info_t *deblock_base;
	h264d_frame_info_t *frame;
	int32_t *top4x4pred_base;
	int32_t *top4x4coef_base;
	prev_mb_t *mb_base;
	h264d_cabac_t *cabac;
	h264d_pps *pps;
	int8_t *num_ref_idx_l0_active_minus1;
	int8_t *num_ref_idx_l1_active_minus1;
	int16_t *qmatc_p;
	int16_t qmaty[16];
	int16_t qmatc[16];
	int offset4x4[16]; /* offset of each 4x4 block in a macroblock. */
	h264d_frame_info_t frame_i;
	h264d_cabac_t cabac_i;
} h264d_mb_current;

typedef struct mb_code {
	int (*mb_dec)(h264d_mb_current *mb, const struct mb_code *mbc, dec_bits *st, int avail);
	int (*mb_pred)(uint8_t *dst, int stride, int avail);
	int8_t cbp;
} mb_code;

typedef struct {
	int src_width, src_height;
	int disp_width, disp_height;
	int frame_num;
	int crop[4];
	int additional_size;
} h264d_info_t;

typedef struct {
	int id;
	h264d_slice_header *slice_header;
	dec_bits *stream;
	dec_bits stream_i;
	h264d_slice_header slice_header_i;
	h264d_mb_current mb_current;
	h264d_pps pps_i[256];
	h264d_sps sps_i[32];
} h264d_context;


int h264d_init(h264d_context *h2d);
int h264d_read_header(h264d_context *h2d, const byte_t *data, size_t len);
int h264d_get_info(h264d_context *h2d, h264d_info_t *info);
int h264d_set_frames(h264d_context *h2d, int num_frame, h264d_frame *frame, uint8_t *second_frame);
int h264d_decode_picture(h264d_context *h2d);
int h246d_get_decoded_frame(h264d_context *h2d, uint8_t **luma, uint8_t **chroma);
void h264d_load_bytes_skip03(dec_bits *ths, intptr_t read_bytes);

#ifdef __cplusplus
}
#endif

#endif /*  __H264_H__ */
