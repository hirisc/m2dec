/** Analyze MPEG-2 sequence layer and below.
 *  Copyright 2008 Takayuki Minegishi
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

#ifndef _MPEG2_H_
#define _MPEG2_H_
  
#include "m2d.h"

#ifdef __cplusplus
extern "C" {
#endif

#define __LIBM2DEC_API


#define ALIGN16(x) (((x) + 15) & ~15)

#ifdef FAST_DECODE
#define MB_LEN 2
#else
#define MB_LEN 16
#endif

enum {
	I_VOP = 1,
	P_VOP,
	B_VOP,
	D_VOP
};

/**Element of VLD tables for DCT coefficients.
 */
typedef struct {
	int16_t run; /**< Number of preceding zero. Minus means ESC. */
	int8_t level; /**< Value of coefficients. Already signed. */
	int8_t length; /**< Negative or zero indicates additional look-up(s). */
} vlc_dct_t;

typedef struct {
	int id;
	int16_t horizontal_size_value;
	int16_t vertical_size_value;
	int16_t display_width;
	int16_t display_height;
	unsigned aspect_ratio_information : 4;
	unsigned frame_rate_code : 4;
	unsigned bit_rate_value : 18;
	unsigned constrained_parameters_flag : 1;
	unsigned load_intra_quantizer_matrix : 1;
	unsigned load_non_intra_quantizer_matrix : 1;
	unsigned vbv_buffer_size;
	unsigned profile_and_level : 8;
	unsigned progressive_sequence : 1;
	unsigned chroma_format : 2;
	unsigned low_delay : 1;
} m2d_seq_header;

typedef struct {
	int time_code;
	int16_t closed_gop;
	int16_t broken_link;
} m2d_gop_header;

typedef struct {
	int16_t temporal_reference;
	uint16_t vbv_delay;
	int8_t picture_coding_type;
	int8_t q_scale;

	unsigned int intra_dc_precision : 2;
	unsigned int picture_structure : 2;
	unsigned int top_field_first : 1;
	unsigned int frame_pred_frame_dct : 1;
	unsigned int concealment_motion_vectors : 1;
	unsigned int q_scale_type : 1;
	unsigned int intra_vlc_format : 1;
	unsigned int alternate_scan : 1;
	unsigned int repeat_first_field : 1;
	unsigned int chroma_420_type : 1;
	unsigned int progressive_frame : 1;
	unsigned int composite_display_flag : 1;

	unsigned int v_axis : 1;
	unsigned int field_sequence : 3;
	unsigned int sub_carrier : 1;
	unsigned int burst_amplitude : 7;
	unsigned int sub_carrier_phase : 8;
} m2d_picture;


enum {
	MAX_FRAME_NUM = 16
};

typedef struct {
	uint8_t *curr_luma;
	uint8_t *curr_chroma;
	ptrdiff_t diff_to_ref[2][2]; /* forward/backward x luma/chroma */
	int idx_of_ref[2];
	int num;
	int index;
	m2d_frame_t frames[MAX_FRAME_NUM];
	int lru[MAX_FRAME_NUM];
} m2d_frames;


enum {
	MV_FIELD,
	MV_FRAME,
	MV_DUALPRIME,
	MV_16_8MC
};

typedef struct {
	int8_t mv_count;
	int8_t format;
	int8_t dmv;
	int8_t type;
} m2d_motion_type_t;

typedef struct {
	int16_t mv[2][2];
} m2d_mv_t;

typedef struct m2d_mb_current_t {
	int16_t frame_width; /**< Width of frame (multiples of 16) */
	int16_t frame_height; /**< Height of frame (multiples of 16 or 32) */
	uint8_t mbmax_x;
	uint8_t mbmax_y;
	int8_t mb_x;
	int8_t mb_y;
	int8_t type;
	int8_t intra_dc_scale;
	int16_t intra_dc_max;
	int16_t q_scale;
	int16_t dc_pred[3];
	int16_t dct_type;
	int8_t r_size[2][2]; /**< f_code - 1 */
	m2d_mv_t mv[2];
	m2d_frames *frames;
	const m2d_motion_type_t *motion_type;
	const int *q_mapping;
	const uint8_t *qmat[4];
	const int8_t *zigzag;
	int frame_mode;
	int (* const *parse_coef)(struct m2d_mb_current_t *mb, dec_bits *stream, int idx);
	void (*skip_mb)(struct m2d_mb_current_t *mb, int mb_increment);
	int (*macroblock_type)(dec_bits *stream); /**< switch according to frame */
	int intra_vlc_format;
	int16_t coef[MB_LEN * MB_LEN / 4];
} m2d_mb_current;

enum {
	MB_FORWARD = 0x01,
	MB_BACKWARD = 0x02,
	MB_INTRA = 0x04,
	MB_PATTERN = 0x08,
	MB_QUANT = 0x10,
	MB_SPACIAL = 0x20,
	MB_MC = MB_BACKWARD | MB_FORWARD,
	MB_COEF = MB_PATTERN | MB_INTRA
};


#define MB_PARAM(type, flg) ((type) & (flg))

typedef struct {
	int id;
	m2d_seq_header *seq_header;
	m2d_picture *picture;
	dec_bits *stream;
	m2d_mb_current *mb_current;
	m2d_gop_header *gop_header;
	int (*header_callback)(void *arg, int seq_id);
	void *header_callback_arg;
	int out_state;
	m2d_seq_header seq_header_i;
	m2d_gop_header gop_header_i;
	dec_bits stream_i;
	m2d_picture picture_i;
	m2d_mb_current mb_current_i;
	m2d_frames frames_i;
	uint8_t qmat[4][64];
} m2d_context;

enum {
	M2D_ERR_USAGE = -1,
	M2D_ERR_NOT_VIDEO = -2,
	M2D_ERR_BIT = -3,
	M2D_ERR_UNKNOWN = -4
};

__LIBM2DEC_API int m2d_init(m2d_context *m2d, int dummy, int (*header_callback)(void *arg, int seq_id), void *arg);
__LIBM2DEC_API int m2d_read_header(m2d_context *m2d, const byte_t *data, size_t len);
__LIBM2DEC_API int m2d_get_info(m2d_context *m2d, m2d_info_t *info);
__LIBM2DEC_API int m2d_set_frames(m2d_context *m2d, int num_mem, m2d_frame_t *mem);
__LIBM2DEC_API int m2d_set_data(m2d_context *m2d, const byte_t *indata, int indata_bytes);
__LIBM2DEC_API int m2d_decode_data(m2d_context *m2d);
__LIBM2DEC_API int m2d_get_decoded_frame(m2d_context *m2d, m2d_frame_t *frame, int is_end);
__LIBM2DEC_API int m2d_skip_frames(m2d_context *m2d, int frame_num);

extern const m2d_func_table_t * const m2d_func;

#ifdef __cplusplus
}
#endif

#endif /* _MPEG2_H_ */

