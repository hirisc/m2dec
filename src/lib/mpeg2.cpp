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

#include <assert.h>
#include <string.h>
#include <setjmp.h>
#include "bitio.h"
#include "mpeg2.h"
#include "vld.h"
#include "motioncomp.h"

#ifdef FAST_DECODE
#define m2d_idct_intra_luma(dst, stride, coef, sum) ((dst)[0] = (((coef)[0] + 4) >> 3))
#define m2d_idct_intra_chroma(dst, stride, coef, sum) m2d_idct_intra_luma(dst, stride, coef, sum)
#define m2d_idct_inter_luma(dst, stride, coef, sum) ((dst)[0] += (((coef)[0] + 4) >> 3))
#define m2d_idct_inter_chroma(dst, stride, coef, sum) m2d_idct_inter_luma(dst, stride, coef, sum)
#define m2d_prefetch_mb(a, b)
#define m2d_clear_coef(a, b) *((a) + (b)) = 0
#else
#include "idct.h"
#endif

#ifdef _M_IX86
#include <crtdbg.h>
#define VC_CHECK assert(_CrtCheckMemory());


#ifdef _WINDLL
#include <Windows.h>

BOOL APIENTRY DllMain (HMODULE hModule, DWORD ul_reason_for_call, LPVOID lpReserved)
{
	switch(ul_reason_for_call) {
	case DLL_PROCESS_ATTACH:
	case DLL_THREAD_ATTACH:
	case DLL_THREAD_DETACH:
	case DLL_PROCESS_DETACH:
	break;
	}
	return TRUE;
}

#endif
#else
#define VC_CHECK
#endif

#ifndef NDEBUG
#include <stdio.h>
#endif

#define SATURATE(x, min, max) (((x) <= (max)) ? (((min) <= (x)) ? (x) : (min)) : (max))
#define SAT255(x) (((uint32_t)(x) <= 255) ? (x) : ((x) < 0 ? 0 : 255))
#define SIGN_EXTEND(x, bitlen) (-((x) & (1 << ((bitlen) - 1))) | (x))
#define ABS_BY_BITLEN(x, len) {int sign = (unsigned)(x) >> ((len) - 1); x = (((((x) ^ (-sign & ((1 << (len)) - 1))) + sign) * 2) | sign);}

#ifdef __RENESAS_VERSION__
extern "C"
#endif
void m2d_copy16xn(const uint8_t *src, uint8_t *dst, int stride, int height);

static void m2d_copy_slice(m2d_mb_current *mb, int mb_row_num_plus1);
static void m2d_skip_mb_P(m2d_mb_current *mb, int mb_increment);
static void m2d_skip_mb_B(m2d_mb_current *mb, int mb_increment);
static int m2d_macroblock_type_I(dec_bits *stream);
static int m2d_macroblock_type_P(dec_bits *stream);
static int m2d_macroblock_type_B(dec_bits *stream);
static int m2d_decode_macroblocks(m2d_context *m2d);
static void m2d_motion_vectors(m2d_mb_current *mb, dec_bits *stream, int s, int *mvxy, int8_t *ref_field);

static int (* const macroblock_type_func[])(dec_bits *stream) = {
	m2d_macroblock_type_I,
	m2d_macroblock_type_P,
	m2d_macroblock_type_B
};

static void (* const skip_mb_func[4])(m2d_mb_current *mb, int mb_increment) = {
	m2d_skip_mb_P, /* dummy */
	m2d_skip_mb_P,
	m2d_skip_mb_B,
};



#define VALID_MEM(x) (((intptr_t)(x) != 0) && (((intptr_t)(x) & 15) == 0))
#define VALID_PTR(x) ((x) && !((intptr_t)(x) & (sizeof(*(x)) - 1)))
#define MARKER_BIT1(stream, err) {err |= !get_onebit(stream);}

__LIBM2DEC_API int m2d_set_frames(m2d_context *m2d, int num_mem, m2d_frame_t *mem)
{
	m2d_frames *frames;
	int i;

	if (!m2d || (MAX_FRAME_NUM < (unsigned)num_mem) || !mem) {
		return -1;
	}
	frames = m2d->mb_current->frames;
	frames->num = num_mem;
	memcpy(frames->frames, mem, sizeof(frames->frames[0]) * num_mem);
	for (i = 0; i < num_mem; ++i) {
		if (!VALID_MEM(mem[i].luma) || !VALID_MEM(mem[i].chroma)) {
			return -1;
		}
	}
	frames->index = -1;
	return 0;
}

static int find_valid_frame(int ref0_idx, int ref1_idx, int *lru, int num_frames)
{
	int i;
	int max_idx = -1;
	int max_val = -1;
	for (i = 0; i < num_frames; ++i) {
		if ((i != ref0_idx) && (i != ref1_idx)) {
			int val = lru[i];
			lru[i] = val + 1;
			if (max_val < val) {
				max_val = val;
				max_idx = i;
			}
		}
	}
	if (max_idx < 0) {
		/* no available frame */
		max_idx = ref0_idx;
	}
	lru[max_idx] = 0;
	return max_idx;
}

static void set_ptrdiff(m2d_frames *frames, int ref, int frm_idx, m2d_frame_t *curr_frame)
{
	frames->diff_to_ref[ref][0] = frames->frames[frm_idx].luma - curr_frame->luma;
	frames->diff_to_ref[ref][1] = frames->frames[frm_idx].chroma - curr_frame->chroma;
}

static void m2d_update_frames(m2d_context *m2d, m2d_frames *frames, int next_coding_type, int temporal_reference)
{
	int curr_idx;
	int ref0_idx, ref1_idx;
	m2d_frame_t *curr_frame;

	curr_idx = frames->index;
	if (curr_idx < 0) {
		/* Just after initialization */
		m2d->out_state = ((next_coding_type == I_VOP) || (next_coding_type == P_VOP)) ? 1 * 2 : 0;
		frames->index = 0;
		return;
	}

	ref0_idx = frames->idx_of_ref[0];
	ref1_idx = frames->idx_of_ref[1];
	curr_idx = find_valid_frame(ref0_idx, ref1_idx, frames->lru, frames->num);
	if ((next_coding_type == I_VOP) || (next_coding_type == P_VOP)) {
		ref0_idx = ref1_idx;
		ref1_idx = curr_idx;
		frames->idx_of_ref[0] = ref0_idx;
		frames->idx_of_ref[1] = ref1_idx;
		if (m2d->out_state < (2 * 2)) {
			m2d->out_state += 2;
		}
	} else {
		/* B_VOP */
		m2d->out_state |= 1;
	}

	frames->index = curr_idx;
	frames->frames[curr_idx].cnt = temporal_reference;
	curr_frame = &frames->frames[curr_idx];
	set_ptrdiff(frames, 0, ref0_idx, curr_frame);
	set_ptrdiff(frames, 1, ref1_idx, curr_frame);
}

static void m2d_frames_inc_mb_x_pos(m2d_frames *frames, int inc_x)
{
	int offset_luma = inc_x * MB_LEN;
	frames->curr_luma = frames->curr_luma + offset_luma;
	frames->curr_chroma = frames->curr_chroma + offset_luma;
}

static void m2d_frames_inc_mb_y_pos(m2d_frames *frames, int inc_y, int width)
{
	int offset_luma = inc_y * width * MB_LEN;
	frames->curr_luma = frames->curr_luma + offset_luma;
	frames->curr_chroma = frames->curr_chroma + (offset_luma >> 1);
}

static void m2d_frames_set_mb_pos(m2d_mb_current *mb, int mb_x, int mb_y, int width)
{
	m2d_frames *frames = mb->frames;
	m2d_frame_t *frame = &frames->frames[frames->index < 0 ? 0 : frames->index];
	mb->mb_x = mb_x;
	mb->mb_y = mb_y;
	frames->curr_luma = frame->luma;
	frames->curr_chroma = frame->chroma;
	m2d_frames_inc_mb_y_pos(frames, mb_y, width);
	m2d_frames_inc_mb_x_pos(frames, mb_x);
}

void m2d_mb_set_default(m2d_mb_current *mb);

__LIBM2DEC_API int m2d_init(m2d_context *m2d, int dummy, int (*header_callback)(void *arg, void *seq_id), void *arg)
{
	memset(m2d, 0, sizeof(*m2d));
	m2d->seq_header = &m2d->seq_header_i;
	m2d->picture = &m2d->picture_i;
	m2d->stream = &m2d->stream_i;
	m2d->mb_current = &m2d->mb_current_i;
	m2d_mb_set_default(m2d->mb_current);
	m2d->gop_header = &m2d->gop_header_i;
	m2d->mb_current->frames = &m2d->frames_i;
	m2d->header_callback = header_callback ? header_callback : header_dummyfunc;
	m2d->header_callback_arg = arg;
	VC_CHECK;
	dec_bits_open(m2d->stream, 0);
	return 0;
}

static dec_bits *m2d_stream_pos(m2d_context *m2d)
{
	return &m2d->stream_i;
}

static void m2d_load_qmat(uint8_t *qmat, const int8_t *zigzag, dec_bits *stream)
{
	int i = 64 / 3;
	do {
		int d = get_bits(stream, 24);
		qmat[*zigzag++] = d >> 16;
		qmat[*zigzag++] = d >> 8;
		qmat[*zigzag++] = d;
	} while (--i);
	qmat[*zigzag] = get_bits(stream, 8);
}

void m2d_set_qmat(m2d_mb_current *mb, dec_bits *stream, uint8_t *qmat, int load_bit, int qmat_idx)
{
	if (load_bit == 0) {
		mb->qmat[qmat_idx] = m2d_qmat_default[qmat_idx];
	} else {
		mb->qmat[qmat_idx] = qmat;
		m2d_load_qmat(qmat, m2d_zigzag[0], stream);
	}
}

static int m2d_parse_coef_mpeg2_intra(m2d_mb_current *mb, dec_bits *stream, int idx);
static int m2d_parse_coef_mpeg2_inter(m2d_mb_current *mb, dec_bits *stream, int idx);
static int m2d_parse_coef_mpeg1_intra(m2d_mb_current *mb, dec_bits *stream, int idx);
static int m2d_parse_coef_mpeg1_inter(m2d_mb_current *mb, dec_bits *stream, int idx);

int (* const m2d_parse_coef[2][2])(m2d_mb_current *mb, dec_bits *stream, int idx) = {
	{
		m2d_parse_coef_mpeg1_intra,
		m2d_parse_coef_mpeg1_inter,
	},
	{
		m2d_parse_coef_mpeg2_intra,
		m2d_parse_coef_mpeg2_inter,
	},
};

/** Switch function pointer to MPEG1 or MPEG2 8x8 coefficients parser
 */
static void m2d_mb_set_mpeg2_mode(m2d_mb_current *mb, int is_mpeg2)
{
	assert((unsigned)is_mpeg2 <= 1);
	mb->parse_coef = m2d_parse_coef[is_mpeg2];
}

void m2d_mb_set_default(m2d_mb_current *mb)
{
	mb->parse_coef = m2d_parse_coef[0];
	mb->intra_vlc_format = 0;
	mb->intra_dc_scale = 3;
	mb->intra_dc_max = (1 << 8) - 1;
	mb->frame_mode = 3;
	mb->zigzag = m2d_zigzag[0];
	m2d_mb_set_mpeg2_mode(mb, 0);
}

static void m2d_mb_set_frame_size(m2d_mb_current *mb, int width, int height)
{
	int mb_x = (width + 15) >> 4;
	int mb_y = (height + 15) >> 4;
	mb->mbmax_x = mb_x;
	mb->mbmax_y = mb->dct_type == 2 ? (mb_y >> 1) : mb_y;
	width = mb_x * MB_LEN;
	height = mb_y * MB_LEN;
#ifdef FAST_DECODE
	mb->frame_width = (mb_x * 2 + 15) & ~15;
	mb->frame_height = (mb_y * 2 + 15) & ~15;
#else
	mb->frame_width = mb_x * 16;
	mb->frame_height = mb_y * 16;
#endif
}

static int m2d_read_seq_header(m2d_context *m2d)
{
	dec_bits *stream;
	m2d_seq_header *header;
	int bit;
	int err = 0;

	stream = m2d->stream;
	header = m2d->seq_header;
	header->horizontal_size_value = get_bits(stream, 12);
	header->vertical_size_value = get_bits(stream, 12);
	header->aspect_ratio_information = get_bits(stream, 4);
	header->frame_rate_code = get_bits(stream, 4);
	header->bit_rate_value = get_bits(stream, 18);
	MARKER_BIT1(stream, err);
	header->vbv_buffer_size = get_bits(stream, 10);
	header->constrained_parameters_flag = get_onebit(stream);

	header->load_intra_quantizer_matrix = bit = get_onebit(stream);
	m2d_set_qmat(m2d->mb_current, stream, m2d->qmat[0], bit, 0);
	header->load_non_intra_quantizer_matrix = bit = get_onebit(stream);
	m2d_set_qmat(m2d->mb_current, stream, m2d->qmat[1], bit, 1);
	m2d_mb_set_frame_size(m2d->mb_current, header->horizontal_size_value, header->vertical_size_value);
	m2d->header_callback(m2d->header_callback_arg, stream->id);
	return err;
}

static int m2d_read_extension_dummy(m2d_context *m2d)
{
	dec_bits *stream = m2d->stream;
	int c;
	byte_align(stream);
	while ((c = show_bits(stream, 24)) != 0x000001) {
		skip_bits(stream, ((c & 0xff) != 0) ? 24 : 8);
	}
	return 0;
}

static int m2d_read_sequence_extension(m2d_context *m2d)
{
	dec_bits *stream;
	m2d_seq_header *header;
	int err = 0;

	stream = m2d->stream;
	header = m2d->seq_header;
	header->profile_and_level = get_bits(stream, 8);
	header->progressive_sequence = get_bits(stream, 1);
	header->chroma_format = get_bits(stream, 2);
	header->horizontal_size_value |= get_bits(stream, 2) << 12;
	header->vertical_size_value |= get_bits(stream, 2) << 12;
	header->bit_rate_value |= get_bits(stream, 12) << 18;
	MARKER_BIT1(stream, err);
	header->vbv_buffer_size |= get_bits(stream, 8) << 10;

	m2d_mb_set_frame_size(m2d->mb_current, header->horizontal_size_value, header->vertical_size_value);
	m2d_mb_set_mpeg2_mode(m2d->mb_current, 1);
	m2d->header_callback(m2d->header_callback_arg, stream->id);
	return err;
}

static int m2d_read_qmatrix_extension(m2d_context *m2d)
{
	dec_bits *stream;
	const uint8_t **mb_qmat;
	const int8_t *zigzag;
	int i;

	stream = m2d->stream;
	mb_qmat = m2d->mb_current->qmat;
	zigzag = m2d->mb_current->zigzag;
	for (i = 0; i < 4; ++i) {
		if (get_onebit(stream) != 0) {
			uint8_t *qmat = m2d->qmat[i];
			mb_qmat[i] = qmat;
			m2d_load_qmat(qmat, zigzag, stream);
		}
	}
	return 0;
}

/**Split bit pattern into bit-field members.
 */
#define BITSET(dst, bits, bitnum, shift_sum) {shift_sum -= (bitnum); dst = (((bits) >> (shift_sum)) & ((1 << (bitnum)) - 1)); }

static void set_coding_extension_param1(m2d_picture *pic, int bits)
{
	int shift_sum = 2 + 2 + 1 * 10;

	BITSET(pic->intra_dc_precision, bits, 2, shift_sum);
	BITSET(pic->picture_structure, bits, 2, shift_sum);
	BITSET(pic->top_field_first, bits, 1, shift_sum);
	BITSET(pic->frame_pred_frame_dct, bits, 1, shift_sum);
	BITSET(pic->concealment_motion_vectors, bits, 1, shift_sum);
	BITSET(pic->q_scale_type, bits, 1, shift_sum);
	BITSET(pic->intra_vlc_format, bits, 1, shift_sum);
	BITSET(pic->alternate_scan, bits, 1, shift_sum);
	BITSET(pic->repeat_first_field, bits, 1, shift_sum);
	BITSET(pic->chroma_420_type, bits, 1, shift_sum);
	BITSET(pic->progressive_frame, bits, 1, shift_sum);
	BITSET(pic->composite_display_flag, bits, 1, shift_sum);
	assert(shift_sum == 0);
}


static void set_coding_extension_param2(m2d_picture *pic, int bits)
{
	int shift_sum = 1 + 3 + 1 + 7 + 8;

	BITSET(pic->v_axis, bits, 1, shift_sum);
	BITSET(pic->field_sequence, bits, 3, shift_sum);
	BITSET(pic->sub_carrier, bits, 1, shift_sum);
	BITSET(pic->burst_amplitude, bits, 7, shift_sum);
	BITSET(pic->sub_carrier_phase, bits, 8, shift_sum);
	assert(shift_sum == 0);
}

static void set_coding_type(m2d_mb_current *mb, int coding_type)
{
	coding_type -= 1;
	mb->macroblock_type = macroblock_type_func[coding_type];
	mb->skip_mb = skip_mb_func[coding_type];
}

static unsigned int guess_picture_coding_type(unsigned int f_codes)
{
	if ((f_codes & 0xff) == 0xff) {
		if ((f_codes & 0xff00) == 0xff00) {
			return I_VOP;
		} else {
			return P_VOP;
		}
	} else {
		return B_VOP;
	}
}

static int m2d_read_coding_extension(m2d_context *m2d)
{
	dec_bits *stream;
	m2d_picture *pic;
	m2d_mb_current *mb;
	int bits;
	int err;
	unsigned int f_codes;

	stream = m2d->stream;
	pic = m2d->picture;
	mb = m2d->mb_current;
	err = 0;

	f_codes = get_bits(stream, 4 * 4);

	mb->r_size[0][0] = (f_codes >> 12) - 1;
	mb->r_size[0][1] = ((f_codes >> 8) & 15) - 1;
	mb->r_size[1][0] = ((f_codes >> 4) & 15) - 1;
	mb->r_size[1][1] = (f_codes & 15) - 1;
	if (!pic->picture_coding_type) {
		/* picture_header() was not present */
		unsigned int coding_type;
		pic->picture_coding_type = coding_type = guess_picture_coding_type(f_codes);
		set_coding_type(mb, coding_type);
	}
	bits = get_bits(stream, 2 + 2 + 1 * 10);
	set_coding_extension_param1(pic, bits);
	mb->intra_vlc_format = (pic->concealment_motion_vectors * 2) | pic->intra_vlc_format;
	mb->intra_dc_scale = 3 - pic->intra_dc_precision;
	mb->intra_dc_max = (1 << (pic->intra_dc_precision + 8)) - 1;
	mb->zigzag = m2d_zigzag[pic->alternate_scan];
	switch (pic->picture_structure) {
	case 1:
		/* FALLTHROUGH */
	case 2:
		mb->frame_mode = 0;
		break;
	case 3:
		mb->frame_mode = pic->frame_pred_frame_dct ? 3 : 1;
		break;
	}
	if (pic->composite_display_flag) {
		bits = get_bits(stream, 1 + 3 + 1 + 7 + 8);
		set_coding_extension_param2(pic, bits);
	}
	return err;
}

static int m2d_read_display_extension(m2d_context *m2d)
{
	dec_bits *stream;
	int format;
	int width_height;

	stream = m2d->stream;
	format = get_bits(stream, 4);
	if (format & 1) {
		skip_bits(stream, 8 * 3);
	}
	width_height = get_bits(stream, 14 + 1 + 14);
	m2d->seq_header->display_height = width_height & 0x3fff;
	m2d->seq_header->display_width = (unsigned)width_height >> 15;
	return 0;
}

static int m2d_read_copyright_extension(m2d_context *m2d)
{
	dec_bits *stream;

	stream = m2d->stream;
	skip_bits(stream, 10 + 8);
	skip_bits(stream, 20 + 1);
	skip_bits(stream, 22 + 1);
	skip_bits(stream, 22);
	return 0;
}

static int (* const m2d_read_extension_func[])(m2d_context *m2d) = {
	m2d_read_extension_dummy,
	m2d_read_sequence_extension,
	m2d_read_display_extension,
	m2d_read_qmatrix_extension,
	m2d_read_copyright_extension,
	m2d_read_extension_dummy,
	m2d_read_extension_dummy,
	m2d_read_extension_dummy,

	m2d_read_coding_extension,
	m2d_read_extension_dummy,
	m2d_read_extension_dummy,
	m2d_read_extension_dummy,
	m2d_read_extension_dummy,
	m2d_read_extension_dummy,
	m2d_read_extension_dummy,
	m2d_read_extension_dummy,
};
	
/**Dispatch apropriate parser for miscellaneous extentions.
 */
static int m2d_read_extension(m2d_context *m2d)
{
	int id = get_bits(m2d->stream, 4);
	return m2d_read_extension_func[id](m2d);
}

int m2d_user_data(m2d_context *m2d)
{
	dec_bits *stream = m2d->stream;
	while (1 < show_bits(stream, 24)) {
		skip_bits(stream, 8);
	}
	return 0;
}

static int m2d_read_gop_header(m2d_context *m2d)
{
	dec_bits *stream;
	m2d_gop_header *header;
	int err;

	stream = m2d->stream;
	header = m2d->gop_header;
	err = 0;

	header->time_code = get_bits(stream, 25);
	header->closed_gop = get_bits(stream, 1);
	header->broken_link = get_bits(stream, 1);
	return err;
}

static void m2d_init_mb_pos(m2d_mb_current *mb);

static int m2d_read_picture_header(m2d_context *m2d)
{
	dec_bits *stream;
	m2d_picture *pic;
	m2d_mb_current *mb;
	int coding_type;
	int err;

	stream = m2d->stream;
	pic = m2d->picture;
	err = 0;

	pic->temporal_reference = get_bits(stream, 10);
	pic->picture_coding_type = coding_type = get_bits(stream, 3);
	pic->vbv_delay = get_bits(stream, 16);
	mb = m2d->mb_current;
	set_coding_type(mb, coding_type);
	m2d_init_mb_pos(m2d->mb_current);
	if ((coding_type == P_VOP) || (coding_type == B_VOP)) {
		int r_size;
		r_size = get_bits(stream, 4) - 1;
		mb->r_size[0][0] = r_size;
		mb->r_size[0][1] = r_size;
		if (coding_type == B_VOP) {
			r_size = get_bits(stream, 4) - 1;
			mb->r_size[1][0] = r_size;
			mb->r_size[1][1] = r_size;
		}
	}
	while (get_bits(stream, 1)) {
		get_bits(stream, 9);
	}
	return err;
}

static int m2d_read_slice(m2d_context *m2d, int code_type)
{
	dec_bits *stream;
	m2d_mb_current *mb;
	m2d_picture *pic;
	int vertical_pos;
	int err;

	stream = m2d->stream;
	err = 0;

	pic = m2d->picture;
	mb = m2d->mb_current;
	mb->q_mapping = m2d_q_mapping[pic->q_scale_type];
	mb->q_scale = mb->q_mapping[get_bits(stream, 5)];
	vertical_pos = (code_type & 255) - 1;
	if (vertical_pos == 0) {
		m2d_update_frames(m2d, mb->frames, pic->picture_coding_type, pic->temporal_reference);
	}
	if (mb->mbmax_y <= vertical_pos) {
		return 0;
	}
	if (1 < vertical_pos - mb->mb_y) {
		m2d_copy_slice(mb, vertical_pos - mb->mb_y);
	}
	m2d_frames_set_mb_pos(mb, -1, vertical_pos, mb->frame_width);
	if (get_onebit(stream) != 0) {
		get_bits(stream, 1 * 2 + 6);
		while (get_onebit(stream) != 0) {
			get_bits(stream, 8);
		}
	}
	err = m2d_decode_macroblocks(m2d);
	return err;
/*	err |= (val == 0); */
}

static int m2d_dispatch_one_nal(m2d_context *m2d, int code_type)
{
	int err;

	if (setjmp(m2d->stream->jmp) != 0) {
		return 0;
	}
	if (code_type < 0xb0) {
		if (code_type == 0) {
			err = m2d_read_picture_header(m2d);
		} else {
			err = m2d_read_slice(m2d, code_type);
		}
	} else {
		switch (code_type) {
		case 0xb2:
			err = m2d_user_data(m2d);
			break;
		case 0xb3:
			err = m2d_read_seq_header(m2d);
			break;
		case 0xb5:
			err = m2d_read_extension(m2d);
			break;
		case 0xb8:
			err = m2d_read_gop_header(m2d);
			break;
		default:
			err = M2D_ERR_NOT_VIDEO;
			break;
		}
	}
	return err;
}

/** Decode macroblock_type for Intra picture.
 * \return bit-fielded information as Table B.2 in the standard.
 */
static int m2d_macroblock_type_I(dec_bits *stream)
{
	int type = show_bits(stream, 2);
	if (type & 2) {
		/* "1" */
		skip_bits(stream, 1);
		return MB_INTRA;
	} else {
		/* "01" */
		assert(type != 0); /* forbidden */
		skip_bits(stream, 2);
		return MB_QUANT | MB_INTRA;
	}
}

static void m2d_copy_slice(m2d_mb_current *mb, int mb_row_num_plus1)
{
	int width, luma_len, mb_row_num;
	m2d_frames *frm;
	uint8_t *dst;
	ptrdiff_t *diff;

	assert(1 < mb_row_num_plus1);
	width = mb->frame_width;
	m2d_frames_set_mb_pos(mb, 0, mb->mb_y + 1, width);
	frm = mb->frames;
	dst = frm->curr_luma;
	diff = frm->diff_to_ref[0];
	mb_row_num = mb_row_num_plus1 - 1;
	luma_len = width * mb_row_num * MB_LEN;
	memcpy(dst, dst + diff[0], luma_len);
	dst = frm->curr_chroma;
	memcpy(dst, dst + diff[1], luma_len >> 1);
}


static void m2d_inc_mb_pos(m2d_mb_current *mb);
static void m2d_mb_reset(m2d_mb_current *mb);
static void dummyfunc(m2d_mb_current *mb) {}

static void m2d_skip_mb_P(m2d_mb_current *mb, int mb_increment)
{
	int skip_num = mb_increment - 1;
	m2d_frames *frm = mb->frames;
	int width = mb->frame_width;
	ptrdiff_t *diff = frm->diff_to_ref[0];
	void (*inc_mb)(m2d_mb_current *mb);

	if (skip_num <= 0) {
		inc_mb = dummyfunc;
		skip_num = 1;
	} else {
		inc_mb = m2d_inc_mb_pos;
	}
	do {
		uint8_t *luma, *chroma;
		inc_mb(mb);
		luma = frm->curr_luma;
		m2d_copy16xn(luma + diff[0], luma, width, MB_LEN);
		chroma = frm->curr_chroma;
		m2d_copy16xn(chroma + diff[1], chroma, width, MB_LEN / 2);
	} while (--skip_num);
	m2d_mb_reset(mb);
}

/** Decode macroblock_type for P picture.
 * \return bit-fielded information as Table B.3 in the standard.
 */
static int m2d_macroblock_type_P(dec_bits *stream)
{
	return m2d_dec_vld_unary(stream, mb_type_p_bit3, 3);
}

static void m2d_skip_mb_B(m2d_mb_current *mb, int mb_increment)
{
	int skip_num = mb_increment - 1;
	m2d_frames *frm = mb->frames;
	int width = mb->frame_width;
	int mb_type = mb->type;
	int dir = MB_PARAM(mb_type, MB_MC);
	ptrdiff_t *diff = frm->diff_to_ref[0];
	void (* const * motion_comp)(const uint8_t *src, uint8_t *dst, int stride, int *mvxy, int height) = m2d_motion_compensation[0];
	int is_bidirectional = (dir == MB_MC);
	int mvxy[4];
	int16_t *mv;
	dir = ~-is_bidirectional & (dir >> 1);
	mv = mb->mv[dir].mv[0];
	mvxy[0] = mv[0];
	mvxy[1] = mv[1];
	if (is_bidirectional) {
		mv = mb->mv[1].mv[0];
		mvxy[2] = mv[0];
		mvxy[3] = mv[1];
	} else {
		diff += dir * 2;
	}
	do {
		uint8_t *luma, *chroma;
		m2d_inc_mb_pos(mb);
		luma = frm->curr_luma;
		motion_comp[0](luma + diff[0], luma, width, mvxy, MB_LEN);
		chroma = frm->curr_chroma;
		motion_comp[1](chroma + diff[1], chroma, width, mvxy, MB_LEN / 2);
		if (is_bidirectional) {
			motion_comp[2](luma + diff[2], luma, width, mvxy + 2, MB_LEN);
			motion_comp[3](chroma + diff[3], chroma, width, mvxy + 2, MB_LEN / 2);
		}
	} while (--skip_num);
}


/** Decode macroblock_type for B picture.
 * \return bit-fielded information as Table B.4 in the standard.
 */
static int m2d_macroblock_type_B(dec_bits *stream)
{
	return m2d_dec_vld_unary(stream, mb_type_b_bit4, 4);
}

static const m2d_motion_type_t m2d_motion_type[2][4] = {
	{
		{2, MV_FIELD, 0, MV_FIELD}, /* dummy */
		{2, MV_FIELD, 0, MV_FIELD},
		{1, MV_FRAME, 0, MV_FRAME},
		{1, MV_FIELD, 1, MV_DUALPRIME},
	},
	{
		{1, MV_FIELD, 0, MV_FIELD}, /* dummy */
		{1, MV_FIELD, 0, MV_FIELD},
		{2, MV_FIELD, 0, MV_16_8MC},
		{1, MV_FIELD, 1, MV_DUALPRIME},
	},
};

static int m2d_decode_macroblock_mode(m2d_mb_current *mb, dec_bits *stream)
{
	int type;
	int frame_mode;
	int idx;

	type = mb->macroblock_type(stream);
	mb->type = type;
	if (MB_PARAM(type, MB_SPACIAL) && 1) {
		assert(0);
	}
	frame_mode = mb->frame_mode;
	if (MB_PARAM(type, MB_MC)) {
		if (frame_mode & 1) {
			if (frame_mode == 1) {
				idx = get_bits(stream, 2);
			} else {
				idx = 2;
			}
			mb->motion_type = &m2d_motion_type[0][idx];
		} else {
			idx = get_bits(stream, 2);
			mb->motion_type = &m2d_motion_type[1][idx];
		}
	} else {
		idx = (frame_mode == 0);
		mb->motion_type = &m2d_motion_type[idx][2 - idx];
	}
	if ((frame_mode == 1) && MB_PARAM(type, MB_COEF)) {
		mb->dct_type = get_onebit(stream);
	} else if (frame_mode != 0) {
		mb->dct_type = 0;
	} else {
		mb->dct_type = 1;
	}
	return type;
}

/**Reset intra DC predictor.
 */
static void m2d_mb_reset_intra(m2d_mb_current *mb)
{
	int dc = (mb->intra_dc_max + 1) >> 1;
	int16_t *dc_pred = mb->dc_pred;
	dc_pred[0] = dc;
	dc_pred[1] = dc;
	dc_pred[2] = dc;
}

/**Reset motion vector predictor.
 */
static void m2d_mb_reset_inter(m2d_mb_current *mb)
{
	memset(mb->mv[0].mv[0], 0, sizeof(mb->mv[0].mv) * 2);
}

/**Reset both intra and inter predictor.
 */
static void m2d_mb_reset(m2d_mb_current *mb)
{
	m2d_mb_reset_intra(mb);
	m2d_mb_reset_inter(mb);
}

typedef struct {
	int bitlen;
	const vlc_t *vld_tab;
} m2d_dc_size_t;

static const m2d_dc_size_t dc_size_table[] = {
	{5, dct_dc_size_luma_bit5},
/*	{5, dct_dc_size_luma_bit5},
	{5, dct_dc_size_luma_bit5},
	{5, dct_dc_size_luma_bit5},*/
	{4, dct_dc_size_chroma_bit4},
	{4, dct_dc_size_chroma_bit4},
	{4, dct_dc_size_chroma_bit4},
	{4, dct_dc_size_chroma_bit4},
};

/** Decode DC coefficient in intra macroblock.
 *\param mb Macroblock layer context.
 *\param stream Input bitstream class
 *\luma_chroma Zero specifies luma DC vlc table. One means chroma DC table.
 *\return Decoded DC coefficient with inverse-quantization.
 */
int m2d_parse_intra_dc(m2d_mb_current *mb, dec_bits *stream, int block_idx)
{
	const m2d_dc_size_t *dc_tbl = &dc_size_table[block_idx];
	int dc_bit_num = m2d_dec_vld_unary(stream, dc_tbl->vld_tab, dc_tbl->bitlen);
	int dc_diff = (dc_bit_num == 0) ? 0 : get_bits(stream, dc_bit_num);
	int16_t *dc_pred = &mb->dc_pred[block_idx];
	int dc = *dc_pred;
	if (dc_bit_num != 0) {
		int half_range = 1 << (dc_bit_num - 1);
		int dc_max;
		if (!(dc_diff & half_range)) {
			dc_diff = dc_diff + 1 - half_range * 2;
		}
		dc += dc_diff;
		dc_max = mb->intra_dc_max;
		*dc_pred = dc;
		dc = SATURATE(dc, 0, dc_max);
	}
	return dc << mb->intra_dc_scale;
}

#define NEG_IF_SIGN(lvl, x) (((lvl) & 1) ? -(x) : (x))

/**A functor for INTRA inverse quantization.
 */
class IntraScaling {
public:
	int operator()(int level, int q) {
		int t = ((level >> 1) * q) >> 4;
		return NEG_IF_SIGN(level, t);
	}
};

/**A functor for INTER inverse quantization.
 */
class InterScaling {
public:
	int operator()(int level, int q) {
		int t = ((level | 1) * q) >> 5;
		return NEG_IF_SIGN(level, t);
	}
};

/**A functor to decode "level" of an MPEG-2 Escape secuence.
 */
class EscapeLevelMpeg2 {
public:
	int operator()(dec_bits *stream) {
		int level = get_bits(stream, 12);
		ABS_BY_BITLEN(level, 12);
		return level;
	}
};

class MismatchMpeg2 {
public:
	int operator()(int mismatch_bit, int idx_sum, int16_t *coef) {
		if (!(mismatch_bit & 1)) {
			coef[63] ^= 1;
			idx_sum |= 1 << 7;
		}
		return idx_sum;
	}
};

/**A functor to decode "level" of an MPEG-1 Escape secuence.
 */
class EscapeLevelMpeg1 {
public:
	int operator()(dec_bits *stream) {
		int level = get_bits(stream, 8);
		if ((level & 0x7f) == 0) {
			level = get_bits(stream, 8) - (level & 0x80) * 2;
		} else {
			level = (int8_t)level;
		}
		return (level < 0) ? ((-level * 2) | 1) : (level * 2);
	}
};

class MismatchMpeg1 {
public:
	int operator()(int mismatch_bit, int idx_sum, int16_t *coef) {
		int i = 64;
		do {
			int c = *coef;
			if (c && !(c & 1)) {
				*coef = (0 < c) ? c - 1 : c + 1;
			}
			coef++;
		} while (--i);
		return 0xff;
	}
};

/**Read DCT coefficients using Table B.14, 15 (Table one.)
 *\param mb Information about current macroblock.
 *\param stream Input bitstream class.
 *\param idx Indicates which 8x8 block is to be decoded.
 */
template <int IsInter, typename IQ, typename Esc, typename Mismatch>
static inline int parse_coef(m2d_mb_current *mb, dec_bits *stream, int dct_type, int idx, IQ _iquant, Esc _esc_level, Mismatch _mismatch)
{
	const vlc_dct_t *vld_tab = m2d_dct_tables[dct_type];
	int16_t *coef = mb->coef;
	const uint8_t *qmat = mb->qmat[IsInter];
	int q_scale = mb->q_scale;
	int mismatch_bit = idx ? coef[0] : 0;
	int idx_sum = 0;
	const int8_t *zigzag = mb->zigzag;

	m2d_clear_coef(coef, idx);
	do {
		int rest_len = VLD_BITLEN;
		int bit = show_bits(stream, rest_len);
		const vlc_dct_t *vld_tab_curr = vld_tab;
		const vlc_dct_t *vlc = &vld_tab_curr[bit];
		int len = vlc->length;
		int level, run, idx_zigzag;

		while (len <= 0) {
			/* additional look-up */
			if (len < 0) {
				/* undefined code */
				longjmp(stream->jmp, 1);
				/* NOTREACHED */
				return 0;
			}
			vld_tab_curr += vlc->run;
			skip_bits(stream, rest_len);
			rest_len = vlc->level < VLD_BITLEN ? vlc->level : VLD_BITLEN;
			bit = show_bits(stream, rest_len);
			vlc = &vld_tab_curr[bit];
			len = vlc->length;
		}
		skip_bits(stream, len);
		level = vlc->level;
		run = vlc->run;
		if (0 <= run) {
			idx += run;
		} else if (level != 0) {
			/* end of block */
			break;
		} else {
			/* ESC */
			idx += get_bits(stream, 6); /* get run */
			level = _esc_level(stream);
		}
		if (64 <= idx) {
			break;
		}
#ifdef FAST_DECODE
		if (idx == 0) {
#endif
		idx_zigzag = zigzag[idx];
		if (idx_zigzag & 7) {
			idx_sum |= 1 << ((unsigned)idx_zigzag >> 3);
		}
		level = _iquant(level, qmat[idx_zigzag] * q_scale);
  		level = SATURATE(level, -2048, 2047);
		mismatch_bit += level;
		coef[idx_zigzag] = level;
#ifdef FAST_DECODE
		}
#endif
	} while (++idx);
#ifndef FAST_DECODE
#if 1
	idx_sum = _mismatch(mismatch_bit, idx_sum, coef);
#else
	if (!(mismatch_bit & 1)) {
		coef[63] ^= 1;
		idx_sum |= 1 << 7;
	}
#endif
#endif
	return idx_sum;
}


static int m2d_parse_coef_mpeg2_intra(m2d_mb_current *mb, dec_bits *stream, int idx)
{
	return parse_coef<0>(mb, stream, mb->intra_vlc_format, idx, IntraScaling(), EscapeLevelMpeg2(), MismatchMpeg2());
}

static int m2d_parse_coef_mpeg2_inter(m2d_mb_current *mb, dec_bits *stream, int idx)
{
	return parse_coef<1>(mb, stream, 0, idx, InterScaling(), EscapeLevelMpeg2(), MismatchMpeg2());
}

static int m2d_parse_coef_mpeg1_intra(m2d_mb_current *mb, dec_bits *stream, int idx)
{
	return parse_coef<0>(mb, stream, mb->intra_vlc_format, idx, IntraScaling(), EscapeLevelMpeg1(), MismatchMpeg1());
}

static int m2d_parse_coef_mpeg1_inter(m2d_mb_current *mb, dec_bits *stream, int idx)
{
	return parse_coef<1>(mb, stream, 0, idx, InterScaling(), EscapeLevelMpeg1(), MismatchMpeg1());
}

#define LUMA_BLOCK_OFFSET(idx, wid, dct_type) ((dct_type == 0) ? ((((idx) & 1) + (-((idx) & 2) & (wid))) * (MB_LEN / 2)) : (((idx) & 1) * (MB_LEN / 2) + (-((idx) & 2) & (wid))))

#ifdef DUMP_COEF
int coefs[64];
void print_coefs()
{
	for (int y = 0; y < 8; ++y) {
		printf("\t");
		for (int x = 0; x < 8 - 1; ++x) {
			printf("%d, ", coefs[y * 8 + x]);
		}
		printf("%d,\n", coefs[y * 8 + 7]);
	}
}
#endif

static void m2d_parse_intra_block_luma(m2d_mb_current *mb, dec_bits *stream, int block_idx)
{
	int sum_coef;
	mb->coef[0] = m2d_parse_intra_dc(mb, stream, 0);
	sum_coef = mb->parse_coef[0](mb, stream, 1);
#ifdef DUMP_COEF
	coefs[sum_coef]++;
#endif
	m2d_idct_intra_luma(mb->frames->curr_luma + LUMA_BLOCK_OFFSET(block_idx, mb->frame_width, mb->dct_type), mb->frame_width << mb->dct_type, mb->coef, sum_coef);
}

static void m2d_parse_intra_block_chroma(m2d_mb_current *mb, dec_bits *stream, int block_idx)
{
	int sum_coef;
	mb->coef[0] = m2d_parse_intra_dc(mb, stream, block_idx + 1);
	sum_coef = mb->parse_coef[0](mb, stream, 1);
	m2d_idct_intra_chroma(mb->frames->curr_chroma + block_idx, mb->frame_width, mb->coef, sum_coef);
}

enum {
	CONCEALMENT_MV = 2
};

/**Read macroblock layer data and store.
 * This procedure doesn't calculate iDCT or motion compensation.
 */
static int m2d_parse_intra_macroblock(m2d_mb_current *mb, dec_bits *stream)
{
	int err = 0;
	int i;

	if (MB_PARAM(mb->type, MB_QUANT)) {
		mb->q_scale = mb->q_mapping[get_bits(stream, 5)];
	}
	if (mb->intra_vlc_format & CONCEALMENT_MV) {
		int mvxy[4];
		int8_t ref_field[2];
		m2d_motion_vectors(mb, stream, 0, mvxy, ref_field);
		MARKER_BIT1(stream, err);
	}
	for (i = 0; i < 4; ++i) {
		m2d_parse_intra_block_luma(mb, stream, i);
	}
	for (i = 0; i < 2; ++i) {
		m2d_parse_intra_block_chroma(mb, stream, i);
	}
	VC_CHECK;
	return err;
}

/**
 *\return Half pel of motion vector;
 */
static int m2d_one_mv(int16_t *mv_p, dec_bits *stream, int r_size, int is_field)
{
	int mv;
 	if (get_onebit(stream) == 0) {
		int limit;
		int residual;
		mv = m2d_dec_vld_unary(stream, motion_code_bit5, 5);
		residual = (0 < r_size) ? 1 + get_bits(stream, r_size) : 1;
		if (0 <= mv) {
			mv = ((mv - 1) << r_size) + residual;
		} else {
			mv = ((mv + 1) << r_size) - residual;
		}
		mv += (*mv_p >> is_field);
		limit = 16 << r_size;
		mv = (-limit <= mv) ? ((mv < limit) ? mv : mv - limit * 2) : mv + limit * 2;
	} else {
		mv = *mv_p >> is_field;
	}
	*mv_p = mv << is_field;
	return mv;
}

static int m2d_one_mv_with_dmv(int16_t *mv_p, dec_bits *stream, int r_size, int is_field)
{
	int mv = m2d_one_mv(mv_p, stream, r_size, is_field);
	int dmvector = get_onebit(stream);
	if (dmvector != 0) {
		dmvector = get_onebit(stream) ? -1 : 1;
	}
	return mv;
}


/**Decode a pair of motion vector(horizontal and vertical).
 *\param mv Destination buffer short[2] in which decoded vector is to be stored.
 *\param stream Input bitstream class.
 *\param r_size Scaling factor of vectors.
 *\return Half pel of vectors.
 * - bit0: Horizontal half pel
 * - bit1: Vertical half pel
 */
static void m2d_mv(int16_t *mv, dec_bits *stream, int8_t *r_size, int is_field, int *dst)
{
	assert((unsigned)is_field <= 1);
	dst[0] = m2d_one_mv(mv, stream, r_size[0], 0);
	dst[1] = m2d_one_mv(mv + 1, stream, r_size[1], is_field);
}

static void m2d_mv_with_dmv(int16_t *mv, dec_bits *stream, int8_t *r_size, int is_field, int *dst)
{
	assert((unsigned)is_field <= 1);
	dst[0] = m2d_one_mv_with_dmv(mv, stream, r_size[0], 0);
	dst[1] = m2d_one_mv_with_dmv(mv + 1, stream, r_size[1], is_field);
}

static void m2d_motion_vectors(m2d_mb_current *mb, dec_bits *stream, int s, int *mvxy, int8_t *ref_field)
{
	const m2d_motion_type_t *mv_type = mb->motion_type;
	m2d_mv_t *mv_s;
	int16_t *mv;
	int vertical_field_select;

	mv_s = &mb->mv[s];
	mv = mv_s->mv[0];
	int8_t *r_size = mb->r_size[s];
	if (mv_type->mv_count == 1) {
		if ((mv_type->format == MV_FIELD) && (mv_type->dmv == 0)) {
			/* motion_vertical_field_select */
			vertical_field_select = get_onebit(stream);
		}
		if (mv_type->dmv) {
			m2d_mv_with_dmv(mv, stream, r_size, (mv_type->format == MV_FIELD), mvxy);
		} else {
			m2d_mv(mv, stream, r_size, (mv_type->format == MV_FIELD), mvxy);
		}
		mv[2] = mv[0];
		mv[3] = mv[1];
	} else {
		for (int i = 0; i < 2; ++i) {
			*ref_field++ = get_onebit(stream);
			m2d_mv(mv, stream, r_size, 1, mvxy);
			mv += 2;
			mvxy+= 2;
		}
	}
}

static void m2d_motion_comp(m2d_mb_current *mb, int s, int is_bidirectional, int *mvxy, int8_t *ref_field)
{
	void (* const * motion_comp)(const uint8_t *, uint8_t *, int, int *, int) = m2d_motion_compensation[is_bidirectional];
	m2d_frames *frames = mb->frames;
	ptrdiff_t *dif = frames->diff_to_ref[s];
	int width = mb->frame_width;
	uint8_t *luma = frames->curr_luma;
	uint8_t *chroma = frames->curr_chroma;
	uint8_t *src_luma = luma + dif[0];
	uint8_t *src_chroma = chroma + dif[1];

	if (mb->motion_type->mv_count == 1) {
		motion_comp[0](src_luma, luma, width, mvxy, MB_LEN);
		motion_comp[1](src_chroma, chroma, width, mvxy, MB_LEN / 2);
	} else {
		for (int i = 0; i < 2; ++i) {
			int src_offset = *ref_field++ ? width : 0;
			motion_comp[0](src_luma + src_offset, luma, width * 2, mvxy, MB_LEN / 2);
#ifdef FAST_DECODE
			if (i == 0) {
				src_offset = 0;
#endif
				motion_comp[1](src_chroma + src_offset, chroma, width * 2, mvxy, MB_LEN / 4);
#ifdef FAST_DECODE
			}
#endif
			mvxy += 2;
			luma += width;
			chroma += width;
		}
	}
}

static inline int m2d_coded_block_pattern(m2d_mb_current *mb, dec_bits *stream)
{
	return m2d_dec_vld_unary(stream, coded_block_pattern_bit5, 5);
}

/** Read inter DC coefficient if VLC were "1s."
 */
static int m2d_parse_inter_dc(dec_bits *stream) 
{
	int bits = show_bits(stream, 2);
	if (bits & 2) {
		skip_bits(stream, 2);
		return bits;
	} else {
		return 0;
	}
}

static uint32_t m2d_parse_inter_block(m2d_mb_current *mb, dec_bits *stream, int block_idx)
{
	int sum_coef;
	int dc = m2d_parse_inter_dc(stream);
	if (dc != 0) {
		mb->coef[0] = InterScaling()(dc, mb->q_scale * (mb->qmat[1][0]));
		dc = 1;
	}
	sum_coef = mb->parse_coef[1](mb, stream, dc);
#ifdef DUMP_COEF
	coefs[sum_coef]++;
#endif
	return sum_coef;
}

static void m2d_parse_inter_block_luma(m2d_mb_current *mb, dec_bits *stream, int block_idx)
{
	int sum_coef = m2d_parse_inter_block(mb, stream, block_idx);
	m2d_idct_inter_luma(mb->frames->curr_luma + LUMA_BLOCK_OFFSET(block_idx, mb->frame_width, mb->dct_type), mb->frame_width << mb->dct_type, mb->coef, sum_coef);
}

static void m2d_parse_inter_block_chroma(m2d_mb_current *mb, dec_bits *stream, int block_idx)
{
	int sum_coef = m2d_parse_inter_block(mb, stream, block_idx);
	m2d_idct_inter_chroma(mb->frames->curr_chroma + block_idx, mb->frame_width, mb->coef, sum_coef);
}

/**Read macroblock layer data and store.
 * This procedure doesn't calculate iDCT or motion compensation.
 */
static int m2d_parse_inter_macroblock(m2d_mb_current *mb, dec_bits *stream)
{
	int mb_type = mb->type;
	int err = 0;

	if (MB_PARAM(mb_type, MB_QUANT)) {
		mb->q_scale = mb->q_mapping[get_bits(stream, 5)];
	}
	if (MB_PARAM(mb_type, MB_MC)) {
		int mvxy[4];
		int8_t ref_field[2];
		int forward = MB_PARAM(mb_type, MB_FORWARD);
		if (forward) {
			m2d_motion_vectors(mb, stream, 0, mvxy, ref_field);
			m2d_motion_comp(mb, 0, 0, mvxy, ref_field);
		}
		if (MB_PARAM(mb_type, MB_BACKWARD)) {
			m2d_motion_vectors(mb, stream, 1, mvxy, ref_field);
			m2d_motion_comp(mb, 1, forward != 0, mvxy, ref_field);
		}
	} else {
		m2d_skip_mb_P(mb, 0);
	}
	if (MB_PARAM(mb_type, MB_PATTERN)) {
		int cbp = m2d_coded_block_pattern(mb, stream);
		int i;
		for (i = 0; i < 4; ++i) {
			if (cbp & (1 << (5 - i))) {
				m2d_parse_inter_block_luma(mb, stream, i);
			}
		}
		for (i = 0; i < 2; ++i) {
			if (cbp & (1 << (1 - i))) {
				m2d_parse_inter_block_chroma(mb, stream, i);
			}
		}
	}
	return err;
}

/**Read macroblock layer data and store.
 * This procedure doesn't calculate iDCT or motion compensation.
 */
static int m2d_parse_macroblock(m2d_mb_current *mb, dec_bits *stream)
{
	int prev_is_intra = MB_PARAM(mb->type, MB_INTRA);
	int mb_type = m2d_decode_macroblock_mode(mb, stream);

	if (MB_PARAM(mb_type, MB_INTRA)) {
		if (!prev_is_intra) {
			m2d_mb_reset_intra(mb);
		}
		return m2d_parse_intra_macroblock(mb, stream);
	} else {
		if (prev_is_intra) {
			m2d_mb_reset_inter(mb);
		}
		return m2d_parse_inter_macroblock(mb, stream);
	}
}

enum {
	M2D_ESCAPE = 0
};


/**Decode macroblock_address_increment.
 * Return maybe "1" in most cases.
 */
static int m2d_macroblock_address_increment(dec_bits *stream)
{
	if (get_onebit(stream) == 0) {
		int val = 0;
		do {
			int t = m2d_dec_vld_unary(stream, mb_inc_bit4, 4);
			val = val + t;
			if (t != M2D_ESCAPE) {
				break;
			} else {
				val += 33;
				if (get_onebit(stream) != 0) {
					val += 1;
					break;
				}
			}
		} while (1);
		return val;
	} else {
		/* short circuit */
		return 1;
	}
}


/**Reset macroblock position.
 * Pointer to destination frame also be updated.
 */
static void m2d_init_mb_pos(m2d_mb_current *mb)
{
	m2d_frames_set_mb_pos(mb, -1, 0, 0);
}

/**Increment macroblock position.
 * Pointer to destination frame also be updated.
 */
static void m2d_inc_mb_pos(m2d_mb_current *mb)
{
	int x = mb->mb_x + 1;
	int mb_width = mb->mbmax_x;
	int inc_x = 1;
	if (mb_width <= x) {
		int inc_y = 0;
		do {
			x -= mb_width;
			inc_x -= mb_width;
			inc_y += 1;
		} while (mb_width < x);
		m2d_frames_inc_mb_y_pos(mb->frames, inc_y, mb->frame_width);
		mb->mb_y += inc_y;
	}
	m2d_frames_inc_mb_x_pos(mb->frames, inc_x);
	mb->mb_x = x;
	if ((x & 1) == 0) {
		m2d_prefetch_mb(mb->frames->curr_luma, mb->frame_width);
	}
}

/** Check whether current position of macroblock have reached
 * boundary of frame or not.
 */
static int m2d_is_last(m2d_mb_current *mb)
{
	int y_max = mb->mbmax_y;
	int y = mb->mb_y;
	int x = mb->mb_x;
	return ((y == (y_max - 1) && (mb->mbmax_x - 1) <= x) || y_max <= y);
}

static void m2d_skip_rest_slices(m2d_mb_current *m2d)
{
}

/**Do macroblocks loop for one slice.
 */
static int m2d_decode_macroblocks(m2d_context *m2d)
{
	int err = 0;
	m2d_mb_current *mb = m2d->mb_current;
	dec_bits *stream = m2d->stream;
	m2d_mb_reset(mb);
	do {
		int mb_inc = m2d_macroblock_address_increment(stream);
		if (1 < mb_inc) {
			mb->skip_mb(mb, mb_inc);
		}
		m2d_inc_mb_pos(mb);
		err = m2d_parse_macroblock(mb, stream);
		if (m2d_is_last(mb)) {
			m2d_init_mb_pos(mb);
			m2d_skip_rest_slices(mb);
			err = 1;
			break;
		}
	} while ((0 <= err) && (show_bits(stream, 23) != 0));
	byte_align(stream);
	return err;
}

static void store_frame_info(m2d_frame_t *frame, const m2d_frames *frames, int idx, const m2d_seq_header *header)
{
	*frame = frames->frames[idx];
	frame->width =
#ifdef FAST_DECODE
		(unsigned)((header->horizontal_size_value + 15) & ~15) >> 3;
#else
		header->horizontal_size_value;
#endif
	frame->height =
#ifdef FAST_DECODE
		(unsigned)((header->vertical_size_value + 15) & ~15) >> 3;
#else
		header->vertical_size_value;
#endif
}

__LIBM2DEC_API int m2d_peek_decoded_frame(m2d_context *m2d, m2d_frame_t *frame, int is_end)
{
	m2d_frames *frames;
	int idx;

	if ((m2d == 0) || ((intptr_t)m2d & 3) || (frame == 0) || ((intptr_t)frame & 3)) {
		return -1;
	}
	frames = m2d->mb_current->frames;
	if (m2d->picture->picture_coding_type == B_VOP) {
		idx = frames->index;
	} else if (is_end && (0 < m2d->out_state) && (m2d->out_state < (2 * 2))) {
		idx = frames->idx_of_ref[1];
	} else {
		idx = frames->idx_of_ref[0];
	}
	store_frame_info(frame, frames, idx, m2d->seq_header);
	if (m2d->picture->picture_coding_type != B_VOP) {
		switch (m2d->out_state >> 1) {
		case 0:
			return 0;
		case 1:
			return (is_end != 0);
		case 2:
			return 1;
		}
	}
	return (m2d->out_state & 1);
}

__LIBM2DEC_API int m2d_get_decoded_frame(m2d_context *m2d, m2d_frame_t *frame, int is_end)
{
	int ret = m2d_peek_decoded_frame(m2d, frame, is_end);
	if (ret < 0) {
		return ret;
	}
	if (ret != 0) {
		if (m2d->picture->picture_coding_type == B_VOP) {
			m2d->out_state &= ~1;
		} else {
			m2d->out_state -= 2;
		}
	}
	return ret;
}

__LIBM2DEC_API int m2d_set_data(m2d_context *m2d, const byte_t *indata, int indata_bytes)
{
	if ((m2d == 0) || (indata == 0) || (indata_bytes == 0)) {
		return -1;
	}
	dec_bits_set_data(m2d->stream, indata, indata_bytes, 0);
	return 0;
}

__LIBM2DEC_API int m2d_decode_data(m2d_context *m2d)
{
	dec_bits *stream;
	int err;

	if (m2d == 0) {
		return -1;
	}
	m2d->picture->picture_coding_type = 0;
	stream = m2d->stream;
	err = 0;
	do {
		if (0 <= (err = m2d_find_mpeg_data(stream))) {
			int code_type = get_bits(stream, 8);
			err = m2d_dispatch_one_nal(m2d, code_type);
		} else {
			break;
		}
		VC_CHECK;
	} while (err != 1);
#ifdef DUMP_COEF
	print_coefs();
#endif
	return err;
}

__LIBM2DEC_API int m2d_read_header(m2d_context *m2d, const byte_t *data, size_t len)
{
	dec_bits *st;
	m2d_mb_current *mb;
	int err;
	int code_type;

	if (!m2d || !data || !m2d->stream) {
		return -1;
	}
	st = m2d->stream;
	err = dec_bits_set_data(st, data, len, 0);
	if (err < 0) {
		return err;
	}
	if (setjmp(st->jmp) != 0) {
		return 0;
	}
	mb = m2d->mb_current;
	bool seq_header_detect = false;
	do {
		err = m2d_find_mpeg_data(st);
		if (err < 0) {
			return err;
		}
		code_type = get_bits(st, 8);
		if ((code_type == 0xb3) || (seq_header_detect && (code_type == 0xb5))) {
			err = m2d_dispatch_one_nal(m2d, code_type);
			if (code_type == 0xb3) {
				seq_header_detect = true;
			} else {
				break;
			}
		} else if (seq_header_detect) {
			break;
		}
	} while (0 <= err);
	return err;
}

__LIBM2DEC_API int m2d_get_info(m2d_context *m2d, m2d_info_t *info)
{
	int src_width, src_height;
	if (!m2d || !info) {
		return -1;
	}
	src_width = m2d->seq_header->horizontal_size_value;
	src_height = m2d->seq_header->vertical_size_value;
	info->src_width = (src_width + 15) & ~15;
	info->src_height = (src_height + 15) & ~15;
	info->disp_width = m2d->seq_header->display_width;
	info->disp_height = m2d->seq_header->display_height;
	info->frame_num = 3;
	info->crop[0] = 0;
	info->crop[1] = info->src_width - src_width;
	info->crop[2] = 0;
	info->crop[3] = info->src_height - src_height;
	info->additional_size = 0;
	return 0;
}

__LIBM2DEC_API int m2d_skip_frames(m2d_context *m2d, int frame_num)
{
	int err;
	dec_bits *stream;

	if ((m2d == 0) || (frame_num <= 0)) {
		return -1;
	}
	if (setjmp(m2d->stream->jmp) != 0) {
		return 0;
	}
	stream = m2d->stream;
	do {
		int code_type;
		do {
			err = m2d_find_mpeg_data(stream);
			if (err < 0) {
				break;
			}
			code_type = get_bits(stream, 8);
			if ((code_type == 0xb3) || (code_type == 0xb5)) {
				m2d_dispatch_one_nal(m2d, code_type);
			}
		} while (code_type);
		if (err < 0) {
			break;
		}
	} while (--frame_num);
	return err;
}

#ifndef NDEBUG

#ifdef FAST_DECODE
int test_parse_coef()
{
	return 0;
}
#else

#include <vector>
#include <algorithm>
#include "txt2bin.h"

#ifndef __RENESAS_VERSION__
using namespace std;
#endif

#define LENGTH(a) (sizeof(a) / sizeof(a[0]))
#define ASSERT_EQUALS(a, b) (assert((a) == (b)))

struct Vld2TextData {
	int pos, level, length;
	const char *input;
};

/**Test pattern for DCT coefficient Table0.
 */
static const Vld2TextData test_vld_dct[] = {
	{1, 1, 5, "11010"},
	{1, -1, 5, "11110"},
	{2, 1, 6, "011010"},
	{2, -1, 6, "011110"},
	{1, 2, 7, "0100010"},
	{1, -2, 7, "0100110"},
	{10, 1, 10, "00001010 10"},
	{10, -1, 10, "00001011 10"},
	{2, 8, 18, "0000 0000 0011 1110 10"},
	{2, -8, 18, "0000 0000 0011 1111 10"},
};

int test_parse_coef()
{
	int err = 0;
//	unsigned char *in_data = new unsigned char[512];
	vector<unsigned char> in_data(512);
	m2d_mb_current mb;
//	int16_t coef[64];
	int16_t *coef;
	uint8_t qmat[2][64];
	fill(qmat[0], qmat[2], 16);
	vector<Vld2TextData> dct;
	for (size_t i = 0; i < LENGTH(test_vld_dct); ++i) {
		dct.push_back(test_vld_dct[i]);
	}
//	for_each(test_vld_dct, test_vld_dct + LENGTH(test_vld_dct), dct.push_back);

	dec_bits stream;
	dec_bits_open(&stream, 0);
	m2d_mb_set_default(&mb);
//	m2d_set_qmat(&mb, 0, qmat[0], 0, 1);
	coef = mb.coef;
	mb.qmat[0] = m2d_qmat_default[1];
	mb.q_mapping = m2d_q_mapping[1];
	mb.q_scale = mb.q_mapping[0];

	for (size_t i = 0; i < LENGTH(test_vld_dct); ++i) {
		const Vld2TextData& vld = test_vld_dct[i];
//		fill(in_data, in_data + 512, 0);
		fill(in_data.begin(), in_data.end(), 0);
		txt2bin(vld.input, &in_data[0]);
		dec_bits_set_data(&stream, &in_data[0], 512, 0);
		fill(coef, coef + 64, 0);
		int ret = m2d_parse_coef[1][0](&mb, &stream, 1);
//		ASSERT_EQUALS(vld.length, dec_bits_sum(&stream));
		printf("%d : %d\n", vld.level, coef[m2d_zigzag[0][vld.pos]]);
		ASSERT_EQUALS(vld.level, coef[m2d_zigzag[0][vld.pos]]);
	}
//	delete[] in_data;
	return err;
}

#endif /* FAST_DECODE */
#endif /* NDEBUG */

const m2d_func_table_t m2d_func_ = {
	sizeof(m2d_context),
	(int (*)(void *, int, int (*)(void *, void *), void *))m2d_init,
	(dec_bits *(*)(void *))m2d_stream_pos,
	(int (*)(void *, m2d_info_t *))m2d_get_info,
	(int (*)(void *, int, m2d_frame_t *, uint8_t *, int))m2d_set_frames,
	(int (*)(void *))m2d_decode_data,
	(int (*)(void *, m2d_frame_t *, int))m2d_peek_decoded_frame,
	(int (*)(void *, m2d_frame_t *, int))m2d_get_decoded_frame
};

const m2d_func_table_t * const m2d_func = &m2d_func_;
