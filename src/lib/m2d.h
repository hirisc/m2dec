/** Yet Another Video decoder
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

#ifndef __M2D_H__
#define __M2D_H__

#include "m2types.h"
#include "bitio.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	uint8_t *luma;
	uint8_t *chroma;
	int16_t width, height;
	int16_t crop[4];
} m2d_frame_t;

typedef struct {
	int16_t pattern;
	int8_t length;
} vlc_t;

typedef struct {
	int16_t src_width, src_height;
	int16_t disp_width, disp_height;
	int16_t frame_num;
	int16_t crop[4];
	int additional_size;
} m2d_info_t;

typedef struct {
	size_t context_size;
	int (*init)(void *, int);
	dec_bits *(*stream_pos)(void *);
	int (*read_header)(void *, const byte_t *, size_t);
	int (*get_info)(void *, m2d_info_t *);
	int (*set_frames)(void *, int, m2d_frame_t *, uint8_t *, int);
	int (*decode_picture)(void *);
	int (*get_decoded_frame)(void *, m2d_frame_t *, int);
} m2d_func_table_t;

int m2d_dec_vld_unary(dec_bits *stream, const vlc_t *vld_tab, int bitlen);
int m2d_find_mpeg_data(dec_bits *stream);

#ifdef __cplusplus
}
#endif

#endif /* __M2D_H__ */
