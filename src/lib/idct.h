/** iDCT for MPEG-2 or WMV8.
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

#ifndef _IDCT_H_
#define _IDCT_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "m2types.h"

#ifdef __RENESAS_VERSION__
void m2d_prefetch_mb(uint8_t *dst, int32_t dst_stride);
void m2d_clear_coef(int16_t *dst, int idx);
#else
#define m2d_prefetch_mb(a, b)
#define m2d_clear_coef(a, b) memset((a) + (b), 0, sizeof(*(a)) * ((MB_LEN * MB_LEN / 4) - b))
#endif

void m2d_idct_intra_luma(uint8_t *dst_base, int32_t dst_stride, int16_t *src_coef, uint32_t coef_exist);
void m2d_idct_intra_chroma(uint8_t *dst_base, int32_t dst_stride, int16_t *src_coef, uint32_t coef_exist);
void m2d_idct_inter_luma(uint8_t *dst_base, int32_t dst_stride, int16_t *src_coef, uint32_t coef_exist);
void m2d_idct_inter_chroma(uint8_t *dst_base, int32_t dst_stride, int16_t *src_coef, uint32_t coef_exist);

#ifdef __cplusplus
}
#endif

#endif /* _IDCT_H_ */
