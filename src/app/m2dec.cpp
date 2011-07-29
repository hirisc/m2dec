/** MPEG-2 decoder application example
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <setjmp.h>
#include <memory>
#include "bitio.h"
#include "mpeg2.h"
#include "frames.h"
#include "mpeg_demux.h"
#include "txt2bin.h"
#include "display.h"
#ifdef ENABLE_AALIB
#include "aadisp.h"
#endif

#ifdef __RENESAS_VERSION__
extern char InputFile[];
extern const int InputSize;
#endif /* __RENESAS_VERSION__ */

const int FRAME_NUM = 6;
const size_t INDATA_MAX = 1 * 1024 * 1024;
const int WIDTH = 720;
const int HEIGHT = 576;
const int LUMA_LEN = WIDTH * HEIGHT;
const int BUF_SIZE = 1920 * 1080 * 3 / 2;
#include "optparser.h"

struct reread_t {
	m2d_context *m2d;
	Options *opt;
	pes_demuxer_t dmx;
};
static int reread_file(void *arg);
static int reread_packet(void *arg);

int error_report(void *arg)
{
	longjmp(*(jmp_buf *)arg, -1);
	return 0;
}


class random_int : public unary_function<int, int> {
public:
	int operator()() {
		int r = rand();
		return (r << 24) | r;
	}
};

static int test_dec_bits()
{
	dec_bits stream0, stream1;
	jmp_buf jb0, jb1;

	vector<int> buf(BUF_SIZE);
	generate(buf.begin(), buf.end(), random_int());
	int err = 0;

	err |= dec_bits_open(&stream0, 0);
	err |= dec_bits_open(&stream1, 0);

	dec_bits_set_callback(&stream0, error_report, &jb0);
	dec_bits_set_callback(&stream1, error_report, &jb1);
	dec_bits_set_data(&stream0, (byte_t *)&buf[0], BUF_SIZE);
	dec_bits_set_data(&stream1, (byte_t *)&buf[0], BUF_SIZE);
	int bit = 0;
	if ((setjmp(jb0) == 0) && (setjmp(jb1) == 0)) {
		for (int i = 0; i < BUF_SIZE / 2; ++i) {
//			bit++;
//			bit = (bit <= 24) ? bit : 1;
			bit = (rand() % 23) + 1;
			cache_t a = show_bits(&stream0, bit);
			//		cache_t b = get_bits(&stream0, bit);
			cache_t b = get_bits(&stream1, bit);
			skip_bits(&stream0, bit);
			if (a != b) {
				printf("%d: %08x != %08x\n", bit, a, b);
				break;
			}
			assert(a == b);
//			printf("%08x\n", a);
		}
	} else {
	}
	return err;
}


extern "C" {

//extern int (* const m2d_parse_coef_mpeg2[2])(m2d_mb_current *mb, dec_bits *stream, int idx, int is_inter);

}

int test_parse_coef();


#define LENGTH(a) (sizeof(a) / sizeof(a[0]))
#define ASSERT_EQUALS(a, b) (assert((a) == (b)))

void m2d_mb_set_default(m2d_mb_current *mb);
void m2d_set_qmat(m2d_mb_current *mb, dec_bits *stream, uint8_t *qmat, int load_bit, int qmat_idx);


#include "vld.h"

struct DcDiffTextData {
	int block_idx; // 0, 1, 2
	int ret;
	const char *input;
};

/**Test pattern for DCT coefficient Table0.
 */
static const DcDiffTextData test_intra_dc_in[] = {
	{0, 0, "100" ""},
	{0, -1, "00" "0"},
	{0, 1, "00" "1"},

	{0, -3, "01" "00"},
	{0, -2, "01" "01"},
	{0, 2, "01" "10"},
	{0, 3, "01" "11"},

	{0, -7, "101" "000"},
	{0, -6, "101" "001"},
	{0, -5, "101" "010"},
	{0, -4, "101" "011"},
	{0, 4, "101" "100"},
	{0, 5, "101" "101"},
	{0, 6, "101" "110"},
	{0, 7, "101" "111"},

	{0, -15, "110" "0000"},
	{0, -14, "110" "0001"},
	{0, -13, "110" "0010"},
	{0, -12, "110" "0011"},
	{0, -11, "110" "0100"},
	{0, -10, "110" "0101"},
	{0, -9, "110" "0110"},
	{0, -8, "110" "0111"},
	{0, 8, "110" "1000"},
	{0, 9, "110" "1001"},
	{0, 10, "110" "1010"},
	{0, 11, "110" "1011"},
	{0, 12, "110" "1100"},
	{0, 13, "110" "1101"},
	{0, 14, "110" "1110"},
	{0, 15, "110" "1111"},

	{0, -2047, "111111111" "00000000000"},
	{0, -1024, "111111111" "01111111111"},
	{0, 1024, "111111111" "10000000000"},
	{0, 2047, "111111111" "11111111111"},
};

int m2d_parse_intra_dc(m2d_mb_current *mb, dec_bits *stream, int block_idx);

static int test_intra_dc()
{
	int err = 0;
	vector<unsigned char> in_data(512);
	m2d_mb_current mb;
	dec_bits stream;
	dec_bits_open(&stream, 0);
	m2d_mb_set_default(&mb);
	for (size_t i = 0; i < LENGTH(test_intra_dc_in); ++i) {
		const DcDiffTextData *in = &test_intra_dc_in[i];
		int len = (txt2bin(in->input, &in_data[0]) + 7) >> 3;
		for (size_t dc_scale = 0; dc_scale <= 3; ++dc_scale) {
			for (size_t dc = 0; dc < 255; ++dc) {
				dec_bits_set_data(&stream, &in_data[0], len);
				mb.intra_dc_scale = (int8_t)dc_scale;
				mb.intra_dc_max = (int16_t)((1 << (8 + 3 - dc_scale)) - 1);
				fill(mb.dc_pred, mb.dc_pred + 3, (int16_t)dc);
				int ret = m2d_parse_intra_dc(&mb, &stream, in->block_idx);
				int ref_ret = in->ret + (int)dc;
				if (ref_ret < 0) {
					ref_ret = 0;
				} else if (mb.intra_dc_max < ref_ret) {
					ref_ret = mb.intra_dc_max;
				}
				ASSERT_EQUALS(ret, ref_ret << dc_scale);
			}
		}
	}

	return err;
}

int test_idct_horizontal();
int test_motioncomp();
int test_buffer();

static int test_all()
{
#ifdef UNIT_TEST

#if defined(__RENESAS_VERSION__)
	return test_idct_horizontal() + test_motioncomp();
#elif !defined(FAST_DECODE)
	return test_buffer() | test_dec_bits() | test_parse_coef() | test_intra_dc();
#endif

#else
	return 0;
#endif
}

#ifdef _M_IX86
#include <crtdbg.h>
#endif

#ifdef ENABLE_DISPLAY
#include <SDL/SDL.h>
#endif

static int set_packet_mode(reread_t *re)
{
	m2d_context *m2d = re->m2d;
	int err;

	pes_demuxer_t *dmx = &re->dmx;
	mpeg_demux_init(dmx, reread_file, re);
	dec_bits_set_callback(m2d->stream, reread_packet, 0);
	re->opt->set_demuxer_mode(true);
	err = m2d_decode_data(m2d);
	return err;
}

static int decode_one_frame(reread_t *re)
{
	m2d_context *m2d = re->m2d;
	int err = m2d_decode_data(m2d);
	return err;
}

#define ALIGN16(x) (((x) + 15) & ~15)
#define ALIGN16_U8(x) (uint8_t *)(((uintptr_t)(x) + 15) & ~15)

static int reread_packet(void *arg)
{
	reread_t *re = (reread_t *)arg;
	pes_demuxer_t *dmx = &re->dmx;
	const byte_t *packet;
	int packet_size;
	packet = mpeg_demux_get_video(dmx, &packet_size);
	if (packet) {
		dec_bits_set_data(re->m2d->stream, packet, (size_t)packet_size);
		return 0;
	} else {
		return -1;
	}
}

static int reread_file(void *arg)
{
	reread_t *re = (reread_t *)arg;
	Options *opt = re->opt;
	dec_bits *stream = re->dmx.stream ? re->dmx.stream : re->m2d->stream;
	const byte_t *indata;
	int size;
	if (opt->next_buffer(&indata, &size) < 0) {
		return -1;
	} else {
		dec_bits_set_data(stream, indata, size);
		return 0;
	}
}

int mpeg_demux_request_set_data(pes_demuxer_t *dmx, const byte_t **indata_p, int *size_p, void *arg)
{
	Options *opt = (Options *)arg;
	if (opt->next_buffer(indata_p, size_p) < 0) {
		return -1;
	} else {
		return 0;
	}
}

#ifdef main
#undef main
#endif


int main(int argc, char **argv)
{
	Options opt;

#ifdef _M_IX86
	_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_WNDW);
//	_CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_WNDW);
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF|_CRTDBG_LEAK_CHECK_DF);
#endif
	assert(test_all() == 0);

	if (opt.parse(argv, argc) < 0) {
		return -1;
	}
	auto_ptr<m2d_context> m2d(new m2d_context);
	Frames frms(WIDTH, HEIGHT, FRAME_NUM, 0);
	int err = 0;
	err |= m2d_init(m2d.get(), 0, 0, 0);
	err |= m2d_set_frames(m2d.get(), FRAME_NUM, frms.aligned());

	reread_t reread_data = {
		m2d.get(),
		&opt
	};
	dec_bits_set_callback(m2d.get()->stream, reread_file, &reread_data);

	int frame_num = opt.frame_num_;
	for (int s = 0; s < 1; ++s) {
		if (opt.is_demuxer_mode()) {
			set_packet_mode(&reread_data);
		}
		m2d_skip_frames(m2d.get(), opt.skip_frames_);
		for (int i = 0; i < frame_num; ++i) {
			m2d_frame_t frm;

			err = decode_one_frame(&reread_data);
			int is_end = (err != 1);
			if (i != 0) {
				if (0 < m2d_get_decoded_frame(m2d.get(), &frm, is_end)) {
					opt.write_file(frm.luma, frm.chroma, frm.width, frm.height);
					if (err != 1) {
						opt.write_file(frm.luma, frm.chroma, frm.width, frm.height); /* output delayed frame */
					}
				}
			}
			if ((err < 0) && (err != M2D_ERR_NOT_VIDEO)) {
				break;
			}
		}
	}

	if (opt.verbose_) {
		m2d_frame_t frm;
		m2d_get_decoded_frame(m2d.get(), &frm, 0);
		printf(
			"width, height: %dx%d\n"
			"sizeof m2d: %d\n"
			"sizeof m2d_mb_current: %d\n",
			frm.width, frm.height,
			sizeof(m2d_context),
			sizeof(m2d_mb_current));
	}

#ifdef _M_IX86
	assert(_CrtCheckMemory());
#endif
	return err;
}

