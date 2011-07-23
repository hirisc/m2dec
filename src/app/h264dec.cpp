#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <vector>
#include <algorithm>
#include <functional>
#include "getopt.h"
#include "bitio.h"
#include "h264.h"
#include "md5.h"

#ifdef __RENESAS_VERSION__

#include <machine.h>

extern "C" {
#pragma section FILES

#define O_RDONLY 1
#define O_WRONLY 2
#define MAX_FILESIZE 8 * 1024 * 1024
#define MAX_MD5SIZE 128 * 1024
char InputFile[MAX_FILESIZE];
char OutputFile[MAX_MD5SIZE];

#pragma section
void abort() {while (1);}
static int infilepos;

void exit(int err) {
	while (err) sleep();
}

int fseek(FILE *fp, long offset, int whence) {
	return 0;
}

int open(char *name, int mode, int flg) {
	if (mode & O_RDONLY) {
		infilepos = 0;
	} else if (mode & O_WRONLY) {
		infilepos = 0;
	}
	return 4;
}
int close(int fineno) {return 0;}
int read(int fileno, char *buf, unsigned int count) {
	if (sizeof(InputFile) < infilepos + count) {
		count = sizeof(InputFile) - infilepos;
		if ((signed)count <= 0) {
			return 0;
		}
	}
	memcpy(buf, InputFile + infilepos, count);
	infilepos += count;
	return count;
}
int lseek() {return 0;}
int write(int fileno, char *buf, unsigned int count) {return count;}
char *getenv(const char *name) {return 0;}

}

#endif /* __RENESAS_VERSION__ */

struct Frames {
	std::vector<m2d_frame_t> frame;
	Frames(int width, int height, int num_mem) : frame(num_mem) {
		int luma_len = ((width + 15) & ~15) * ((height + 15) & ~15) + 15;
		if (num_mem <= 0 || luma_len <= 0) {
			return;
		}
		std::for_each(frame.begin(), frame.end(), std::bind2nd(Create(), luma_len));
	}
	~Frames() {
		if (frame.empty()) {
			return;
		}
		std::for_each(frame.begin(), frame.end(), Delete());
		frame.clear();
	}
private:
	struct Create : std::binary_function<m2d_frame_t, int, void> {
		void operator()(m2d_frame_t& frm, int luma_len) const {
			frm.luma = new uint8_t[luma_len];
			frm.chroma = new uint8_t[luma_len >> 1];
		}
	};
	struct Delete {
		void operator()(m2d_frame_t& frm) const {
			delete[] frm.luma;
			delete[] frm.chroma;
		}
	};
};

static bool outfilename(char *infilename, char *outfilename, size_t size)
{
	char *start;
	char *end;
	if (!(start = strrchr(infilename, '/')) && !(start = strrchr(infilename, '\\'))) {
		start = infilename;
	} else {
		start += 1;
	}
	end = strrchr(start, '.');
	if (end) {
		*end = '\0';
	}
	if (size <= strlen(start)) {
		return false;
	}
	sprintf(outfilename, "%s.out", start);
	return true;
}

struct input_data_t {
	uint8_t *data_;
	size_t len_;
	size_t pos_;
	void *var_;
	FILE *fo_;
	int dpb_;
	bool raw_;
	bool md5_;
	bool force_exec_;
	input_data_t(int argc, char *argv[], void *v) : pos_(0), var_(v), fo_(0), dpb_(-1), raw_(false), md5_(false), force_exec_(false) {
		FILE *fi;
		int opt;
		while ((opt = getopt(argc, argv, "bd:oOx")) != -1) {
			switch (opt) {
			case 'b':
				dpb_ = 1;
				break;
			case 'd':
				dpb_ = strtol(optarg, 0, 0);
				if (32 < (unsigned)dpb_) {
					BlameUser();
					/* NOTREACHED */
				}
				break;
			case 'O':
				md5_ = true;
				break;
			case 'o':
				raw_ = true;
				break;
			case 'x':
				force_exec_ = true;
				break;
			default:
				BlameUser();
				/* NOTREACHED */
			}
		}
		if ((argc <= optind) ||	!(fi = fopen(argv[optind], "rb"))) {
			BlameUser();
			/* NOTREACHED */
		}
		if (md5_ || raw_) {
			char outfile[256];
			if (outfilename(argv[optind], outfile, sizeof(outfile))) {
				fo_ = fopen(outfile, "wb");
			}
		}
#ifdef __RENESAS_VERSION__
		data_ = (uint8_t *)InputFile;
		len_ = sizeof(InputFile);
#else
		fseek(fi, 0, SEEK_END);
		len_ = ftell(fi);
		fseek(fi, 0, SEEK_SET);
		data_ = new uint8_t[len_];
		fread(data_, 1, len_, fi);
#endif
		fclose(fi);
	}
	~input_data_t() {
#ifndef __RENESAS_VERSION__
		delete[] data_;
#endif
		if (fo_) {
			fclose(fo_);
		}
	}
	static size_t write_md5(void *dst, const uint8_t *src, size_t size) {
		md5_append((md5_state_t *)dst, (const md5_byte_t *)src, (int)size);
		return size;
	}
	static size_t write_raw(void *dst, const uint8_t *src, size_t size) {
		return fwrite((const void *)src, 1, size, (FILE *)dst);
	}
	static void write_cropping(const m2d_frame_t *frame, size_t (*writef)(void *dst, const uint8_t *src, size_t size), void *dst) {
		const uint8_t *src;
		int stride, height, width;

		stride = frame->width;
		src = frame->luma + stride * frame->crop[2] + frame->crop[0];
		height = frame->height - frame->crop[2] - frame->crop[3];
		width = stride - frame->crop[0] - frame->crop[1];
		for (int y = 0; y < height; ++y) {
			writef(dst, src, width);
			src += stride;
		}
		src = frame->chroma + stride * (frame->crop[2] >> 1) + frame->crop[0];
		height >>= 1;
		for (int y = 0; y < height; ++y) {
			writef(dst, src, width);
			src += stride;
		}
	}
	static void md5txt(const unsigned char *md5, char *txt)
	{
		for (int i = 0; i < 16; ++i) {
			sprintf(txt + i * 2, "%02xd", *md5++);
		}
		txt[32] = 0x0d;
		txt[33] = 0x0a;
	}
	size_t writeframe(const m2d_frame_t *frame) {
		if (fo_) {
			int luma_size = frame->width * frame->height;
			if (md5_) {
				md5_state_t md5;
				md5_byte_t digest[16];
				char txt[16 * 2 + 2];

				md5_init(&md5);
				write_cropping(frame, write_md5, (void *)&md5);
				md5_finish(&md5, digest);
				md5txt(digest, txt);
				fwrite(txt, 1, sizeof(txt), fo_);
			} else {
				write_cropping(frame, write_raw, (void *)fo_);
			}
			fflush(fo_);
			return luma_size;
		}
		return 0;
	}
	void BlameUser() {
		fprintf(stderr,
			"Usage:\n"
			"\th264dec [-b] [-d <dpb_size>] [-o|O ] <infile>\n"
			"\t\t-b: Bypass DPB\n"
			"\t\t-d <dpb_size>: Specify number of DPB frames -1, 1..16 (default: -1(auto))\n"
			"\t\t-o: RAW output\n"
			"\t\t-O: MD5 output\n"
			"\t\t-x: Mask SIGABRT on error."
			);
		exit(1);
	}
};

static int reread_file(void *var)
{
	input_data_t *d = (input_data_t *)var;
	if (d->pos_ < d->len_) {
		dec_bits_set_data(((h264d_context *)d->var_)->stream, d->data_ + d->pos_, d->len_);
		d->pos_ += d->len_;
		return 0;
	} else {
		return -1;
	}
}

#ifdef _M_IX86
#include <crtdbg.h>
#elif defined(__linux__)
#include <pthread.h>
#include <signal.h>
static void trap(int no)
{
	fprintf(stderr, "trap %d\n", no);
	exit(0);
}
#endif

int main(int argc, char *argv[])
{
	h264d_context *h2d;
	int err;

#ifdef _M_IX86
	_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_WNDW);
//	_CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_WNDW);
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF|_CRTDBG_LEAK_CHECK_DF);
#endif
	h2d = new h264d_context;
	input_data_t data(argc, argv, (void *)h2d);
	if (data.len_ <= 0) {
		delete h2d;
		return -1;
	}
#ifdef __linux__
	if (data.force_exec_) {
		struct sigaction sa;
		memset(&sa, 0, sizeof(sa));
		sa.sa_handler = trap;
		sigaction(SIGABRT, &sa, 0);
		sigaction(SIGSEGV, &sa, 0);
	}
#endif
	err = h264d_init(h2d, data.dpb_, 0, 0);
	if (err) {
		return err;
	}
	err = h264d_read_header(h2d, data.data_, data.len_);
	data.pos_ = data.len_;
	if (err) {
		return err;
	}
	m2d_info_t info;
	h264d_get_info(h2d, &info);
	fprintf(stderr,
		"length: %ld\n"
		"width x height x num: %d x %d x %d\n"
		"Context Size: %ld\n",
		data.len_, info.src_width, info.src_height, info.frame_num, sizeof(h264d_context));
	info.frame_num = (info.frame_num < 3 ? 3 : info.frame_num) + (data.dpb_ < 0 ? 16 : data.dpb_);
	Frames frm(info.src_width, info.src_height, info.frame_num);
	std::vector<uint8_t> second_frame(info.additional_size);
	err = h264d_set_frames(h2d, info.frame_num, &frm.frame[0], &second_frame[0], info.additional_size);
	if (err) {
		return err;
	}
	dec_bits_set_callback(h2d->stream, reread_file, &data);
	m2d_frame_t frame;
	for (int i = 0; i < INT_MAX; ++i) {
		err = h264d_decode_picture(h2d);
		if (err == -1) {
			break;
		}
		while (h264d_get_decoded_frame(h2d, &frame, (data.dpb_ == 1))) {
			data.writeframe(&frame);
		}
		if (err < 0) {
			break;
		}
	}
	while (h264d_get_decoded_frame(h2d, &frame, 1)) {
		data.writeframe(&frame);
	}
	delete h2d;

#ifdef _M_IX86
	assert(_CrtCheckMemory());
#endif
	return (data.force_exec_) ? 0 : ((err == -2) ? 0 : err);
}

