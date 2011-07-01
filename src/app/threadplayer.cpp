#ifdef ENABLE_DISPLAY
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <list>
#include <map>
#include "getopt.h"
#include "md5.h"
#include "bitio.h"
#include "mpeg2.h"
#include "mpeg_demux.h"
#include "h264.h"
#ifdef __GNUC__
#include <stdint.h>
#endif


#include "unithread.h"

struct Buffer {
	uint8_t *data;
	int len;
};

template <typename T>
class Queue {
	T *data_;
	int max_;
	int head_;
	int tail_;
	UniMutex *mutex_;
	UniCond *cond_;
	bool terminated_;
	int next(int idx) const {
		++idx;
		return idx < max_ ? idx : 0;
	}
public:
	Queue(T *data, int buf_num)
		: data_(data), max_(buf_num), head_(0), tail_(0),
		  mutex_(UniCreateMutex()), cond_(UniCreateCond()),
		  terminated_(false) {}
	~Queue() {
		UniDestroyCond(cond_);
		UniDestroyMutex(mutex_);
	}
	bool full() const {
		return next(head_) == tail_;
	}
	bool empty() const {
		return head_ == tail_;
	}
	int size() const {
		return max_;
	}
	bool terminated() const {
		return terminated_;
	}
	void terminate() {
		UniLockMutex(mutex_);
		terminated_ = true;
		UniCondSignal(cond_);
		UniUnlockMutex(mutex_);
	}
	T& emptybuf() {
		UniLockMutex(mutex_);
		while (full()) {
			UniCondWait(cond_, mutex_);
		}
		T& d = data_[head_];
		UniUnlockMutex(mutex_);
		return d;
	}
	void setfilled(T& dat) {
		assert(!full());
		UniLockMutex(mutex_);
		data_[head_] = dat;
		head_ = next(head_);
		UniCondSignal(cond_);
		UniUnlockMutex(mutex_);
	}
	T& getfilled() {
		UniLockMutex(mutex_);
		while (empty() && !terminated()) {
			UniCondWait(cond_, mutex_);
		}
		int tail = tail_;
		tail_ = next(tail_);
		T& d = data_[tail];
		UniCondSignal(cond_);
		UniUnlockMutex(mutex_);
		return d;
	}
};

typedef Queue<Buffer> FileQueueType;

class FileReader {
	FileQueueType outqueue_;
	std::list<char *> infiles_;
public:
	FILE *fd_;
	int insize_;
	FileReader(Buffer *src_p, int buf_num, FILE *fd, int insize, std::list<char *> &infiles)
		: outqueue_(FileQueueType(src_p, buf_num)), infiles_(infiles), fd_(fd), insize_(insize) {}
	FileQueueType& outqueue() {
		return outqueue_;
	}
	std::list<char *> &infiles() {
		return infiles_;
	}

	static int run(void *data) {
		RecordTime(1);
		FileReader *fr = (FileReader *)data;
		while (1) {
			Buffer& buf = fr->outqueue().emptybuf();
			int read_size = fread(buf.data, 1, fr->insize_, fr->fd_);
			buf.len = read_size;
			if (read_size == 0) {
				fclose(fr->fd_);
				if (!fr->infiles().empty() && (fr->fd_ = fopen(fr->infiles().front(), "rb"))) {
					fr->infiles().pop_front();
					buf.len = fread(buf.data, 1, fr->insize_, fr->fd_);
				} else {
					fr->outqueue().setfilled(buf);
					fr->outqueue().terminate();
					fprintf(stderr, "File terminate.\n");
					RecordTime(0);
					return 0;
				}
			}
			fr->outqueue().setfilled(buf);
		}
		return -1;
	}
};

struct Frame {
	uint8_t *data, *data2;
	int width, height;
	int crop_right, crop_bottom;
};

typedef Queue<Frame> FrameQueueType;

static int reread_file(void *arg);
static int reread_packet(void *arg);

enum {
	MODE_MPEG2,
	MODE_MPEG2PS,
	MODE_H264
};

class M2Decoder {
	uint8_t **frame_org_;
	uint8_t *second_frame_;
	int frame_num_;
	FileQueueType& inqueue_;
	FrameQueueType outqueue_;
	byte_t *context_;
	const m2d_func_table_t *func_;
	pes_demuxer_t demux_;
	const byte_t *curr_data_;
	int curr_data_len_;
public:
	M2Decoder(FileQueueType& inqueue, Frame *dst, int dstnum, int codec_mode)
		: inqueue_(inqueue), outqueue_(FrameQueueType(dst, dstnum)), context_(0), curr_data_(0) {
		switch (codec_mode) {
		case MODE_MPEG2:
			/* FALLTHROUGH */
		case MODE_MPEG2PS:
			func_ = m2d_func;
			break;
		case MODE_H264:
			func_ = h264d_func;
			break;
		}
		context_ = new byte_t[func_->context_size];
		func_->init(context_, -1);
		memset(&demux_, 0, sizeof(demux_));
		if (codec_mode == MODE_MPEG2PS) {
			mpeg_demux_init(&demux_, reread_file, this);
			dec_bits_set_callback(func_->stream_pos(context_), reread_packet, this);
			reread_packet((void *)this);
		} else {
			dec_bits_set_callback(func_->stream_pos(context_), reread_file, this);
			reread_file((void *)this);
		}
		func_->read_header(context_, curr_data_, curr_data_len_);
		m2d_info_t info;
		func_->get_info(context_, &info);
		int len = ((info.src_width + 15) & ~15) * ((info.src_height + 15) & ~15);
		int bufnum = info.frame_num + ((codec_mode == MODE_H264) ? 16 : 0);
		frame_num_ = bufnum;
		std::vector<m2d_frame_t> frame(bufnum);
		frame_org_ = new uint8_t *[bufnum];
		for (int i = 0; i < bufnum; ++i) {
			frame_org_[i] = new uint8_t[((len * 3) >> 1) + 31];
			frame[i].luma = (uint8_t *)((((intptr_t)frame_org_[i]) + 15) & ~15);
			frame[i].chroma = frame[i].luma + len;
		}
		second_frame_ = new byte_t[info.additional_size];
		func_->set_frames(context_, bufnum, &frame[0], second_frame_, info.additional_size);
	}
	~M2Decoder() {
		delete[] second_frame_;
		for (int i = 0; i < frame_num_; ++i) {
			delete[] frame_org_[i];
		}
		delete[] frame_org_;
		delete[] context_;
	}
	void memory_current(const byte_t *curr_data, int curr_data_len) {
		curr_data_ = curr_data;
		curr_data_len_ = curr_data_len;
	}
	FileQueueType& inqueue() {
		return inqueue_;
	}
	FrameQueueType& outqueue() {
		return outqueue_;
	}
	void *context() {
		return context_;
	}
	const m2d_func_table_t *func() {
		return func_;
	}
	dec_bits *stream() {
		return func_->stream_pos(context_);
	}
	pes_demuxer_t *demuxer() {
		return &demux_;
	}
	void post_dst(m2d_frame_t& frm) {
		Frame& dst = outqueue().emptybuf();
		dst.width = (frm.width + 15) & ~15;
		dst.height = (frm.height + 15) & ~15;
		dst.crop_right = frm.crop[0] + frm.crop[1];
		dst.crop_bottom = frm.crop[2] + frm.crop[3];
		dst.data = frm.luma + frm.crop[2] * frm.width + frm.crop[0];
		dst.data2 = frm.chroma + ((frm.crop[2] * frm.width) >> 1) + frm.crop[0];
		outqueue().setfilled(dst);
	}
	static int run(void *data) {
		RecordTime(1);
		M2Decoder *module = (M2Decoder *)data;
		while (1) {
			m2d_frame_t frm;
			while (module->func()->get_decoded_frame(module->context(), &frm, 0) <= 0) {
				int err = module->func()->decode_picture(module->context());
				if (err < 0) {
					while (module->func()->get_decoded_frame(module->context(), &frm, 1)) {
						module->post_dst(frm);
					}
					module->outqueue().terminate();
					RecordTime(0);
					return 0;
				}
			}
			module->post_dst(frm);
		}
		return -1;
	}
};

static int reread_file(void *arg)
{
	M2Decoder *module = (M2Decoder *)arg;
	Buffer& src = module->inqueue().getfilled();
	if (src.len <= 0) {
		module->memory_current(0, 0);
		return -1;
	} else {
		dec_bits *stream = module->demuxer()->stream;
		stream = stream ? stream : module->stream();
		module->memory_current(src.data, src.len);
		dec_bits_set_data(stream, src.data, src.len);
		return 0;
	}
}

static int reread_packet(void *arg)
{
	M2Decoder *module = (M2Decoder *)arg;
	pes_demuxer_t *dmx = module->demuxer();
	const byte_t *packet;
	int packet_size;
	packet = mpeg_demux_get_video(dmx, &packet_size);
	module->memory_current(packet, packet_size);
	if (packet) {
		dec_bits_set_data(module->stream(), packet, (size_t)packet_size);
		return 0;
	} else {
		return -1;
	}
}

void display_write(uint8_t **dst, const uint8_t *src_luma, const uint8_t *src_chroma, int src_stride,
		   uint16_t *pitches, int width, int height)
{
	uint8_t *dst0 = dst[0];
	for (int i = 0; i < height; ++i) {
		memcpy(dst0, src_luma, width);
		src_luma += src_stride;
		dst0 += pitches[0];
	}
	width >>= 1;
	height >>= 1;
	uint8_t *dst1 = dst[1];
	uint8_t *dst2 = dst[2];
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			dst1[j] = src_chroma[j * 2];
			dst2[j] = src_chroma[j * 2 + 1];
		}
		src_chroma += src_stride;
		dst1 += pitches[1];
		dst2 += pitches[2];
	}
}

const int MAX_WIDTH = 1920;
const int MAX_HEIGHT = 1088;
const int MAX_LEN = (MAX_WIDTH * MAX_HEIGHT * 3) >> 1;
const int FILE_READ_SIZE = 65536 * 7;
const int BUFNUM = 5;

static void BlameUser() {
	fprintf(stderr,
		"Usage: srview [-s] [-r] [-t interval] [-m outfile(MD5)] [-o outfile(Raw)] infile [infile ...]\n"
		"\t-h : H.264 Elementary Data\n"
		"\t-s : MPEG-2 Program Stream (PS)\n"
		"\t-r : repeat\n"
		"\t-l : log dump\n"
		"\t-t interval : specify interval of each frame in ms unit\n");
	exit(-1);
}

struct Options {
	int interval_;
	std::list<char *> infile_list_;
	FILE *infile_, *outfile_, *outfilemd5_;
	int codec_mode_;
	bool repeat_;
	bool logdump_;

	Options(int argc, char **argv)
		: interval_(0), infile_(0), outfile_(0), outfilemd5_(0),
		  codec_mode_(MODE_MPEG2), repeat_(false), logdump_(false) {
		int opt;
		while ((opt = getopt(argc, argv, "hlm:o:rst:")) != -1) {
			switch (opt) {
			case 'h':
				codec_mode_ = MODE_H264;
				break;
			case 'l':
				logdump_ = true;
				break;
			case 'm':
				outfilemd5_ = fopen(optarg, "wb");
				break;
			case 'o':
				outfile_ = fopen(optarg, "wb");
				break;
			case 'r':
				repeat_ = true;
				break;
			case 's':
				codec_mode_ = MODE_MPEG2PS;
				break;
			case 't':
				interval_ = strtoul(optarg, 0, 0);
				if (interval_ == 0) {
					interval_ = 1;
				}
				break;
			default:
				BlameUser();
				/* NOTREACHED */
			}
		}
		if ((argc <= optind) ||	!(infile_ = fopen(argv[optind], "rb"))) {
			BlameUser();
			/* NOTREACHED */
		}
		while (++optind < argc) {
			infile_list_.push_back(argv[optind]);
		}
	}
};

static void waitevents_with_timer(SDL_Overlay *yuv, SDL_Rect *rect, int &suspended)
{
	UniEvent ev;
	UniWaitEvent(&ev);
	do {
		switch (ev.type) {
		case SDL_USEREVENT:
			SDL_DisplayYUVOverlay(yuv, rect);
			break;
		case SDL_MOUSEBUTTONDOWN:
			suspended ^= 1;
			break;
		case SDL_QUIT:
			exit(0);
			/* NOTREACHED */
			break;
		}
	} while (SDL_PollEvent(&ev));
}

static void waitevents(SDL_Overlay *yuv, SDL_Rect *rect,  int &suspended)
{
	UniEvent ev;
	SDL_DisplayYUVOverlay(yuv, rect);
	while (SDL_PollEvent(&ev)) {
		switch (ev.type) {
		case SDL_MOUSEBUTTONDOWN:
			suspended ^= 1;
			break;
		case SDL_QUIT:
			exit(0);
			/* NOTREACHED */
			break;
		}
	}
}

static void fwrites(void *fout, const byte_t *src, int size)
{
	fwrite(src, 1, size, (FILE *)fout);
}

static int file_write(const byte_t *luma, const byte_t *chroma, int src_stride, int width, int height, void *out, void (*Write)(void *out, const byte_t *src, int size))
{
	for (int y = 0; y < height; ++y) {
		Write(out, luma, width);
		luma += src_stride;
	}
	height >>= 1;
	for (int y = 0; y < height; ++y) {
		Write(out, chroma, width);
		chroma += src_stride;
	}
	return ((width * height) * 3) >> 1;
}

static void bin2hex(const md5_byte_t *md5, char *hex)
{
	for (int i = 0; i < 16; ++i) {
		sprintf(hex + i * 2, "%02xd", *md5++);
	}
	hex[32] = 0x0d;
	hex[33] = 0x0a;
}

void run_loop(Options& opt) {
	int width = 0;
	int height = 0;
	SDL_Surface *surface;
	SDL_Overlay *yuv = 0;
	SDL_Rect rect = {
		0, 0, width, height
	};

	LogTags.insert(std::pair<int, const char *> (UniThreadID(), "Main"));
	RecordTime(1);

	/* Run File Loader */
	Buffer src[BUFNUM];
	for (int i = 0; i < BUFNUM; ++i) {
		src[i].data = new unsigned char[FILE_READ_SIZE];
		src[i].len = FILE_READ_SIZE;
	}
	FileReader fr(src, BUFNUM, opt.infile_, FILE_READ_SIZE, opt.infile_list_);
	UniThread *thr_file = UniCreateThread(FileReader::run, (void *)&fr);
	LogTags.insert(std::pair<int, const char *> (UniGetThreadID(thr_file), "FileLoader"));

	/* Run Video Decoder */
	Frame dst_align[BUFNUM];
 	M2Decoder m2dec(fr.outqueue(), dst_align, BUFNUM, opt.codec_mode_);
	UniThread *thr_m2d = UniCreateThread(M2Decoder::run, (void *)&m2dec);
	FrameQueueType &outqueue = m2dec.outqueue();
	LogTags.insert(std::pair<int, const char *> (UniGetThreadID(thr_m2d), "Decoder"));

	UniEvent ev;
	while (UniPollEvent(&ev)) ;
	int suspended = 0;
	while (1) {
		if (!suspended) {
			if (outqueue.terminated() && outqueue.empty()) {
				goto endloop;
			}
			Frame& out = outqueue.getfilled();
			if ((width != (out.width - out.crop_right)) || (height != (out.height - out.crop_bottom))) {
				width = out.width - out.crop_right;
				height = out.height - out.crop_bottom;
				rect.w = width;
				rect.h = height;
				fprintf(stderr, "size: %d x %d\n", width, height);
				surface = SDL_SetVideoMode(width, height, 0, SDL_SWSURFACE);
				if (yuv) {
					SDL_FreeYUVOverlay(yuv);
				}
				yuv = SDL_CreateYUVOverlay(width, height, SDL_IYUV_OVERLAY, surface);
			}
			SDL_LockYUVOverlay(yuv);
			display_write(yuv->pixels, out.data, out.data2, out.width, yuv->pitches, yuv->w, yuv->h);
			SDL_UnlockYUVOverlay(yuv);
			if (opt.outfile_) {
				file_write(out.data, out.data2, out.width, width, height, opt.outfile_, fwrites);
			}
			if (opt.outfilemd5_) {
				md5_state_t md5;
				md5_byte_t digest[16];
				char hex[32 + 2];

				md5_init(&md5);
				file_write(out.data, out.data2, out.width, width, height, (void *)&md5, (void (*)(void *out, const byte_t *src, int size))md5_append);
				md5_finish(&md5, digest);
				bin2hex(digest, hex);
				fwrite(hex, 1, sizeof(hex), opt.outfilemd5_);
			}
		} else {
			UniDelay(100);
		}
		if (opt.interval_) {
			waitevents_with_timer(yuv, &rect, suspended);
		} else {
			waitevents(yuv, &rect, suspended);
		}
	}
endloop:
	int status;

	UniWaitThread(thr_m2d, &status);
	RecordTime(0);
	if (yuv) {
		SDL_FreeYUVOverlay(yuv);
		yuv = 0;
	}
	for (int i = 0; i < BUFNUM; ++i) {
		delete[] src[i].data;
	}
}

Uint32 DispTimer(Uint32 interval, void *param)
{
	SDL_Event ev;
	SDL_UserEvent user_ev;
	user_ev.type = SDL_USEREVENT;
	user_ev.code = 0;
	user_ev.data1 = 0;
	user_ev.data2 = 0;
	ev.type = SDL_USEREVENT;
	ev.user = user_ev;
	UniPushEvent(&ev);
	return interval;
}

#ifdef main
#undef main
#endif
#ifdef _M_IX86
#include <crtdbg.h>
#endif
#endif /* ENABLE_DISPLAY */

int main(int argc, char **argv)
{
#ifdef ENABLE_DISPLAY
	SDL_TimerID timer;

#ifdef _M_IX86
	_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_WNDW);
//	_CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_WNDW);
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF|_CRTDBG_LEAK_CHECK_DF);
#endif
	if (SDL_Init(SDL_INIT_TIMER | SDL_INIT_VIDEO) < 0) {
		return -1;
	}
	atexit(SDL_Quit);
	Options opt(argc, argv);
	LogInit();
	atexit(LogFin);
	if (opt.logdump_) {
		atexit(LogDump);
	}
#ifdef _M_IX86
	atexit((void (*)(void))_CrtCheckMemory);
#endif
	SDL_EventState(SDL_ACTIVEEVENT, SDL_IGNORE);
	SDL_EventState(SDL_MOUSEMOTION, SDL_IGNORE);
	if (opt.interval_) {
		timer = SDL_AddTimer(opt.interval_, DispTimer, 0);
	}
	do {
		run_loop(opt);
	} while (opt.repeat_);

	if (opt.interval_) {
		SDL_RemoveTimer(timer);
	}
	if (opt.outfile_) {
		fclose(opt.outfile_);
	}
	if (opt.outfilemd5_) {
		fclose(opt.outfilemd5_);
	}
#endif /* ENABLE_DISPLAY */
	return 0;
}

