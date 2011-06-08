#ifdef ENABLE_DISPLAY
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <list>
#include <map>
#include "getopt.h"
#include "bitio.h"
#include "mpeg2.h"
#include "mpeg_demux.h"
#ifdef __GNUC__
#include <stdint.h>
#endif

#ifndef __RENESAS_VERSION__
using namespace std;
#endif

#include "unithread.h"

struct Buffer {
	uint8_t *data, *data2;
	int width, height;
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
	UniMutex *mutex() const {
		return mutex_;
	}
	UniCond *cond() const {
		return cond_;
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
		while (!terminated() && empty()) {
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

typedef Queue<Buffer> QueueType;

class FileReader {
	QueueType outqueue_;
	list<char *> infiles_;
public:
	FILE *fd_;
	int insize_;
	FileReader(Buffer *src_p, int buf_num, FILE *fd, int insize, list<char *> &infiles)
		: outqueue_(QueueType(src_p, buf_num)), infiles_(infiles), fd_(fd), insize_(insize) {}
	QueueType& outqueue() {
		return outqueue_;
	}
	list<char *> &infiles() {
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

static int reread_file_ps(void *arg);

class PesDemuxer {
	QueueType outqueue_;
	list<char *> infiles_;
	pes_demuxer_t demux_;
public:
	FILE *fd_;
	PesDemuxer(Buffer *src_p, int buf_num, FILE *fd, list<char *> &infiles)
		: outqueue_(QueueType(src_p, buf_num)), infiles_(infiles), fd_(fd) {
		//		mpeg_demux_init(&demux_, reread_file_ps, this, jmp_buf *jmp);
	}
	QueueType& outqueue() {
		return outqueue_;
	}
	list<char *> &infiles() {
		return infiles_;
	}
	pes_demuxer_t *demuxer() {
		return &demux_;
	}
	static int run(void *data) {
		RecordTime(1);
		PesDemuxer *fr = (PesDemuxer *)data;
		while (1) {
			int read_size;
			int ret;
			Buffer& buf = fr->outqueue().emptybuf();
//			int read_size = fread(buf.data, 1, fr->insize_, fr->fd_);
			ret = 0;
			if ((fread(&read_size, 1, 4, fr->fd_) == 4)
			    && ((ret = fread(buf.data, 1, read_size, fr->fd_)) == read_size)) {
				buf.len = read_size;
			} else {
#if 0
				fclose(fr->fd_);
				if (!fr->infiles().empty() && (fr->fd_ = fopen(fr->infiles().front(), "rb"))) {
					fr->infiles().pop_front();
					buf.len = fread(buf.data, 1, fr->insize_, fr->fd_);
				} else {
					fr->outqueue().terminate();
					RecordTime(0);
					return 0;
				}
#endif
			}
			fr->outqueue().setfilled(buf);
		}
		return -1;
	}
};

static int reread_file_ps(void *arg)
{
	PesDemuxer *module = (PesDemuxer *)arg;
	Buffer& dst = module->outqueue().emptybuf();
	dec_bits *stream = module->demuxer()->stream;
	#if 0
	if (src.len <= 0) {
		return -1;
	} else {
		dec_bits_set_data(stream, src.data, src.len);
		return 0;
	}
	#endif
}

static int reread_file(void *arg);
static int reread_packet(void *arg);

class M2Decoder {
	m2d_frame *frame_;
	uint8_t **frame_org;
	QueueType& inqueue_;
	QueueType outqueue_;
	jmp_buf jmp_;
	jmp_buf jmp2_;
	m2d_context m2d_;
	pes_demuxer_t demux_;
public:
	M2Decoder(QueueType& inqueue, Buffer *dst, int dstnum, int bufnum,
		int width, int height, int pes_mode)
		: inqueue_(inqueue), outqueue_(QueueType(dst, dstnum)) {
		int len = ((width + 15) & ~15) * ((height + 15) & ~15);
		frame_ = new m2d_frame[bufnum];
		frame_org = new uint8_t *[bufnum];
		for (int i = 0; i < bufnum; ++i) {
			frame_org[i] = new uint8_t[((len * 3) >> 1) + 31];
			frame_[i].luma = (uint8_t *)((((intptr_t)frame_org[i]) + 15) & ~15);
			frame_[i].chroma = frame_[i].luma + len;
		}
		m2d_init(&m2d_, bufnum, frame_);
		if (pes_mode) {
			mpeg_demux_init(&demux_, reread_file, this, &jmp_);
			dec_bits_set_callback(m2d_.stream, reread_packet, this, &jmp2_);
		} else {
			memset(&demux_, 0, sizeof(demux_));
			dec_bits_set_callback(m2d_.stream, reread_file, this, &jmp_);
		}
	}
	~M2Decoder() {
		for (int i = 0; i < outqueue_.size(); ++i) {
			delete[] frame_org[i];
		}
		delete[] frame_org;
		delete[] frame_;
	}
	QueueType& inqueue() {
		return inqueue_;
	}
	QueueType& outqueue() {
		return outqueue_;
	}
	m2d_context *context() {
		return &m2d_;
	}
	dec_bits *stream() {
		return m2d_.stream;
	}
	pes_demuxer_t *demuxer() {
		return &demux_;
	}
	static int run(void *data) {
		RecordTime(1);
		M2Decoder *module = (M2Decoder *)data;
		m2d_context *m2d = module->context();
		int first = 1;
		while (1) {
			int width, height;
			const m2d_frame *frm;
			if (module->inqueue().terminated() && module->inqueue().empty()) {
				module->outqueue().terminate();
				RecordTime(0);
				return 0;
			}
			m2d_decode_data(m2d);
			if (first) {
				first = 0;
				continue;
			}
			frm = m2d_get_decoded_frame(m2d, &width, &height, 0);
			Buffer& dst = module->outqueue().emptybuf();
			dst.width = (width + 15) & ~15;
			dst.height = (height + 15) & ~15;
			dst.data = frm->luma;
			dst.data2 = frm->chroma;
			module->outqueue().setfilled(dst);
		}
		return -1;
	}
};

static int reread_file(void *arg)
{
	M2Decoder *module = (M2Decoder *)arg;
	Buffer& src = module->inqueue().getfilled();
	dec_bits *stream = module->demuxer()->stream;
	stream = stream ? stream : module->stream();
	if (src.len <= 0) {
		return -1;
	} else {
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
	if (packet) {
		dec_bits_set_data(module->context()->stream, packet, (size_t)packet_size);
		return 0;
	} else {
		return -1;
	}
}

void display_write(uint8_t **dst, const uint8_t *src_luma, const uint8_t *src_chroma,
		   uint16_t *pitches, int width, int height)
{
	uint8_t *dst0 = dst[0];
	for (int i = 0; i < height; ++i) {
		memcpy(dst0, src_luma, width);
		src_luma += width;
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
		src_chroma += width * 2;
		dst1 += pitches[1];
		dst2 += pitches[2];
	}
}

const int MAX_WIDTH = 1920;
const int MAX_HEIGHT = 1088;
const int MAX_LEN = (MAX_WIDTH * MAX_HEIGHT * 3) >> 1;
const int FILE_READ_SIZE = 65536 * 7;
const int BUFNUM = 3;

static void BlameUser() {
	fprintf(stderr,
		"Usage: srview [-s] [-r] [-t interval] [-o outfile] infile [infile ...]\n"
		"\t-s : Program Stream (PS)\n"
		"\t-r : repeat\n"
		"\t-l : log dump\n"
		"\t-t interval : specify interval of each frame in ms unit\n");
	exit(-1);
}

struct Options {
	int interval_;
	list<char *> infile_list_;
	FILE *infile_, *outfile_;
	int pes_mode_;
	bool repeat_;
	bool logdump_;

	Options(int argc, char **argv)
		: interval_(0), infile_(0), outfile_(0),
		  pes_mode_(0), repeat_(false), logdump_(false) {
		int opt;
		while ((opt = getopt(argc, argv, "lo:rst:")) != -1) {
			switch (opt) {
			case 'l':
				logdump_ = true;
				break;
			case 'o':
				outfile_ = fopen(optarg, "wb");
				break;
			case 'r':
				repeat_ = true;
				break;
			case 's':
				pes_mode_ = 1;
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

void run_loop(Options& opt) {
	int width = 0;
	int height = 0;
	SDL_Surface *surface;
	SDL_Overlay *yuv = 0;
	SDL_Rect rect = {
		0, 0, width, height
	};

	LogTags.insert(pair<int, const char *> (UniThreadID(), "Main"));
	RecordTime(1);

	/* Run File Loader */
	Buffer src[BUFNUM];
	for (int i = 0; i < BUFNUM; ++i) {
		src[i].data = new unsigned char[FILE_READ_SIZE];
		src[i].len = FILE_READ_SIZE;
	}
	FileReader fr(src, BUFNUM, opt.infile_, FILE_READ_SIZE, opt.infile_list_);
	UniThread *thr_file = UniCreateThread(FileReader::run, (void *)&fr);
	LogTags.insert(pair<int, const char *> (UniGetThreadID(thr_file), "FileLoader"));

	/* Run Mpeg-2 Decoder */
	Buffer dst_align[BUFNUM];
	M2Decoder m2dec(fr.outqueue(), dst_align, BUFNUM, BUFNUM + 2, MAX_WIDTH, MAX_HEIGHT, opt.pes_mode_);
	UniThread *thr_m2d = UniCreateThread(M2Decoder::run, (void *)&m2dec);
	QueueType &outqueue = m2dec.outqueue();
	LogTags.insert(pair<int, const char *> (UniGetThreadID(thr_m2d), "Mpeg2Dec"));

	UniEvent ev;
	while (UniPollEvent(&ev)) ;
	int suspended = 0;
	while (1) {
		if (!suspended) {
			if (outqueue.terminated() && outqueue.empty()) {
				goto endloop;
			}
			Buffer& out = outqueue.getfilled();
			if ((width != out.width) || (height != out.height)) {
				width = out.width;
				height = out.height;
				rect.w = width;
				rect.h = height;
				fprintf(stderr, "size: %d x %d\n", width, height);
				surface = SDL_SetVideoMode(width, height, 0, SDL_SWSURFACE);
				if (yuv) {
					SDL_FreeYUVOverlay(yuv);
				}
				yuv = SDL_CreateYUVOverlay(width, height, SDL_IYUV_OVERLAY, surface);
				fprintf(stderr, "stride: %d\n", yuv->pitches[0]);
			}
			SDL_LockYUVOverlay(yuv);
			display_write(yuv->pixels, out.data, out.data2, yuv->pitches, yuv->w, yuv->h);
			SDL_UnlockYUVOverlay(yuv);
			int lumasize = width * height;
			if (opt.outfile_) {
				fwrite(out.data, 1, lumasize, opt.outfile_);
				fwrite(out.data2, 1, lumasize >> 1, opt.outfile_);
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
	for (int i = 0; i < BUFNUM; ++i) {
		delete[] src[i].data;
	}
	fseek(opt.infile_, 0, SEEK_SET);
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
	if (SDL_Init(SDL_INIT_EVERYTHING) < 0) {
		return -1;
	}
	atexit(SDL_Quit);
	Options opt(argc, argv);
	LogInit();
	atexit(LogFin);
	if (opt.logdump_) {
		atexit(LogDump);
	}
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
#endif /* ENABLE_DISPLAY */
	return 0;
}

