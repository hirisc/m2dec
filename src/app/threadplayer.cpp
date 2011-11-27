#ifdef ENABLE_DISPLAY
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <functional>
#include <vector>
#include <list>
#include <map>
#include "frames.h"
#include "m2decoder.h"
#include "filewrite.h"
#include "getopt.h"
#include "md5.h"
#ifdef __GNUC__
#include <stdint.h>
#endif

#include "unithread.h"

template <typename T>
class Queue {
	T *data_;
	int max_;
	int size_;
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
		: data_(data), max_(buf_num), size_(1), head_(0), tail_(buf_num - 1),
		  mutex_(UniCreateMutex()),
		  cond_(UniCreateCond()),
		  terminated_(false) {}
	~Queue() {
		UniDestroyCond(cond_);
		UniDestroyMutex(mutex_);
	}
	bool full() const {
		return (head_ == tail_) && (size_ != 0);
	}
	bool empty() const {
		return (head_ == tail_) && (size_ == 0);
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
	T& front() {
		UniLockMutex(mutex_);
		while (full()) {
			UniCondWait(cond_, mutex_);
		}
		UniUnlockMutex(mutex_);
		return data_[head_];
	}
	void push_front(const T& dat) {
		UniLockMutex(mutex_);
		while (full()) {
			UniCondWait(cond_, mutex_);
		}
		size_++;
		data_[head_] = dat;
		head_ = next(head_);
		UniCondSignal(cond_);
		UniUnlockMutex(mutex_);
	}
	T& back() {
		UniLockMutex(mutex_);
		while (empty() && !terminated()) {
			UniCondWait(cond_, mutex_);
		}
		UniUnlockMutex(mutex_);
		return data_[tail_];
	}
	bool pop_back() {
		UniLockMutex(mutex_);
		if (empty() && terminated()) {
			UniUnlockMutex(mutex_);
			return false;
		}
		assert(!empty());
		size_--;
		tail_ = next(tail_);
		UniCondSignal(cond_);
		UniUnlockMutex(mutex_);
		return true;
	}
};

struct Buffer {
	uint8_t *data;
	int len;
	int type;
};

class FileReader {
	std::list<char *> infiles_;
	M2Decoder::type_t codec_;
	FILE *fd_;
	int insize_;
	int file_open() {
		while (!infiles_.empty()) {
			const char *filename = infiles_.front();
			codec_ = detect_file(filename);
			fd_ = fopen(filename, "rb");
			infiles_.pop_front();
			if (fd_) {
				return 0;
			} else {
				fprintf(stderr, "Error on %s.\n", filename);
			}
		}
		fd_ = 0;
		return -1;
	}
public:
	FileReader(std::list<char *> &infiles, int insize)
		: infiles_(infiles), codec_(M2Decoder::MODE_NONE), fd_(0), insize_(insize) {
		file_open();
	}
	~FileReader() {
		if (fd_) {
			fclose(fd_);
		}
	}
	int read_block(Buffer& dst) {
		if ((fd_ == 0) && (file_open() < 0)) {
			return -1;
		}
		int read_size = fread(dst.data, 1, insize_, fd_);
		dst.len = read_size;
		if (read_size == 0) {
			fclose(fd_);
			if (file_open() < 0) {
				return -1;
			}
			dst.len = read_size = fread(dst.data, 1, insize_, fd_);
		}
		dst.type = (int)codec_;
		return read_size;
	}
	M2Decoder::type_t next_codec() {
		return codec_;
	}
	static void to_lower_ext(const char *src, char *dst, int len) {
		for (int i = 0; i < len; ++i) {
			const char s = src[i];
			if (s == '\0') {
				dst[i] = s;
				break;
			}
			dst[i] = tolower(s);
		}
	}
	static M2Decoder::type_t detect_file(const char *filename) {
		static const struct {
			M2Decoder::type_t type;
			const char *ext;
		} ext_map[] = {
			{ M2Decoder::MODE_MPEG2, "m2v"},
			{ M2Decoder::MODE_MPEG2PS, "vob"},
			{ M2Decoder::MODE_H264, "264"},
			{ M2Decoder::MODE_H264, "jsv"},
			{ M2Decoder::MODE_NONE, ""}
		};
		char ext[16];
		const char *ext_p = strrchr(filename, '.');
		if (ext_p++ != 0) {
			to_lower_ext(ext_p, ext, sizeof(ext));
			int i = -1;
			while (ext_map[++i].type != M2Decoder::MODE_NONE) {
				if (strcmp(ext_map[i].ext, ext) == 0) {
					return ext_map[i].type;
				}
			}
		}
		return M2Decoder::MODE_MPEG2;
	}
};

class FileReaderUnit {
public:
	typedef Queue<Buffer> QueueType;
	FileReaderUnit(int insize, std::list<char *> &infiles, FileReaderUnit::QueueType& outqueue)
		: fr_(infiles, insize), outqueue_(outqueue) {}
	FileReaderUnit::QueueType& outqueue() {
		return outqueue_;
	}
	M2Decoder::type_t next_codec() {
		return fr_.next_codec();
	}
	static int run(void *data) {
		return ((FileReaderUnit *)data)->run_impl();
	}
private:
	FileReader fr_;
	QueueType& outqueue_;
	int run_impl() {
		int err;
		RecordTime(1);
		for (;;) {
			Buffer& buf = outqueue_.front();
			err = fr_.read_block(buf);
			if (err < 0) {
				break;
			}
			outqueue_.push_front(buf);
		}
		outqueue_.terminate();
		fprintf(stderr, "File terminate.\n");
		RecordTime(0);
		return 0;
	}
};

class M2DecoderUnit {
public:
	typedef Queue<Frame> QueueType;
	M2DecoderUnit(M2Decoder::type_t codec, FileReaderUnit::QueueType& inqueue, Frame *dst, int dstnum, bool dpb_emptify)
		: m2dec_(codec, dstnum, reread_file, this), dstnum_(dstnum), inqueue_(inqueue), outqueue_(QueueType(dst, dstnum)),
		  src_next_(0),
		  dpb_emptify_(dpb_emptify) {}
	~M2DecoderUnit() {}
	FileReaderUnit::QueueType& inqueue() {
		return inqueue_;
	}
	QueueType& outqueue() {
		return outqueue_;
	}
	M2Decoder& dec() {
		return m2dec_;
	}
	static int run(void *data) {
		return ((M2DecoderUnit *)data)->run_impl();
	}
private:
	M2Decoder m2dec_;
	int dstnum_;
	FileReaderUnit::QueueType& inqueue_;
	M2DecoderUnit::QueueType outqueue_;
	Buffer *src_next_;
	bool dpb_emptify_;
	static void post_dst(void *obj, Frame& frm) {
		M2DecoderUnit *ths = (M2DecoderUnit *)obj;
		ths->outqueue().push_front(frm);
	}
	static int reread_file(void *arg) {
		return ((M2DecoderUnit *)arg)->reread_file_impl();
	}
	int reread_file_impl() {
		Buffer *src;
		if (!src_next_) {
			if (!inqueue().pop_back()) {
				return -1;
			}
			src = &inqueue().back();
			M2Decoder::type_t codec = static_cast<M2Decoder::type_t>(src->type);
			if (codec != m2dec_.codec_mode()) {
				m2dec_.change_codec(codec);
				src_next_ = src;
				return -1;
			}
		} else {
			src = src_next_;
			src_next_ = 0;
		}
		if (src->len <= 0) {
			return -1;
		} else {
			dec_bits *stream = dec().demuxer()->stream;
			stream = stream ? stream : dec().stream();
			dec_bits_set_data(stream, src->data, src->len);
			return 0;
		}
	}
	int run_impl() {
		RecordTime(1);
		do {
			while (0 <= dec().decode(this, post_dst, dpb_emptify_)) {}
			dec().decode_residual(this, post_dst);
		} while (src_next_);
		static Frame nullframe = {0,};
		post_dst(this, nullframe);
		outqueue().terminate();
		RecordTime(0);
		return 0;
	}
};

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
const int FILE_READ_SIZE = 65536 * 2;
const int IBUFNUM = 3;

static void BlameUser() {
	fprintf(stderr,
		"Usage: srview [-r] [-t interval] [-m outfile(MD5)] [-o outfile(Raw)] infile [infile ...]\n"
		"\t-r : repeat\n"
		"\t-l : log dump\n"
		"\t-f frame_num(3-256) : specify number of frames before display.\n"
		"\t-e : DPB emptify mode\n"
		"\t-t interval : specify interval of each frame in ms unit\n");
	exit(-1);
}

struct Options {
	int interval_;
	std::list<char *> infile_list_;
	std::list<FileWriter *> fw_;
	int outbuf_;
	bool repeat_;
	bool logdump_;
	bool dpb_mode_;

	Options(int argc, char **argv)
		: interval_(0),
		  outbuf_(3), repeat_(false), logdump_(false), dpb_mode_(false) {
		int opt;
		while ((opt = getopt(argc, argv, "ef:hlm:o:rst:")) != -1) {
			FILE *fo;
			switch (opt) {
			case 'e':
				dpb_mode_ = true;
				break;
			case 'f':
				outbuf_ = strtoul(optarg, 0, 0);
				if (253U < (unsigned)(outbuf_ - 3)) {
					BlameUser();
				}
				break;
			case 'l':
				logdump_ = true;
				break;
			case 'm':
				fo = fopen(optarg, "wb");
				if (fo) {
					fw_.push_back(new FileWriterMd5(fo));
				}
				break;
			case 'o':
				fo = fopen(optarg, "wb");
				if (fo) {
					fw_.push_back(new FileWriterRaw(fo));
				}
				break;
			case 'r':
				repeat_ = true;
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
		if (argc <= optind) {
			BlameUser();
			/* NOTREACHED */
		}
		do {
			infile_list_.push_back(argv[optind]);
		} while (++optind < argc);
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

struct WriteFrame : std::binary_function<FileWriter *, Frame *, void> {
	void operator()(FileWriter *fw, const Frame *out) const {
		fw->writeframe(out);
	}
};

struct DeleteFw {
	void operator()(FileWriter *fw) const {
		delete fw;
	}
};

void run_loop(Options& opt, int outbuf) {
	int width = 0;
	int height = 0;
	SDL_Surface *surface;
	SDL_Overlay *yuv = 0;
	SDL_Rect rect = {
		0, 0, width, height
	};

	LogTags.insert(std::pair<int, const char *> (UniThreadID(), "Main"));
	RecordTime(1);

	Buffer src[IBUFNUM];
	for (int i = 0; i < IBUFNUM; ++i) {
		src[i].data = new unsigned char[FILE_READ_SIZE];
		src[i].len = FILE_READ_SIZE;
	}
	FileReaderUnit::QueueType fr_queue(src, IBUFNUM);
	FileReaderUnit fr(FILE_READ_SIZE, opt.infile_list_, fr_queue);

	std::vector<Frame> dst_align(outbuf);
 	M2DecoderUnit m2dec(fr.next_codec(), fr.outqueue(), &dst_align[0], outbuf, opt.dpb_mode_);

	UniThread *thr_file = UniCreateThread(FileReaderUnit::run, (void *)&fr);
	LogTags.insert(std::pair<int, const char *> (UniGetThreadID(thr_file), "FileLoader"));
	UniThread *thr_m2d = UniCreateThread(M2DecoderUnit::run, (void *)&m2dec);
	LogTags.insert(std::pair<int, const char *> (UniGetThreadID(thr_m2d), "Decoder"));

	UniEvent ev;
	while (UniPollEvent(&ev)) ;
	int suspended = 0;
	M2DecoderUnit::QueueType &outqueue = m2dec.outqueue();
	while (1) {
		if (!suspended) {
			if (!outqueue.pop_back()) {
				goto endloop;
			}
			Frame& out = outqueue.back();
			if (out.luma == 0) {
				goto endloop;
			}
			if ((width != (out.width - out.crop[1])) || (height != (out.height - out.crop[3]))) {
				width = out.width - out.crop[1];
				height = out.height - out.crop[3];
				rect.w = width;
				rect.h = height;
				surface = SDL_SetVideoMode(width, height, 0, SDL_SWSURFACE);
				if (yuv) {
					SDL_FreeYUVOverlay(yuv);
				}
				yuv = SDL_CreateYUVOverlay(width, height, SDL_IYUV_OVERLAY, surface);
			}
			SDL_LockYUVOverlay(yuv);
			display_write(yuv->pixels, out.luma, out.chroma, out.width, yuv->pitches, yuv->w, yuv->h);
			SDL_UnlockYUVOverlay(yuv);
			for_each(opt.fw_.begin(), opt.fw_.end(), std::bind2nd(WriteFrame(), &out));
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
	for_each(opt.fw_.begin(), opt.fw_.end(), DeleteFw());
	if (yuv) {
		SDL_FreeYUVOverlay(yuv);
		yuv = 0;
	}
	for (int i = 0; i < IBUFNUM; ++i) {
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

#if defined(_M_IX86) && !defined(NDEBUG)
	_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_WNDW);
//	_CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_WNDW);
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF|_CRTDBG_LEAK_CHECK_DF);
	atexit((void (*)(void))_CrtCheckMemory);
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
	SDL_EventState(SDL_ACTIVEEVENT, SDL_IGNORE);
	SDL_EventState(SDL_MOUSEMOTION, SDL_IGNORE);
	if (opt.interval_) {
		timer = SDL_AddTimer(opt.interval_, DispTimer, 0);
	}
	do {
		run_loop(opt, opt.outbuf_);
	} while (opt.repeat_);

	if (opt.interval_) {
		SDL_RemoveTimer(timer);
	}
#endif /* ENABLE_DISPLAY */
	return 0;
}

