#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <functional>
#include <vector>
#include <queue>
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

#define NELEM(a) (sizeof(a) / sizeof(a[0]))

class Uncopyable {
protected:
	Uncopyable() {}
	~Uncopyable() {}
private:
	Uncopyable(const Uncopyable&);
	Uncopyable& operator=(const Uncopyable&);
};

class Lock: private Uncopyable {
public:
	explicit Lock(UniMutex m) : mutex_(m) {
		UniLockMutex(m);
	}
	~Lock() {
		UniUnlockMutex(mutex_);
	}
private:
	UniMutex mutex_;
};

template <typename T, typename CONTAINER>
class AsyncContainer: private Uncopyable {
	typedef CONTAINER container_type;
	const typename container_type::size_type max_;
	UniMutex mutex_;
	UniCond cond_;
	container_type container_;
public:
	AsyncContainer(int max)
		: max_(static_cast<typename container_type::size_type>(max)),
		  mutex_(UniCreateMutex()),
		  cond_(UniCreateCond()) {}
	virtual ~AsyncContainer() {
		UniDestroyCond(cond_);
		UniDestroyMutex(mutex_);
	}
	void push_nolock(const T& data) {
		container_.push(data);
	}
	void push(const T& data) {
		Lock lk(mutex_);
		while (max_ <= container_.size()) {
			UniCondWait(cond_, mutex_);
		}
		container_.push(data);
		UniCondSignal(cond_);
	}
	T pop() {
		Lock lk(mutex_);
		while (container_.empty()) {
			UniCondWait(cond_, mutex_);
		}
		T result = container_.front();
		container_.pop();
		UniCondSignal(cond_);
		return result;
	}
	size_t size() const {
		return container_.size();
	}
	void clear() {
		while (!container_.empty()) {
			container_.pop();
		}
	}
};

template <typename T>
class AsyncQueue: public AsyncContainer<T, std::queue<T> > {
public:
	AsyncQueue(int max) : AsyncContainer<T, std::queue<T> >(max) {}
};

template <typename T> class WritableFrame;
template <typename T> class ReadableFrame;

template <typename T>
class QueuedBuffer : private Uncopyable {
	friend class WritableFrame<T>;
	friend class ReadableFrame<T>;
	int num_;
	int filled_num_;
	std::vector<T> frames_;
	AsyncQueue<T*> empty_frames_;
	AsyncQueue<T*> stuffed_frames_;
	T* reserveFrame() {
		return empty_frames_.pop();
	}
	void pushFrame(T* frm) {
		stuffed_frames_.push(frm);
	}
public:
	QueuedBuffer(int num)
		: num_(num), filled_num_(0),
		  frames_(num),
		  empty_frames_(num),
		  stuffed_frames_(num)
	{
		for (int i = 0; i < num_; ++i) {
			empty_frames_.push(&frames_[i]);
		}
	}
	~QueuedBuffer() {}
	T* popFrame() {
		return stuffed_frames_.pop();
	}
	void discardFrame(T* frm) {
		empty_frames_.push(frm);
	}
	void assign(std::vector<T> dat) {
		frames_.assign(dat.begin(), dat.end());
		drop();
		for (int i = 0; i < num_; ++i) {
			empty_frames_.push(&frames_[i]);
		}
	}
	void drop() {
		stuffed_frames_.clear();
		empty_frames_.clear();
	}
};

template <typename T>
class ReadableFrame: private Uncopyable {
	QueuedBuffer<T>& buffer_;
	T* frame_;
public:
	ReadableFrame(QueuedBuffer<T>& buf) : buffer_(buf), frame_(buffer_.popFrame()) {}
	~ReadableFrame() {
		buffer_.discardFrame(frame_);
	}
	const T* get() {
		return frame_;
	}
};

template <typename T>
class WritableFrame: private Uncopyable {
	QueuedBuffer<T>& buffer_;
	T* frame_;
public:
	WritableFrame(QueuedBuffer<T>& buf) : buffer_(buf), frame_(buffer_.reserveFrame()) {}
	~WritableFrame() {
		buffer_.pushFrame(frame_);
	}
	T* get() {
		return frame_;
	}
};

struct Buffer {
	uint8_t *data;
	int len;
	int type;
	void *id;
};

class BufferArray {
	std::vector<Buffer> buf_;
public:
	BufferArray(int num, int insize) : buf_(num) {
		if (insize <= 0) {
			return;
		}
		for (int i = 0; i < num; ++i) {
			buf_[i].data = new unsigned char[insize];
			buf_[i].len = insize;
		}
	}
	~BufferArray() {
		for (int i = 0; i < buf_.size(); ++i) {
			delete[] buf_[i].data;
		}
	}
	std::vector<Buffer>& get() {
		return buf_;
	}
};

class FileReader {
	std::list<const char *> infiles_;
	M2Decoder::type_t codec_;
	const char *filename_;
	FILE *fd_;
	int insize_;
	int file_open() {
		while (!infiles_.empty()) {
			filename_ = infiles_.front();
			codec_ = detect_file(filename_);
			fd_ = fopen(filename_, "rb");
			infiles_.pop_front();
			if (fd_) {
				return 0;
			} else {
				fprintf(stderr, "Error on %s.\n", filename_);
			}
		}
		fd_ = 0;
		codec_ = M2Decoder::MODE_NONE;
		return -1;
	}
public:
	FileReader(std::list<const char *> &infiles, int insize)
		: infiles_(infiles), codec_(M2Decoder::MODE_NONE), filename_(0), fd_(0), insize_(insize) {
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
			fd_ = 0;
		}
		dst.id = (void *)filename_;
		dst.type = (int)codec_;
		return read_size;
	}
	M2Decoder::type_t next_codec() {
		return codec_;
	}
};

class FileReaderUnit {
public:
	typedef QueuedBuffer<Buffer> QueueType;
	FileReaderUnit(int insize, int inbufnum, std::list<const char *> &infiles)
		: fr_(infiles, insize), dst_(inbufnum, insize), outqueue_(inbufnum) {
		outqueue_.assign(dst_.get());
	}
	~FileReaderUnit() {}
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
	BufferArray dst_;
	QueueType outqueue_;
	int run_impl() {
		int err;
		RecordTime(1);
		for (;;) {
			WritableFrame<Buffer> wr(outqueue_);
			Buffer* buf = wr.get();
			err = fr_.read_block(*buf);
			if (err < 0) {
				break;
			}
		}
		fprintf(stderr, "File terminate.\n");
		RecordTime(0);
		WritableFrame<Buffer> wr(outqueue_);
		wr.get()->len = -1;
		return 0;
	}
};

class M2DecoderUnit {
public:
	typedef QueuedBuffer<Frame> QueueType;
	M2DecoderUnit(M2Decoder::type_t codec, FileReaderUnit::QueueType& inqueue, int dstnum, bool dpb_emptify)
		: m2dec_(codec, dstnum, reread_file, this),
		  dst_align_(dstnum),
		  inqueue_(inqueue), outqueue_(dstnum),
		  src_(0), dpb_emptify_(dpb_emptify), terminated_(false) {
	}
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
	std::vector<Frame> dst_align_;
	FileReaderUnit::QueueType& inqueue_;
	M2DecoderUnit::QueueType outqueue_;
	Buffer* src_;
	bool dpb_emptify_;
	bool terminated_;
	static void post_dst(void *obj, Frame& frm) {
		((M2DecoderUnit *)obj)->post_dst_impl(frm);
	}
	static int reread_file(void *arg) {
		return ((M2DecoderUnit *)arg)->reread_file_impl();
	}
	void post_dst_impl(Frame& frm) {
		WritableFrame<Frame> wr(outqueue_);
		*wr.get() = frm;
	}
	int reread_file_impl() {
		if (src_) {
			inqueue_.discardFrame(src_);
		}
		src_ = inqueue_.popFrame();
		if (src_->len <= 0) {
			terminated_ = true;
			return -1;
		}
		M2Decoder::type_t codec = static_cast<M2Decoder::type_t>(src_->type);
		if (codec != m2dec_.codec_mode()) {
			m2dec_.change_codec(codec);
		}
		dec_bits *stream = dec().demuxer()->stream;
		stream = stream ? stream : dec().stream();
		dec_bits_set_data(stream, src_->data, src_->len, src_->id);
		return 0;
	}
	int run_impl() {
		RecordTime(1);
		do {
			while (0 <= dec().decode(this, post_dst, dpb_emptify_)) {}
			dec().decode_residual(this, post_dst);
		} while (!terminated_);
		static Frame nullframe = {0,};
		post_dst(this, nullframe);
		RecordTime(0);
		return 0;
	}
};

#include <emmintrin.h>

void deinterleave(const uint8_t *src, uint8_t *dst, int stride, int height) {
	uint8_t* dst1 = dst + (stride >> 1);
	__m128i msk = _mm_shuffle_epi32(_mm_cvtsi32_si128(0x00ff00ff), 0);
	int width = (stride + 1) >> 1;
	do {
		for (int x = 0; x < width; x += 16) {
			__m128i d0 = _mm_load_si128((__m128i const *)&src[x * 2]);
			__m128i d1 = _mm_load_si128((__m128i const *)&src[x * 2 + 16]);
			__m128i d0l = _mm_and_si128(d0, msk);
			__m128i d1l = _mm_and_si128(d1, msk);
			d0 = _mm_srli_epi16(d0, 8);
			d1 = _mm_srli_epi16(d1, 8);
			d0l = _mm_packus_epi16(d0l, d1l);
			d0 = _mm_packus_epi16(d0, d1);
			_mm_storeu_si128((__m128i *)&dst[x], d0l);
			_mm_storeu_si128((__m128i *)&dst1[x], d0);
		}
		src += stride;
		dst += stride;
		dst1 += stride;
	} while (--height);
}

const int FILE_READ_SIZE = 65536 * 2;
const int IBUFNUM = 3;

static void BlameUser() {
	fprintf(stderr,
		"Usage: srview [-m] [-o] [-r] [-t interval] infile [infile ...]\n"
		"\t-m : outfile(MD5)\n"
		"\t-o : outfile(Raw)\n"
		"\t-r : repeat\n"
		"\t-l : log dump\n"
		"\t-f frame_num(3-256) : specify number of frames before display.\n"
		"\t-e : DPB emptify mode\n"
		"\t-t interval : specify interval of each frame in ms unit\n");
	exit(1);
}

struct DeleteFw {
	void operator()(FileWriter *fw) const {
		delete fw;
	}
};

struct WriteFrame : std::binary_function<FileWriter *, const Frame *, void> {
	void operator()(FileWriter *fw, const Frame *out) const {
		fw->writeframe(out);
	}
};

class FileWriterUnit {
	typedef std::list<FileWriter *> outfiles_t;
	outfiles_t fw_;
	std::list<const char *> infiles_;
	void *id_;
	void change() {
		if (infiles_.empty()) {
			fw_.clear();
		} else {
			for (outfiles_t::iterator p = fw_.begin(); p != fw_.end(); ++p) {
				(*p)->set_file(infiles_.front(), true);
			}
			infiles_.pop_front();
		}
	}
public:
	FileWriterUnit() : id_(0) {}
	~FileWriterUnit() {
		for_each(fw_.begin(), fw_.end(), DeleteFw());
	}
	void set_mode(FileWriter::write_type t) {
		FileWriter *fw;
		switch (t) {
		case FileWriter::WRITE_RAW:
			fw = new FileWriterRaw();
			break;
		case FileWriter::WRITE_MD5:
			fw = new FileWriterMd5();
			break;
		}
		fw_.push_back(fw);
	}
	void set_filelist(const std::list<const char *>& infiles) {
		infiles_ = infiles;
	}
	bool empty() const {
		return fw_.empty();
	}
	void write(const Frame& out) {
		if (out.id != id_) {
			id_ = out.id;
			change();
		}
		for_each(fw_.begin(), fw_.end(), std::bind2nd(WriteFrame(), &out));
	}
};

struct Options {
	int interval_;
	std::list<const char *> infile_list_;
	FileWriterUnit fw_;
	int outbuf_;
	int repeat_;
	bool logdump_;
	bool dpb_mode_;
	~Options() {}
	Options(int argc, char **argv)
		: interval_(0),
		  outbuf_(3), repeat_(1), logdump_(false), dpb_mode_(false) {
		int opt;
		while ((opt = getopt(argc, argv, "ef:hlmor:st:")) != -1) {
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
				fw_.set_mode(FileWriter::WRITE_MD5);
				break;
			case 'o':
				fw_.set_mode(FileWriter::WRITE_RAW);
				break;
			case 'r':
				repeat_ = strtoul(optarg, 0, 0);
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
		if (!fw_.empty() && !infile_list_.empty()) {
			fw_.set_filelist(infile_list_);
		}
	}
};

struct UniSurface {
#ifdef ENABLE_DISPLAY
	SDL_Window *screen_;
	SDL_Renderer *renderer_;
	SDL_Texture *texture_;
	std::vector<uint8_t> cbcr_base_;
	uint8_t* cbcr_;
	int width_, height_;
#endif /* ENABLE_DISPLAY */
	int suspended_;
	void *id_;
	UniSurface() :
#ifdef ENABLE_DISPLAY
		screen_(0), renderer_(0), texture_(0),
		width_(0), height_(0),
#endif /* ENABLE_DISPLAY */
		suspended_(0), id_(0) {
#ifdef ENABLE_DISPLAY
		UniEvent ev;
		while (UniPollEvent(&ev)) ;
#endif /* ENABLE_DISPLAY */
	}
	~UniSurface() {
		teardown();
	}
	int suspended() const {
		return suspended_;
	}
	void display(const Frame& out) {
#ifdef ENABLE_DISPLAY
		if (out.id != id_) {
			id_ = out.id;
			change(out);
		}
		deinterleave(out.chroma, &cbcr_[0], out.width, (out.height - out.crop[3]) >> 1);
		SDL_UpdateYUVTexture(texture_, 0, out.luma, out.width, &cbcr_[0], out.width, &cbcr_[out.width >> 1], out.width);
		SDL_RenderClear(renderer_);
		SDL_RenderCopy(renderer_, texture_, NULL, NULL);
#endif /* ENABLE_DISPLAY */
	}
	void waitevents(int interval) {
#ifdef ENABLE_DISPLAY
		if (interval) {
			waitevents_with_timer();
		} else {
			waitevents_without_timer();
		}
#endif /* ENABLE_DISPLAY */
	}
private:
	void teardown() {
#ifdef ENABLE_DISPLAY
		if (texture_) {
			SDL_DestroyTexture(texture_);
		}
		if (renderer_) {
			SDL_DestroyRenderer(renderer_);
		}
		if (screen_) {
			SDL_DestroyWindow(screen_);
		}
#endif /* ENABLE_DISPLAY */
	}
#ifdef ENABLE_DISPLAY
	void waitevents_with_timer() {
		UniEvent ev;
		UniWaitEvent(&ev);
		do {
			switch (ev.type) {
			case SDL_USEREVENT:
				SDL_RenderPresent(renderer_);
				break;
			case SDL_MOUSEBUTTONDOWN:
				suspended_ ^= 1;
				break;
			case SDL_QUIT:
				exit(0);
				/* NOTREACHED */
				break;
			}
		} while (SDL_PollEvent(&ev));
	}
	void waitevents_without_timer() {
		UniEvent ev;
		SDL_RenderPresent(renderer_);
		while (SDL_PollEvent(&ev)) {
			switch (ev.type) {
			case SDL_MOUSEBUTTONDOWN:
				suspended_ ^= 1;
				break;
			case SDL_QUIT:
				exit(0);
				/* NOTREACHED */
				break;
			}
		}
	}
	static void adjust_windowsize(int width, int height, int& dwidth, int& dheight) {
		while ((1440 < width) || (720 < height)) {
			width >>= 1;
			height >>= 1;
		}
		dwidth = width;
		dheight = height;
	}
	void change(const Frame& out) {
		if ((width_ != (out.width - out.crop[1])) || (height_ != (out.height - out.crop[3]))) {
			static Uint32 winflg_cand[] = {
				SDL_WINDOW_OPENGL, 0
			};
			teardown();
			width_ = out.width - out.crop[1];
			height_ = out.height - out.crop[3];
			cbcr_base_.resize(((width_ * height_) >> 1) + 15);
			cbcr_ = reinterpret_cast<uint8_t*>(ALIGN16(reinterpret_cast<uintptr_t>(&cbcr_base_[0])));
			int dwidth, dheight;
			adjust_windowsize(width_, height_, dwidth, dheight);
			for (int i = 0; i < NELEM(winflg_cand); ++i) {
				screen_ = SDL_CreateWindow("disp", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, dwidth, dheight, winflg_cand[i]);
				if (screen_) {
					renderer_ = SDL_CreateRenderer(screen_, -1, 0);
					texture_ = SDL_CreateTexture(renderer_, SDL_PIXELFORMAT_IYUV, SDL_TEXTUREACCESS_STREAMING, width_, height_);
					break;
				}
			}
		}
	}
#endif /* ENABLE_DISPLAY */
};

int run_loop(Options& opt, int outbuf) {
	LogTags.insert(std::pair<int, const char *> (UniThreadID(), "Main"));
	RecordTime(1);

	FileReaderUnit fr(FILE_READ_SIZE, IBUFNUM, opt.infile_list_);
	if (fr.next_codec() == M2Decoder::MODE_NONE) {
		return -1;
	}
 	M2DecoderUnit m2dec(fr.next_codec(), fr.outqueue(), outbuf, opt.dpb_mode_);

	UniThread thr_file = UniCreateThread(FileReaderUnit::run, (void *)&fr, "FileLoader");
	UniThread thr_m2d = UniCreateThread(M2DecoderUnit::run, (void *)&m2dec, "Decoder");

	UniSurface surface;
	while (1) {
		if (!surface.suspended()) {
			ReadableFrame<Frame> rd(m2dec.outqueue());
			const Frame* out = rd.get();
			if (!out->luma) {
				break;
			}
			surface.display(*out);
			opt.fw_.write(*out);
		} else {
			UniDelay(100);
		}
		surface.waitevents(opt.interval_);
	}
	int status;
	UniWaitThread(thr_m2d, &status);
	RecordTime(0);
	return 0;
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

void logOut(void* p, int category, SDL_LogPriority priority, const char* message) {
	fprintf(stderr, "%s", message);
}

#ifdef main
#undef main
#endif
#ifdef _M_IX86
#include <crtdbg.h>
#endif

int main(int argc, char **argv)
{
#if defined(_MSC_VER) && !defined(NDEBUG)
	_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_WNDW);
//	_CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_WNDW);
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF|_CRTDBG_LEAK_CHECK_DF);
	atexit((void (*)(void))_CrtCheckMemory);
#endif
	SDL_TimerID timer;
	if (SDL_Init(SDL_INIT_TIMER | SDL_INIT_EVENTS
#ifdef ENABLE_DISPLAY
			| SDL_INIT_VIDEO
#endif /* ENABLE_DISPLAY */
			) < 0) {
		return -1;
	}
	atexit(SDL_Quit);
	Options opt(argc, argv);
	LogInit();
	atexit(LogFin);
	if (opt.logdump_) {
		atexit(LogDump);
	}
#ifdef ENABLE_DISPLAY
	SDL_EventState(SDL_WINDOWEVENT, SDL_IGNORE);
	SDL_EventState(SDL_MOUSEMOTION, SDL_IGNORE);
	SDL_LogSetOutputFunction(logOut, 0);
	if (opt.interval_) {
		timer = SDL_AddTimer(opt.interval_, DispTimer, 0);
	}
#endif /* ENABLE_DISPLAY */
	opt.repeat_ = std::max(1, opt.repeat_);
	do {
		if (run_loop(opt, opt.outbuf_) < 0) {
			break;
		}
	} while (--opt.repeat_);
#ifdef ENABLE_DISPLAY
	if (opt.interval_) {
		SDL_RemoveTimer(timer);
	}
#endif /* ENABLE_DISPLAY */
	return 0;
}

