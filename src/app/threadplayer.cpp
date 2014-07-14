#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
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
	explicit Lock(UniMutex *m) : mutex_(m) {
		UniLockMutex(m);
	}
	~Lock() {
		UniUnlockMutex(mutex_);
	}
private:
	UniMutex *mutex_;
};

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
		Lock lk(mutex_);
		terminated_ = true;
		UniCondSignal(cond_);
	}
	T& front() {
		Lock lk(mutex_);
		while (full()) {
			UniCondWait(cond_, mutex_);
		}
		return data_[head_];
	}
	void push_front(const T& dat) {
		Lock lk(mutex_);
		while (full()) {
			UniCondWait(cond_, mutex_);
		}
		size_++;
		data_[head_] = dat;
		head_ = next(head_);
		UniCondSignal(cond_);
	}
	T& back() {
		Lock lk(mutex_);
		while (empty() && !terminated()) {
			UniCondWait(cond_, mutex_);
		}
		return data_[tail_];
	}
	bool pop_back() {
		Lock lk(mutex_);
		if (empty() && terminated()) {
			UniUnlockMutex(mutex_);
			return false;
		}
		assert(!empty());
		size_--;
		tail_ = next(tail_);
		UniCondSignal(cond_);
		return true;
	}
};

struct Buffer {
	uint8_t *data;
	int len;
	int type;
	void *id;
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
			if (file_open() < 0) {
				return -1;
			}
			dst.len = read_size = fread(dst.data, 1, insize_, fd_);
		}
		dst.id = (void *)filename_;
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
			{ M2Decoder::MODE_H265, "265"},
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
	FileReaderUnit(int insize, int inbufnum, std::list<const char *> &infiles)
		: fr_(infiles, insize), dst_(inbufnum), outqueue_(&dst_[0], inbufnum) {
		for (int i = 0; i < inbufnum; ++i) {
			dst_[i].data = new unsigned char[insize];
			dst_[i].len = insize;
		}
	}
	~FileReaderUnit() {
		for (size_t i = 0; i < dst_.size(); ++i) {
			delete[] dst_[i].data;
		}
	}
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
	std::vector<Buffer> dst_;
	QueueType outqueue_;
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
	M2DecoderUnit(M2Decoder::type_t codec, FileReaderUnit::QueueType& inqueue, int dstnum, bool dpb_emptify)
		: m2dec_(codec, dstnum, reread_file, this),
		  dst_align_(dstnum),
		  inqueue_(inqueue), outqueue_(QueueType(&dst_align_[0], dstnum)),
		  src_next_(0), dpb_emptify_(dpb_emptify) {}
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
			if (!inqueue().pop_back() || inqueue().empty()) {
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
			dec_bits_set_data(stream, src->data, src->len, src->id);
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

#include <emmintrin.h>

void deinterleave(const uint8_t *src, uint8_t *dst, int stride, int height) {
	uint8_t* dst1 = dst + (stride >> 1);
	__m128i msk = _mm_cvtsi32_si128(0x00ff00ff);
	msk = _mm_shuffle_epi32(msk, 0);
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

#if 0
void display_write_halve(uint8_t **dst, const uint8_t *src_luma, const uint8_t *src_chroma, int src_stride,
		   uint16_t *pitches, int width, int height)
{
	const uint8_t *src = src_luma;
	uint8_t *dst0 = dst[0];
	int pitch0 = pitches[0];
	int y = height - 1;
	__m128i msk = _mm_cvtsi32_si128(0x00ff00ff);
	msk = _mm_shuffle_epi32(msk, 0);
	do {
		for (int x = 0; x < width; x += 16) {
			__m128i d0 = _mm_load_si128((__m128i const *)&src[x * 2]);
			__m128i d1 = _mm_load_si128((__m128i const *)&src[x * 2 + 16]);
			__m128i d0l = _mm_and_si128(d0, msk);
			__m128i d1l = _mm_and_si128(d1, msk);
			d0 = _mm_srli_epi16(d0, 8);
			d1 = _mm_srli_epi16(d1, 8);
			d0 = _mm_add_epi16(d0, d0l);
			d1 = _mm_add_epi16(d1, d1l);
			__m128i d2 = _mm_load_si128((__m128i const *)&src[x * 2 + src_stride]);
			__m128i d3 = _mm_load_si128((__m128i const *)&src[x * 2 + 16 + src_stride]);
			__m128i d2l = _mm_and_si128(d2, msk);
			__m128i d3l = _mm_and_si128(d3, msk);
			d0 = _mm_add_epi16(d0, d2l);
			d1 = _mm_add_epi16(d1, d3l);
			d2 = _mm_srli_epi16(d2, 8);
			d3 = _mm_srli_epi16(d3, 8);
			d0 = _mm_add_epi16(d0, d2);
			d1 = _mm_add_epi16(d1, d3);
			d0 = _mm_srli_epi16(d0, 2);
			d1 = _mm_srli_epi16(d1, 2);
			d0 = _mm_packus_epi16(d0, d1);
			_mm_storeu_si128((__m128i *)&dst0[x], d0);
		}
		src += src_stride * 2;
		dst0 += pitch0;
	} while (--y);

	src = src_chroma;
	dst0 = dst[1];
	uint8_t *dst1 = dst[2];
	width = width >> 1;
	height = height >> 1;
	src_stride *= 2;
	pitch0 = pitches[1];
	int pitch1 = pitches[2];
	do {
		for (int x = 0; x < width; x += 8) {
			__m128i d0 = _mm_load_si128((__m128i const *)&src[x * 4]);
			__m128i d1 = _mm_load_si128((__m128i const *)&src[x * 4 + 16]);
			__m128i d0l = _mm_and_si128(d0, msk);
			__m128i d1l = _mm_and_si128(d1, msk);
			d0 = _mm_srli_epi16(d0, 8);
			d1 = _mm_srli_epi16(d1, 8);
			d0l = _mm_packus_epi16(d0l, d1l);
			d0 = _mm_packus_epi16(d0, d1);
			d0l = _mm_and_si128(d0l, msk);
			d0 = _mm_and_si128(d0, msk);
			d0l = _mm_packus_epi16(d0l, d0l);
			d0 = _mm_packus_epi16(d0, d0);
			_mm_storel_epi64((__m128i *)&dst0[x], d0l);
			_mm_storel_epi64((__m128i *)&dst1[x], d0);
		}
		src += src_stride;
		dst0 += pitch0;
		dst1 += pitch1;
	} while (--height);
}

void display_write_normal(uint8_t **dst, const uint8_t *src_luma, const uint8_t *src_chroma, int src_stride,
		   uint16_t *pitches, int width, int height)
{
	const uint8_t *src = src_luma;
	uint8_t *dst0 = dst[0];
	int pitch0 = pitches[0];
	int y = height - 1;
	do {
		for (int x = 0; x < width; x += 32) {
			__m128i d0 = _mm_load_si128((__m128i const *)&src[x]);
			__m128i d1 = _mm_load_si128((__m128i const *)&src[x + 16]);
			_mm_storeu_si128((__m128i *)&dst0[x], d0);
			_mm_storeu_si128((__m128i *)&dst0[x + 16], d1);
		}
		src += src_stride;
		dst0 += pitch0;
	} while (--y);
	memcpy(dst0, src, width);

	src = src_chroma;
	dst0 = dst[1];
	uint8_t *dst1 = dst[2];
	height = (height >> 1) - 1;
	pitch0 = pitches[1];
	int pitch1 = pitches[2];
	__m128i msk = _mm_cvtsi32_si128(0x00ff00ff);
	msk = _mm_shuffle_epi32(msk, 0);
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
			_mm_storeu_si128((__m128i *)&dst0[x], d0l);
			_mm_storeu_si128((__m128i *)&dst1[x], d0);
		}
		src += src_stride;
		dst0 += pitch0;
		dst1 += pitch1;
	} while (--height);
	for (int j = 0; j < width; j += 2) {
		dst0[j] = src[j];
		dst1[j] = src[j + 1];
	}
}

void display_write_small(uint8_t **dst, const uint8_t *src_luma, const uint8_t *src_chroma, int src_stride,
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
#endif

const int MAX_WIDTH = 1920;
const int MAX_HEIGHT = 1088;
const int MAX_LEN = (MAX_WIDTH * MAX_HEIGHT * 3) >> 1;
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
	bool repeat_;
	bool logdump_;
	bool dpb_mode_;
	~Options() {}
	Options(int argc, char **argv)
		: interval_(0),
		  outbuf_(3), repeat_(false), logdump_(false), dpb_mode_(false) {
		int opt;
		while ((opt = getopt(argc, argv, "ef:hlmorst:")) != -1) {
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
	int scale_;
//	void (*display_write)(uint8_t **dst, const uint8_t *src_luma, const uint8_t *src_chroma, int src_stride, uint16_t *pitches, int width, int height);
#endif /* ENABLE_DISPLAY */
	int suspended_;
	void *id_;
	UniSurface() :
#ifdef ENABLE_DISPLAY
		screen_(0), renderer_(0), texture_(0),
		width_(0), height_(0),
		scale_(1),
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
//		SDL_LockYUVOverlay(yuv);
//		display_write(yuv->pixels, out.luma, out.chroma, out.width, yuv->pitches, yuv->w, yuv->h);
//		SDL_UnlockYUVOverlay(yuv);
		deinterleave(out.chroma, &cbcr_[0], out.width, out.height >> 1);
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
	void change(const Frame& out) {
		int width = width_ * scale_;
		int height = height_ * scale_;
		if ((width != (out.width - out.crop[1])) || (height != (out.height - out.crop[3]))) {
			teardown();
			width = out.width - out.crop[1];
			height = out.height - out.crop[3];
			cbcr_base_.resize(((width * height) >> 1) + 15);
			cbcr_ = reinterpret_cast<uint8_t*>(ALIGN16(reinterpret_cast<uintptr_t>(&cbcr_base_[0])));
			if ((1440 < width) || (720 < height)) {
				scale_ = 2;
//				width >>= 1;
//				height >>= 1;
//				display_write = display_write_halve;
//			} else if ((32 <= width) && (2 < height)) {
//				display_write = display_write_normal;
//			} else {
//				display_write = display_write_small;
			}
			width_ = width;
			height_ = height;
			screen_ = SDL_CreateWindow("disp", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, width, height, SDL_WINDOW_OPENGL);
			renderer_ = SDL_CreateRenderer(screen_, -1, 0);
			texture_ = SDL_CreateTexture(renderer_, SDL_PIXELFORMAT_IYUV, SDL_TEXTUREACCESS_STREAMING, width, height);
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

	UniThread *thr_file = UniCreateThread(FileReaderUnit::run, (void *)&fr, "FileLoader");
	UniThread *thr_m2d = UniCreateThread(M2DecoderUnit::run, (void *)&m2dec, "Decoder");

	UniSurface surface;
	while (1) {
		if (!surface.suspended()) {
			if (!m2dec.outqueue().pop_back()) {
				break;
			}
			const Frame& out = m2dec.outqueue().back();
			if (out.luma == 0) {
				break;
			}
			surface.display(out);
			opt.fw_.write(out);
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
#if defined(_M_IX86) && !defined(NDEBUG)
	_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_WNDW);
//	_CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_WNDW);
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF|_CRTDBG_LEAK_CHECK_DF);
	atexit((void (*)(void))_CrtCheckMemory);
#endif
#ifdef ENABLE_DISPLAY
	SDL_TimerID timer;
	if (SDL_Init(SDL_INIT_TIMER | SDL_INIT_VIDEO) < 0) {
		return -1;
	}
	atexit(SDL_Quit);
#endif /* ENABLE_DISPLAY */
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
	do {
		if (run_loop(opt, opt.outbuf_) < 0) {
			break;
		}
	} while (opt.repeat_);
#ifdef ENABLE_DISPLAY
	if (opt.interval_) {
		SDL_RemoveTimer(timer);
	}
#endif /* ENABLE_DISPLAY */
	return 0;
}

