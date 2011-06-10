#ifndef _OPTPARSER_H_
#define _OPTPARSER_H_

#include <stdio.h>
#include <limits.h>
#include <vector>
#include <list>
#include <functional>
#include "getopt.h"

using namespace std;

class FileWriter {
protected:
	FILE *out_file_;
public:
	FileWriter() {
		out_file_ = 0;
	}
	virtual ~FileWriter() {
		if (out_file_) {
			printf("closed\n");
			fclose(out_file_);
		}
	}
	virtual FILE *Open(const char *filename) = 0;
	virtual size_t Write(const void *src, size_t size) = 0;
	virtual void Sync() = 0;
};

class RawFileWriter : public FileWriter {
public:
	virtual FILE *Open(const char *filename) {
		out_file_ = fopen(filename, "wb");
		return out_file_;
	}
	virtual size_t Write(const void *src, size_t size) {
		return fwrite(src, 1, size, out_file_);
	}
	virtual void Sync() {
	}
};

#include "md5.h"

class Md5FileWriter : public FileWriter {
	md5_state_t md5state_;
	void md5textize(md5_byte_t *md5, char *md5text) {
		for (int i = 0; i < 16; ++i) {
			sprintf(&md5text[i * 2], "%02x", md5[i]);
		}
		md5text[32] = 0x0d;
		md5text[33] = 0x0a;
	}
public:
	virtual FILE *Open(const char *filename) {
		md5_init(&md5state_);
		out_file_ = fopen(filename, "wb");
		return out_file_;
	}
	virtual size_t Write(const void *src, size_t size) {
		md5_append(&md5state_, (const md5_byte_t *)src, (int)size);
		return size;
	}
	virtual void Sync() {
		md5_byte_t md5[16];
		char md5text[16 * 2 + 2];
		md5_finish(&md5state_, md5);
		md5textize(md5, md5text);
		md5_init(&md5state_);
		fwrite(md5text, 1, 16 * 2 + 2, out_file_);
	}
};

template <typename _T>
struct Buffer {
	_T *buffer_;
	int size_;
};

template <typename _T, int _LEN>
struct BufferAlloc {
	void operator()(Buffer<_T>& elem) {
		elem.buffer_ = new _T[_LEN];
	}
};

template <typename _T>
struct BufferFree {
	void operator()(Buffer<_T>& elem) {
		delete[] elem.buffer_;
	}
};

template <typename _T, int _LEN, int _NUM>
class RingBuffer {
	vector <Buffer<_T> > buffer_;
	int front_, back_;
public:
	RingBuffer() : front_(0), back_(0), buffer_(_NUM + 1) {
		for_each(buffer_.begin(), buffer_.end(), BufferAlloc<_T, _LEN>());
	}
	~RingBuffer() {
		for_each(buffer_.begin(), buffer_.end(), BufferFree<_T>());
	}
	bool empty() {
		return (front_ == back_);
	}
	bool full() {
		return (front_ == next(back_));
	}
	int size() {
		int len = back_ - front_;
		return 0 <= len ? len : len + _LEN;
	}
	Buffer<_T> *back() {
		if (full()) {
			return 0;
		} else {
			return &buffer_[back_];
		}
	}
	Buffer<_T> *front() {
		if (empty()) {
			return 0;
		} else {
			return &buffer_[front_];
		}
	}
	void shift_back() {
		if (!full()) {
			back_ = next(back_);
		}
	}
	void pop_front() {
		if (!empty()) {
			front_ = next(front_);
		}
	}
private:
	static int next(int pos) {
		return _NUM < ++pos ? 0 : pos;
	}
};


struct Print {
	void operator()(const char *d) {
		printf("%s\n", d);
	}
};

typedef list<const char *> fname_t;

template <int _BUFF_LEN, int _BUFF_NUM>
class FileLoader {
	FILE *in_file_;
	fname_t filenames_;
	RingBuffer<byte_t, _BUFF_LEN, _BUFF_NUM> buffer_;
	bool finished_;
public:
	FileLoader(fname_t filenames) : in_file_(0), finished_(false) {
		for_each(filenames.begin(), filenames.end(), Print());
		if (filenames.empty()) {
			return;
		}
		filenames_ = fname_t(filenames);
		next_file();
		preload_buffers();
	}
	~FileLoader() {
		if (in_file_) {
			fclose(in_file_);
			in_file_ = 0;
		}
	}
#if 0
	void test_out(void *buffer, int size)
	{
		FILE *fout = fopen("binout.vob", "ab");
		fwrite(buffer, 1, size, fout);
		fclose(fout);
	}
#endif
	int next_buffer(const byte_t **buffer_p, int *size_p) {
		if (!buffer_p || !size_p || buffer_.empty()) {
			return -1;
		}
		Buffer<byte_t> *curr = buffer_.front();
		*buffer_p = (const byte_t *)curr->buffer_;
		*size_p = curr->size_;
		buffer_.pop_front();
		preload_buffers();
		return 0;
	}
private:
	int load_buffer() {
		Buffer<byte_t> *back = buffer_.back();
		if (back && (in_file_ || next_file())) {
			size_t len = fread(back->buffer_, 1, _BUFF_LEN, in_file_);
			back->size_ = len;
			buffer_.shift_back();
			if (len < _BUFF_LEN || feof(in_file_)) {
				fclose(in_file_);
				in_file_ = 0;
			}
			return 0;
		} else {
			return -1;
		}
	}
	int preload_buffers() {
		while (!buffer_.full()) {
			if (load_buffer() < 0) {
				return -1;
			}
		}
		return 0;
	}
	bool next_file() {
		if (filenames_.empty()) {
			in_file_ = 0;
			return false;
		} else if (in_file_) {
			fclose(in_file_);
		}
		in_file_ = fopen(filenames_.front(), "rb");
		filenames_.pop_front();
		return in_file_ ? true : false;
	}
};

struct Options {
	vector<byte_t> indata_[2];
	FileLoader<1024 * 1024, 2> *infile_;
	FileWriter *out_file_;
	int frame_num_;
	int skip_frames_;
	Display *disp_;
#ifdef ENABLE_AALIB
	aadisplay *aadisp_;
#endif
	bool demuxer_;
	bool verbose_;

	Options() {
		infile_ = 0;
		out_file_ = 0;
		frame_num_ = INT_MAX;
		skip_frames_ = 0;
		disp_ = 0;
#ifdef ENABLE_AALIB
		aadisp_ = 0; /* TODO: Must be changed to use polymorphing */
#endif
		demuxer_ = false;
		verbose_ = false;
	}
	~Options() {
		indata_[0].clear();
		indata_[1].clear();
		if (infile_) {
			delete infile_;
		}
		if (out_file_) {
			delete out_file_;
			out_file_ = 0;
		}
		if (disp_) {
			delete disp_;
		}
#ifdef ENABLE_AALIB
		if (aadisp_) {
			delete aadisp_;
		}
#endif
	}

	int next_buffer(const byte_t **indata_p, int *size_p) {
		return infile_->next_buffer(indata_p, size_p);
	}

	void write_file(uint8_t *luma, uint8_t *chroma, size_t width, size_t height) const {
		if (disp_) {
			disp_->display((void *)luma, (void *)chroma, (int)width, (int)height);
		}
#ifdef ENABLE_AALIB
		if (aadisp_) {
			aadisp_->writeback(luma, (size_t)width, (size_t)height);
			aadisp_->update();
		}
#endif
		if (out_file_) {
			size_t luma_size = width * height;
			out_file_->Write((const void *)luma, luma_size);
			out_file_->Write((const void *)chroma, luma_size >> 1);
			out_file_->Sync();
		}
	}

	static void blame_user() {
		fprintf(stderr,
			"Usage:\tm2dec [-o <outfile>] [-n <frame_num>] [-p] [-b <bytes>] <infile1> [<infile2> ....]\n"
			" -o outfile\tName of output file (NV12 Raw)\n"
			" -O outfile\tName of output file (NV12 MD5)\n"
			" -s \t\tPacketized stream (.vob files)\n"
			" -n N\t\tNumber of frame(s) to be decoded (default: Unlimited)\n"
			" -p \t\tDisplay output using SDL\n"
#ifdef ENABLE_AALIB
			" -a \t\tDisplay output using aalib\n"
#endif
			" -b Skip specified bytes\n"
			" -v \t\tShow verbose info\n"
			" infile(s)\tName of MPEG-2 input file\n"
			); 
	}

	int parse(char **argv, int argc) {
		int opt;
		if (argc <= 1) {
			blame_user();
			return -1;
		}
		while ((opt = getopt(argc, argv, "af:n:o:O:psv")) != -1) {
 			switch (opt) {
			case 'O':
				out_file_ = new Md5FileWriter();
				out_file_->Open(optarg);
				break;
			case 'a':
#ifdef ENABLE_AALIB
				aadisp_ = new aadisplay();
#endif
				break;
			case 'f':
				skip_frames_ = strtol(optarg, 0, 0);
				break;
			case 'o':
				out_file_ = new RawFileWriter();
				out_file_->Open(optarg);
				break;
			case 'p':
				disp_ = new Display;
				break;
			case 'n':
				frame_num_ = strtoul(optarg, 0, 0);
				break;
			case 's':
				demuxer_ = true;
				break;
			case 'v':
				verbose_ = true;
				break;
			default:
				blame_user();
				return -1;
			}
		}
		if (argc <= optind) {
#ifndef __RENESAS_VERSION__
			fprintf(stderr, "Invalid input\n");
			blame_user();
			return -1;
#endif
		}
		fname_t infiles;
		do {
			infiles.push_back(argv[optind++]);
		} while (optind < argc);
		infile_ = new FileLoader<1024 * 1024, 2>(infiles);
		return 0;
	}
	void set_demuxer_mode(bool mode) {
		demuxer_ = mode;
	}

	bool is_demuxer_mode() {
		return demuxer_;
	}
};

#endif /* _OPTPARSER_H_ */
