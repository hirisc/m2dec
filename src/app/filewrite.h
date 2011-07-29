#ifndef _FILEWRITE_H_
#define _FILEWRITE_H_

#include "md5.h"
#include "m2d.h"

class FileWriter {
protected:
	FILE *fo_;
	template <typename F>
	static void write_cropping(const m2d_frame_t *frame, void *dst, F writef) {
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
public:
	enum {
		WRITE_NONE,
		WRITE_MD5,
		WRITE_RAW
	};
	FileWriter(FILE *fo) : fo_(fo) {}
	~FileWriter() {
		if (fo_) {
			fclose(fo_);
		}
	}
	virtual size_t writeframe(const m2d_frame_t *frame) = 0;
};

class FileWriterRaw : public FileWriter {
	struct write_raw {
		size_t operator()(void *dst, const uint8_t *src, size_t size) {
			return fwrite((const void *)src, 1, size, (FILE *)dst);
		}
	};
public:
	FileWriterRaw(FILE *fo) : FileWriter(fo) {}
	size_t writeframe(const m2d_frame_t *frame) {
		if (!fo_) {
			return 0;
		}
		int luma_size = frame->width * frame->height;
		write_cropping(frame, (void *)fo_, write_raw());
		fflush(fo_);
		return luma_size;
	}
};

class FileWriterMd5 : public FileWriter {
protected:
	struct write_md5 {
		size_t operator()(void *dst, const uint8_t *src, size_t size) {
			md5_append((md5_state_t *)dst, (const md5_byte_t *)src, (int)size);
			return size;
		}
	};
	static void md5txt(const unsigned char *md5, char *txt) {
		for (int i = 0; i < 16; ++i) {
			sprintf(txt + i * 2, "%02xd", *md5++);
		}
		txt[32] = 0x0d;
		txt[33] = 0x0a;
	}
public:
	FileWriterMd5(FILE *fo) : FileWriter(fo) {}
	size_t writeframe(const m2d_frame_t *frame) {
		if (!fo_) {
			return 0;
		}
		int luma_size = frame->width * frame->height;
		md5_state_t md5;
		md5_byte_t digest[16];
		char txt[16 * 2 + 2];

		md5_init(&md5);
		write_cropping(frame, (void *)&md5, write_md5());
		md5_finish(&md5, digest);
		md5txt(digest, txt);
		fwrite(txt, 1, sizeof(txt), fo_);
		fflush(fo_);
		return luma_size;
	}
};

#endif /* _FILEWRITE_H_ */
