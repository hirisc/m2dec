#ifdef ENABLE_AALIB

#ifndef _AADISP_H_
#define _AADISP_H_

#include <aalib.h>

class aadisplay {
	aa_context *context_;
	aa_renderparams *params_;
	aa_palette palette_;
	unsigned char *bitmap_;
	size_t width_, height_;
	void initpalette();
public:
	aadisplay();
	aadisplay(size_t width, size_t height) {
		set_size(width, height);
	}
	~aadisplay();
	void set_size(size_t width, size_t height) {
		width_ = width;
		height_ = height;
	}
	size_t get_width() {
		return aa_scrwidth(context_);
	}
	size_t get_height() {
		return aa_scrheight(context_);
	}
	unsigned char *get_bitmap() {
		return bitmap_;
	}
	unsigned char *get_txt() {
		return (unsigned char *)aa_text(context_);
	}
	void update() {
		aa_renderpalette(context_, palette_, params_, 0, 0, width_, height_);
		aa_flush(context_);
	}
	void writeback(const unsigned char *src, size_t width, size_t height);
	int display(void *luma, void *chroma, int width, int height) {return 0;}
};

#endif /* _AADISP_H_ */
#endif /* ENABLE_AALIB */
