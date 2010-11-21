
#ifdef ENABLE_AALIB

#include <stdint.h>
#include <assert.h>
#include "aadisp.h"

void aadisplay::initpalette()
{
	for (int i = 0; i < 255; ++i) {
		aa_setpalette(palette_, i, i, i, i);
	}
}

aadisplay::aadisplay()
{
	aa_hardware_params prm;
	prm = aa_defparams;
	context_ = aa_autoinit(&prm);
	initpalette();
	aa_autoinitkbd(context_, 0);
	params_ = aa_getrenderparams();
	bitmap_ = aa_image(context_);
	set_size(aa_imgwidth(context_), aa_imgheight(context_));
}

aadisplay::~aadisplay()
{
	if (context_) {
		aa_close(context_);
	}
}

void aadisplay::writeback(const unsigned char *src, size_t width, size_t height)
{
	assert((0 < (intptr_t)width) && (0 < (intptr_t)height));
	size_t y_max = height_;
	size_t x_max = width_;
	int y_delta = (height / y_max) * width;
	int x_delta = width / x_max;
	unsigned char *dst = bitmap_;
	for (size_t y = 0; y < y_max; ++y) {
		const unsigned char *s = src;
		for (size_t x = 0; x < x_max; x += 2) {
			dst[x] = s[0];
			dst[x + 1] = s[x_delta];
			s += x_delta * 2;
		}
		dst += x_max;
		src += y_delta;
	}
}

#endif /* ENABLE_AALIB */
