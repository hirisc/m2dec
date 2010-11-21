#ifndef _DISPLAY_H_
#define _DISPLAY_H_

#ifdef ENABLE_DISPLAY

#if defined(HAVE_SDL_H)
#include <SDL.h>
#include <SDL_thread.h>
#else
#include <SDL/SDL.h>
#include <SDL/SDL_thread.h>
#endif

#endif /* ENABLE_DISPLAY */

typedef enum {
	None = 0,
	Timer = 1
} DisplayMode;

class Display {
	int width_;
	int height_;
	int width16_;
	int height16_;
	int right_crop_;
	int bottom_crop_;
#ifdef ENABLE_DISPLAY
	SDL_Surface *disp_;
	void (*yuvtorgb)( unsigned char *yAdr, unsigned char *uvAdr, int x, int y, unsigned char *rgbAdrOrg, int rgbWidth );
#endif
	DisplayMode mode_;
	int timer_;
	int draw(void *luma, void *chroma, int width, int height);
public:
	Display();
	~Display();
	int set_size(int width, int height);
	int display(void *luma, void *chroma, int width, int height);
	void update();
	void set_mode(int mode, bool on, void *arg);
};

#endif /* _DISPLAY_H_ */
