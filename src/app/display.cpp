#ifdef __GNUC__
#include "config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include "display.h"

//#include "usr_sdltimer.h"

#ifdef ENABLE_DISPLAY

static int TBL_UV[ 256 * 4 ];

/* function pointer */
static void yuvtorgb565( unsigned char *yAdr, unsigned char *uvAdr, int x, int y, unsigned char *rgbAdrOrg, int rgbWidth );
static void yuvtorgb565_2( unsigned char *yAdr, unsigned char *uvAdr, int x, int y, unsigned char *rgbAdrOrg, int rgbWidth );
static void yuvtorgb888( unsigned char *yAdr, unsigned char *uvAdr, int x, int y, unsigned char *rgbAdrOrg, int rgbWidth );

#define clip255char( a ) ( !( (a) & ~255 ) ? (a) : ~( (unsigned)(a) >> 16 ) )
#define clip255int( a ) ( !( (a) & ~255 ) ? (a) : ( (unsigned)~(a) >> 24 ) )



#endif

Display::~Display()
{
#ifdef ENABLE_DISPLAY
	SDL_FreeSurface(disp_);
	timer_ = 0;
#endif
}

int Display::set_size(int width, int height)
{
#ifdef ENABLE_DISPLAY
	const SDL_VideoInfo *video_info;
	int width16, height16;
	int bpp;

	if (width <= 0 || height <= 0) {
		return -1;
	} else if ((width_ == width) && (height_ == height)) {
		return 0;
	}
	if ((0 < width_) || (0 < height_)) {
		SDL_FreeSurface(disp_);
	}
	width_ = width;
	height_ = height;
	width16_ = width16 = (width + 15) & ~15;
	height16_ = height16 = (height + 15) & ~15;
	right_crop_ = width16 - width;
	bottom_crop_ = height16 - height;

	video_info = SDL_GetVideoInfo();
	SDL_PixelFormat *fmt = video_info->vfmt;
	if ((fmt->BitsPerPixel == 32) && (fmt->Rshift == 16) && (fmt->Gshift == 8) && (fmt->Bshift == 0)) {
		yuvtorgb = yuvtorgb888;
		bpp = 32;
	} else {
		yuvtorgb = yuvtorgb565;
		bpp = 16;
	}
	disp_ = SDL_SetVideoMode(width16_, height16_, bpp,
//				      SDL_FULLSCREEN );
				      SDL_SWSURFACE );
//				      SDL_ANYFORMAT );
	if (disp_ == NULL) {
		return -1;
	}
#endif /* ENABLE_DISPLAY */
	return 0;
}

#ifdef ENABLE_DISPLAY
#if 0
static void usr_print_table( char *name, int *table )
{
	/* just for debug. */
	printf( "static int %s[] = {\xd\xa", name );
	for ( int i = 0; i < 512 / 8; i++ ) {
		int j;
		printf( "\t" );
		for ( j = 0; j < 8 - 1; j++ ) {
			printf( "%d, ", table[ j ] );
		}
		printf( "%d,\xd\xa", table[ j ] );
		table += 8;
	}
	printf( "};\xd\xa\xd\xa" );
}
#endif

#define BT709

static void table_generate()
{
	for ( int i = 0; i < 256; i++ ) {
#ifdef BT709
		/* BT.709 YCbCr==>RGB */
		TBL_UV[ i * 4 ] = -12275 * ( i - 128 ) + 0x8000; /* g */
		TBL_UV[ i * 4 + 2 ] = -30677 * ( i - 128 ); /* g */
		TBL_UV[ i * 4 + 1 ] = ( 121609 * ( i - 128 ) + 0x8000 ) >> 16;
		TBL_UV[ i * 4 + 3 ] = ( 103206 * ( i - 128 ) + 0x8000 ) >> 16;
#else
		/* BT.601(extended) YCbCr==>RGB */
		TBL_UV[ i * 4 ] = -22553 * ( i - 128 ) + 0x8000; /* g */
		TBL_UV[ i * 4 + 2 ] = -46802 * ( i - 128 ); /* g */
		TBL_UV[ i * 4 + 1 ] = ( 116130 * ( i - 128 ) + 0x8000 ) >> 16;
		TBL_UV[ i * 4 + 3 ] = ( 91882 * ( i - 128 ) + 0x8000 ) >> 16;
#endif
	}
}
#endif /* ENABLE_DISPLAY */

#ifdef ENABLE_DISPLAY
static void scan_event()
{
	SDL_Event ev;
	while (SDL_PollEvent(&ev)) {
		switch (ev.type) {
		case SDL_QUIT:
			exit(0);
			/* NOTREACHED */
		default:
			break;
		}
	}
}
#endif /* ENABLE_DISPLAY */

Display::Display()
{
#ifdef ENABLE_DISPLAY
	char Caption[ 32 ];
        /* Init Display */
	width_ = 0;
	height_ = 0;
	disp_ = 0;
	if (SDL_WasInit(SDL_INIT_VIDEO) == 0) {
		if (SDL_Init(SDL_INIT_VIDEO
#if 0
			       | SDL_INIT_EVENTTHREAD
#endif
			       ) < 0 ) {
			exit(1);
		}
		atexit(SDL_Quit);
//		SDL_EventState( SDL_KEYDOWN, SDL_IGNORE );
//		SDL_EventState( SDL_KEYUP, SDL_IGNORE );
///		play_mutex = SDL_CreateMutex();
///		if ( play_mutex == NULL ) {
///			return -1;
///		}
	}
	SDL_ShowCursor(SDL_DISABLE);
	sprintf(Caption, "MPEG-2 Decoder");
	SDL_WM_SetCaption(Caption, NULL);
///	usr_timer = new CSdlTimer();
	table_generate();
#endif /* ENABLE_DISPLAY */
}

#ifdef ENABLE_DISPLAY
#if 1
#define clipbit( val, bit ) ( !( ( (val) >> ( 8 - (bit) ) ) & ~( ( 1 << (bit) ) - 1 ) ) ? ( (val) >> ( 8 - (bit) ) ) : ( 0 < ( (val) >> ( 8 - (bit) ) ) ) ? ( ( 1 << (bit) ) - 1 ) : 0 )
#define pack_rgb( r, g, b, rbit, gbit, bbit ) ( ( clipbit( r, rbit ) << ( gbit + bbit ) ) | ( clipbit( g, gbit ) << bbit ) | clipbit( b, bbit ) )
#ifdef _M_IX86
#define construct_rgb( y, vr, uvg, ub, rbit, gbit, bbit ) ( r = ((y) + (vr)) >> 5, g = ((y) + (uvg)) >> 5, b = ((y) + (ub)) >> 5, pack_rgb( r, g, b, rbit, gbit, bbit ) )
#else
#define construct_rgb( y, vr, uvg, ub, rbit, gbit, bbit ) ( r = ((y) + (vr)), g = ((y) + (uvg)), b = ((y) + (ub)), pack_rgb( r, g, b, rbit, gbit, bbit ) )
#endif
#else

static inline int clipbit( int val, const int bit )
{
	val = val >> ( 8 - bit );
	return !( val & ~( ( 1 << bit ) - 1 ) ) ? val : ( 0 < val ) ? ( ( 1 << bit ) - 1 ) : 0;
}

static inline int pack_rgb( int r, int g, int b, const int rbit, const int gbit, const int bbit )
{
	return ( clipbit( r, rbit ) << ( gbit + bbit ) ) | ( clipbit( g, gbit ) << bbit ) | clipbit( b, bbit );
}

static inline int construct_rgb( int y, int vr, int uvg, int ub, const int rbit, const int gbit, const int bbit )
{
	int r, g, b;
	r = (y + vr) >> 5,
	g = (y + uvg) >> 5,
	b = (y + ub) >> 5;
	return pack_rgb( r, g, b, rbit, gbit, bbit );
}

#endif

static void yuvtorgb565( unsigned char *yAdr, unsigned char *uvAdr, int width, int height, unsigned char *rgbAdrOrg, int rgbWidth )
{
	height >>= 1;
	rgbWidth = (unsigned )rgbWidth >> 1;
	do {
		short *rgbAdr = (short *)rgbAdrOrg;
		for ( int j = 0; j < width; j += 2 ) {
			int u, v;
			int uvg, ub, vr;
			int r, g, b;
			u = uvAdr[ j ];
			v = uvAdr[ j + 1 ];
			uvg = ( TBL_UV[ u * 4 ] + TBL_UV[ v * 4 + 2 ] ) >> 16;
			ub = TBL_UV[ u * 4 + 1 ];
			vr = TBL_UV[ v * 4 + 3 ];
			rgbAdr[ j ] = construct_rgb( yAdr[ j + 0 ], vr, uvg, ub, 5, 6, 5 );
#ifdef NO_SKIP
			rgbAdr[ j + 1 ] = construct_rgb( yAdr[ j + 1 ], vr, uvg, ub, 5, 6, 5 );
			rgbAdr[ j + rgbWidth ] = construct_rgb( yAdr[ j + width ], vr, uvg, ub, 5, 6, 5 );
			rgbAdr[ j + rgbWidth + 1 ] = construct_rgb( yAdr[ j + width + 1 ], vr, uvg, ub, 5, 6, 5 );
#endif
		}
		rgbAdrOrg += rgbWidth * 2 * sizeof(short);
		yAdr += width * 2;
		uvAdr += width;
	} while ( --height );
}
 
/**Scale 1/2.
 */
static void yuvtorgb565_2( unsigned char *yAdr, unsigned char *uvAdr, int width, int height, unsigned char *rgbAdrOrg, int rgbWidth )
{
	height >>= 1;
	width >>= 1;
	rgbWidth = (unsigned )rgbWidth >> 2;
	do {
		short *rgbAdr = (short *)rgbAdrOrg;
		for ( int j = 0; j < width; j++ ) {
			int y, u, v;
			int uvg, ub, vr;
			int r, g, b;
			u = uvAdr[ j * 2 ];
			y = yAdr[ j * 2 + 0 ];
			v = uvAdr[ j * 2 + 1 ];
			uvg = ( ( TBL_UV[ u * 4 ] + TBL_UV[ v * 4 + 2 ] ) >> 16 );
			ub = TBL_UV[ u * 4 + 1 ];
			vr = TBL_UV[ v * 4 + 3 ];
			rgbAdr[ j ] = construct_rgb( y, vr, uvg, ub, 5, 6, 5 );
		}
		rgbAdrOrg += rgbWidth * 2 * sizeof(short);
		yAdr += width * 2 * 2;
		uvAdr += width * 2;
	} while ( --height );
}

#if defined(_M_IX86) || (defined(__GNUC__) && defined(__i386__))
#if defined(_M_IX86) && !defined(__attribute__)
#define __attribute__(x)
#endif
extern "C" const unsigned long long
#if defined(__GNUC__) && !defined(__CYGWIN__) && !defined(__MINGW32__)
_rgbScale[]
#else
rgbScale[]
#endif
__attribute__((used)) = {
	0x8080808080808080ULL, /* de-bias chroma */
#ifdef BT709
	0xffdefff3ffdefff3ULL, /* C(G): -0.2134 * 64, -0.5330 * 64 */
	0x0073008700730087ULL, /* C(R,B): 2.11218 * 64, 1.79287 * 64 */
	0x004b004b004b004bULL, /* Y: 1.1644 * 64 */
#else
	0xffccffe7ffccffe7ULL, /* C(G): -0.392 * 64, -0.813 * 64 */
	0x0066008100660081ULL, /* C(R,B): 2.017 * 64, 1.596 * 64 */
	0x004b004b004b004bULL, /* Y: 1.164 * 64 */
#endif
	0x0000ffff0000ffffULL, /* mask alpha */
	0x0000000000000000ULL,
};
#include <malloc.h>
#endif


static void yuvtorgb888( unsigned char *yAdr, unsigned char *uvAdr, int width, int height, unsigned char *rgbAdrOrg, int rgbWidth )
{
#if defined(__GNUC__) && defined(__i386__)
	asm volatile ("\n\t"
		"push		%%ebp\n\t"
		"push		%%ebx\n\t"
		"push		%%edi\n\t"
		"push		%%esi\n\t"
		"movl		%3, %%eax\n\t"
		"shr		$1, %%eax\n\t"
		"push		%%eax\n\t"
		"movl		%2, %%ebx\n\t"
		"movl		%0, %%eax\n\t"
		"movl		%1, %%edx\n\t"
		"movl		%%eax, %%esi\n\t"
		"add		%%ebx, %%esi\n\t"
		"add		%%ebx, %%ebx\n\t"
		"add		%%ebx, %%ebx\n\t"
		"movl		%4, %%ebp\n\t"
		"movl		%%ebx, %%edi\n\t"
		"add		%%ebp, %%edi\n\t"
		"shr		$5, %%ebx\n\t"
	"h00_yuvtorgb888:\n\t"
		"xor		%%ecx,%%ecx\n\t"
	"h01_yuvtorgb888:\n\t"
		"movq		(%%edx, %%ecx, 8), %%mm4\n\t"	// cbcr
		"pxor		%%mm3, %%mm3\n\t"
		"paddb		_rgbScale, %%mm4\n\t"		// -128s

		"pcmpgtb	%%mm4, %%mm3\n\t"
		"punpcklbw	%%mm3, %%mm4\n\t"
		"movq		%%mm4, %%mm5\n\t"
		"pmaddwd	_rgbScale + 8, %%mm4\n\t"	// 00 00 00 g1 00 00 00 g0
		"pmullw		_rgbScale + 16, %%mm5\n\t"	// 00 r1 00 b1 00 r0 00 b0
		"movq		(%%eax, %%ecx, 8), %%mm6\n\t"	// y
		"punpcklbw	_rgbScale + 40, %%mm6\n\t"
		"movq		(%%esi, %%ecx, 8), %%mm7\n\t"	// y(next line)
		"punpcklbw	_rgbScale + 40, %%mm7\n\t"
		"punpckldq	%%mm4, %%mm4\n\t"		// 00 00 00 g0 00 00 00 g0
		"punpckldq	%%mm5, %%mm5\n\t"		// 00 r0 00 b0 00 r0 00 b0
		"pmullw		_rgbScale + 24, %%mm6\n\t"	// 00 y3 00 y2 00 y1 00 y0
		"movq		%%mm6, %%mm0\n\t"
		"punpcklwd	%%mm6, %%mm6\n\t"		// 00 y1 00 y1 00 y0 00 y0
		"pmullw		_rgbScale + 24, %%mm7\n\t"
		"movq		%%mm7, %%mm1\n\t"
		"punpcklwd	%%mm7, %%mm7\n\t"
		"punpckhwd	%%mm0, %%mm0\n\t"		// 00 y3 00 y3 00 y2 00 y2
		"punpckhwd	%%mm1, %%mm1\n\t"
		"movq		%%mm6, %%mm2\n\t"
		"paddsw		%%mm4, %%mm6\n\t"
		"paddsw		%%mm5, %%mm2\n\t"
		"psraw		$6, %%mm6\n\t"
		"psraw		$6, %%mm2\n\t"
		"pand		_rgbScale + 32, %%mm6\n\t"
		"movq		%%mm2, %%mm3\n\t"
		"punpckhwd	%%mm6, %%mm2\n\t"		// 00 00 00 r0 00 g0 00 b0
		"punpcklwd	%%mm6, %%mm3\n\t"		// 00 00 00 r1 00 g1 00 b1
		"packuswb	%%mm2, %%mm3\n\t"
		"movntq		%%mm3, (%%ebp)\n\t"
		"movq		%%mm0, %%mm2\n\t"
		"paddsw		%%mm4, %%mm0\n\t"
		"paddsw		%%mm5, %%mm2\n\t"
		"psraw		$6, %%mm0\n\t"
		"psraw		$6, %%mm2\n\t"
		"pand		_rgbScale + 32, %%mm0\n\t"
		"movq		%%mm2, %%mm3\n\t"
		"punpcklwd	%%mm0, %%mm2\n\t"
		"punpckhwd	%%mm0, %%mm3\n\t"
		"packuswb	%%mm3, %%mm2\n\t"
		"movntq		%%mm2, 8(%%ebp)\n\t"

		"movq		%%mm7, %%mm2\n\t"
		"paddsw		%%mm4, %%mm7\n\t"
		"paddsw		%%mm5, %%mm2\n\t"
		"paddsw		%%mm1, %%mm4\n\t"
		"paddsw		%%mm1, %%mm5\n\t"
		"psraw		$6, %%mm7\n\t"
		"psraw		$6, %%mm4\n\t"
		"pand		_rgbScale + 32, %%mm7\n\t"
		"psraw		$6, %%mm2\n\t"
		"psraw		$6, %%mm5\n\t"
		"pand		_rgbScale + 32, %%mm4\n\t"
		"movq		%%mm2, %%mm3\n\t"
		"movq		%%mm5, %%mm1\n\t"
		"punpckhwd	%%mm7, %%mm2\n\t"		// 00 00 00 r0 00 g0 00 b0
		"punpcklwd	%%mm4, %%mm5\n\t"
		"punpcklwd	%%mm7, %%mm3\n\t"		// 00 00 00 r1 00 g1 00 b1
		"punpckhwd	%%mm4, %%mm1\n\t"
		"packuswb	%%mm2, %%mm3\n\t"
		"packuswb	%%mm1, %%mm5\n\t"
		"movntq		%%mm3, (%%edi)\n\t"
		"movntq		%%mm5, 8(%%edi)\n\t"

		"movq		(%%edx, %%ecx, 8), %%mm4\n\t"
		"pxor		%%mm3, %%mm3\n\t"
		"paddb		_rgbScale, %%mm4\n\t"		// -128s
		"pcmpgtb	%%mm4, %%mm3\n\t"
		"punpckhbw	%%mm3, %%mm4\n\t"
		"movq		%%mm4, %%mm5\n\t"
		"pmaddwd	_rgbScale + 8, %%mm4\n\t"	// 00 00 00 g1 00 00 00 g0
		"pmullw		_rgbScale + 16, %%mm5\n\t"	// 00 r1 00 b1 00 r0 00 b0
		"punpckldq	%%mm4, %%mm4\n\t"		// 00 00 00 g0 00 00 00 g0
		"punpckldq	%%mm5, %%mm5\n\t"		// 00 r0 00 b0 00 r0 00 b0

		"movq		(%%eax, %%ecx, 8), %%mm6\n\t"	// y
		"movq		(%%esi, %%ecx, 8), %%mm7\n\t"	// y(next line)
		"punpckhbw	_rgbScale + 40, %%mm6\n\t"
		"punpckhbw	_rgbScale + 40, %%mm7\n\t"
		"pmullw		_rgbScale + 24, %%mm6\n\t"
		"pmullw		_rgbScale + 24, %%mm7\n\t"
		"movq		%%mm6, %%mm0\n\t"
		"movq		%%mm7, %%mm1\n\t"
		"punpcklwd	%%mm6, %%mm6\n\t"
		"punpcklwd	%%mm7, %%mm7\n\t"
		"punpckhwd	%%mm0, %%mm0\n\t"
		"punpckhwd	%%mm1, %%mm1\n\t"
		"movq		%%mm6, %%mm3\n\t"
		"paddsw		%%mm4, %%mm6\n\t"
		"paddsw		%%mm5, %%mm3\n\t"
		"psraw		$6, %%mm6\n\t"
		"psraw		$6, %%mm3\n\t"
		"pand		_rgbScale + 32, %%mm6\n\t"
		"movq		%%mm3, %%mm2\n\t"
		"punpcklwd	%%mm6, %%mm3\n\t"
		"punpckhwd	%%mm6, %%mm2\n\t"
		"packuswb	%%mm2, %%mm3\n\t"
		"movntq		%%mm3, 16(%%ebp)\n\t"
		"movq		%%mm0, %%mm3\n\t"
		"paddsw		%%mm4, %%mm0\n\t"
		"paddsw		%%mm5, %%mm3\n\t"
		"psraw		$6, %%mm0\n\t"
		"psraw		$6, %%mm3\n\t"
		"pand		_rgbScale + 32, %%mm0\n\t"
		"movq		%%mm3, %%mm6\n\t"
		"punpcklwd	%%mm0, %%mm3\n\t"
		"punpckhwd	%%mm0, %%mm6\n\t"
		"packuswb	%%mm6, %%mm3\n\t"
		"movntq		%%mm3, 24(%%ebp)\n\t"

		"movq		%%mm7, %%mm2\n\t"
		"paddsw		%%mm4, %%mm7\n\t"
		"paddsw		%%mm5, %%mm2\n\t"
		"paddsw		%%mm1, %%mm4\n\t"
		"paddsw		%%mm1, %%mm5\n\t"
		"psraw		$6, %%mm7\n\t"
		"psraw		$6, %%mm2\n\t"
		"psraw		$6, %%mm4\n\t"
		"psraw		$6, %%mm5\n\t"
		"pand		_rgbScale + 32, %%mm7\n\t"
		"movq		%%mm2, %%mm6\n\t"
		"pand		_rgbScale + 32, %%mm4\n\t"
		"movq		%%mm5, %%mm1\n\t"
		"punpcklwd	%%mm7, %%mm2\n\t"
		"punpckhwd	%%mm7, %%mm6\n\t"
		"punpcklwd	%%mm4, %%mm5\n\t"
		"punpckhwd	%%mm4, %%mm1\n\t"
		"packuswb	%%mm6, %%mm2\n\t"
		"packuswb	%%mm1, %%mm5\n\t"
		"movntq		%%mm2, 16(%%edi)\n\t"
		"movntq		%%mm5, 24(%%edi)\n\t"

		"add	$1, %%ecx\n\t"
		"add	$32, %%ebp\n\t"
		"add	$32, %%edi\n\t"
		"cmp	%%ebx, %%ecx\n\t"
		"jne	h01_yuvtorgb888\n\t"
		"push	%%ebx\n\t"
		"shl	$3, %%ebx\n\t"
		"add	%%ebx, %%edx\n\t"
		"add	%%ebx, %%ebx\n\t"
		"add	%%ebx, %%eax\n\t"			// src += width * 2
		"add	%%ebx, %%esi\n\t"			// src2 += width * 2
		"add	%%ebx, %%ebx\n\t"
		"add	%%ebx, %%ebp\n\t"			// dst += width * 4
		"add	%%ebx, %%edi\n\t"
		"pop	%%ebx\n\t"
		"pop	%%ecx\n\t"
		"add	$-1, %%ecx\n\t"
		"push	%%ecx\n\t"
		"jnz	h00_yuvtorgb888\n\t"

		"pop	%%eax\n\t"
		"pop	%%esi\n\t"
		"pop	%%edi\n\t"
		"pop	%%ebx\n\t"
		"pop	%%ebp\n\t"
		"emms"
		:
		: "m"(yAdr), "m"(uvAdr), "m"(width), "m"(height),
		"m"(rgbAdrOrg), "m"(rgbWidth));
#elif defined(_M_IX86)
	__asm {
		push	ebp
		push	ebx
		push	edi
		push	esi
		mov	eax,height
		shr	eax,1
		push	eax		; height
		mov	ebx,width
		mov	eax,yAdr
		mov	edx,uvAdr
		mov	esi,eax
		add	esi,ebx
		add	ebx,ebx
		add	ebx,ebx
		mov	ebp,rgbAdrOrg
		mov	edi,ebx
		add	edi,ebp
		shr	ebx,5
h00_yuvtorgb888:
		xor	ecx,ecx
h01_yuvtorgb888:
		movq	mm4,[edx][ecx * 8]	; cbcr
		pxor	mm3,mm3
		paddb	mm4,rgbScale		; -128s

		pcmpgtb	mm3,mm4
		punpcklbw	mm4,mm3
		movq	mm5,mm4
		pmaddwd	mm4,rgbScale + 8	; 00 00 00 g1 00 00 00 g0
		movq	mm6,[eax][ecx * 8]	; y
		movq	mm7,[esi][ecx * 8]	; y(next line)
		pmullw	mm5,rgbScale + 16	; 00 r1 00 b1 00 r0 00 b0
		punpckldq	mm4,mm4		; 00 00 00 g0 00 00 00 g0
		punpckldq	mm5,mm5		; 00 r0 00 b0 00 r0 00 b0
		punpcklbw	mm6,rgbScale + 40
		pmullw	mm6,rgbScale + 24	; 00 y3 00 y2 00 y1 00 y0
		punpcklbw	mm7,rgbScale + 40
		pmullw	mm7,rgbScale + 24
		movq	mm0,mm6
		movq	mm1,mm7
		punpcklwd	mm6,mm6		; 00 y1 00 y1 00 y0 00 y0
		punpcklwd	mm7,mm7
		punpckhwd	mm0,mm0		; 00 y3 00 y3 00 y2 00 y2
		punpckhwd	mm1,mm1
		movq	mm2,mm6
		paddsw	mm6,mm4
		paddsw	mm2,mm5
		psraw	mm6,6
		psraw	mm2,6
		pand	mm6,rgbScale + 32
		movq	mm3,mm2
		punpckhwd	mm2,mm6		; 00 00 00 r0 00 g0 00 b0
		punpcklwd	mm3,mm6		; 00 00 00 r1 00 g1 00 b1
		packuswb	mm3,mm2
		movntq	[ebp],mm3
		movq	mm2,mm0
		paddsw	mm0,mm4
		paddsw	mm2,mm5
		psraw	mm0,6
		psraw	mm2,6
		pand	mm0,rgbScale + 32
		movq	mm3,mm2
		punpcklwd	mm2,mm0
		punpckhwd	mm3,mm0
		packuswb	mm2,mm3
		movntq	[ebp + 8],mm2

		movq	mm2,mm7
		paddsw	mm7,mm4
		paddsw	mm2,mm5
		paddsw	mm4,mm1
		paddsw	mm5,mm1
		psraw	mm7,6
		psraw	mm4,6
		psraw	mm2,6
		psraw	mm5,6
		pand	mm7,rgbScale + 32
		pand	mm4,rgbScale + 32
		movq	mm3,mm2
		movq	mm1,mm5
		punpckhwd	mm2,mm7		; 00 00 00 r0 00 g0 00 b0
		punpcklwd	mm5,mm4
		punpcklwd	mm3,mm7		; 00 00 00 r1 00 g1 00 b1
		punpckhwd	mm1,mm4
		packuswb	mm3,mm2
		packuswb	mm5,mm1
		movntq	[edi],mm3
		movntq	[edi + 8],mm5


		movq	mm4,[edx][ecx * 8]
		pxor	mm3,mm3
		paddb	mm4,rgbScale		; -128s
		pcmpgtb	mm3,mm4
		punpckhbw	mm4,mm3
		movq	mm5,mm4
		pmaddwd	mm4,rgbScale + 8	; 00 00 00 g1 00 00 00 g0
		pmullw	mm5,rgbScale + 16	; 00 r1 00 b1 00 r0 00 b0
		punpckldq	mm4,mm4		; 00 00 00 g0 00 00 00 g0
		punpckldq	mm5,mm5		; 00 r0 00 b0 00 r0 00 b0

		movq	mm6,[eax][ecx * 8]	; y
		movq	mm7,[esi][ecx * 8]	; y(next line)
		punpckhbw	mm6,rgbScale + 40
		punpckhbw	mm7,rgbScale + 40
		pmullw	mm6,rgbScale + 24
		pmullw	mm7,rgbScale + 24
		movq	mm0,mm6
		movq	mm1,mm7
		punpcklwd	mm6,mm6
		punpcklwd	mm7,mm7
		punpckhwd	mm0,mm0
		punpckhwd	mm1,mm1
		movq	mm3,mm6
		paddsw	mm6,mm4
		paddsw	mm3,mm5
		psraw	mm6,6
		psraw	mm3,6
		pand	mm6,rgbScale + 32
		movq	mm2,mm3
		punpcklwd	mm3,mm6
		punpckhwd	mm2,mm6
		packuswb	mm3,mm2
		movntq	[ebp + 16],mm3
		movq	mm3,mm0
		paddsw	mm0,mm4
		paddsw	mm3,mm5
		psraw	mm0,6
		psraw	mm3,6
		pand	mm0,rgbScale + 32
		movq	mm6,mm3
		punpcklwd	mm3,mm0
		punpckhwd	mm6,mm0
		packuswb	mm3,mm6
		movntq	[ebp + 24],mm3
;;;;
		movq	mm2,mm7
		paddsw	mm7,mm4
		paddsw	mm2,mm5
		paddsw	mm4,mm1
		paddsw	mm5,mm1
		psraw	mm7,6
		psraw	mm2,6
		psraw	mm4,6
		psraw	mm5,6
		pand	mm7,rgbScale + 32
		movq	mm6,mm2
		pand	mm4,rgbScale + 32
		movq	mm1,mm5
		punpcklwd	mm2,mm7
		punpckhwd	mm6,mm7
		punpcklwd	mm5,mm4
		punpckhwd	mm1,mm4
		packuswb	mm2,mm6
		packuswb	mm5,mm1
		movntq	[edi + 16],mm2
		movntq	[edi + 24],mm5
;;;;

		add	ecx,1
		add	ebp,32
		add	edi,32
		cmp	ecx,ebx
		jne	h01_yuvtorgb888

		push	ebx
		shl	ebx,3
		add	edx,ebx
		add	ebx,ebx
		add	eax,ebx			; src += width * 2
		add	esi,ebx			; src2 += width * 2
		add	ebx,ebx
		add	ebp,ebx			; dst += width * 4
		add	edi,ebx
		pop	ebx
		pop	ecx
		add	ecx,-1
		push	ecx
		jnz	h00_yuvtorgb888

		pop	eax			; dummy
		pop	esi
		pop	edi
		pop	ebx
		pop	ebp
	}
	__asm EMMS;
	/*
		movq	mm3,mm7
		movq	mm2,mm1
		paddsw	mm7,mm4
		paddsw	mm3,mm5
		paddsw	mm1,mm4
		paddsw	mm2,mm5
		psraw	mm7,6
		psraw	mm3,6
		psraw	mm1,6
		psraw	mm2,6
		pand	mm7,rgbScale + 32
		movq	mm6,mm3
		pand	mm1,rgbScale + 32
		movq	mm0,mm2
		punpcklwd	mm3,mm7
		punpckhwd	mm6,mm7
		punpcklwd	mm2,mm1
		punpckhwd	mm0,mm1
		packuswb	mm3,mm6
		packuswb	mm2,mm0
		movntq	[edi + 16],mm3
		movntq	[edi + 24],mm2
	*/
#else
	height >>= 1;
	rgbWidth = (unsigned )rgbWidth >> 2;
	do {
		unsigned int *rgbAdr = (unsigned int *)rgbAdrOrg;
		for ( int j = 0; j < width; j += 2 ) {
			int y, u, v;
			int uvg, ub, vr;
			int r, g, b;
			u = uvAdr[ j ];
			v = uvAdr[ j + 1 ];
			y = yAdr[ j + 0 ];
			uvg = ( ( TBL_UV[ u * 4 ] + TBL_UV[ v * 4 + 2 ] ) >> 16 );
			ub = TBL_UV[ u * 4 + 1 ];
			vr = TBL_UV[ v * 4 + 3 ];
			rgbAdr[ j ] = construct_rgb( y, vr, uvg, ub, 8, 8, 8 );
#ifdef NO_SKIP
			y = yAdr[ j + 1 ];
			rgbAdr[ j + 1 ] = construct_rgb( y, vr, uvg, ub, 8, 8, 8 );
			y = yAdr[ j + width ];
			rgbAdr[ j + rgbWidth ] = construct_rgb( y, vr, uvg, ub, 8, 8, 8 );
			y = yAdr[ j + width + 1 ];
			rgbAdr[ j + rgbWidth + 1 ] = construct_rgb( y, vr, uvg, ub, 8, 8, 8 );
#endif
		}
		rgbAdrOrg += rgbWidth * 2 * 4;
		yAdr += width * 2;
		uvAdr += width;
	} while ( --height );
#endif
}
#endif /* ENABLE_DISPLAY */

int Display::display(void *luma, void *chroma, int width, int height)
{
#ifdef ENABLE_DISPLAY
	int ret;

	if (!luma || !chroma || width <= 0 || height <= 0) {
		return -1;
	}
///	SDL_mutexP( play_mutex );
	ret = draw(luma, chroma, width, height);
///	SDL_mutexV( play_mutex );
	return ret;
#else /* ENABLE_DISPLAY */
	return 0;
#endif
}

int Display::draw(void *luma, void *chroma, int width, int height)
{
#ifdef ENABLE_DISPLAY
	if (set_size(width, height) < 0 ) {
		return -3;
	}

	if (SDL_MUSTLOCK(disp_) && (SDL_LockSurface(disp_) < 0)) {
		return -2;
	}
	yuvtorgb((unsigned char *)luma, (unsigned char *)chroma, width16_, height16_, (unsigned char *)disp_->pixels, disp_->pitch);
	if (SDL_MUSTLOCK(disp_)) {
		SDL_UnlockSurface(disp_);
	}
//	if ( DisplayInfo.displayMode & ( 1 << ELLAPSE_MODE ) ) {
//		usr_timer->EllapseTime();
//	}
	SDL_UpdateRect(disp_, 0, 0, width, height);
	scan_event();
#endif /* ENABLE_DISPLAY */
	return 0;
}

void Display::update()
{
#ifdef ENABLE_DISPLAY
#if 1
#else
	SDL_DisplayYUVOverlay( YuvOverlay, &Rect );
#endif
#endif
}

void Display::set_mode(int mode, bool on, void *arg)
{
#ifdef ENABLE_DISPLAY
#if 0
	if ( 31 < iCommand ) {
		return;
	}
	int iMask = 1 << iCommand;
	if ( iSet ) {
		DisplayInfo.displayMode |= iMask;
	} else {
		DisplayInfo.displayMode &= ~iMask;
	}
#endif
#endif
}

