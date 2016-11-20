if FSTDEC
DEFFST = -DFAST_DECODE
else
DEFFST =
endif

INCLUDES = -I$(includedir) -I../include
DEFS = $(DEFFST) -DUNIT_TEST #-fno-inline-functions
noinst_LTLIBRARIES = libm2dec.la libmpegdemux.la
libm2dec_la_SOURCES = config.h m2d.h m2d.cpp idct.cpp idct.h motioncomp.cpp motioncomp.h bitio.c bitio.h mpeg2.cpp mpeg2.h vld.h h265.cpp h265.h h265modules.h h264.cpp h264.h h264vld.h m2types.h m2d_macro.h txt2bin.c txt2bin.h h265_x86.cpp
libmpegdemux_la_SOURCES = mpeg_demux.cpp mpeg_demux.h bitio.h mpeg2.h
EXTRA_DIST = idct_asm.src motioncomp_asm.src
