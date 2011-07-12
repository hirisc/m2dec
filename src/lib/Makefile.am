if FSTDEC
DEFFST = -DFAST_DECODE
else
DEFFST =
endif

INCLUDES = -I$(includedir) -I../include
DEFS = $(DEFFST) -DUNIT_TEST #-fno-inline-functions
noinst_LIBRARIES = libm2dec.a libmpegdemux.a
libm2dec_a_SOURCES = config.h m2d.h m2d.cpp idct.cpp idct.h motioncomp.cpp motioncomp.h bitio.c bitio.h mpeg2.cpp mpeg2.h vld.h h264.cpp h264.h h264vld.h m2types.h txt2bin.c txt2bin.h
libmpegdemux_a_SOURCES = mpeg_demux.cpp mpeg_demux.h bitio.h mpeg2.h
EXTRA_DIST = idct_asm.src motioncomp_asm.src
