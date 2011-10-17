if DISP
DEFDISP = -DENABLE_DISPLAY
LDDISP = -L/usr/local/lib -lSDL
else
DEFDISP =
LDDISP =
endif

if AALIB
DEFAA = -DENABLE_AALIB
LDAA = -laa
else
DEFAA =
LDAA =
endif

INCLUDES = -I$(includedir) -I../lib
DEFS = $(DEFDISP) $(DEFAA) -DUNIT_TEST -D'_ANSI_ARGS_(x)=x' #-fno-inline-functions
noinst_PROGRAMS = m2dec thrplay h264dec
m2dec_LDADD = ../lib/libm2dec.a ../lib/libmpegdemux.a $(LDDISP) $(LDAA)
m2dec_SOURCES = m2dec.cpp buffer.cpp buffer.h display.h display.cpp aadisp.cpp aadisp.h lowlevel.c md5.c md5.h my_getopt.c my_getopt.h getopt.h optparser.h
thrplay_LDADD = ../lib/libm2dec.a ../lib/libmpegdemux.a $(LDDISP)
thrplay_SOURCES = threadplayer.cpp unithread.h md5.c md5.h my_getopt.c my_getopt.h getopt.h frames.h filewrite.h m2decoder.h
h264dec_LDADD = ../lib/libm2dec.a ../lib/libmpegdemux.a $(LDDISP)
h264dec_SOURCES = h264dec.cpp lowlevel.c md5.c md5.h my_getopt.c my_getopt.h getopt.h frames.h filewrite.h m2decoder.h
