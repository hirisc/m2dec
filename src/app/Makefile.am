if DISP
DEFDISP = -DENABLE_DISPLAY
else
DEFDISP =
endif

LDM2LIB = ../lib/libm2dec.la ../lib/libmpegdemux.la
LDDISP = -L/usr/local/lib -lSDL2

INCLUDES = -I$(includedir) -I../lib
DEFS = $(DEFDISP) -DUNIT_TEST -D'_ANSI_ARGS_(x)=x' #-fno-inline-functions
noinst_PROGRAMS = thrplay h264dec
thrplay_LDADD = $(LDM2LIB) $(LDDISP)
thrplay_SOURCES = threadplayer.cpp unithread.h md5.c md5.h my_getopt.c my_getopt.h getopt.h frames.h filewrite.h m2decoder.h
h264dec_LDADD = $(LDM2LIB) $(LDDISP)
h264dec_SOURCES = h264dec.cpp lowlevel.c md5.c md5.h my_getopt.c my_getopt.h getopt.h frames.h filewrite.h m2decoder.h
