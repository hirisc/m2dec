INCLUDES = -I$(includedir) -I../lib
DEFS = #-fno-inline-functions
noinst_PROGRAMS = omxtest
omxtest_LDADD = ../lib/libm2dec.a ../lib/libmpegdemux.a -L/usr/local/lib -lSDL
omxtest_SOURCES = main.cpp fmx_core.cpp fmx_componentbase.h fmx_componentbase.cpp fmx_port.h fmx_port.cpp fmx_common.h OMX_Audio.h OMX_Component.h OMX_ContentPipe.h OMX_Core.h OMX_IVCommon.h OMX_Image.h OMX_Index.h OMX_Other.h OMX_Types.h OMX_Video.h
