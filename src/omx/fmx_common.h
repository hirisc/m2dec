#ifndef _OMX_COMMON_H_
#define _OMX_COMMON_H_
#include <deque>
#include "OMX_Audio.h"
#include "OMX_Component.h"
#include "OMX_ContentPipe.h"
#include "OMX_Core.h"
#include "OMX_IVCommon.h"
#include "OMX_Image.h"
#include "OMX_Index.h"
#include "OMX_Other.h"
#include "OMX_Types.h"
#include "OMX_Video.h"

#ifndef __RENESAS_VERSION__
using namespace std;
#endif

#ifdef _M_IX86
#define snprintf sprintf_s
#define strncpy(a, b, c) strcpy_s((a), (c), (b))
#elif defined(__GNUC__) && !defined(__STRICT_ANSI__)
#else
/* TODO: define snprintf */
#define snprintf(a, b, c, d) sprintf((a), (c), (d))
#endif

#define NUM_OF_ELEM(a) (sizeof(a) / sizeof((a)[0]))

#define OMX_VERSION_MAJOR 1
#define OMX_VERSION_MINOR 1
#define OMX_VERSION_REVISION 0
#define OMX_VERSION_STEP 0
#define SET_VERSION(x) {(x)->nVersion.s.nVersionMajor = OMX_VERSION_MAJOR; (x)->nVersion.s.nVersionMinor = OMX_VERSION_MINOR; (x)->nVersion.s.nRevision = OMX_VERSION_REVISION;	(x)->nVersion.s.nStep = OMX_VERSION_STEP;}
#define SET_SIZE_VERSION(x) {(x)->nSize = sizeof(*(x)); SET_VERSION(x);}

#endif /* _OMX_COMMON_H_ */
