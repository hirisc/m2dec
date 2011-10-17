#ifndef __UNITHREAD_H__
#define __UNITHREAD_H__

#ifndef __RENESAS_VERSION__
#include <SDL/SDL.h>

static void RecordTime(int mark);

#ifdef __cplusplus
extern "C" {
#endif

typedef SDL_Thread UniThread;
typedef SDL_mutex UniMutex;
typedef SDL_cond UniCond;
typedef SDL_Event UniEvent;

#define UniCreateThread(a, b) SDL_CreateThread((a), (b))
#define UniThreadID SDL_ThreadID
#define UniGetThreadID(x) SDL_GetThreadID(x)
#define UniWaitThread(x, y) SDL_WaitThread((x), (y))
#define UniKillThread(x) SDL_KillThread(x)

#define UniCreateMutex SDL_CreateMutex
#define UniLockMutex(x) {RecordTime(0); SDL_LockMutex(x); RecordTime(1);}
#define UniUnlockMutex(x) SDL_UnlockMutex(x)
#define UniDestroyMutex(x) SDL_DestroyMutex(x)
#define UniCreateCond SDL_CreateCond
#define UniDestroyCond(x) SDL_DestroyCond(x)
#define UniCondSignal(x) SDL_CondSignal(x)
#define UniCondBroadcast(x) SDL_CondBroadcast(x)
#define UniCondWait(x, y) {RecordTime(0); SDL_CondWait((x), (y)); RecordTime(1);}
#define UniCondWaitTimeout(x, y, z) {RecordTime(0); SDL_CondWaitTimeout((x), (y), (z)); RecordTime(1);}
#define UniPollEvent(x) SDL_PollEvent(x)
#define UniWaitEvent(x) {RecordTime(0); SDL_WaitEvent(x); RecordTime(1);}
#define UniPushEvent(x) SDL_PushEvent(x)
#define UniDelay(x) {RecordTime(0); SDL_Delay(x); RecordTime(1);}

#ifdef __cplusplus
}
#endif

#endif /* __RENESAS_VERSION__ */


#ifndef __RENESAS_VERSION__

typedef std::map<int, const char *> Logmap;
static Logmap LogTags;

static UniMutex *LogMutex;
struct LogRecord {
	unsigned long long time;
	int threadid;
	int mark;
	LogRecord(unsigned long long t, int id, int m) : time(t), threadid(id), mark(m) {}
};

typedef std::vector<LogRecord> LogListType;
static LogListType LogList;

static void LogInit()
{
	LogMutex = UniCreateMutex();
	LogList.reserve(65536);
}

static void LogFin()
{
	LogList.clear();
	LogTags.clear();
	UniDestroyMutex(LogMutex);
	LogMutex = 0;
}

#ifdef __GNUC__
static inline unsigned long long gettimer() {
	unsigned long long ret;
	__asm__ volatile ("rdtsc" : "=A" (ret));
	return ret;
}

#elif defined(_M_IX86)

inline __int64 __fastcall gettimer() {
	__asm{
		cpuid
		rdtsc
	}
}
#endif

#endif /* __RENESAS_VERSION__ */

static void RecordTime(int mark)
{
#ifndef __RENESAS_VERSION__
	SDL_LockMutex(LogMutex);
	LogList.push_back(LogRecord(gettimer(), UniThreadID(), mark));
	SDL_UnlockMutex(LogMutex);
#endif
}

#ifndef __RENESAS_VERSION__

struct PrintLabel {
	void operator() (std::pair<int, const char *> el) const {
		printf("\"%s\",", el.second);
	}
};

struct PrintMark {
	void operator() (int mark) const {
		printf("%d,", mark);
	}
};

#endif

static void LogDump()
{
#ifndef __RENESAS_VERSION__
	printf("\"time\",");
	std::for_each(LogTags.begin(), LogTags.end(), PrintLabel());
	printf("\n");
	std::map<int, int> LogCnv;
	int i = 0;
	for (Logmap::iterator p = LogTags.begin(); p != LogTags.end(); ++p) {
		LogCnv[p->first] = i++;
	}
	std::vector<int> marks(LogTags.size());
	std::fill(marks.begin(), marks.end(), 0);
	unsigned long long start = LogList.front().time;
	for (LogListType::iterator p = LogList.begin(); p != LogList.end(); ++p) {
		marks[LogCnv[p->threadid]] = p->mark;
		printf("%lld,", p->time - start);
		std::for_each(marks.begin(), marks.end(), PrintMark());
		printf("\n");
	}
#endif
}

#endif /* __UNITHREAD_H__ */

