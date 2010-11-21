#include <memory>
#include <assert.h>
#include <stdio.h>
using namespace std;

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
#include "bitio.h"
#include "mpeg2.h"
#include "mpeg_demux.h"
#ifdef _M_IX86
#include <crtdbg.h>
#endif

#define pico(EXP) ( !!(EXP) || (printf("%s(%d) : error : " #EXP "\n", __FILE__,__LINE__),0))

OMX_ERRORTYPE null_handler(OMX_IN OMX_HANDLETYPE hComponent, OMX_IN OMX_PTR pAppData, OMX_IN OMX_EVENTTYPE eEvent, OMX_IN OMX_U32 nData1, OMX_IN OMX_U32 nData2, OMX_IN OMX_PTR pEventData)
{
	return OMX_ErrorNone;
}

OMX_ERRORTYPE dummyfunc(OMX_OUT OMX_HANDLETYPE hComponent, OMX_OUT OMX_PTR pAppData, OMX_OUT OMX_BUFFERHEADERTYPE* pBuffer)
{
	return OMX_ErrorNone;
}

OMX_CALLBACKTYPE null_callback = {
	null_handler,
	dummyfunc,
	dummyfunc
};

static void StateCheck(OMX_HANDLETYPE module, OMX_STATETYPE ref_state)
{
	OMX_ERRORTYPE err;
	OMX_STATETYPE state;

	err = OMX_GetState(module, &state);
	assert(err == OMX_ErrorNone);
	assert(state == ref_state);
}


int main(int argc, char **argv)
{
	OMX_ERRORTYPE err;
	OMX_HANDLETYPE module0 = 0;
#ifdef _M_IX86
	_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_WNDW);
//	_CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_WNDW);
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF|_CRTDBG_LEAK_CHECK_DF);
#endif
	auto_ptr<m2d_context> m2d(new m2d_context);

	err = OMX_APIENTRY OMX_Init();

	char component_name[127];
//	pico(argc == 3);
	int i;
	for (i = 0; i < 255; ++i) {
		err = OMX_ComponentNameEnum(component_name, sizeof(component_name), i);
		if (err != OMX_ErrorNone) {
			break;
		}
		printf("%s\n", component_name);
	}
	printf("number of components: %d\n", i);

	err = OMX_GetHandle(&module0, "OMX.MINE.VIDEO.DEC.MPEG2", (OMX_PTR)m2d.get(), &null_callback);
	if (err != OMX_ErrorNone) {
		if (module0) {
			OMX_FreeHandle(module0);
			module0 = 0;
		}
		return err;
	}
	StateCheck(module0, OMX_StateLoaded);

	char comp_name[128];
	OMX_UUIDTYPE uuid;
	OMX_VERSIONTYPE comp_ver, spec_ver;
	OMX_GetComponentVersion(module0, comp_name, &comp_ver, &spec_ver, &uuid);
	OMX_BUFFERHEADERTYPE *buf[2];
	OMX_STATETYPE state;

	err = OMX_SendCommand(module0, OMX_CommandStateSet, (OMX_U32)OMX_StateLoaded, 0);
	assert(err == OMX_ErrorSameState);
	StateCheck(module0, OMX_StateLoaded);

	err = OMX_SendCommand(module0, OMX_CommandStateSet, (OMX_U32)OMX_StateExecuting, 0);
	assert(err == OMX_ErrorIncorrectStateTransition);
	StateCheck(module0, OMX_StateLoaded);

	err = OMX_SendCommand(module0, OMX_CommandStateSet, (OMX_U32)OMX_StatePause, 0);
	assert(err == OMX_ErrorIncorrectStateTransition);
	StateCheck(module0, OMX_StateLoaded);

	err = OMX_SendCommand(module0, OMX_CommandStateSet, (OMX_U32)OMX_StateInvalid, 0);
	assert(err == OMX_ErrorNone);
	StateCheck(module0, OMX_StateInvalid);

	err = OMX_SendCommand(module0, OMX_CommandStateSet, (OMX_U32)OMX_StateLoaded, 0);
	assert(err == OMX_ErrorInvalidState);
	StateCheck(module0, OMX_StateInvalid);

	err = OMX_FreeHandle(module0);
	assert(err == OMX_ErrorNone);

// ----------------------------------------

	err = OMX_GetHandle(&module0, "OMX.MINE.VIDEO.DEC.MPEG2", (OMX_PTR)m2d.get(), &null_callback);
	for (i = 0; i < 2; ++i) {
		err = OMX_AllocateBuffer(module0, &buf[i], 0, 0, (640 * 480 * 3) / 2);
	}
	assert(err == OMX_ErrorNone);

	err = OMX_SendCommand(module0, OMX_CommandStateSet, (OMX_U32)OMX_StateIdle, 0);
	assert(err == OMX_ErrorNone);
	StateCheck(module0, OMX_StateIdle);

	err = OMX_SendCommand(module0, OMX_CommandStateSet, (OMX_U32)OMX_StateExecuting, 0);
	assert(err == OMX_ErrorNone);
	StateCheck(module0, OMX_StateExecuting);

	err = OMX_FillThisBuffer(module0, buf[0]);
	err = OMX_FillThisBuffer(module0, buf[0]);
	err = OMX_EmptyThisBuffer(module0, buf[1]);

	err = OMX_SendCommand(module0, OMX_CommandFlush, (OMX_U32)0, 0);
	err = OMX_SendCommand(module0, OMX_CommandStateSet, (OMX_U32)OMX_StateIdle, 0);
	assert(err == OMX_ErrorNone);
	do {
		err = OMX_GetState(module0, &state);
	} while (state != OMX_StateIdle);
	err = OMX_FreeHandle(module0);

#ifdef _M_IX86
	assert(_CrtCheckMemory());
#endif
	for (i = 0; i < 2; ++i) {
		if (buf[i]) {
			err = OMX_FreeBuffer(module0, 0, buf[i]);
		}
	}
	err = OMX_Deinit();

#ifdef _M_IX86
	assert(_CrtCheckMemory());
#endif
	return (int)err;
}
