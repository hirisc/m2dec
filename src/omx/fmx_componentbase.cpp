#include <string.h>
#include <vector>
#include <functional>
#include <algorithm>
#include "fmx_common.h"
#include "fmx_port.h"
#include "fmx_componentbase.h"

OMX_ERRORTYPE dummy_init(OMX_COMPONENTTYPE *comp)
{
	return OMX_ErrorNone;
}

static const ComponentList Components[] = {
	{"OMX.MINE.TIMER", "651a54ac-3806-44eb-b4ba-0ac902fff139", {dummy_init}},
	{"OMX.MINE.FILE.READER", "6cbf060b-97e1-45c8-998d-f07c5e60024c", {dummy_init}},
	{"OMX.MINE.FILE.WRITER", "179186de-7d12-4ce9-bd7c-02671a81f150", {dummy_init}},
	{"OMX.MINE.VIDEO.DEC.MPEG2", "e0055041-d543-4df7-bea2-d17ca227cfe0", {dummy_init}},
	{"OMX.MINE.IMAGE.DISPLAY", "0ceaba6e-06a6-4e84-92d0-9be33ddbe2ea", {dummy_init}},
};

OMX_API OMX_ERRORTYPE OMX_APIENTRY OMX_ComponentNameEnum(OMX_OUT OMX_STRING cComponentName, OMX_IN  OMX_U32 nNameLength, OMX_IN  OMX_U32 nIndex)
{
	if (!cComponentName || (128 <= nNameLength) || (NUM_OF_ELEM(Components) <= nIndex)) {
		return OMX_ErrorBadParameter;
	}
	snprintf(cComponentName, nNameLength, "%s", Components[nIndex].name);
	return OMX_ErrorNone;
}

const ModuleFunction *find_component_function(OMX_STRING cComponentName)
{
	for (const ComponentList *lst = Components; lst != Components + NUM_OF_ELEM(Components); ++lst) {
		if (strcmp(lst->name, cComponentName)) {
			return &lst->func;
		}
	}
	return 0;
}

static OMX_ERRORTYPE null_handler(OMX_IN OMX_HANDLETYPE hComponent, OMX_IN OMX_PTR pAppData, OMX_IN OMX_EVENTTYPE eEvent, OMX_IN OMX_U32 nData1, OMX_IN OMX_U32 nData2, OMX_IN OMX_PTR pEventData)
{
	return OMX_ErrorNone;
}

static OMX_ERRORTYPE dummyfunc(OMX_OUT OMX_HANDLETYPE hComponent, OMX_OUT OMX_PTR pAppData, OMX_OUT OMX_BUFFERHEADERTYPE* pBuffer)
{
	return OMX_ErrorNone;
}


static const OMX_CALLBACKTYPE dummycallback = {
	null_handler,
	dummyfunc,
	dummyfunc
};

OMX_ERRORTYPE omx_AllocateBuffer(OMX_IN OMX_HANDLETYPE hComponent, OMX_INOUT OMX_BUFFERHEADERTYPE **ppBuffer, OMX_IN OMX_U32 nPortIndex, OMX_IN OMX_PTR pAppPrivate, OMX_IN OMX_U32 nSizeBytes);
OMX_ERRORTYPE omx_FreeBuffer(OMX_IN OMX_HANDLETYPE hComponent, OMX_IN OMX_U32 nPortIndex, OMX_IN OMX_BUFFERHEADERTYPE *pBuffer);
OMX_ERRORTYPE omx_SetCallbacks(OMX_IN OMX_HANDLETYPE hComponent, OMX_IN OMX_CALLBACKTYPE *pCallbacks, OMX_IN OMX_PTR pAppData);
OMX_ERRORTYPE omx_GetState(OMX_IN OMX_HANDLETYPE hComponent, OMX_OUT OMX_STATETYPE* pState);
OMX_ERRORTYPE omx_SendCommand(OMX_IN OMX_HANDLETYPE hComponent, OMX_IN OMX_COMMANDTYPE Cmd, OMX_IN OMX_U32 nParam1, OMX_IN OMX_PTR pCmdData);
OMX_ERRORTYPE omx_GetComponentVersion(OMX_IN OMX_HANDLETYPE hComponent, OMX_OUT OMX_STRING pComponentName, OMX_OUT OMX_VERSIONTYPE* pComponentVersion, OMX_OUT OMX_VERSIONTYPE* pSpecVersion, OMX_OUT OMX_UUIDTYPE* pComponentUUID);
OMX_ERRORTYPE omx_EmptyThisBuffer(OMX_IN OMX_HANDLETYPE hComponent, OMX_IN OMX_BUFFERHEADERTYPE* pBuffer);
OMX_ERRORTYPE omx_FillThisBuffer(OMX_IN OMX_HANDLETYPE hComponent, OMX_IN OMX_BUFFERHEADERTYPE* pBuffer);

void ComponentBase::PushQueue(Port& port, OMX_BUFFERHEADERTYPE *buf) {
#ifdef ENABLE_THREAD
	Lock();
	while (MAX_QUEUE <= port.queue().size()) {
		SDL_CondWait(cond_, mutex_);
	}
	Unlock();
#endif
	port.queue().push_back(buf);
#ifdef ENABLE_THREAD
	SDL_CondSignal(cond_);
#endif
}

#if 0
OMX_BUFFERHEADERTYPE *ComponentBase::PopQueue(Port &port) {
	OMX_BUFFERHEADERTYPE *buf
		Lock();
	if (port.queue.empty()) {
		Unlock();
		return 0;
	} else {
		buf = port.queue.front();
		queue.pop_front();
		Unlock();
#ifdef ENABLE_THREAD
		SDL_CondSignal(cond_);
#endif
		return buf;
	}
}
#endif

OMX_ERRORTYPE ComponentBase::FlushQueue(OMX_U32 portnum) {
	if (portnum == OMX_ALL) {
	}
	return OMX_ErrorNone;
}

OMX_ERRORTYPE ComponentBase::MarkBuffer(OMX_U32 portnum, OMX_MARKTYPE *mark) {
	return OMX_ErrorNone;
}

void ComponentBase::TransState(OMX_STATETYPE state) {
	Lock();
	state_ = state;
	Unlock();
}

/**State machine. Fall through down to error if invalid state transision.
 */
OMX_ERRORTYPE ComponentBase::SetState(OMX_STATETYPE state) {
	if (state == state_) {
		/* Not affected */
		return OMX_ErrorSameState;
	}
	if (state == OMX_StateInvalid) {
		/* Invalidation are always accept */
		TransState(state);
		return OMX_ErrorNone;
	}
	for_each(ports_.begin(), ports_.end(), bind2nd(mem_fun_ref(&Port::SetState), state));
	switch (state_) {
	case OMX_StateLoaded:
		switch (state) {
		case OMX_StateIdle:
			TransState(state);
			return OMX_ErrorNone;
			break;
		default:
			break;
		}
		break;
	case OMX_StateIdle:
		switch (state) {
		case OMX_StateExecuting:
		case OMX_StatePause:
			TransState(state);
			GoExecuting();
			return OMX_ErrorNone;
			break;
		case OMX_StateLoaded:
			TransState(state);
			return OMX_ErrorNone;
			break;
		default:
			break;
		}
		break;
	case OMX_StateExecuting:
		switch (state) {
		case OMX_StatePause:
			TransState(state);
			return OMX_ErrorNone;
			break;
		case OMX_StateIdle:
			TransState(state);
			return OMX_ErrorNone;
			break;
		default:
			break;
		}
		break;
	case OMX_StatePause:
		switch (state) {
		case OMX_StateExecuting:
			TransState(state);
			return OMX_ErrorNone;
			break;
		case OMX_StateIdle:
			TransState(state);
			return OMX_ErrorNone;
			break;
		default:
			break;
		}
		break;
	case OMX_StateInvalid:
		return OMX_ErrorInvalidState;
		break;
	default:
		break;
	}
	return OMX_ErrorIncorrectStateTransition;
}

#ifdef ENABLE_THREAD
static int BufferHandler(void *data);
#endif

void ComponentBase::CreateThread() {
#ifdef ENABLE_THREAD
	thread_ = SDL_CreateThread(BufferHandler, (void *)this);
	while (!loaded_);
#endif
}

void ComponentBase::CondWait() {
#ifdef ENABLE_THREAD
	SDL_CondWait(cond_, mutex_);
#endif
}

ComponentBase::ComponentBase(OMX_COMPONENTTYPE **host, OMX_PTR pAppData, OMX_CALLBACKTYPE* pCallBacks)
     : state_(OMX_StateLoaded), callback_(dummycallback), loaded_(false), id_(&Components[3])
{
	OMX_COMPONENTTYPE *comp = new OMX_COMPONENTTYPE;
	*host = host_ = comp;
	memset((void *)comp, 0, sizeof(*comp));
	comp->pComponentPrivate = (OMX_PTR)this;
	comp->nSize = sizeof(OMX_COMPONENTTYPE) + sizeof(ComponentBase);
	SET_VERSION(comp);
	comp->pApplicationPrivate = pAppData;
	comp->SetCallbacks = omx_SetCallbacks;
	comp->AllocateBuffer = omx_AllocateBuffer;
	comp->FreeBuffer = omx_FreeBuffer;
	comp->GetComponentVersion = omx_GetComponentVersion;
	comp->EmptyThisBuffer = omx_EmptyThisBuffer;
	comp->FillThisBuffer = omx_FillThisBuffer;
	comp->GetState = omx_GetState;
	comp->SendCommand = omx_SendCommand;
	comp->SetCallbacks(comp, pCallBacks, comp);

#ifdef ENABLE_THREAD
	thread_ = 0;
	if (!SDL_WasInit(SDL_INIT_VIDEO)) {
		if (SDL_Init(SDL_INIT_VIDEO
#if 0
			     | SDL_INIT_EVENTTHREAD
#endif
			     ) < 0 ) {
			exit(1);
		}
		atexit(SDL_Quit);
	}
	mutex_ = SDL_CreateMutex();
	cond_ = SDL_CreateCond();
#endif
}

ComponentBase::~ComponentBase() {
	delete host_;
#ifdef ENABLE_THREAD
	loaded_ = false;
	SDL_CondSignal(cond_);
	if (thread_) {
		SDL_WaitThread(thread_, 0);
		thread_ = 0;
	}
	SDL_DestroyCond(cond_);
	SDL_DestroyMutex(mutex_);
#endif
}

/** Allocate initial resources.
 */
OMX_ERRORTYPE ComponentBase::GoExecuting() {
	Lock();
	CreateThread();
	Unlock();
	return OMX_ErrorNone;
}

/** Disllocate initial resources.
 */
OMX_ERRORTYPE ComponentBase::Unload() {
	loaded_ = false;
	return OMX_ErrorNone;
}

OMX_ERRORTYPE ComponentBase::GetState(OMX_STATETYPE *state) {
	*state = state_;
	return OMX_ErrorNone;
}

OMX_ERRORTYPE ComponentBase::SendCommand(OMX_IN OMX_COMMANDTYPE Cmd, OMX_IN OMX_U32 nParam1, OMX_IN OMX_PTR pCmdData) {
	OMX_ERRORTYPE err;
	if ((Cmd != OMX_CommandStateSet) && (ports_.size() <= nParam1)) {
		return OMX_ErrorBadPortIndex;
	}
	err = OMX_ErrorNone;
	Lock();
	switch (Cmd) {
	case OMX_CommandStateSet:
		err = SetState((OMX_STATETYPE)nParam1);
		break;
	case OMX_CommandFlush:
		err = FlushQueue(nParam1);
		break;
	case OMX_CommandPortDisable:
		ports_[nParam1].enable(false);
		break;
	case OMX_CommandPortEnable:
		ports_[nParam1].enable(true);
		break;
	case OMX_CommandMarkBuffer:
		err = MarkBuffer(nParam1, (OMX_MARKTYPE *)pCmdData);
		break;
	default:
		err = OMX_ErrorBadParameter;
		break;
	}
	Unlock();
	callback_.EventHandler(host_, ((OMX_COMPONENTTYPE *)host_)->pApplicationPrivate, OMX_EventCmdComplete, Cmd, 0, 0);
	return err;
}

OMX_ERRORTYPE ComponentBase::EmptyThisBuffer(OMX_BUFFERHEADERTYPE* buf_head) {
	//		PushQueue(in_queue_, cmd);
	return OMX_ErrorNone;
}

OMX_ERRORTYPE ComponentBase::FillThisBuffer(OMX_BUFFERHEADERTYPE* buf_head) {
	OMX_U32 portnum;
	Port *port;

	portnum = buf_head->nOutputPortIndex;
	if (ports_.size() <= portnum) {
		return OMX_ErrorBadPortIndex;
	}
	port = &ports_[portnum];
	if (port->port_def().eDir != OMX_DirOutput) {
		return OMX_ErrorBadPortIndex;
	}
	PushQueue(*port, buf_head);
	return OMX_ErrorNone;
}

int ComponentBase::Handler() {
	int err = 0;

	loaded_ = true;
	Lock();
	do {
		vector<Port>::iterator port;

		CondWait();
		for (port = ports_.begin(); port != ports_.end(); ++port) {
			if (port[0].queue().empty()) {
				break;
			}
		}
		if (port != ports_.end()) {
			continue;
		}
		err = Execute();
		Lock();
	} while (err || (state_ == OMX_StateExecuting) || (state_ == OMX_StatePause));
	Unlock();
	return 0;
}

int ComponentBase::Execute() {
#if 0
	OMX_BUFFERHEADERTYPE *in, *out;
	in = PopQueue(ports_in_queue_.front();
	out = out_queue_.front();
	in_queue_.pop_front();
	out_queue_.pop_front();
#endif
	Unlock();
	return 0;
}

void ComponentBase::SetCallback(OMX_CALLBACKTYPE *callback)
{
	assert(callback);
	callback_ = *callback;
}

#ifdef ENABLE_THREAD
static int BufferHandler(void *data) {
	ComponentBase *comp = (ComponentBase *)data;
	return comp->Handler();
}
#endif

