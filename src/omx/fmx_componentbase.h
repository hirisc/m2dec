#ifndef _FMX_COMPONENTBASE_H_
#define _FMX_COMPONENTBASE_H_

struct ModuleFunction {
	OMX_ERRORTYPE (*init)(OMX_COMPONENTTYPE *comp);
};

struct ComponentList {
	const char *name;
	const char *uuid;
	ModuleFunction func;
};

class ComponentBase {
private:
	OMX_COMPONENTTYPE *host_;
	OMX_STATETYPE state_;
	OMX_CALLBACKTYPE callback_;
	vector<Port> ports_;
	bool loaded_;
	const ComponentList *id_;
#ifdef ENABLE_THREAD
	SDL_mutex *mutex_;
	SDL_cond *cond_;
	SDL_Thread *thread_;
#endif
	void Lock() {
#ifdef ENABLE_THREAD
		SDL_mutexP(mutex_);
#endif
	}
	void Unlock() {
#ifdef ENABLE_THREAD
		SDL_mutexV(mutex_);
#endif
	}

	void PushQueue(Port& port, OMX_BUFFERHEADERTYPE *buf);
#if 0
	OMX_BUFFERHEADERTYPE *PopQueue(Port &port);
#endif
	OMX_ERRORTYPE FlushQueue(OMX_U32 portnum);
	OMX_ERRORTYPE MarkBuffer(OMX_U32 portnum, OMX_MARKTYPE *mark);
	void TransState(OMX_STATETYPE state);
	/**State machine. Fall through down to error if invalid state transision.
	 */
	OMX_ERRORTYPE SetState(OMX_STATETYPE state);
	void CreateThread();
	void CondWait();
public:
	void SetCallback(OMX_CALLBACKTYPE *callback);
	ComponentBase(OMX_COMPONENTTYPE **host, OMX_PTR pAppData, OMX_CALLBACKTYPE* pCallBacks);
	~ComponentBase();

	/** Allocate initial resources.
	 */
	OMX_ERRORTYPE GoExecuting();
	/** Disllocate initial resources.
	 */
	OMX_ERRORTYPE Unload();
	OMX_ERRORTYPE GetState(OMX_STATETYPE *state);
	OMX_ERRORTYPE SendCommand(OMX_IN OMX_COMMANDTYPE Cmd, OMX_IN OMX_U32 nParam1, OMX_IN OMX_PTR pCmdData);
	OMX_ERRORTYPE EmptyThisBuffer(OMX_BUFFERHEADERTYPE* buf_head);
	OMX_ERRORTYPE FillThisBuffer(OMX_BUFFERHEADERTYPE* buf_head);
	int Handler();
	const ComponentList *Id() {
		return id_;
	}
	int Execute();
};

#endif /*  _FMX_COMPONENTBASE_H_ */
