/**OpenMAX(TM) Integration Layer Core Module.
 *  OpenMAX(TM) implementation
 *  Copyright 2008 Takayuki Minegishi
 *
 *  Permission is hereby granted, free of charge, to any person
 *  obtaining a copy of this software and associated documentation
 *  files (the "Software"), to deal in the Software without
 *  restriction, including without limitation the rights to use, copy,
 *  modify, merge, publish, distribute, sublicense, and/or sell copies
 *  of the Software, and to permit persons to whom the Software is
 *  furnished to do so, subject to the following conditions:
 *  
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the Software.
 *  
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *  HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *  WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 *  DEALINGS IN THE SOFTWARE.
 */

#include <stdio.h>
#include <assert.h>
#include <list>
#include <vector>
#include <algorithm>
#include <functional>
#include "fmx_common.h"
#include "fmx_port.h"
#include "fmx_componentbase.h"

#define ENABLE_THREAD

#ifdef ENABLE_THREAD
#include <SDL/SDL.h>
#include <SDL/SDL_thread.h>
#endif

struct Tunnel {
	OMX_HANDLETYPE out_;
	OMX_U32 out_port_;
	OMX_HANDLETYPE in_;
	OMX_U32 in_port_;
	Tunnel(OMX_HANDLETYPE out, OMX_U32 out_port, OMX_HANDLETYPE in, OMX_U32 in_port) : out_(out), out_port_(out_port), in_(in), in_port_(in_port) {
	}
};

template <typename _T>
struct DeleteObject {
	void operator()(const _T* obj) const {
		delete obj;
	}
};

class OpenmaxCore {
	list<Tunnel *> Tunnels_;
public:
	OpenmaxCore() {
	}
	~OpenmaxCore() {
		if (!Tunnels_.empty()) {
			for_each(Tunnels_.begin(), Tunnels_.end(), DeleteObject<Tunnel>());
		}
	}
	void PushTunnel(OMX_HANDLETYPE out, OMX_U32 out_port, OMX_HANDLETYPE in, OMX_U32 in_port) {
		Tunnels_.push_front(new Tunnel(out, out_port, in, in_port));
	}
};

static OpenmaxCore *Core = 0;

OMX_API OMX_ERRORTYPE OMX_APIENTRY OMX_Init(void)
{
	Core = new OpenmaxCore;
	return OMX_ErrorNone;
}

OMX_API OMX_ERRORTYPE OMX_APIENTRY OMX_Deinit(void)
{
	if (Core) {
		delete Core;
		Core = 0;
	}
	return OMX_ErrorNone;
}

const ModuleFunction *find_component_function(OMX_STRING cComponentName);

OMX_API OMX_ERRORTYPE OMX_APIENTRY OMX_GetHandle(OMX_OUT OMX_HANDLETYPE* pHandle, OMX_IN  OMX_STRING cComponentName, OMX_IN  OMX_PTR pAppData, OMX_IN  OMX_CALLBACKTYPE* pCallBacks)
{
	ComponentBase *base;
	OMX_COMPONENTTYPE *comp;
	const ModuleFunction *specific_func;

	if (!pHandle || !cComponentName || !pAppData || !pCallBacks || !(specific_func = find_component_function(cComponentName))) {
		return OMX_ErrorBadParameter;
	}
	base = new ComponentBase(&comp, pAppData, pCallBacks);
	*pHandle = comp;
	return specific_func->init(comp);
}

OMX_API OMX_ERRORTYPE OMX_APIENTRY OMX_FreeHandle(OMX_IN  OMX_HANDLETYPE hComponent)
{
	if (!hComponent) {
		return OMX_ErrorBadParameter;
	}
	delete (ComponentBase *)(((OMX_COMPONENTTYPE *)hComponent)->pComponentPrivate);
	return OMX_ErrorNone;
}

OMX_ERRORTYPE omx_AllocateBuffer(OMX_IN OMX_HANDLETYPE hComponent, OMX_INOUT OMX_BUFFERHEADERTYPE **ppBuffer, OMX_IN OMX_U32 nPortIndex, OMX_IN OMX_PTR pAppPrivate, OMX_IN OMX_U32 nSizeBytes)
{
	OMX_BUFFERHEADERTYPE *buf;

	if (!hComponent || !ppBuffer|| ((OMX_S32)nPortIndex < 0) || ((OMX_S32)nSizeBytes <= 0)) {
		return OMX_ErrorBadParameter;
	}
	*ppBuffer = buf = (OMX_BUFFERHEADERTYPE *)malloc(sizeof(OMX_BUFFERHEADERTYPE));
	if (!buf) {
		return OMX_ErrorInsufficientResources;
	}
	memset(buf, 0, sizeof(OMX_BUFFERHEADERTYPE));
	SET_SIZE_VERSION(buf);
	if (!(buf->pBuffer = (OMX_U8 *)malloc(nSizeBytes))) {
		free(buf);
		*ppBuffer = 0;
		return OMX_ErrorInsufficientResources;
	}
	buf->nAllocLen = nSizeBytes;
	buf->pAppPrivate = pAppPrivate;
	buf->nInputPortIndex = nPortIndex;
	return OMX_ErrorNone;
}

OMX_ERRORTYPE omx_FreeBuffer(OMX_IN OMX_HANDLETYPE hComponent, OMX_IN OMX_U32 nPortIndex, OMX_IN OMX_BUFFERHEADERTYPE *pBuffer)
{
	if (!hComponent || ((OMX_S32)nPortIndex < 0) || !pBuffer || !pBuffer->pBuffer) {
		return OMX_ErrorBadParameter;
	}
	free(pBuffer->pBuffer);
	memset(pBuffer, 0, sizeof(*pBuffer));
	free(pBuffer);
	return OMX_ErrorNone;
}

OMX_ERRORTYPE omx_SetCallbacks(OMX_IN OMX_HANDLETYPE hComponent, OMX_IN OMX_CALLBACKTYPE *pCallbacks, OMX_IN OMX_PTR pAppData)
{
	ComponentBase *comp;
	if (!hComponent || !pCallbacks) {
		return OMX_ErrorBadParameter;
	}
	comp = (ComponentBase *)(((OMX_COMPONENTTYPE *)hComponent)->pComponentPrivate);
	comp->SetCallback(pCallbacks);
	return OMX_ErrorNone;
}


OMX_ERRORTYPE omx_GetState(OMX_IN OMX_HANDLETYPE hComponent, OMX_OUT OMX_STATETYPE* pState)
{
	ComponentBase *comp;
	if (!hComponent || !pState) {
		return OMX_ErrorBadParameter;
	}
	comp = (ComponentBase *)(((OMX_COMPONENTTYPE *)hComponent)->pComponentPrivate);
	return comp->GetState(pState);
}

OMX_ERRORTYPE omx_SendCommand(OMX_IN OMX_HANDLETYPE hComponent, OMX_IN OMX_COMMANDTYPE Cmd, OMX_IN OMX_U32 nParam1, OMX_IN OMX_PTR pCmdData)
{
	ComponentBase *comp;
	if (!hComponent) {
		return OMX_ErrorBadParameter;
	}
	comp = (ComponentBase *)(((OMX_COMPONENTTYPE *)hComponent)->pComponentPrivate);
	return comp->SendCommand(Cmd, nParam1, pCmdData);
}

OMX_ERRORTYPE omx_EmptyThisBuffer(OMX_IN OMX_HANDLETYPE hComponent, OMX_IN OMX_BUFFERHEADERTYPE* pBuffer)
{
	ComponentBase *comp;
	if (!hComponent || !pBuffer) {
		return OMX_ErrorBadParameter;
	}
	comp = (ComponentBase *)(((OMX_COMPONENTTYPE *)hComponent)->pComponentPrivate);
	return comp->EmptyThisBuffer(pBuffer);
}

OMX_ERRORTYPE omx_FillThisBuffer(OMX_IN OMX_HANDLETYPE hComponent, OMX_IN OMX_BUFFERHEADERTYPE* pBuffer)
{
	ComponentBase *comp;
	if (!hComponent || !pBuffer) {
		return OMX_ErrorBadParameter;
	}
	comp = (ComponentBase *)(((OMX_COMPONENTTYPE *)hComponent)->pComponentPrivate);
	return comp->FillThisBuffer(pBuffer);
}

OMX_ERRORTYPE omx_GetComponentVersion(OMX_IN OMX_HANDLETYPE hComponent, OMX_OUT OMX_STRING pComponentName, OMX_OUT OMX_VERSIONTYPE* pComponentVersion, OMX_OUT OMX_VERSIONTYPE* pSpecVersion, OMX_OUT OMX_UUIDTYPE* pComponentUUID)
{
	ComponentBase *comp;
	OMX_COMPONENTTYPE *cp;
	const ComponentList *id;

	if (!hComponent) {
		return OMX_ErrorBadParameter;
	}
	cp = (OMX_COMPONENTTYPE *)hComponent;
	comp = (ComponentBase *)(cp->pComponentPrivate);
	id = comp->Id();
	if (pComponentName) {
		strncpy(pComponentName, id->name, 128);
	}
	if (pComponentVersion) {
		*pComponentVersion = cp->nVersion;
	}
	if (pSpecVersion) {
		*pComponentVersion = cp->nVersion;
	}
	if (pComponentUUID) {
		strncpy((char *)pComponentUUID, id->uuid, sizeof(*pComponentUUID));
	}
	return OMX_ErrorNone;
}

OMX_API OMX_ERRORTYPE OMX_APIENTRY OMX_SetupTunnel(OMX_IN  OMX_HANDLETYPE hOutput, OMX_IN  OMX_U32 nPortOutput, OMX_IN  OMX_HANDLETYPE hInput, OMX_IN  OMX_U32 nPortInput)
{
	OMX_COMPONENTTYPE *out, *in;
	OMX_ERRORTYPE err;
	/* FIXME */
	OMX_TUNNELSETUPTYPE tunnel_in = {
		OMX_PORTTUNNELFLAG_READONLY,
		OMX_BufferSupplyOutput
	};
	OMX_TUNNELSETUPTYPE tunnel_out = {
		0,
		OMX_BufferSupplyOutput
	};

	if (!hOutput || !hInput) {
		return OMX_ErrorBadParameter;
	}
	Core->PushTunnel(hOutput, nPortOutput, hInput, nPortInput);
	out = (OMX_COMPONENTTYPE *)hOutput;
	in = (OMX_COMPONENTTYPE *)hInput;
	err = out->ComponentTunnelRequest(out, nPortOutput, in, nPortInput, &tunnel_out);
	if (err != OMX_ErrorNone) {
		return err; /* OMX_ErrorTunnelingUnsupported */
	}
	if (out) {
	}
	if (tunnel_in.eSupplier != OMX_BufferSupplyInput) {
		tunnel_out.eSupplier = OMX_BufferSupplyInput;
	}
	if (tunnel_out.nTunnelFlags & OMX_PORTTUNNELFLAG_READONLY) {
	}
	err = in->ComponentTunnelRequest(in, nPortInput, out, nPortOutput, &tunnel_in);
	return OMX_ErrorNone;
}

