#ifndef _OMX_PORT_H_
#define _OMX_PORT_H_

#include <stdio.h>
#include <assert.h>

const unsigned int MAX_QUEUE = 5;

class Port {
	OMX_STATETYPE state_;
	deque<OMX_BUFFERHEADERTYPE *> queue_;
	OMX_PARAM_PORTDEFINITIONTYPE port_def_;
	int SetFormat(OMX_PORTDOMAINTYPE domain, OMX_PTR domain_def);
public:
	Port() : state_(OMX_StateInvalid) {
	}
	Port(int portnum, OMX_DIRTYPE is_out, OMX_PORTDOMAINTYPE domain, OMX_PTR domain_def) : state_(OMX_StateInvalid) {
		Init(portnum, is_out, domain, domain_def);
	}
	int Init(int portnum, OMX_DIRTYPE is_out, OMX_PORTDOMAINTYPE domain, OMX_PTR domain_def);
	void enable(bool on);
	void SetState(OMX_STATETYPE state);
	void populate(bool on);
	OMX_PARAM_PORTDEFINITIONTYPE &port_def() {
		return port_def_;
	}
	deque<OMX_BUFFERHEADERTYPE *> &queue() {
		return queue_;
	}
};

#endif /* _OMX_PORT_H_ */

