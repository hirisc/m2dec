#include "fmx_common.h"
#include "fmx_port.h"

void Port::enable(bool on) {
	port_def_.bEnabled = (OMX_BOOL)on;
	if (!on) {
		populate(false);
	}
}

void Port::SetState(OMX_STATETYPE state) {
	if ((state_ == OMX_StateIdle) && (state == OMX_StateLoaded)) {
		if (port_def_.bEnabled) {
			populate(false);
		}
	} else if ((state_ == OMX_StateLoaded) && (state == OMX_StateIdle)) {
		if (port_def_.bEnabled) {
			populate(true);
		}
	}
	state_ = state;
}

void Port::populate(bool on) {
	port_def_.bPopulated = (OMX_BOOL)on;
}

static const OMX_AUDIO_PORTDEFINITIONTYPE port_template_audio = {
	"audio/PCMA",
	0, /* pnativeRender */
       OMX_FALSE, /* bFlagErrorConcealment */
       OMX_AUDIO_CodingPCM
};

static const OMX_PARAM_PORTDEFINITIONTYPE port_template = {
	sizeof(port_template),
	{OMX_VERSION_MAJOR, OMX_VERSION_MINOR, OMX_VERSION_REVISION, OMX_VERSION_STEP},
	0, /* nPortIndex */
	OMX_DirInput,
	0, /* nBufferCountActual */
	0, /* nBufferCountMin */
	0, /* nBufferSize */
	OMX_FALSE, /* bEnabled */
	OMX_FALSE, /* bPopulated */
	OMX_PortDomainAudio,
	port_template_audio
};

int Port::SetFormat(OMX_PORTDOMAINTYPE domain, OMX_PTR domain_def) {
	switch (domain) {
	case OMX_PortDomainAudio:
		port_def_.format.audio = *(OMX_AUDIO_PORTDEFINITIONTYPE *)domain_def;
		break;
	case OMX_PortDomainVideo:
		port_def_.format.video = *(OMX_VIDEO_PORTDEFINITIONTYPE *)domain_def;
		break;
	case OMX_PortDomainImage:
		port_def_.format.image = *(OMX_IMAGE_PORTDEFINITIONTYPE *)domain_def;
		break;
	case OMX_PortDomainOther:
		port_def_.format.other = *(OMX_OTHER_PORTDEFINITIONTYPE *)domain_def;
		break;
	default:
		return -1;
		break;
	}
	return 0;
}

int Port::Init(int portnum, OMX_DIRTYPE is_out, OMX_PORTDOMAINTYPE domain, OMX_PTR domain_def) {
	port_def_ = port_template;
	port_def_.nPortIndex = portnum;
	port_def_.eDir = is_out;
	port_def_.eDomain = domain;
	return SetFormat(domain, domain_def);
}

