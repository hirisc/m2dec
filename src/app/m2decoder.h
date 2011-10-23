#ifndef _M2DECODER_H_
#define _M2DECODER_H_

#include "frames.h"
#include "mpeg2.h"
#include "h264.h"
#include "mpeg_demux.h"

typedef m2d_frame_t Frame;

class M2Decoder {
	Frames *frames_;
	int codec_mode_;
	int outbuf_;
	byte_t *context_;
	const m2d_func_table_t *func_;
	pes_demuxer_t demux_;
	int reread_packet_impl() {
		pes_demuxer_t *dmx = demuxer();
		const byte_t *packet;
		int packet_size;
		packet = mpeg_demux_get_video(dmx, &packet_size);
		if (packet) {
			dec_bits_set_data(stream(), packet, (size_t)packet_size);
			return 0;
		} else {
			return -1;
		}
	}
	static int reread_packet(void *arg) {
		return ((M2Decoder *)arg)->reread_packet_impl();
	}
	static int header_callback(void *arg, int id) {
		((M2Decoder *)arg)->SetFrames(id);
		return 0;
	}
public:
	enum {
		MODE_MPEG2,
		MODE_MPEG2PS,
		MODE_H264
	};
	M2Decoder(int codec_mode, int outbuf, int (*reread_file)(void *arg), void *reread_arg)
		: frames_(0), codec_mode_(codec_mode), outbuf_(outbuf), context_(0) {
		switch (codec_mode) {
		case MODE_MPEG2:
			/* FALLTHROUGH */
		case MODE_MPEG2PS:
			func_ = m2d_func;
			break;
		case MODE_H264:
			func_ = h264d_func;
			break;
		}
		context_ = new byte_t[func_->context_size];
		func_->init(context_, -1, header_callback, this);
		memset(&demux_, 0, sizeof(demux_));
		if (codec_mode == MODE_MPEG2PS) {
			mpeg_demux_init(&demux_, reread_file, reread_arg);
			dec_bits_set_callback(func_->stream_pos(context_), reread_packet, this);
		} else {
			dec_bits_set_callback(func_->stream_pos(context_), reread_file, reread_arg);
		}
	}
	~M2Decoder() {
		if (frames_) {
			delete frames_;
		}
		delete[] context_;
	}
	void SetFrames(int id) {
		m2d_info_t info;
		func_->get_info(context_, &info);
		int width = (info.src_width + 15) & ~15;
		int height = (info.src_height + 15) & ~15;
		int luma_len = width * height;
		int bufnum = outbuf_ + info.frame_num + ((codec_mode_ == MODE_H264) ? 16 : 0);
		if (codec_mode_ == MODE_H264) {
			if (H264D_MAX_FRAME_NUM < bufnum) {
				bufnum = H264D_MAX_FRAME_NUM;
			}
		} else if (MAX_FRAME_NUM < bufnum) {
			bufnum = MAX_FRAME_NUM;
		}
		if (frames_) {
			if (frames_->sufficient(bufnum, luma_len, info.additional_size)) {
				return;
			}
			delete frames_;
		}
		fprintf(stderr, "%d x %d x %d\n", info.src_width, info.src_height, info.frame_num);
		frames_ = new Frames(width, height, bufnum, info.additional_size);
		func_->set_frames(context_, bufnum, frames_->aligned(), frames_->second(), info.additional_size);
	}
	void *context() {
		return context_;
	}
	const m2d_func_table_t *func() {
		return func_;
	}
	dec_bits *stream() {
		return func_->stream_pos(context_);
	}
	pes_demuxer_t *demuxer() {
		return &demux_;
	}
	int decode(void *obj, void (*post_dst)(void *, m2d_frame_t&), bool emptify_mode) {
		m2d_frame_t frm;
		int err = -1;
		while (func()->peek_decoded_frame(context(), &frm, 0) <= 0) {
			err = func()->decode_picture(context());
			if (err < 0) {
				while (func()->peek_decoded_frame(context(), &frm, 1)) {
					post_dst(obj, frm);
					func()->get_decoded_frame(context(), &frm, 1);
				}
				return err;
			}
		}
		do {
			func()->get_decoded_frame(context(), &frm, 0);
			post_dst(obj, frm);
		} while (emptify_mode && (0 < func()->peek_decoded_frame(context(), &frm, 0)));
		return func()->decode_picture(context());
	}
	void decode_residual(void *obj, void (*post_dst)(void *, m2d_frame_t&)) {
		m2d_frame_t frm;
		while (0 < func()->peek_decoded_frame(context(), &frm, 1)) {
			post_dst(obj, frm);
			func()->get_decoded_frame(context(), &frm, 1);
		}
	}
};

#endif /* _M2DECODER_H_*/
