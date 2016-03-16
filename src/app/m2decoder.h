#ifndef _M2DECODER_H_
#define _M2DECODER_H_

#include <ctype.h>
#include "frames.h"
#include "m2d.h"
#include "mpeg2.h"
#include "h264.h"
#ifdef __cplusplus
extern "C" {
#endif
extern const m2d_func_table_t * const h265d_func;
#ifdef __cplusplus
}
#endif

#include "mpeg_demux.h"
#include <deque>

typedef m2d_frame_t Frame;

class M2Decoder {
public:
	typedef enum {
		MODE_MPEG2,
		MODE_MPEG2PS,
		MODE_H264,
		MODE_H265,
		MODE_NONE
	} type_t;
	typedef std::pair<const uint8_t *, int> header_data_t;
	typedef std::deque<header_data_t> header_data_list_t;
	M2Decoder(type_t codec_mode, int outbuf, int (*reread_file)(void *arg), void *reread_arg)
		: frames_(0),
		outbuf_(outbuf), reread_file_(reread_file), reread_arg_(reread_arg),
		context_(0), func_(0) {
		set_codec(codec_mode);
	}
	~M2Decoder() {
		if (frames_) {
			delete frames_;
		}
		if (context_) {
			delete[] context_;
		}
	}
	void change_codec(type_t codec_mode) {
		if (context_) {
			delete[] context_;
			context_ = 0;
		}
		set_codec(codec_mode);
	}
	void SetFrames(void *id) {
		m2d_info_t info;
		func_->get_info(context_, &info);
		int width = (info.src_width + 15) & ~15;
		int height = (info.src_height + 15) & ~15;
		int luma_len = width * height;
		int bufnum = outbuf_ + info.frame_num + (((codec_mode_ == MODE_H264) || (codec_mode_ == MODE_H265)) ? 16 : 0);
		if (codec_mode_ == MODE_H264) {
			if (H264D_MAX_FRAME_NUM < bufnum) {
				bufnum = H264D_MAX_FRAME_NUM;
			}
		} else if (MAX_FRAME_NUM < bufnum) {
			bufnum = MAX_FRAME_NUM;
		}
		if (frames_) {
			if (frames_->sufficient(bufnum, luma_len, info.additional_size)) {
				frames_->set_id(id);
				return;
			}
			delete frames_;
		}
		fprintf(stderr, "%d x %d x %d\n", info.src_width - info.crop[0] - info.crop[1], info.src_height - info.crop[2] - info.crop[3], info.frame_num);
		frames_ = new Frames(width, height, bufnum, info.additional_size, id);
		if (func_->set_frames(context_, bufnum, frames_->aligned(), frames_->second(), info.additional_size) < 0) {
			throw std::bad_exception();
		}
	}
	type_t codec_mode() {
		return codec_mode_;
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
	int skip_frames(const uint8_t *indata, int indata_bytes, int skip_frm, int& skipped_bytes, header_data_list_t& headers) {
		int pos = 0;
		int skipped_frm = 0;
		int skipped_frm_key = 0;
		const uint8_t *indata_key = 0;
		while (pos < indata_bytes) {
			int read_bytes = m2d_next_start_code(indata + pos, indata_bytes - pos);
			if (read_bytes < 0) {
				break;
			}
			pos += read_bytes;
			bool is_keyframe, is_header;
			if (is_h264frame_head(indata + pos, indata_bytes - pos, is_keyframe, is_header)) {
				if (is_keyframe) {
					indata_key = indata + pos - 3;
					skipped_frm_key = skipped_frm;
				}
				if (skip_frm < ++skipped_frm) {
					break;
				}
			} else if (is_header) {
				int size = m2d_next_start_code(indata + pos, indata_bytes - pos);
				headers.push_back(header_data_t(indata + pos - 3, size));
			}
		}
		while (!headers.empty()) {
			headers.push_back(header_data_t(reinterpret_cast<const uint8_t *>(0), 0));
			func()->decode_picture(context());
		}
		if (indata_key) {
			skipped_bytes = indata_key - indata;
			return skipped_frm_key;
		} else {
			return -1;
		}
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
private:
	Frames *frames_;
	type_t codec_mode_;
	int outbuf_;
	int (*reread_file_)(void *arg);
	void *reread_arg_;
	byte_t *context_;
	const m2d_func_table_t *func_;
	pes_demuxer_t demux_;
	void set_codec(type_t codec_mode) {
		codec_mode_ = codec_mode;
		switch (codec_mode) {
		case MODE_NONE:
			return;
		case MODE_MPEG2:
			/* FALLTHROUGH */
		case MODE_MPEG2PS:
			func_ = m2d_func;
			break;
		case MODE_H264:
			func_ = h264d_func;
			break;
		case MODE_H265:
			func_ = h265d_func;
			break;
		}
		context_ = new byte_t[func_->context_size];
		func_->init(context_, -1, header_callback, this);
		memset(&demux_, 0, sizeof(demux_));
		if (codec_mode == MODE_MPEG2PS) {
			mpeg_demux_init(&demux_, reread_file_, reread_arg_);
			dec_bits_set_callback(func_->stream_pos(context_), reread_packet, this);
		} else {
			dec_bits_set_callback(func_->stream_pos(context_), reread_file_, reread_arg_);
		}
	}
	int reread_packet_impl() {
		pes_demuxer_t *dmx = demuxer();
		const byte_t *packet;
		int packet_size;
		void *id;
		packet = mpeg_demux_get_video(dmx, &packet_size, &id);
		if (packet) {
			dec_bits_set_data(stream(), packet, (size_t)packet_size, id);
			return 0;
		} else {
			return -1;
		}
	}
	static int reread_packet(void *arg) {
		return ((M2Decoder *)arg)->reread_packet_impl();
	}
	static int header_callback(void *arg, void *id) {
		((M2Decoder *)arg)->SetFrames(id);
		return 0;
	}
	bool is_h264frame_head(const uint8_t *indata, int indata_bytes, bool& is_keyframe, bool& is_header) {
		if (indata_bytes < 2) {
			return false;
		}
		int nal_type = indata[0] & 31;
		is_keyframe = (nal_type == SLICE_IDR_NAL);
		is_header = (nal_type == SPS_NAL) || (nal_type == PPS_NAL);
		return (indata[1] & 128) && ((nal_type == SLICE_IDR_NAL) || (nal_type == SLICE_NONIDR_NAL));
	}
};

static void to_lower_ext(const char *src, char *dst, int len) {
	for (int i = 0; i < len; ++i) {
		const char s = src[i];
		if (s == '\0') {
			dst[i] = s;
			break;
		}
		dst[i] = tolower(s);
	}
}

static M2Decoder::type_t detect_file(const char *filename) {
	static const struct {
		M2Decoder::type_t type;
		const char *ext;
	} ext_map[] = {
		{ M2Decoder::MODE_MPEG2, "m2v"},
		{ M2Decoder::MODE_MPEG2PS, "vob"},
		{ M2Decoder::MODE_H264, "264"},
		{ M2Decoder::MODE_H264, "jsv"},
		{ M2Decoder::MODE_H265, "265"},
		{ M2Decoder::MODE_NONE, ""}
	};
	char ext[16];
	const char *ext_p = strrchr(filename, '.');
	if (ext_p++ != 0) {
		to_lower_ext(ext_p, ext, sizeof(ext));
		int i = -1;
		while (ext_map[++i].type != M2Decoder::MODE_NONE) {
			if (strcmp(ext_map[i].ext, ext) == 0) {
				return ext_map[i].type;
			}
		}
	}
	return M2Decoder::MODE_MPEG2;
}

#endif /* _M2DECODER_H_*/
