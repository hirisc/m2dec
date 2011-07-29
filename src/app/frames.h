#ifndef _FRAMES_H_
#define _FRAMES_H_

#include <vector>
#include <functional>
#include <algorithm>
#include "m2d.h"

class Frames {
	size_t luma_len_;
	typedef std::vector<m2d_frame_t> frames_type;
	frames_type frame_, aligned_;
	std::vector<uint8_t> second_;
	struct Create : std::binary_function<m2d_frame_t, int, void> {
		void operator()(m2d_frame_t& frm, int luma_len) const {
			frm.luma = new uint8_t[luma_len + 15];
			frm.chroma = new uint8_t[(luma_len >> 1) + 15];
		}
	};
	struct Delete {
		void operator()(m2d_frame_t& frm) const {
			delete[] frm.luma;
			delete[] frm.chroma;
		}
	};
	static uint8_t *align16(uint8_t *src) {
		return (uint8_t *)(((uintptr_t)src + 15) & ~15);
	}
	void align_frame() {
		for (int i = 0; i < frame_.size(); ++i) {
			aligned_[i].luma = align16(frame_[i].luma);
			aligned_[i].chroma = align16(frame_[i].chroma);
		}
	}
public:
	Frames(int width, int height, int num_mem, int second_len)
		: luma_len_(((width + 15) & ~15) * ((height + 15) & ~15)),
		frame_(num_mem),
		aligned_(num_mem),
		second_(second_len) {
		if (num_mem <= 0 || luma_len_ <= 0) {
			return;
		}
		std::for_each(frame_.begin(), frame_.end(), std::bind2nd(Create(), luma_len_));
		align_frame();
	}
	~Frames() {
		if (frame_.empty()) {
			return;
		}
		std::for_each(frame_.begin(), frame_.end(), Delete());
		frame_.clear();
	}
	m2d_frame_t *aligned() {
		return &aligned_[0];
	}
	uint8_t *second() {
		return &second_[0];
	}
	bool sufficient(size_t bufnum, int luma_len, size_t additional_size) {
		return (bufnum <= frame_.size())
			&& (luma_len <= luma_len_)
			&& (additional_size <= second_.size());
	}
};

#endif /* _FRAMES_H_ */
