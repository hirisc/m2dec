#!/usr/bin/env ruby
require './printarray.rb'

# Base offset of CABAC sig_coeff_flag decoding except for 4x4 block.

def siginit(px, py, prev_sbf)
  case prev_sbf
  when 0
    sum = px + py
    (sum == 0) ? 2 : (sum < 3) ? 1 : 0
  when 1
    (py == 0) ? 2 : (py == 1) ? 1 : 0
  when 2
    (px == 0) ? 2 : (px == 1) ? 1 : 0
  else
    2
  end
end

array = PrintableArray.new
(0..3).each { |py|
  (0..3).each { |px|
    (0..3).each { |sbf|
      array << siginit(px, py, sbf)
    }
  }
}

print "static const int8_t sig_coeff_inc_init[", array.length, "] = "
array.dump(8, 0)
print ";\n"
