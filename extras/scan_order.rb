#!/usr/bin/env ruby
require './printarray.rb'

# Output three zigzag-scan LUT and its inverse scan.
# For horizontal/vertical scan, inverse LUT aren't necessary.
# Horizontal LUT for all size are same as largest one.

def mux_xy(x, y, log2size)
  (y << log2size) | x
end

class DiagonalArray < PrintableArray
  def initialize(log2size)
    size = 1 << log2size
    max_length = size * size
    x = 0
    y = 0
    loop do
      if ((x < size) && (y < size)) then
        self << mux_xy(x, y, log2size)
        break if (max_length <= self.length)
      end
      y -= 1
      x += 1
      if (y < 0) then
        y = x
        x = 0
      end
    end
  end
end

class GeneratedArray < PrintableArray
  def initialize(log2size, generator)
    size = 1 << log2size
    (0...size).each do |y|
      (0...size).each do |x|
        self << generator.call(x, y, log2size)
      end
    end
  end
end

def print_lut(suffix, size, array)
  print "static const int8_t h265d_scan_order", size, "x", size, suffix, "[", size, " * ", size, "] = "
  array.dump(8, 0)
  print ";\n\n"
end

(1..3).each do |log2size|
  diag = DiagonalArray.new(log2size)
  size = 1 << log2size
  print_lut("diag_inverse", size, GeneratedArray.new(log2size, lambda {|x, y, log2size| diag.index(mux_xy(x, y, log2size))}))
  print_lut("diag", size, diag)
  print_lut("vertical", size, GeneratedArray.new(log2size, lambda {|x, y, log2size| mux_xy(y, x, log2size)}))
end
print_lut("horizontal", 8, PrintableArray.new(8 * 8) do |i| i end)
