#! /bin/env ruby

#require 'strscan'
require 'csv'

#
# code consists of: value, bitlen
#   bitlen > 0:  valid code, complete
#   bitlen < 0:  valid code, extra look-up is needed.
#     rest bitlen = -bitlen
#     offset to next LUT = value
#   bitlen == 0: invalid code(ERROR)
#

class Code
  attr_accessor :len, :pat, :level
  def initialize(pattern, len, level)
    @pat = pattern
    @len = len
    @level = level
  end
end

class Vlc
  attr_accessor :level, :length
  def initialize(level, length)
    @level = level
    @length = length
  end
  def dump
    printf("{%d, %d},", @level, @length)
  end
end

# VLD String to integer code and its length.
def pattern_to_code(pattern)
  pat = pattern.split(//)
  code = 0
  len = pattern.size
  pat.each { |e|
    code = (code << 1) | e.to_i
  }
  return code, len
end

# read VLD codes from file and convert to array
def read_codes(codes, h)
  codes.each { |line|
    t = line[0].gsub(/\s/, "")
    if line[1] == "B"
      val = -1
    elsif line[1] == "E"
      val = 0
    elsif line[1] =~ /\d/
      val = line[1].to_i
    else
      p "Syntax error."
    end

    if t =~ /\A[01]+s\z/
      p "not supported."
    elsif t =~ /\A[01]+/
      pat, len = pattern_to_code(t)
      h.push(Code.new(pat, len, val))
    else
      p "Syntax error."
    end
  }
end

def extract_table(h, tab, bit_num)
  next_bound = 1 << bit_num
  ex_tables = Array.new(bit_num)
  h.each { |elem|
    shortbit = bit_num - elem.len
    if 0 <= shortbit
      # full code can be stored in this table
      idx = elem.pat << shortbit
      level = elem.level
      len = elem.len
      (1 << shortbit).times { |i|
        tab[idx + i] = Vlc.new(level, len)
      }
    else
      # excess bits, extra table is needed
      idx = elem.pat >> (-shortbit)
      elem.len = -shortbit
      elem.pat &= ~(((1 << bit_num) - 1) << -shortbit)
      ex_tables[idx] = Array.new if !ex_tables[idx]
      ex_tables[idx].push(elem)
    end
  }

  # create and append extra tables recursively
  ex_tables.each_with_index { |elem, i|
    if elem
      maxlen = 0
      elem.each { |code|
        maxlen = code.len if maxlen < code.len
      }
      limited_maxlen = bit_num < maxlen ? bit_num : maxlen
      ex_tab = Array.new(1 << limited_maxlen)
      extract_table(elem, ex_tab, limited_maxlen)
      tab[i] = Vlc.new(tab.size, -maxlen)
      tab.concat(ex_tab)
    else
    end
  }
end

def print_table(tab, bit_num)
  printf("/** This is generated from standard document using vldbuild.rb. */\n")
  printf("static const vlc_t bit%d[%d] = {\n\t", bit_num, tab.size)
  dummy = Vlc.new(0, 0)
  tab.each_with_index { |elem, i|
    (elem ? elem : dummy).dump
    printf((i & 3) == 3 ? "\n\t" : " ")
  }
  printf("};\n")
end

def devide_tables(filename, bit_num)
  a = CSV.readlines(filename)
  h = Array.new
  read_codes(a, h)
  tab = Array.new(1 << bit_num)
  extract_table(h, tab, bit_num)
  print_table(tab, bit_num)
end

if !ARGV[0]
    printf("Usage: makevld.rb infile [bit-length]\n")
    printf("\tDevide DCT variable length code specified by infile.\n")
    printf("\tinfile shall be copy&pasted from MPEG-2 standard document.\n")
else
  devide_tables(ARGV[0], ARGV[1] ? ARGV[1].to_i : 6)
end

