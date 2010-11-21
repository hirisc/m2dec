#! /bin/env ruby

require 'strscan'

class Code
  attr_accessor :len, :pat, :run, :level
  def initialize(pattern, len, run, level)
    @pat = pattern
    @len = len
    @run = run
    @level = level
  end
end

class Vlc
  attr_accessor :run, :level, :length
  def initialize(run, level, length)
    @run = run
    if 0 < length
      level = level * 2
      level = -level + 1 if level < 0
    end
    @level = level
    @length = length
  end
  def dump
    printf("{%d, %d, %d},", @run, @level, @length)
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
def read_codes(lines, h)
  #    lines = file.read
  s = StringScanner.new(lines)
  until s.eos?
    t = s.scan(/^[10][10 ]*s/)
    if t
      t.gsub!(/\s/, "")
      run = s.scan_until(/\d+/).to_i
      level = s.scan_until(/\d+/).to_i
      pat, len = pattern_to_code(t.gsub(/s/, "0"))
      h.push(Code.new(pat, len, run, level))
      pat, len = pattern_to_code(t.gsub(/s/, "1"))
      h.push(Code.new(pat, len, run, -level))
    else
      t = s.scan(/^[10][10 ]*[BE]/)
      if t
        t.gsub!(/\s/, "")
        #          p t[t.size - 1].chr
        pat, len = pattern_to_code(t.gsub(/[BE]/, ""))
        level = t[t.size - 1].chr == "B" ? -1 : 0
        h.push(Code.new(pat, len, -1, level))
        #          p len
      else
        s.scan_until(/[.\s]/)
      end
    end
  end
end

def extract_table(h, tab, bit_num)
  next_bound = 1 << bit_num
  ex_tables = Array.new(bit_num)
  h.each { |elem|
    shortbit = bit_num - elem.len
    if 0 <= shortbit
      # full code can be stored in this table
      idx = elem.pat << shortbit
      run = elem.run
      level = elem.level
      len = elem.len
      (1 << shortbit).times { |i|
        tab[idx + i] = Vlc.new(run, level, len)
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
      tab[i] = Vlc.new(tab.size, maxlen, 0)
      tab.concat(ex_tab)
    else
    end
  }
end

def print_table(tab, bit_num)
  printf("/** This is generated from standard document using makevld.rb. */\n")
  printf("static const vlc_dct_t bit%d[%d] = {\n\t", bit_num, tab.size)
  dummy = Vlc.new(-1, -1, -1)
  tab.each_with_index { |elem, i|
    (elem ? elem : dummy).dump
    printf((i & 3) == 3 ? "\n\t" : " ")
  }
  printf("};\n")
end

def devide_tables(filename, bit_num)
  h = Array.new
  File.open(filename) { |file|
    read_codes(file.read, h)
  }
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

