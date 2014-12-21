#!/usr/bin/env ruby

require './printarray.rb'

# generate LUT to detect filter mode for intra prediction.


def min_dist(dir)
  [(dir - 26).abs, (dir - 10).abs].min
end

class PrintableArray
  def elem(v)
    sprintf("%d", v)
  end
end

THRESHOLD = [7, 1, 0].reverse
=begin
PrintableArray.new((2..34).map {|dir|
  dist = min_dist(dir)
  THRESHOLD.reduce(0) {|i, thr|
    i * 2 + ((thr < dist) ? 1 : 0)
  }
}).dump(16, 0)
print ";\n"
#PrintableArray.new(
(2..34).each {|dir|
  dist = min_dist(dir)
  p THRESHOLD.map {|thr|
    thr < dist
  }
}
=end

ANGLE = [32, 26, 21, 17, 13, 9, 5, 2, 0, -2, -5, -9, -13, -17, -21, -26, -32, -26, -21, -17, -13, -9, -5, -2, 0, 2, 5, 9, 13, 17, 21, 26, 32]

INV_ANGLE = {11 => -4096, 12 => -1638, 13 => -910, 14 => -630, 15 => -482, 16 => -390, 17 => -315, 18 => -256, 19 => -315, 20 => -390, 21 => -482, 22 => -630, 23 => -910, 24 => -1638, 25 => -4096}

def build_ref18_34(dir, size)
  ref = {}
  ang = ANGLE[dir - 2]
  (0..size).each {|x|
    ref[x] = [-1, x - 1]
  }
  if ang < 0 then
    lower = (size * ang) >> 5
    if lower < -1 then
      inv_ang = INV_ANGLE[dir]
      (lower..-1).each {|x|
        ref[x] = [((x * inv_ang + 128) >> 8) - 1, -1]
      }
    end
  else
    ((size + 1)..(size * 2)).each {|x|
      ref[x] = [-1, x - 1]
    }
  end
  ref
end

def build_ref2_17(dir, size)
  ref = {}
  (0..size).each {|x|
    ref[x] = [x - 1, -1]
  }
  ang = ANGLE[dir - 2]
  if ang < 0 then
    lower = (size * ang) >> 5
    if lower < -1 then
      inv_ang = INV_ANGLE[dir]
      (lower..-1).each {|x|
        ref[x] = [-1, ((x * inv_ang + 128) >> 8) - 1]
      }
    end
  else
    ((size + 1)..(size * 2)).each {|x|
      ref[x] = [x - 1, -1]
    }
  end
  ref
end

def fact_idx(x, ang)
    base = (x + 1) * ang
    return (base & 31), (base >> 5)
end

def position_table(pos, xy)
  pos.uniq!
  extra = pos.find_all {|e|
    0 <= e[xy]
  }
  run = (pos - extra).map {|e|
    e[xy ^ 1]
  }
  return [extra.length, extra.map {|e|
    e[xy]
  }.reverse, run[0], run.length].flatten
end

def extract(pos, xy, fact, inc)
  pos_tbl = position_table(pos, xy)
  inc[0] = pos_tbl[0]
  return [pos_tbl, fact, inc]
end

def lookup_top(dir, size)
  ref = build_ref18_34(dir, size)
  ang = ANGLE[dir - 2]
  fact = []
  pos = []
  inc = []
  prev_idx = ang >> 5
  (0...size).each {|y|
    fct, idx = fact_idx(y, ang)
    fact << fct
    inc << idx - prev_idx
    prev_idx = idx
    (0...size).each {|x|
      pos << ref[x + idx + 1]
      if fct != 0 then
        pos << ref[x + idx + 2]
      end
    }
  }
  extract(pos, 0, fact, inc)
end

def lookup_left(dir, size)
  ref = build_ref2_17(dir, size)
  ang = ANGLE[dir - 2]
  fact = []
  inc = []
  idx = []
  (0...size).each {|x|
    f, i = fact_idx(x, ang)
    fact << f
    inc << (idx.empty? ? 0 : i - idx.last)
    idx << i
  }
  pos = []
  (0...size).each {|y|
    (0...size).each {|x|
      pos << ref[y + idx[x] + 1]
      if fact[x] != 0 then
        pos << ref[y + idx[x] + 2]
      end
    }
  }
  extract(pos, 1, fact, inc)
end

def dump_row(elm, num, br = true)
  eol = elm.length - 1
  print "\t" * num, br ? "{" : ""
  elm.each_with_index {|e, x|
    if e.kind_of?(Array) then
      print "{"
      e.each {|v|
        print v, ","
      }
      print "}"
    else
      print e
    end
    if x == eol then
      print (br ? "}" : ""), (0 < num ? "," : ""), "\n"
    else
      print ", "
    end
  }
end

def dump_tbl(tbl)
  tail = tbl.length - 1
  tbl.each_with_index {|row, y|
    dump_row(row)
  }
end

def strtorange(str, default)
  if str then
    num = str.to_i
    if num != 0 then
      (num..num)
    else
      eval(str)
    end
  else
    default
  end
end

log2_range = strtorange(ARGV[0], (2..5))
dir_range = strtorange(ARGV[1], (2..34))

pos = []
coef = []
inc = []
dir_range.each {|dir|
  log2_range.each {|lg2|
    size = 1 << lg2
    if (18 <= dir) then
      tbl = lookup_top(dir, size)
    else
      tbl = lookup_left(dir, size)
    end
    pos << tbl[0]
    if (lg2 == 5) then
      coef << tbl[1]
      inc << tbl[2]
    end
  }
}

dir_len = dir_range.max - dir_range.min + 1
log2_len = log2_range.max - log2_range.min + 1

names = []

dir_range.each {|dir|
  log2_range.each {|lg2|
    size = (1 << lg2).to_s
    name = "intra_pred_pos_dir" + dir.to_s + "_" + size + "x" + size
    names << name
    print "static const int8_t ", name, "[] = {\n"
    dump_row(pos.shift, 1, false)
    print "};\n\n"
  }
}

print "static const int8_t* intra_pred_pos[", dir_len, "][", log2_len, "] = {\n"
dir_range.each {|dir|
  print "\t{\n"
  log2_range.each {|lg2|
    print "\t\t", names.shift, ",\n"
  }
  print "\t},\n"
}
print "};\n\n"

print "static const int8_t intra_pred_coef[", dir_len, "][2][32] = {\n"
dir_range.each {|dir|
  print "\t{\n"
  dump_row(coef.shift, 2)
  dump_row(inc.shift, 2)
  print "\t},\n"
}
print "};\n\n"
