#!/usr/bin/ruby

include Math

N = 8

def c(x)
  x == 0 ? 1.0 / sqrt(2) : 1.0
end

def f(x, y)
  ((x == N - 1) && (y == N - 1)) ? 1 : ((x == 0) && (y == 0)) ? 0 : 0
end

def ccos(c, u)
  cos(c * u)
end
#  cos((x * 2 + 1) * u * PI / (N * 2))

def idct_mismatch(x, y)
  sum = 0
  cx = (x * 2.0 + 1) * PI / (N * 2)
  cy = (y * 2.0 + 1) * PI / (N * 2)
  (0...N).each {|v|
    ccosy = cos(cy * v)
    ccosy = ccosy / sqrt(2.0) if v == 0
    (0...N).each {|u|
      sum += c(u) * f(u, v) * cos(cx * u) * ccosy
    }
  }
  sum * 2 / N
end

(0...N).each {|y|
  (0...N).each {|x|
    printf(" %+f,", idct_mismatch(x, y))
#    printf(" %+f,", f(x, y))
  }
  puts
}

