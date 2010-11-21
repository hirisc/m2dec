#!/usr/bin/ruby

#MB_SPACIAL = 0x01  ==> 0x20
#MB_INTRA = 0x02 ==> 0x04
#MB_PATTERN = 0x04 ==> 0x08
#MB_BACKWARD = 0x08 ==> 0x02
#MB_FORWARD = 0x10 ==> 0x01
#MB_QUANT = 0x20 ==> 0x10

def bit(val, b)
  (val & (1 << b)) != 0 ? 1 : 0
end

def convert_mbtype(mbtype)
  sp = bit(mbtype, 0)
  intra = bit(mbtype, 1)
  pat = bit(mbtype, 2)
  backw = bit(mbtype, 3)
  forw = bit(mbtype, 4)
  quant = bit(mbtype, 5)

  (sp << 5) | (quant << 4) | (pat << 3) | (intra << 2) | (backw << 1) | forw
end

def convert_line(line)
  val = line.split(',')
  printf "%s, %d\r\n", val[0], convert_mbtype(val[1].to_i)
end


while (line = gets) do
  convert_line(line)
end

