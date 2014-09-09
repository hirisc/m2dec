class PrintableArray < Array
  def elem(v)
    sprintf("0x%02x", v)
  end
  def dump(line_max, indent)
    print("\t" * indent, "{\n")
    max = self.length - 1
    prefix = "\t" * (indent + 1)
    self.each_with_index do |v, i|
      e = elem(v)
      mod = i % line_max
      if mod == line_max - 1 then
        print(" ", e, ",\n")
      elsif mod == 0 then
        print(prefix, e, ",")
      elsif i == max then
        print(" ", e, "\n")
      else
        print(" ", e, ",")
      end
    end
    print("\n") if (line_max < self.length) && ((self.length % line_max) != 0)
    print("\t" * indent, "}")
    self
  end
end
