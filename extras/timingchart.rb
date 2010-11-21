#! /bin/ruby

require 'csv'
require "tk"

class TimingChart
  attr_accessor :width, :height, :scale, :offset

  def initialize(list, width, height)
    @list = list
    @width = width
    @height = height
    @offset = 0
    @scale = 500000000
  end

  def array_to_be_displayed(block, offset, scale)
    block.find_all {|v|
      ((offset <= v[0]) && (v[0] <= offset + scale)) || ((offset <= (v[0] + v[1])) && ((v[0] + v[1]) <= offset + scale))
    }
  end

  def reset(canvas)
    @offset = 0
    @scale = 5000000000
    draw(canvas)
  end

  def draw(canvas)
    filled_rect(canvas, 0, 0, @width - 1, @height - 1, 'white')
    raw_height = (@height - 10) /  @list.size
    pos_y = (raw_height - 20) / 2
    pos_x = 60
    @list.each { |el|
      draw_string(canvas, 0, pos_y, el[0], 'black')
      array_to_be_displayed(el[1], @offset, @scale).each { |slot|
        left = (slot[0] - @offset) * (@width - pos_x) / @scale;
        right = (slot[0] + slot[1] - @offset) * (@width - pos_x) / @scale;
        if @width <= right then
          right = @width - 1
        end
        if left < 0 then
          left = 0
        end
        filled_rect(canvas, pos_x + left, pos_y, pos_x + right, pos_y + 10, 'red', 'black')
      }
      pos_y += raw_height
    }
  end

  def line(gd, x1, y1, x2, y2, col)
    TkcLine.new(gd, x1, y1, x2, y2) {fill col}
  end

  def rect(gd, x1, y1, x2, y2, col)
    TkcRectangle.new(gd, x1, y1, x2, y2, "outline" => col)
  end

  def filled_rect(gd, x1, y1, x2, y2, col, line_col = col)
    TkcRectangle.new(gd, x1, y1, x2, y2, "fill" => col, "outline" => line_col)
  end

  def draw_string(gd, offset_x, offset_y, string, color, anc = 'nw')
    TkcText.new(gd, offset_x, offset_y) {|t|
      text string
      anchor anc
      font TkFont.new(:size => 8)
    }
  end
end

def readTiming(filename)
  timelist = {}
  starttime = {}
  CSV.open(filename, "r") do |line|
    key = line[1].to_s
    if timelist[key] == nil
      timelist[key] = []
      starttime[key] = 0
    end
    edge = line[2].to_i
    if edge == 0 then
      timelist[key] << [starttime[key], line[0].to_i - starttime[key]]
    else
      starttime[key] = line[0].to_i
    end
  end
  return timelist
end

if 1 <= ARGV.size
  filename = ARGV[0]
else
  filename = Tk.getOpenFile
end
$graph = TimingChart.new(readTiming(filename), 640, 480)

TkFrame.new {
  TkMenubutton.new(self) {
    text "File"
    menu TkMenu.new(self) {
      add "command", {"label"=>"Open",
        "command"=>proc{read_file(Tk.getOpenFile); $graph.reset($c)},
        "accelerator"=>"Ctrl-O"
      }
      add "command", {"label"=>"Exit",
        "command"=>proc{exit},
        "accelerator"=>"Ctrl-Q"
      }
    }
  }.pack("side"=>"left", 'expand' => false)
}.pack

TkFrame.new { |f|
  TkButton.new(f, 'text'=>'Quit', 'command'=>proc{exit}).pack('side' => 'left')
  TkButton.new(f, 'text'=>'Reset', 'command'=>proc{$graph.reset($c)}).pack('side' => 'left')
  TkButton.new(f, 'text'=>'Zoom In', 'command'=>proc{
                 $graph.scale *= 0.5
                 $graph.draw $c
               }).pack('side' => 'left')
  TkButton.new(f, 'text'=>'Zoom Out', 'command'=>proc{
                 $graph.scale *= 2
                 $graph.draw($c)
               }).pack('side' => 'left')
  TkButton.new(f, 'text'=>'Right', 'command'=>proc{
                 $graph.offset += $graph.scale / 8
                 $graph.draw($c)
               }).pack('side' => 'left')
}.pack

TkLabel.new(nil, 'text'=>'Timing chart').pack
$c = TkCanvas.new { |c|
  width $graph.width
  height $graph.height
  bg 'white'
  bind "ButtonPress-1", proc { |x, y|
    $mouse = [x, y]}, "%x %y"
  bind 'Motion', proc { |x, y|
    next unless $mouse
    $graph.offset = ($mouse[0] - x) * $graph.scale / 128
    $graph.offset = 0 if $graph.offset < 0
    $mouse = [x, y]
    $graph.draw($c)}, "%x %y"
  bind 'ButtonRelease-1', proc {$mouse = nil}
}.pack('side' => 'bottom', 'expand' => true)

$graph.reset($c)

Tk.mainloop

