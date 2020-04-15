require 'pp'

tool_dirs = Dir.glob("tools/*/")

tools = tool_dirs.map{|dir| [dir.split("/")[-1], dir] }

tools.each do |tool, dir|
#  pp Dir.pwd+dir
#  pp tool
  `ln -s #{Dir.pwd}/#{dir}/#{tool} bin/#{tool}`
end
