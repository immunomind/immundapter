function RawBlock(block)
  if block.format == "html" and block.text:match("^<pre>[^<]*</pre>%s*$") then
    local code = pandoc.read(block.text, "html").blocks[1]
    code.text = code.text:gsub("\\%$", "$")
    code.classes = {"r"}
    return code
  end
end
