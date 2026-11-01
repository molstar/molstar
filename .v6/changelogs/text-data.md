## Graphics: memory

- [Breaking] Text keeps its data once per glyph, background and tether, `textData` (read by the shaders as `tText`,
  64 bytes), not per vertex with an index buffer of its own (about 170 bytes)
