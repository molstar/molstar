## Graphics: memory

- [Breaking] Lines keep their data once per line, `lineData` (read by the shaders as `tLine`, 32 bytes), not per
  vertex with an index buffer of their own (about 170 bytes; strips too): 4V6X as lines, CPU 110 → 21 MB, GPU vertex
  attributes 94 → 0 MB (21 MB of data instead); merged lines 2× faster on WebGPU on an integrated GPU (31.2 → 15.3 ms)
