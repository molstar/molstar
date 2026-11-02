## Graphics: memory

- [Breaking] Cylinders keep their data once per cylinder, `cylinderData` (read by the shaders as `tCylinder`, 48 bytes),
  not per vertex with an index buffer of their own (about 360 bytes): 4V6X ball-and-stick bonds (650 k cylinders),
  CPU 235.7 → 31.4 MB, GPU vertex attributes 204 → 0 MB (31 MB of data instead); the opaque pass up to 2× faster on
  WebGPU (integrated GPU 29.5 → 14.4 ms), building the representation 40% faster
