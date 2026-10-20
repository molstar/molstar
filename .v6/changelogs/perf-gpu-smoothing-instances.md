## Graphics: faster surface recoloring

- GPU gaussian surfaces (`tryUseGpu`, WebGL2 and WebGPU) on instanced structures smooth overpaint, transparency,
  emissive and substance only over the instances the layer reaches and their neighbours, as the mesh path does, with
  the same grid (identical in the channels the shaders read). One chain, headless Chrome on an Apple M4 (median of
  three runs): 1tup ×32 on WebGL2 (resolution 1) draws 13 of 32 instances, overpaint 180 → 95 ms, transparency 188 →
  95 ms; 1tup ×192 on WebGPU (`low`) draws 18 of 192, overpaint 102 → 59 ms, transparency 63 → 45 ms. At ×864 (18 of
  864 drawn) the full-size grid dominates and the gain is within the noise.
