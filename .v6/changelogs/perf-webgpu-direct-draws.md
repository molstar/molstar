## Graphics

- WebGPU: GPU surfaces are drawn directly once their vertex count is read back, indirectly only until then. Indirect
  draws cost far more GPU time in Chrome: 3j3q (1356 Gaussian surfaces) renders and hovers at 25 ms a frame instead of
  105–540 ms, as fast as on WebGL
