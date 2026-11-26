## Headless rendering

- `HeadlessScreenshotHelperOptions.graphics`: a headless plugin or screenshot helper draws with an existing WebGPU
  context (its drawing buffer resized to the canvas size) and leaves it when disposed: several headless renderers in
  turn on one device, its pipelines compiled once (the headless WebGPU tests share one a worker: their summed time
  350 s → about 215 s)
