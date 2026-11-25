## Graphics: illumination input culled on the GPU

- Illumination's tracing input is culled on the GPU with occlusion in two phases, as the opaque pass (WebGPU, objects
  culled on the GPU; WebGL unchanged; the back depth of the auto thickness not by occlusion): a moving frame on
  mesoscale models 1.4-1.5× faster (MG_1184 zoomed, 2560×1440: 72.7 → 51.0 ms; PetWorld zoomed, 1280×720: 8.8 →
  5.7 ms).
