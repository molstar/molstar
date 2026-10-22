## Graphics: labels

- Labels, measurements and MVS labels share one GPU texture per font atlas instead of one per text object (WebGL and
  WebGPU), and upload it only when glyphs were added: 10 distance measurements 983.7 → 613.8 MiB of textures, an MVS
  scene with labels on 6 structures 819.3 → 613.8 MiB (41 MiB per text object saved). `TextureSpecOptions.shared` opts a
  texture of a renderable schema in; `TextureImage.version` marks in-place writes.
