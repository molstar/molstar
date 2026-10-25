## Graphics

- WebGPU: passes are replayed from render bundles by default (`RendererProps.drawLists`), objects changing every frame
  (a hover, a recolor, a rotation) drawn directly. 1000 separate units: 2.4 → 0.6 ms a frame still, hover and rotation
  33-66% faster; merged objects, small scenes and PetWorld unchanged
