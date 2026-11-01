## Graphics: culling, viewports and compiling ahead

- Two-phase occlusion culling of objects culled on the GPU, the sphere raster included (`hiZ.twoPhase`, on by default):
  what the previous frame's pyramid left out is tested against the frame's depth so far and drawn in the same frame, so
  nothing is missing while the camera moves
- A Canvas3D frame extending past the canvas is cropped by Canvas3D, the camera projecting the part drawn; WebGPU viewports
  outside the target are clipped where the device rejects them
- Representations compile what they draw while they build; frames after an update wait for the programs they need
- Several Canvas3D instances sharing a context on WebGPU: a frame of one keeps the image of the others
