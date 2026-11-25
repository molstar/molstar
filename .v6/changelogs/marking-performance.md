## Graphics

- Marking draws only the spheres a pass needs: unmarked ones dropped in the vertex stage of the mask, marked ones in
  that of the depth (PetWorld zoomed in, one particle highlighted: the mask 1.76 → 0.92 ms on WebGPU, 4.14 → 2.05 ms on
  WebGL2)
- Marking reads what hides marked objects from the scene's depth instead of drawing the scene again: its opaque depth
  and, where unmarked transparent objects are less transparent than the new `marking.transparentThreshold` (0.4, as
  SSAO's), its front transparent layer; more transparent objects no longer hide marked ones. The depth pass stays for
  direct volumes, path tracing and without depth textures (PetWorld zoomed in, a hover redrawn: 3.6 → 1.3 ms on WebGPU,
  9.9 → 8.1 ms on WebGL2; rotating with a hover 4.7 → 2.7 ms on WebGPU, 16.5 → 10.3 ms on WebGL2)
- Marking draws the marked instances alone where they are few, at the levels of detail they are seen at: the markers
  keep a count of marked entries per instance (`MarkerData.markerInstanceCounts`). A hover's cost follows what it marks,
  not the object it is merged into (PetWorld zoomed in, a hover redrawn: 1.4 → 0.6 ms on WebGPU, 7.1 → 0.5 ms on
  WebGL2, where the redraw no longer culls the scene again)
