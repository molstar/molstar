## Graphics: ghosts of clipped parts

- `Canvas3DParams.clipGhost` (off by default; also `ImageParams.clipGhost`): what clip objects remove is drawn as a
  faint ghost, a tint and its edges, where it is nearer than the scene (behind a transparent surface where that is
  fainter than `marking.transparentThreshold`). Its color, opacity, edge opacity and resolution (half by default) are
  set there. Drawn as depth at that resolution with the clip test inverted: meshes by fragment, spheres, cylinders,
  points and lines whole by their centers. Over the illuminated image too.
- `clipGhost` geometry param (on by default): an object opted out has no ghost (a membrane whose ghost would cover the
  view, solvent).
- Prism: Ghost in the clipping section's Settings (`clipping.settings` `ghost`), what the cuts remove drawn as a faint
  ghost; its color, opacity, edges and resolution below it while on (`ghostColor`, `ghostOpacity`, `ghostEdgeOpacity`,
  `ghostResolution`).
