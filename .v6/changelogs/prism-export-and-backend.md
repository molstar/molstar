## Prism: graphics backend switch, scale line in exports

- The Canvas card's Graphics backend section and the `canvas.backend` operation (and command) switch between the
  viewer's backends (WebGPU, WebGL 2, WebGL 1) with the scene kept; shown where there is more than one. `canvas.get`
  reports the backend.
- Image and video exports take the scale line as they take the axes: a `scale` option of their own (Prism's
  `export.image`, the Image section, the video export; the classic screenshot helper), off by default, no longer the
  canvas's.
