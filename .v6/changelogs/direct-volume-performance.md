## Graphics: faster direct volume rendering

- [Breaking] The volume direct-volume representation no longer has `jumpLength`: volumes skip transparent bricks
  instead; the Gaussian volume keeps it
- Direct volume rendering marches only through the volume's box (from its entry to its exit or the nearest opaque
  surface), skips the work of fully transparent samples and skips bricks of 8³ cells the transfer function makes
  transparent: sparse maps and tomograms 7-9x faster (256³ map at 2560x1440: 8.5 → 1.0 ms on WebGL2, 8.5 → 1.1 ms
  on WebGPU, float data 15.9 → 1.6 ms; 512x512x96 slab: 8.9 → 1.2 ms)
- Direct volumes offset each pixel's samples by a fraction of a step (varied per multi-sample sample), and take 2
  samples per cell by default (was 3): the same look, a quarter faster
- `volumeResolutionScale` Canvas3D and image param: direct volumes at that resolution relative to the drawing buffer,
  composited with depth-aware upsampling; 0.5 by default (frames with a volume filling the view 2.5x faster); in
  Prism, the Canvas card's Advanced section (Volume resolution)
- Lit direct volumes render up to a quarter faster on WebGL (the normal matrix computed per instance, not per pixel)
- Direct volumes are created faster: their 3d texture is built 2-3x faster (256³ map: 0.9 → 0.4 s)
- Direct volumes look the same at any `stepsPerCell` also for dense transfer functions: the opacity per step is
  corrected exactly, not linearly
- Highlighted or selected direct volumes (with the marking pass off) no longer fill their whole box with a faint haze;
  only the visible density is highlighted
- Highlighted or selected direct volumes show the marking pass's edges and fill (the default; were not marked at all),
  drawn where the volume looks solid
- `surfaceOpacity` param of direct volumes (0.3 by default): the opacity along a ray at which a volume counts as
  solid, where marking draws its edges and picking finds it
- Direct volumes can be picked (were not): hovering or clicking where a volume looks solid gives its cell there
- Fix: in the `dpoit` transparency mode, draws before the first depth-peeling pass of a canvas no longer read stale
  texture units (a WebGL error now and then, e.g. after a direct volume was drawn)
- Direct volumes are drawn under transparent objects (were over them, also over those in front): transparent surfaces
  in front of a volume cover it, as do those behind it (volumes are not split at transparent surfaces)
- Fix: direct volumes in the `wboit` transparency mode (the default) no longer show holes (white speckles) where their
  rays turn fully opaque, e.g. a transfer function reaching 1
- Fix: direct volumes of grids with non-cubic cells (e.g. 1x1x4 Å) are no longer cut off at their far side
- Fix: direct volumes colored per cell no longer show black pixels on grid-aligned samples
