## Graphics: merged units and GPU-culled views

- [Breaking] The geometry param `renderMerged` is removed: units representations always merge their groups (one group
  too), at small counts measured no worse than unmerged (2 to 30 units, WebGL2 and WebGPU) and cheaper to update at
  30; the particle targets merge two or more objects of a type, or one culled on the GPU.
- Overpaint, transparency, emissive and substance layers on merged surfaces are color-smoothed where they reach (they
  were left unsmoothed on merged ones): the reached units drawn as a merged object of their own, smoothed per vertex on
  the CPU. 60 copies of 7qpd (420 units), a layer set on one unit 18.8 ms against 12.9 ms for per-unit grids, on all
  2.75 s against 2.04 s (WebGPU), the same memory (35 MB). Smoothing colors per vertex is faster: 3.67 to 2.75 s of
  that.
