## Graphics: merged units and GPU-culled views

- [Breaking] The geometry param `renderMerged` is removed: units representations always merge their groups (one group
  too), at small counts measured no worse than unmerged (2 to 30 units, WebGL2 and WebGPU) and cheaper to update at
  30; the particle targets merge two or more objects of a type, or one culled on the GPU.
