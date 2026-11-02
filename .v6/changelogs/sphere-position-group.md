## Graphics: memory

- [Breaking] Spheres and points keep center and group in one array, `positionGroup` (read by the shaders as
  `tPositionGroup`), not also in `centerBuffer` and `groupBuffer`: half their CPU bytes (16 instead of 32 per sphere;
  3J3Q as spacefill, 2.44 M spheres: 39 MB less)
