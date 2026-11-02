## Core

- `Clip.classifySphere` and `Clip.classifySphereObject`: how a sphere lies against clip objects (no point clipped,
  some, every point), from distance bounds valid for every shape; `Clip.getSignedDistance`, the shaders' signed
  distance on the CPU. `Clip.canIntersectSphere` on them: exact for planes, cubes, cylinders and round spheres,
  bounded (no longer always true) for ellipsoids and infinite cones

## Graphics

- What clip objects (without a clipping layer) remove whole is culled: objects, and the cells and instances of
  instanced objects on the CPU and the GPU (PetWorld zoomed and cut through its center by a plane: 5.12 → 1.76 ms
  culled on the GPU, 8.24 → 4.47 on the CPU)
- Sphere and cylinder impostors, points and lines clipped by pixel: those clip objects remove whole are dropped in the vertex stage, and
  those they do not cut skip the per-fragment clip tests (300k atoms as spacefill at 2560×1440 by three clip objects
  cutting nothing: 25.4 → 15.6 ms on WebGPU, 21.4 → 12.5 on WebGL2)
