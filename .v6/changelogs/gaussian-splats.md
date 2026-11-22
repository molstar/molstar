## Graphics: Gaussian splats

- Gaussian splats (3D Gaussian splatting) as a geometry kind on WebGL2 and WebGPU: anisotropic Gaussians projected per
  splat (EWA), blended back to front (sorted when the view turns: on the GPU on WebGPU, 1 M splats 0.4 ms; on the CPU on
  WebGL2, 8-10 ms) in a pass of their own after the other transparency, in every transparency mode and with
  illumination; params `sizeFactor`, `kernelSize`, `alphaCutoff`. Not on WebGL1.
- PLY files of Gaussian splats load as splats: the standard layout of 3DGS training and the compressed layout of
  PlayCanvas and SuperSplat, turned upright by default (`upright`); view-dependent color (spherical harmonics, degree 1
  to 3) with `shDegree` (0 by default).
- `.splat` (antimatter15), `.spz` (Niantic, versions 1 to 3) and `.sog` (PlayCanvas, versions 2 and 3; its WebP images
  decoded by the browser, not in Node) files load as gaussian splats.
- Binary PLY tables of one property type read as one typed view, not value by value.
- `countingSortIndices` in `@molstar/core` (`data/util/sort`): a stable counting sort of small integer keys.
