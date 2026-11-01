# Plan: shader compile time

## Status

- **Current**: analysis in progress. Solid-interior cylinders fixed (their cap functions at one call site each: fresh
  shader text, WebGPU pipeline 955-1327 -> 522 ms, WebGL2 case 2.8 -> 1.0 s). An empty cache takes about 6 s for scenes
  of solid interiors with pixel clipping (WebGPU, 3 s on WebGL2), which no single pipeline explains.
- **Why**: the compilers inline every call site with what it calls (FXC on D3D12, ANGLE's on WebGL), so a large function
  reached from several exclusive branches is compiled that many times.

## How to measure

Every unique program of `shaders/coverage` compiled alone, one after the other, with a unique comment in its text, on the
RTX 4070 Ti. Timing pipelines while they compile together misleads: what is queued in one batch reads alike. Nothing else
may run alongside. A pipeline's state counts: unique by text and constants only left out pipelines that differ by state.

## Measured

- **Clip object count does not matter**: slowest pipeline at 1, 2, 4, 8 and 16 objects, ms: cylinders interior 872, 851,
  875, 972, 988 (WebGPU), 822, 797, 830, 821, 818 (WebGL2); spheres interior 539, 531, 543, 512, 535 (WebGPU).
- **By state**: solid-interior cylinders, the `colorBlended` pipeline 238 ms, the `color` pipeline 872 ms. Why: open.
- **WebGL2**: 950 programs, 84 s in all, median 64 ms, slowest 882 ms, 87 over 200 ms. Direct volumes are the most, 25 s
  of 99 programs (median 243 ms, 86 ms on WebGPU). Each program once by its text: only the `colorBlended` variants of
  what was timed by clip count.
- **WebGPU**: 1179 pipelines, 71 s in all, median 55 ms: too low, without the pipelines that differ by state only.

## Dropped

- `clipExit`'s intervals computed once per object, not within its loops: no gain (cylinders 955 -> 1251 ms, spheres 328
  -> 324 ms, WebGL2 unchanged).
- Noise (`fbm`), wiggle, FXAA and marching cubes: far from the 200 ms worth pursuing.

## Open

1. WebGPU again with pipeline state, and why `color` costs more than `colorBlended`.
2. Spheres' cap functions at one call site each, if their branches are exclusive.
3. Direct volumes on WebGL2: what costs, central against forward differences.
4. A real scene's total with an empty cache: many variants adding up.

## What a change must show

Only branches merged so that one call site remains, where one branch runs anyway; work moved ahead of an early exit needs
a runtime measurement. Compile time clearly better, repeated with fresh text; frame time unchanged on the matching
benchmark; images identical before and after.

## Also found

`shaders/coverage` fails on WebGPU: `prepareAllKernels` does not compile the sphere raster's four kernels. On demand only,
so not in CI.
