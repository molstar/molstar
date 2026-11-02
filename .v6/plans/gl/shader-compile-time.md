# Plan: shader compile time

## Status

- **Current**: analysis in progress. Solid-interior cylinders fixed (their cap functions at one call site each: fresh
  shader text, WebGPU pipeline 955-1327 -> 522 ms, WebGL2 case 2.8 -> 1.0 s). No pipeline is slow enough to matter now.
- **Why**: the compilers inline every call site with what it calls (FXC on D3D12, ANGLE's on WebGL), so a large function
  reached from several exclusive branches is compiled that many times.

## How to measure

Every unique program of `shaders/coverage` compiled alone, one after the other, with a unique comment in its text, on the
RTX 4070 Ti. Timing pipelines while they compile together misleads: what is queued in one batch reads alike. Nothing else
may run alongside. A WebGPU case's reference page is WebGL2's: both write to what logs the times.

## Measured

- **Clip object count matters little**: slowest pipeline at 1, 2, 4, 8 and 16 objects, ms: WebGPU, cylinders interior 238,
  248, 256, 294, 347, spheres interior 172, 182, 186, 185, 193, without interiors flat; WebGL2, cylinders interior 872,
  851, 875, 972, 988 and 822, 797, 830, 821, 818 in two runs, spheres interior 539, 531, 543, 512, 535.
- **WebGPU**: 2019 pipelines (unique with their state), 103 s in all, median 45 ms, slowest 260 ms, one over 200 ms:
  cylinders with solid interiors and pixel clipping, 260 color, 200 marking, 186 pick; spheres alike 168, 148, 136.
- **WebGL2**: 950 programs, 84 s in all, median 64 ms, slowest 882 ms, 87 over 200 ms. Direct volumes are the most, 25 s
  of 99 programs (median 243 ms, 83 ms on WebGPU).

## Dropped

- `clipExit`'s intervals computed once per object, not within its loops: no gain (cylinders 955 -> 1251 ms, spheres 328
  -> 324 ms, WebGL2 unchanged).
- Noise (`fbm`), wiggle, FXAA and marching cubes: far from the 200 ms worth pursuing.
- Spheres' cap functions: one call site each already in the color, pick and marking variants, the back cap behind
  `dRenderVariant_depth`.
- Cylinders' back cap compiled into the depth variant only: three runs each, WebGPU 252 -> 224 ms, WebGL2 862 -> 865 ms.
  Too little for another branch of the preprocessor.
- Direct volumes by forward differences instead of central: three runs each, WebGL2 231 -> 232 ms, WebGPU 96 -> 93 ms.
  The doubled reads cost no compile time.

## Not reproduced

About 6 s (WebGPU, 3 s on WebGL2) with an empty cache for a scene of solid interiors with pixel clipping, seen once; no
pipeline takes close to that, and it did not come back. In the same period `renderer/clipping/impostor-interior` stalled
in the GPU process for minutes and `renderer/marking/direct-volume` drew wrongly on WebGPU, repeatably. All three are gone
on the same commits and the same Chrome, without a reboot: the machine's state, by inference. If it comes back, before a
reboot: `chrome://gpu`, the GPU process's stack, a run at a known commit.

## Open

- Why WebGL2 takes 2.5-3.5 times WebGPU on the slowest programs (solid-interior cylinders 862 against 252 ms, direct
  volumes 231 against 92 ms), though about the same at the median: ANGLE's translation, not inlining the shaders avoid.

## What a change must show

Only branches merged so that one call site remains, where one branch runs anyway; work moved ahead of an early exit needs
a runtime measurement. Compile time clearly better, repeated with fresh text; frame time unchanged on the matching
benchmark; images identical before and after.
