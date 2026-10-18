# Plan: Request the adapter's limits for WebGPU devices

## Status

- Done: devices request the adapter's maxima by default (`requestWebGPUDevice`, `createWebGPUBackend({ limits })`); the
  `webgpu` visual profile runs with them (16384 textures and 2 GB storage bindings on the development machine), the new
  `webgpu-core` profile with the core limits (8192, 128 MiB); both pass all 118 cases, the Viewer smoke test passes. A
  step of [backend-selection.md](backend-selection.md), Default backend.
- Before, `requestWebGPUDevice` requested no limits, so devices had the core defaults of the spec: 2d textures of 8192
  (WebGL usually allows 16384), storage bindings of 128 MiB, buffers of 256 MiB, 8 vertex buffers. Code sizes by
  `device.limits` (texture size, storage buffer size and slabs, dispatch limits), so it would take larger limits as they
  are; WebGL-sized criteria applied to WebGPU at half the texture size caused silent fallbacks before
  ([webgpu-compute.md](webgpu-compute.md), Limits).

## Design

- `requestWebGPUDevice(powerPreference, gpu, limits)`, `limits`: `'adapter'` (default) requests every `max*` limit at
  the adapter's value, `'core'` requests none (the spec defaults, the minimum spec of
  [minimum-capabilities.md](minimum-capabilities.md)). Alignment limits (`min*OffsetAlignment`) keep their defaults:
  smaller is "better" there and the uniform ring assumes 256 bytes.
- `createWebGPUBackend({ strict, limits })` passes it on, also to the device that replaces a lost one; headless takes
  the default.
- Visual tests: the `webgpu` profile runs with the adapter's limits (as users will), a `webgpu-core` profile with the
  core limits, so paths sized by the limits (slabs, texture sizes) stay covered at the minimum spec.
- A unit test of the requested limits against a stand-in adapter (limits as prototype getters, as in browsers).

## Decisions

1. All `max*` limits rather than a curated list: higher maxima cost nothing unless used, and every size computed from
   `device.limits` benefits; a new limit of the spec is picked up without a change.
2. A core profile in the visual tests rather than only the adapter one: the minimum spec stays the tested baseline.
