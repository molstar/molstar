# Plan: Synchronous readback on WebGPU via canvas copies

## Status

- Planning. Follows [webgpu-viewer.md](webgpu-viewer.md) milestone 2 (picking with `readTextureAsync` and async `PickBuffers`). Related: [backend-split.md](backend-split.md) (readback and sync), [neutral-conventions.md](neutral-conventions.md).

## Goal

`identify` (clicks, the debug occlusion click, API callers) finds the picked element on WebGPU in the same call, as on WebGL, so that the click path needs no async special case and `Canvas3D.identify` behaves the same on both backends. Hover keeps the async read.

## Background

- WebGPU has no blocking readback: `mapAsync` resolves on a later task, the main thread can't block (`Atomics.wait` is not allowed there) and a blocking wait would deadlock, because the map callback needs the event loop.
- Copying the contents of a WebGPU canvas is synchronous: `drawImage` of a canvas (or `OffscreenCanvas`) with a `webgpu` context into a 2D canvas, then `getImageData`. The WebGPU spec defines this as "get a copy of the image contents of a context", which includes the current texture with work that hasn't been submitted yet; the browser flushes and waits internally (a stall like `readPixels` on WebGL).
- Today:
  - WebGPU pick buffers (`AsyncPickBuffers` in `canvas3d/passes/pick.ts`) have `synchronous = false`;
  - `PickHelper.identify` starts a read and finds nothing until it has finished;
  - `Canvas3dInteractionHelper` routes clicks through `canvasAsyncIdentify` when `synchronousIdentify` is false.
- Pick targets: object, instance and group ids are rgb with alpha 1, cleared to white with alpha 1 (`renderer.clear(false)`). The depth target packs depth into all four channels (`packDepthToRGBA`), so its alpha is data.
- Other synchronous reads (`readPixels` in histogram pyramid sum, grid3d, Gaussian density, texture mesh) are float or integer data of WebGL-only GPU compute. They are out of scope: rgba8 canvas copies can't carry them, and they become backend compute services ([backend-selection.md](backend-selection.md) step 3).

## Design

### Neutral API

- `GraphicsContext.readTextureSync?: (texture: Texture, x: number, y: number, width: number, height: number, out: Uint8Array) => boolean`: an rgba8 2d texture, a rectangle with a bottom-left origin, rows top to bottom (as `readTextureAsync`). Returns false if the backend can't read synchronously (the caller falls back to the async read). Optional, or present on both backends with WebGL implemented by `readPixels`.
- Capability flag `capabilities.syncReadback` (WebGL true; WebGPU true when the canvas copy works, see "Exactness check"), so that `createPickBuffers` can choose without trying.

### WebGPU implementation (`gl/webgpu/sync-readback.ts`)

- A lazily created `OffscreenCanvas` with a `webgpu` context, configured `format: 'rgba8unorm'`, `alphaMode: 'opaque'`, `usage: COPY_DST | RENDER_ATTACHMENT`, sized to the requested rectangle (resized on demand). A second `OffscreenCanvas` with a `2d` context, `willReadFrequently: true`, `colorSpace: 'srgb'`.
- Read: end the open pass, write the source rectangle (flipped to the top-left origin) into `getCurrentTexture()` (a `copyTextureToTexture` for rgb data, or the repack draw below), submit, then `drawImage` (with `globalCompositeOperation = 'copy'`, image smoothing off) and `getImageData` into `out`.
- With `alphaMode: 'opaque'` the alpha channel reads back as 255, so only rgb is exact: the read covers rgb data. Rgba data goes through a repack draw (below).
- Several textures per call (the four pick targets): one canvas sized `width × 4·height` with the rectangles stacked, one `drawImage`/`getImageData` round, instead of four stalls.

### Depth with its alpha channel

- The depth pick target uses all four channels. Options:
  1. **Repack to rgb** (proposed): a small WGSL quad program reads the packed rgba depth and writes `packDepthToRGB` (24 bits, like `packIntToRGB` precision) into the canvas texture instead of a copy. 24 bits are enough for the picked position (unprojected depth); the async path keeps full precision.
  2. Two reads (rgb, then alpha moved to rgb), doubling the cost.
  3. A pick depth target written as rgb by the pick variant on WebGPU, which would split the neutral pick output by backend; rejected, as the neutral pick layout should not change for a readback detail.
- The repack program also lets the id targets go through one code path (a quad program copying the rgb of each source into its slot of the stacked canvas texture), so no copy-specific format checks are needed.

### Pick buffers

- `AsyncPickBuffers` gets a synchronous `read()` when `readTextureSync` is available (one stacked read of the four targets); `synchronous` follows the capability. `asyncRead` stays for hover.
- With synchronous reads on WebGPU, `Canvas3dInteractionHelper` gets `synchronousIdentify = true` and the pending-click path is unused on both backends. Keep it for backends without synchronous reads, or remove it if both backends always have them.

## Exactness check

Canvas copies go through the browser's color and alpha handling, so values must be checked, not assumed:

- Unit-like browser test (visual harness or a Playwright check): fill a texture with all 2^24 rgb values (4096 × 4096, or chunks of it) and with random packed depths, read back synchronously, compare bitwise with the async read.
- Run on Chrome (D3D12, Vulkan, Metal where available), Firefox and Safari when WebGPU is enabled there; disable `syncReadback` on browsers that don't match (feature detection at context creation with a small known pattern, e.g. 256 values, and fall back to the async path).
- Known risks: color management between the canvases (both sRGB, no conversion expected), `getImageData` rounding of opaque pixels (none expected), canvas readback disabled by privacy settings (fingerprinting protection randomizes `getImageData` in some browsers, e.g. Firefox `privacy.resistFingerprinting`, Brave): the detection pattern catches it and falls back.

## Cost

- One GPU stall per synchronous read, as `readPixels` on WebGL; the pick viewport at pick scale 0.25 is small (a 1200 × 800 canvas gives 300 × 200 × 4 texels).
- Only `identify` (clicks) reads synchronously; hover keeps `asyncIdentify`. Measure the click latency against WebGL with the timing mode.

## Steps

1. **Prototype and exactness check**: `readTextureSync` for rgb data on WebGPU with the stacked canvas; the browser check over all rgb values; detection pattern and `capabilities.syncReadback`.
2. **Depth repack**: WGSL quad program for packed depth to 24-bit rgb (and the id copy into the stacked slots); compare picked positions with WebGL (unprojected position within a depth tolerance).
3. **Pick buffers**: synchronous `read()` in the WebGPU pick buffers, `synchronous` from the capability; `identify` finds elements in the same call; Viewer check (click selects without the pending-click path) and the hover/click map against WebGL.
4. **Cleanup**: decide whether the pending-click path in `Canvas3dInteractionHelper` stays (backends without sync readback) or goes; document the fallback.

## Validation

- WebGL visual profiles byte-identical (no WebGL change beyond the optional `readTextureSync` wrapper around `readPixels`).
- The exactness check passes or disables the capability per browser.
- Viewer: clicks select on the first click after a camera change, `canvas3d.identify` returns the element in the same call; the hover map matches WebGL (81 of 81 grid points as in milestone 2).

## Open questions

- Is `readTextureSync` worth exposing beyond picking (e.g. screenshots, which already use `readRenderTarget` asynchronously)? Proposed: picking only, until another synchronous need appears.
- Worker rendering (`OffscreenCanvas` transferred to a worker) keeps working: the copy canvases are `OffscreenCanvas` too.
