# Plan: Viewports outside the target, cropped by Canvas3D

## Status

- **Current**: done: the audit, the camera's crop, Canvas3D cropping its frame, the WebGPU guard; checked on the iOS
  device (Safari 26, `?debug-log=1`).
- **The problem**: the neutral `setViewport` takes a viewport partly outside the pass target and crops it, as WebGL
  does and as the WebGPU spec allows (`setViewport`: `x, y ≥ -2·maxTextureDimension2D`,
  `x + width, y + height ≤ 2·maxTextureDimension2D − 1`; the width and height at most `maxTextureDimension2D`). Chrome
  (Dawn) follows the spec. Safari 26 on iOS still validates the old rule (the viewport inside the attachment): the
  `setViewport` is a validation error that invalidates the pass encoder, so the whole pass draws nothing, not just the
  part outside.
- **Who sets such viewports**: the SSAO and volume passes did (their scaled viewports rounded past the target), and are
  fixed. The `gpu-paths/viewport-cropping` visual case does it on purpose. Canvas3D's `viewport` prop (`static-frame`,
  `relative-frame`) is user input: a frame extending past the canvas reaches every pass through `camera.viewport`.
  Nothing else (the audit, step 1): the other `setViewport` callers set a camera's viewport, a pass target's full size
  or tiles inside one (Gaussian density, color smoothing, bloom's blur), and the image pass, the ray helper and the XR
  eyes have cameras of their own whose viewport is their target.

## Goal

A Canvas3D frame extending past the canvas draws the same on every backend and implementation: the part inside the
canvas, as if the whole frame were drawn and cropped. No pass loses its draws on an implementation that rejects the
viewport.

## Approach

Canvas3D crops the frame itself, so no viewport outside the target reaches the backend: the passes draw into the frame's
part inside the canvas, and the camera projects only that part of the frame's frustum (a sub-frustum, as an off-axis
projection). Every pass, shader and readback then sees one consistent viewport (the cropped one): `uViewport`,
`uDrawingBufferSize`, SSAO, outlines, Hi-Z, picking and `Camera.project`/`unproject` agree without knowing about the
frame. What measures the view rather than draws it keeps the whole frame (`Camera.getFrame`, below): rotation, pan and
zoom speeds, sizes in pixels, as today.

- **The readers of `camera.viewport`** (the audit, step 1), by what they want:
  - the cropped viewport, as drawn: every pass (draw, multi-sample, jitter, marking, postprocessing, FXAA, SMAA, CAS,
    DOF, bloom, SSAO, shadow, tracing, illumination, volume, background, Hi-Z), the pick helper and pick pass,
    `Camera.project`/`unproject`/`getRay` (canvas pixels through the cropped projection), the renderer's `uCssHeight`
    (with `projection[5]`, a pixel's size: the two change inversely under the crop);
  - the frame, the extent the field of view spans: the controls (`object.ts` `panDirection` with `mouseOnScreen`, the
    trackball's pan; the trackball's own `viewport` is `controls.viewport`, the frame already), the pointer helper's
    `getSizeForPixels`, the focus manager's distance check (`focus.ts`), volume streaming's box aspect
    (`volume-streaming/behavior.ts`), the camera helper's corner;
  - neither: `drag-frame.ts` (gizmo drags) takes the camera viewport for the canvas's CSS size, and is off for any frame
    smaller than the canvas, crop or not: the canvas's size there, a fix of its own; stereo builds the eyes'
    projections from `camera.viewport` and the field of view (no view offset either): a frame past the canvas in stereo
    is not cropped, but passed as given (cropped by WebGL and Dawn, the guard's squashing where rejected).

- **The camera's crop**: a `Camera.Crop` (`enabled`, the frame's `fullWidth`/`fullHeight`, the cropped part's
  `offsetX`/`offsetY` in it), applied in `updatePers` and `updateOrtho` before the view offset: the aspect and the
  ortho extent from the full frame, the frustum narrowed to the cropped part. It is separate from `viewOffset`, which
  stays the jitter: `setJitter` overwrites `viewOffset` with the viewport as the full size (`jitter.ts`), and the
  renderer passes its offset to the shaders as the jitter (`uViewOffset`, `renderer.ts`: the dithering in
  `material.wgsl`), so a crop folded into it would be lost or read as jitter. Jitter composes with the crop as today:
  sub-pixel offsets of the cropped viewport. `Camera.getFrame(out, camera)` gives the frame (the viewport where the
  crop is off) for the readers that want it.
- **Canvas3D**: `updateViewport` intersects the frame with the drawing buffer, sets the passes' viewport
  (`camera.viewport`, the renderer, the pick helper, Hi-Z) to the intersection and the camera's crop to its place in
  the frame; `controls.viewport` is the frame. A frame with nothing inside the canvas, and any frame in stereo, is passed
  as given.
- **The backends** keep `setViewport` as given (the neutral API follows WebGPU, `neutral-conventions.md`). Only the WebGPU
  backend guards against the implementations that reject it: a probe at context creation (a pass with a viewport at
  `-1` in an error scope, beside the `fragDepthInBundles` probe) decides; where rejected, a viewport outside the target
  is clipped to it (squashed into the target: wrong, but its pass is drawn) and logged once in debug mode. Until the
  probe answers, viewports outside the target are clipped.

## The gizmo layers

The `gizmo-layer` branch ([gizmo-layer.md](gizmo-layer.md)), to be merged first or alongside, adds cameras derived from
the main one; how each meets the crop:

- **`OverlayCamera`** (`overlay-helper.ts`) copies the main camera's projection and changes only its depth terms
  (`setProjectionPlanes`): the crop's narrowed frustum is carried as is, nothing to do but copy the crop with the view
  offset.
- **`ScreenCamera`** (`screen-helper.ts`) builds its own orthographic projection in CSS pixels of `source.viewport` and
  anchors the corner axes, the view cube and the scale line to that viewport's corners and edges. It already takes
  `viewOffset` as a tile of a larger image ("keeps the place of the full image in a tile of it", `screen-camera.test.ts`),
  but nothing produces such a view offset: the only producer is the jitter, whose full size is the viewport. The crop
  is that tile: `ScreenCamera` anchors to the crop's full frame and projects the cropped part of it, the jitter's view
  offset applied on top as sub-pixels; the tile reading of `viewOffset` goes (`setScreenProjection`, the test rewritten
  against the crop).
- **Screen-sized objects and the scale line** take a pixel's size from `projection[5]` and `uCssHeight` (the viewport's
  height in CSS pixels). Under the crop both change by the frame's height over the cropped height, inversely, so their
  product, and the size, stays: nothing to do; a quick test of a screen-sized mesh's size cropped and not.
- **Hit areas, occlusion by the scene, depth of field**: from `camera.viewport`, the target's pixels and unprojection
  with the cropped camera, consistent as for the passes.

## Steps

1. **Audit** (done): the callers and readers above. The camera helper (the axes gizmo, `camera-helper.ts`) is placed
   in a corner of `camera.viewport`: it is placed from the frame instead and its camera takes the crop as it takes the
   jitter's view offset, so the gizmo stays in the frame's corner (cropped with the rest where that corner is outside
   the canvas), as today. After the gizmo layers merge, the same holds for the screen layer
   ([The gizmo layers](#the-gizmo-layers)).
2. **The camera's crop**: `Camera.Crop` and `Camera.getFrame`, in `updatePers` and `updateOrtho`; quick tests that the
   projection of a cropped camera equals the full frame's projection followed by the crop (a point in the frame lands
   on the same pixel), perspective and ortho, with and without jitter.
3. **Canvas3D**: `updateViewport` crops the frame (not in stereo); the frame's readers on `Camera.getFrame`; a Dawn
   test (`webgpu-viewport-crop.test.ts`) drawing a scene through a frame off two sides of the canvas (left and top,
   right and bottom; perspective and ortho) against the frame drawn whole on a canvas of its size: the pixels inside
   equal; picking at an atom in the cropped part finds it.
4. **The WebGPU guard**: the probe (`testViewportOutsideTarget`, `viewport-outside-target.ts`; the context option
   `viewportOutsideTarget` forces it, `WebGPUContext.viewportOutsideTarget` tells it), the clipping
   (`toTargetViewport` in `getPassState` and the bundle clear's scissor check; also on devices that take a viewport
   outside, where it is beyond the spec's bounds: larger than `maxTextureDimension2D`, as a frame passed as given may be)
   and the debug log, once; a Dawn test (`webgpu-viewport-outside-target.test.ts`: the probe true, the viewport as
   given cropped, forced clipped squashed and drawn, an oversized one clipped, no errors); the Canvas3D test covers a
   frame larger than the canvas, off every side; the `viewport cropping` visual case checks the cropping where the device takes
   the viewport and, where it does not, that the pass is drawn. Checked on the iOS device: the Viewer with
   `?debug-log=1` draws and logs the probe's line.

## Not doing

- **Mapping clip positions in the WebGPU shaders** (the backend clips the viewport and every vertex shader maps its
  position from the full viewport to the clipped one, a scale and offset per draw state): exact for any caller, but
  about twenty `position` writes, a uniform the neutral global schema (shared with WebGL) or a fourth bind group (the
  last of the default `maxBindGroups`) would carry, draw-list offsets that shift when the viewport changes inside a pass,
  and passes that read a viewport the backend silently remaps. The fallback if a caller other than Canvas3D turns out
  to need viewports outside the target on Safari.
- **Requiring viewports inside the target in the neutral API**: the spec allows them and WebGL crops them; the guard
  covers the implementations that do not.
- **Rendering through an intermediate target the size of the frame**: a copy per pass and the memory of a second
  target.
