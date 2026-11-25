# Plan: follow-ups from testing the overlay layer

## Status

- **Current**: done. 1: the scale line an option of the image and video exports, off by default (the breaking changes
  say how an app keeps it out of its canvas). 2: `Representation.State.pickAtAnyAlpha`
  (applied by the canvas, as `layer`) sets the material uniform `uPickingAlphaMin` to 1, a floor of the alpha the pick
  check tests (GLSL and WGSL, the image shader's own check too); the corner axes set it on their render objects, picked
  at opacity 0.2 (`webgpu-screen.test.ts`, failing without it); gizmo bodies and fills keep the threshold. With it,
  `occludedAlpha` moved from the geometry params to the state too: not in every representation's params. 3:
  `GizmoPart.lit` (mesh parts only, unlit by default) draws the part with the canvas's lights (`ignoreLight` off); no
  visual case of its own, the overlay's occlusion cases draw lit meshes in that layer. From a test of the overlay and
  screen layers (gizmo-layer.md) on Linux with a GTX 980, Viewer and Prism, WebGL2 and WebGPU: no blockers, three small
  things.

## 1. The scale line in exported images

- **Now**: the image pass has a screen layer of its own (`ImagePass.init`, `cameraHelper`), but the plugin's screenshot
  helper passes the canvas's scale line (`viewport-screenshot.ts`, three places:
  `scale: c.props.camera.helper.scale`). The axes have an export option (`axes`, off by default); the scale line has
  none. `CameraHelperParams.scale` defaults to `on`, so exported images include it, and an app setting
  `camera.helper: { axes }` only gets the scale line too (the default fills it in). Headless screenshots turn both off
  (`defaultImagePassParams`).
- **Change**: a `scale` option in `ViewportScreenshotHelperParams` next to `axes`, of `CameraHelperParams.scale`, off
  by default as the axes are; the three places pass it. The screenshot UI lists it under the axes.
- **Breaking changes**: an entry for `camera.helper.scale` (on by default; `camera.helper` with only `axes` gets the
  scale line; `scale: { name: 'off', params: {} }` turns it off).
- **Open**: whether the canvas's default should be `off` instead (the Viewer turning it on), so that apps without it
  see no change. Recommended: keep `on`; it is useful and the entry says how to turn it off.

## 2. Translucent corner axes not picked

- **Now**: the screen layer is picked as the scene is (`PickPass`, `renderPick` per screen scene), against the global
  `pickingAlphaThreshold` (0.5; the shaders' `check_picking_alpha` discards below it). The corner axes' default
  opacity of 0.51 is just above it; below 0.5 they cannot be hovered or clicked.
- **Not a global change**: the overlay and screen layers rely on the threshold. A gizmo's `body` is translucent so that
  it is not picked, and a faint fill below 0.5 is left out on purpose (`GizmoPart.alpha`). Picking those layers at
  any alpha would break both.
- **Change**: a per-object opt-in, picked at any alpha: a render value (e.g. `uPickingAlphaThreshold` overridable per
  object, or a `dPickAnyAlpha` define) that the check reads in GLSL and WGSL. The corner axes set it; widgets that
  want it can too. The fog check stays.
- **Test**: the corner axes at opacity 0.3 identify their loci (a quick test with the pick helper), a gizmo body
  still does not.

## 3. Lit gizmo meshes

- **Now**: the gizmo controller draws every mesh part with `ignoreLight: true` (`controller.ts`, the mesh case), so
  handles are flat; lines and points have no lighting anyway.
- **Change**: an opt-in on `GizmoPart` (`lit?: boolean`, default unlit) or on the layout for all of its meshes; a lit
  part is drawn with `ignoreLight: false` and the renderer's light, its hover glow (`emissive`) unchanged. Defaults stay
  as they are.
- **Open**: on the part or on the layout. Recommended: on the part (a handle lit, its rim not), the layout's as a
  default for its parts if that is wanted later.

## Order

1, then 2, then 3: each is small and independent. 1 touches the plugin and the breaking changes, 2 the shaders of both
backends (the shader parity test covers the new value), 3 the gizmo controller only. Visual cases: the corner axes at
low opacity picked (2) and a lit gizmo (3), on both backends.
