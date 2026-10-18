# Plan: XR on WebGPU

## Status

- **Current**: not started. The neutral interface is in place: `GraphicsXR` (`gl/context.ts`) with `isPresenting()` and
  `getView(view, out)`, used by `StereoCamera` and the XR manager; WebGL implements it with `XRWebGLLayer`, WebGPU's
  is a stub (`set` throws).

## Where things stand

- The camera needs per view of a pose: the viewport in the drawing buffer (top-left origin) and the projection with
  GL's clip depth (`[-1, 1]`, the backends convert tagged projections, see `DepthRangeTag`). The views share the
  drawing buffer, side by side, each draw sets the eye's viewport.
- WebGL: `XRWebGLLayer` as `baseLayer`, its framebuffer is the drawing buffer; its viewports have GL's bottom-left
  origin and are flipped in `getView`; projections have GL's clip depth.

## Decisions

1. **WebGPU through the WebXR/WebGPU binding**
   ([explainer](https://github.com/immersive-web/WebXR-WebGPU-Binding/blob/main/explainer.md)): sessions requested
   with the `'webgpu'` feature (WebGPU sessions cannot use `XRWebGLLayer`), the adapter with `xrCompatible: true`,
   `new XRGPUBinding(session, device)`, a projection layer (`createProjectionLayer` with the preferred color format
   and the drawing buffer's depth-stencil format), set with `updateRenderState({ layers: [layer] })`.
2. **Drawing buffer**: per frame the layer's sub-images (`getViewSubImage(layer, view)`) are the drawing buffer; with
   one texture array per eye (texture layers) the views do not share a target, so the draw pass renders the eyes into
   their sub-images, or into one target of both views copied to the sub-images.
3. **getView**: the sub-image viewport as given (WebGPU's top-left origin); `view.projectionMatrix` has a `[0, 1]`
   clip depth in WebGPU sessions, converted to GL's (`z' = 2z - w`) in double precision.

## Steps

1. Session setup and teardown in the WebGPU context (`set`, `end`, `isPresenting`), the adapter request with
   `xrCompatible` where XR is wanted (an option of `WebGPUContextOptions`).
2. The drawing buffer per frame from the sub-images (decision 2), `getView` (decision 3).
3. Checks: the WebXR emulator (or a headset) on WebGPU and WebGL, the same scene; the XR manager's pointer rays.

## Open questions

- Whether the session's requested features (`'webgpu'`) are chosen by the Viewer from the context's backend, or by the
  XR manager through the context (`GraphicsXR.getSessionInit()`).
