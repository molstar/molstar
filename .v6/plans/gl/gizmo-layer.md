# Plan: Gizmo layers, overlay scenes replacing `HandleHelper` and `CameraHelper`

## Status

- **Current**: proposed, nothing done. Replaces Canvas3D's two helper scenes with two layers that any representation
  can be drawn in:
  - `overlay`, anchored in the world, takes over from `HandleHelper`;
  - `screen`, anchored to the viewport and turning with the camera, takes over from `CameraHelper` (the corner axes).

  Moves the gizmo library (`@molstar/graphics/util/gizmo`) and Prism's gizmos onto the overlay, and makes the corner
  axes the first user of the screen layer. Answers Prism's TODO "Gizmo handles behind what is drawn cannot be
  grabbed".

## Today

- **`HandleHelper`** (`canvas3d/helper/handle-helper.ts`): one built-in translate handle (`Canvas3DParams.handle`) in
  its own blended `Scene`. It is drawn by `DrawPass._renderHelpers` after the scene, into the postprocessed target
  (before antialiasing, DOF and the marking's presentation), and picked after the scene (`PickPass.renderVariant`). It
  has its own loci (`isHandleLoci`) and highlight-only marking.
  - Nothing positions it: `update` has no caller.
  - Its depth depends on the path: with the opaque depth as a texture it is detached (no depth: no self-occlusion);
    with packed depth it is the scene's renderbuffer (depth tested). So it is neither reliably on top nor
    self-occluding.
- **`CameraHelper`** (`canvas3d/helper/camera-helper.ts`, `Canvas3DParams.camera.helper.axes`): the corner axes,
  three render objects (axes, planes, labels) in their own blended `Scene`. Its own camera is orthographic in device
  pixels, its view the main camera's rotation, placed in a corner (`location`, offsets; `viewOffset` for tiled
  screenshots). It is drawn last in `_renderHelpers` and picked last.
  - It has its own loci: `CameraAxesLoci` with `CameraHelperAxis` group ids. The plugin's `CameraAxisHelper` behavior
    turns a click on them into a camera alignment.
  - Its own copy in `ImagePass` (`cameraHelper` props) draws it into screenshots.
  - Canvas3D turns it off in XR (no screen to anchor to) and back on after.
  - Its scene's depth is whatever the target holds at that point. Its parts sort by object only, so it splits its
    translucent parts into separate objects (axes, planes) to hide its own insides.
  - The helper's TODO: a scale line or grid.
- **Gizmos** (`GizmoController` over a `VisualSet`; Prism's clipping and structure-move gizmos) are shape
  representations in the main scene. So they:
  - are hidden by, and lose picks to, what is in front (the pick takes the nearest object);
  - are counted in the scene's bounds (camera reset and focus, near and far planes);
  - are clipped by clip objects, fogged, shaded by SSAO, outlines and illumination, and blurred by DOF;
  - are in screenshots and geometry exports;
  - are sized in world units (the item's reach): tiny when zoomed out, screen-filling when zoomed in;
  - turn to face the camera a frame late (`facingCamera` on `didDraw`, then one more frame).

## Goal

Gizmos are drawn and picked above the scene, and widgets anchored to the screen (the corner axes, a scale bar, a
view cube) above both, all as representations. Per object there is control over:

- how scene depth applies;
- an optional size that stays constant on screen;
- their look and blending.

They are left out of everything that is about the scene (bounds, post-processing, exports, screenshots).

## Design

1. **A layer as representation state**: `Representation.State.layer: 'scene' | 'overlay' | 'screen'` (`scene` by
   default). `Canvas3D.add`, `update` and `remove` put a representation's render objects into its layer's `Scene`
   (`Helper.overlay` and `Helper.screen`, in place of `Helper.handle` and `Helper.camera`). `reprRenderObjects` holds
   every layer, so `getLoci`, `mark`, hover, click owners and backend switching work unchanged. `VisualSet` passes the
   layer (`add`/`create`/`addShape` state).
   - Why state, not a Canvas3D flag: the state survives a backend switch (the canvas re-creates the set's
     representations), and later layers (a highlight layer, one with its own post-processing) are more values of the
     same field.
   - The layers are drawn in order (scene, overlay, screen) and picked in reverse. Each has its own depth (cleared).
2. **The screen layer**: drawn with the camera `CameraHelper` builds today (orthographic in device pixels, its view the
   main camera's rotation, `viewOffset` kept for tiled screenshots), generalised from one corner to any.
   - An object's place is its anchor (a corner or the center, plus an offset in CSS px) and its unit is a CSS pixel,
     both scaled by the pixel ratio when the layer's camera is built (`CameraHelper.update` does it for the corner
     today). So it is screen-sized by construction and needs neither item 5 (scene depth) nor item 6 (screen size).
   - A per-object `followRotation` (default on) turns it with the camera (axes, a view cube) or keeps it upright (a
     scale bar, a compass label).
   - Picked first, depth cleared: a corner widget always takes its own pixels.
   - Off in XR, as the axes are now; the overlay layer stays (gizmos are useful in XR).
3. **The corner axes as a representation**: Canvas3D adds a shape representation of the axes, planes and labels to the
   screen layer while `camera.helper.axes` is on. Its params stay the same. Its parts can be one object again, since
   the layer's own depth hides its insides.
   - Its loci stay `CameraAxesLoci`, with the same `CameraHelperAxis` group ids: the representation's `getLoci` maps
     its shape's groups to them. So the plugin's `CameraAxisHelper` and other code that reads them keep working.
   - `ImagePass` draws the screen layer when its `cameraHelper` props say so (as now), and the overlay only when asked
     (item 4).
   - Later on the same layer: the helper's scale line, and a view cube (a click on a face, edge or corner aligns the
     camera; a drag orbits) built with the gizmo library's handles.
4. **Drawing**: the overlay comes after the scene's image is finished (post-processing, illumination's composite, the
   multi-sample accumulation's per-sample draw) and before antialiasing. It is drawn as a blended scene
   (`renderBlended`) with its own depth buffer (a pass target, cleared), so a gizmo occludes itself. The screen layer
   follows, with its own camera and the same depth buffer cleared again; neither is post-processed.
   - The scene's opaque depth is bound as a texture (`tDepth`, already a global texture) for the occlusion test of
     item 5.
   - Its own near and far planes come from its bounds (as the pointer and camera helpers have their own camera), so a
     gizmo outside the scene is not cut by the scene's planes. Its bounds stay out of `boundingSphere` and
     `boundingSphereVisible`.
   - The same place in `IlluminationPass` and `MultiSamplePass`.
   - `ImagePass` draws it only when asked (`overlay`, off by default: screenshots without gizmos). Geometry export
     skips it.
   - Out of merging (`renderMerged`), GPU culling and Hi-Z (both layers): few objects, always in view. Draw lists
     and the precompile recorder take them as any scene (a gizmo shown the first time compiles nothing).
5. **Scene depth, per object**: a material uniform `uOccludedAlpha` (0 to 1). A fragment behind the scene's opaque
   depth has its alpha multiplied by it:
   - 1: on top;
   - 0: depth tested as in the scene;
   - in between: x-ray (the part behind the scene drawn faint).

   One code path for every mode on both backends: a fragment chunk in GLSL and WGSL behind a define (`dOverlay`), so
   scene objects pay nothing. An optional `uOccludedSaturation` greys the occluded part, as 3D tools draw hidden
   lines.
   - Packed depth (WebGL1 without depth textures): no per-fragment test; 1 draws on top, below 1 is tested against the
     scene's renderbuffer.
   - Picking: the overlay is picked after the scene. Objects with `uOccludedAlpha` ≥ 0.5 (the alpha a click needs,
     as the gizmo bodies use) are drawn after clearing the pick depth, so they always win. The rest are tested against
     the scene's pick depth.
6. **Screen-constant size, per object**: an object uniform `uScreenSize` (CSS pixels per unit of the mesh; 0: world
   units). The vertex stage scales about the instance's origin by the world size of a pixel there:
   - perspective: `-viewZ * 2 / (P[1][1] * height)`;
   - orthographic: `2 / (P[1][1] * height)`;
   - divided by the pixel ratio.

   Mesh first (all the gizmo library draws); lines and text later.
   - Why the shader, not the controller each frame: every camera of a frame gets it right (the pick pass at its pick
     scale, both stereo eyes, XR, multi-sample jitter), with no frame of lag.
   - Bounds for the near and far planes: the object's sphere at the largest scale the frame's camera gives it.
   - Facing the camera can move to the same stage later (a billboard flag), lag-free as well.
7. **Look and blending**, the overlay's and the screen layer's own:
   - no fog and no clip objects (`dClipObjects` off for overlay objects): a clipping gizmo is not cut by its own cut;
   - unlit or flat with a rim (`ignoreLight` plus a rim term), so it reads as 3D on any background;
   - hover and the active handle by the marker's highlight (as now). During a drag the other handles are dimmed (the
     controller's transparency over their groups);
   - facing fade (the controller): a slide or plane handle fades as its axis turns toward the viewer, and is hidden and
     not picked below a threshold (per group transparency, updated when the camera turns);
   - the screen layer's alpha as a whole (the axes' `alpha` today) and dimmed while the pointer is elsewhere, brighter
     when hovered (as navigation widgets do);
   - glow: bloom does not reach the overlay. A halo, if wanted, is an overlay effect of its own (later), not the
     scene's bloom.
8. **The gizmo library on it**:
   - `GizmoLayout` gets `screenSize?` and per-part depth (`occludedAlpha` for body and handles; defaults: body 0.25,
     handles 1).
   - `GizmoController` adds through `VisualSet` with `layer: 'overlay'`.
   - Prism's bindings choose per gizmo: a clip volume's walls x-ray, handles on top, the structure-move handles at a
     fixed screen size.
9. **The helpers removed**:
   - `HandleHelper`: its handle is never positioned, and the gizmo library's transform handles
     (`util/gizmo/transform.ts`) do what it was for. This removes `Canvas3DParams.handle`, `HandleHelper`,
     `isHandleLoci` and `HandleGroup`, recorded in `breaking-v6-changes.md`.
   - `CameraHelper`: its camera moves to the screen layer and its meshes to the axes representation. Its params
     (`CameraHelperParams`), `CameraAxesLoci`, `isCameraAxesLoci` and `CameraHelperAxis` stay, so nothing outside
     Canvas3D changes. `Helper.camera` goes, `ImagePass`'s own copy of it with it.

## Steps

1. Layer state, Canvas3D routing, `Helper.overlay`, drawn and picked by every pass where the handle is now. Own depth
   target, own near and far planes, out of the bounds. Tests: routing, loci and marking through the overlay, pick
   priority, bounds unchanged by a gizmo.
2. Occluded alpha (GLSL, WGSL, define, uniform) and the pick split. Visual cases: a gizmo inside a structure on top,
   x-ray and depth tested, on every profile.
3. Screen size (mesh vertex stage). Tests: the same size at two zooms, the pick at another pick scale, stereo.
4. `HandleHelper` removed (ledger, v6 changelog entry).
5. The screen layer, with the corner axes on it and `CameraHelper` removed. Tests: the axes' loci and the plugin's
   alignment by a click, every corner and offset, tiled screenshots (`viewOffset`), off and back on around XR. Visual
   cases: the axes over a structure, in a screenshot with and without them.
6. Gizmo library and Prism on the overlay, facing fade and drag dimming. Prism's TODO item closed.
7. Later: hit areas larger than what is drawn (pick-only proxies), billboards in the vertex stage, more layers
   (highlight, overlay effects such as a halo), thin handles as screen-space lines; on the screen layer, the scale
   line and a view cube.

## Open questions

- DOF runs after antialiasing today. The overlay belongs after DOF (a gizmo is never blurred) but before antialiasing.
  Swap them (antialiasing after DOF, the usual order), or draw the overlay after both, unsmoothed, when DOF is on?
- Occlusion against the opaque depth only (a gizmo inside a translucent surface reads as in front), or the transparent
  depth too (`depthTextureTransparent`)?
- Two controllers' gizmos overlapping: drawn in the order added, or a priority on the state?
- Handles lit by the scene's lights, or always unlit?
- The screen layer's place among the layers of other apps (a host's own corner widgets): one screen layer, objects in
  the order added, or anchors that stack (two widgets in one corner side by side)?
