# Plan: Gizmo layers, overlay scenes replacing `HandleHelper` and `CameraHelper`

## Status

- **Current**: steps 1 to 9 done (the overlay layer; occluded alpha and saturation; screen-constant size of meshes;
  `HandleHelper` removed; the screen layer with the corner axes; the gizmo library and Prism on the overlay; DOF before
  the overlay, hit areas, billboards, unlit handles with a rim, the axes' hover, a scale line, line handles; the clipping gizmos' hit areas; the lines' pick width, already a pick pixel; the clipping gizmos' shape at
  the scene's scale, their move handles at a fixed size). Replaces Canvas3D's two helper scenes with two layers that any representation can be drawn in:
  - `overlay`, anchored in the world, takes over from `HandleHelper`;
  - `screen`, anchored to the viewport and turning with the camera, takes over from `CameraHelper` (the corner axes).

  Moves the gizmo library (`@molstar/graphics/util/gizmo`) and Prism's gizmos onto the overlay, and makes the corner
  axes the first user of the screen layer. Answers Prism's TODO "Gizmo handles behind what is drawn cannot be
  grabbed".

- **Step 1, as built** (`canvas3d/helper/overlay-helper.ts`, `canvas3d/passes/overlay.ts`; tests:
  `canvas3d/_test/overlay-camera.test.ts`, `canvas3d/passes/_test/overlay.test.ts`, the headless
  `webgpu-overlay.test.ts`):
  - `Layer` is `'scene' | 'overlay'` until step 5 adds `screen`. Canvas3D keeps each representation's scene
    (`reprScenes`); `update` of a representation whose layer changed (`setState`) moves it.
  - The overlay scene is committed after the scene, and only once the scene's queue is empty: a representation moved
    between the layers loses its scene renderable (disposed, with resources it may share) before the overlay creates
    one.
  - `OverlayPass` draws into a framebuffer of its own: the target's color texture with the pass's depth renderbuffer
    (cleared). On the drawing buffer (no attachments to swap) it clears and uses the buffer's own depth, which nothing
    reads after.
  - `OverlayCamera`: the camera's view and projection with near and far planes fitting the overlay's visible bounds
    (only the depth terms of the projection change, so offsets and a sample's jitter stay), fog off. One per source
    camera (stereo eyes, an image pass's, the ray helper's).
  - Draw and pick run without the renderer's occlusion test (Hi-Z tests the scene's depth): `Renderer.occlusionTest`.
  - Pick positions on overlay objects are unprojected with the overlay camera (its depth is of its own planes), in
    `PickHelper` and `RayHelper`; the ray helper also casts at the overlay's bounds.
  - A mark change on an overlay representation re-renders (the overlay is in the scene color, as the helpers were);
    its highlight is the color marker (the overlay scene's `dColorMarker`), not the marking pass.
  - `ImageProps.overlay` (off by default) draws it into images.

- **Step 2, as built** (`apply_occlusion` in GLSL and `apply_occlusion_wgsl`, between the material and the transparency
  check of mesh, spheres, cylinders, lines, points and text; tests: the headless `webgpu-overlay.test.ts`, the visual
  cases `overlay/**`):
  - `dOverlay` is a define of the object (`BaseSchema`), set by Canvas3D when it routes a representation, not a scene
    global: no global define to thread through every renderable, and precompile signatures carry it.
  - `occludedAlpha` and `occludedSaturation` are hidden geometry params (`BaseGeometry.Params`), the material uniforms
    `uOccludedAlpha`, `uOccludedSaturation`; default 1 (on top, as step 1).
  - The occlusion compares view z, not depth: the overlay's depth is of its own planes, the scene's opaque depth
    (`tDepth`) of the scene's (the globals `uOccluderNear`, `uOccluderFar`, set by `Renderer.renderOverlay` and
    `renderPick` from a `Renderer.Occluder`).
  - Occluded fragments of an opaque object fail the opaque pass (alpha below 1), so the transparent pass takes objects
    with `DrawFact.Occludable` (`dOverlay` and `uOccludedAlpha` < 1) as well; `renderOverlay` draws both passes.
  - Packed depth needs no path of its own: `tDepth` is the packed opaque depth target there, which `getDepth` reads.
  - While an overlay is shown, the draw pass doesn't draw straight to the drawing buffer (that path writes no opaque
    depth texture); the image goes through the color target and a copy, as with volumes.
  - A transparent object occludes where it is opaque enough to be picked: the front transparent depth with its alpha
    (`tDepthAlpha`) tested at `uPickingAlphaThreshold`, so what can be clicked hides what is behind it, in the image and
    the pick. Only the front surface counts: a fainter one in front of a stronger one lets the overlay through. The
    frame writes that depth while an overlay is over transparent objects (`DrawPass.isTransparentDepthRequired`,
    `IlluminationPass.renderInput`), one more pass only then. Test: the headless `webgpu-overlay.test.ts` (a surface
    above and below the threshold, drawn and picked), visual cases `overlay/occlusion/transparent/**`.
  - Picking: one per-fragment test instead of the design's split by `uOccludedAlpha` ≥ 0.5 (the pick depth buffer is
    of the scene's planes, so a depth test against it can't work). An occluded fragment is dropped where
    `uAlpha * uOccludedAlpha` is below the picking threshold, so the scene keeps the pick there. The scene depth is
    that of the last draw (`OverlayHelper.occluderDepth`), read at the fragment's place relative to the target (the
    pick targets are smaller). The ray helper's picks (XR) have no occlusion: their camera is not the drawn one.

- **Step 3, as built** (`assign_position` in GLSL, the mesh vertex stage in WGSL; tests: `overlay-camera.test.ts`, the
  headless `webgpu-overlay.test.ts`, the visual cases `overlay/screen-size/**`):
  - `screenSize` is a hidden mesh param (`Mesh.Params`), the object uniform `uScreenSize`, any layer. The instance's
    own scale multiplies (CSS pixels per unit of the mesh at instance scale 1).
  - The CSS height is the global `uCssHeight`, from the camera's viewport and the context's pixel ratio (as the levels
    of detail take it), not a pass target's: the pick at its scale, the stereo eyes and images of any size get the
    same size on screen. (The design's "divided by the pixel ratio" is this: CSS, not device, pixels.)
  - Culling: `DrawFact.ScreenSized` objects skip the frustum and lod tests (their world bounds depend on the camera).
  - The overlay camera's planes fit each screen-sized object per instance at the size the camera gives it
    (`getScreenSizeScale`, the shader's formula on the CPU); the overlay scene's bounding sphere is not used.

- **Step 5, as built** (`canvas3d/helper/screen-helper.ts`, `CameraAxes` in `canvas3d/helper/camera-helper.ts`; tests:
  `canvas3d/_test/screen-camera.test.ts`, the headless `webgpu-screen.test.ts`, the visual cases `screen/axes/**`):
  - `Layer` adds `'screen'`, placed by `State.screen` (`ScreenPlacement`: `anchor`, the center, a corner or an edge's
    middle; `offsetX`, `offsetY`; `followRotation`; `inset`, by default the bounding radius). An anchor is its sides
    (`screenAnchorSides`: x left, middle or right, y bottom, middle or top): inward from an anchored edge by the inset
    and the offset, the offset as is along a middle. The widgets' `location` takes the corners and the edges' middles
    (`ScreenEdgeAnchors`); the scale line is centered on a middle.
  - `ScreenHelper` keeps a `Scene` per owner (a representation, the corner axes), each drawn with a `ScreenCamera` of
    its own: orthographic in CSS pixels (the design's device pixels times the pixel ratio, as `CameraHelper` scaled its
    meshes), `viewOffset` kept, near and far planes fitting the object about its origin, so culling sees where it is.
  - `OverlayPass` draws it after the overlay, its depth cleared, without the occlusion test. With the screen layer
    only, on the drawing buffer it uses the buffer's own depth: the default viewer (axes on) still draws straight to
    the drawing buffer, without a copy.
  - Picked last, after clearing the pick depth, which wins the pixel as the design's "picked first" would.
  - The corner axes are not a shape representation but `CameraAxes` (was `CameraHelper`): its three render objects
    in a screen scene of its own. Still three: the transparent pass writes no depth, so the layer's own depth does not
    hide the insides of one object (design item 3). Params, `CameraAxesLoci` and `CameraHelperAxis` are unchanged;
    `Helper.camera` is `Helper.axes` and `Helper.screen`. `ImagePass` has its own screen layer and axes
    (`cameraHelper`).
  - Off in XR by `ScreenHelper.isEnabled` (while the context has a session), not by setting the axes' props off and
    back: nothing to restore after a session.
  - `reprCount` counts the scene layer's representations only (the overlay's and the screen's are not what the scene
    shows).
  - Not tested end to end: XR (no session headless); `isEnabled` is.

- **Step 6, as built** (`util/gizmo/controller.ts`; Prism's `sources/move/gizmos.ts`; tests:
  `util/gizmo/_test/controller.test.ts`, the headless `webgpu-gizmo.test.ts`):
  - `GizmoController` adds through its `VisualSet` with `layer: 'overlay'`. `GizmoLayout` gets `screenSize` and
    `occludedAlpha` per part (`DefaultOccludedAlpha`: body 0.25, handles 1).
  - Facing fade and drag dimming are one per-group transparency of the handles (`State.transparency`, a layer per
    value), set after the frame moves, the camera moves (`didDraw`) or a drag starts or ends. A slide fades as its
    axis turns toward the view (|cos| 0.9 to 0.98), a plane as it turns edge-on (0.25 to 0.08); hidden at 1, and not
    picked once the alpha is below the pick threshold (the shader's per-group transparency test). During a drag every
    other handle is at 0.75 (not picked either). The dragged handle is never faded: it stays under the cursor. Kept by
    its group, not its object: a layout made again during a drag has handles of its own. Steps of 0.05, so a turning
    camera updates the state only when a step changes.
  - Rings are not faded: edge-on a ring is still grabbed by its tangent (`rotateAbout`).
  - Prism: the clipping gizmos take the defaults (walls x-ray, handles on top), sized with the scene; the
    structure-move handles are a unit's at 90 CSS pixels (`screenSize`), so `StructurePlacement` loses its `radius`.
  - Unlit handles with a rim and the screen layer's hover came after (step 7, items 4 and 5).

## Today

- **`HandleHelper`** (removed by step 4; was `canvas3d/helper/handle-helper.ts`): one built-in translate handle (`Canvas3DParams.handle`) in
  its own blended `Scene`. It is drawn by `DrawPass._renderHelpers` after the scene, into the postprocessed target
  (before antialiasing, DOF and the marking's presentation), and picked after the scene (`PickPass.renderVariant`). It
  has its own loci (`isHandleLoci`) and highlight-only marking.
  - Nothing positions it: `update` has no caller.
  - Its depth depends on the path: with the opaque depth as a texture it is detached (no depth: no self-occlusion);
    with packed depth it is the scene's renderbuffer (depth tested). So it is neither reliably on top nor
    self-occluding.
- **`CameraHelper`** (replaced by step 5 with `CameraAxes` in the screen layer; `canvas3d/helper/camera-helper.ts`,
  `Canvas3DParams.camera.helper.axes`): the corner axes,
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

Gizmos are drawn and picked above the scene, and widgets anchored to the screen (the corner axes, a scale line) above both, all as representations. Per object there is control over:

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
     representations), and a later layer would be one more value of the same field. None is planned: highlight and
     selection are the marking passes' (their own depth and mask, composed onto the finished image), not a layer.
   - The layers are drawn in order (scene, overlay, screen) and picked in reverse. Each has its own depth (cleared).
2. **The screen layer**: drawn with the camera `CameraHelper` builds today (orthographic in device pixels, its view the
   main camera's rotation, `viewOffset` kept for tiled screenshots), generalised from one corner to any.
   - An object's place is its anchor (a corner or the center, plus an offset in CSS px) and its unit is a CSS pixel,
     both scaled by the pixel ratio when the layer's camera is built (`CameraHelper.update` does it for the corner
     today). So it is screen-sized by construction and needs neither item 5 (scene depth) nor item 6 (screen size).
   - A per-object `followRotation` (default on) turns it with the camera (axes) or keeps it upright (a
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
   - Later on the same layer: the helper's scale line.
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
   - glow: bloom does not reach the overlay (the handles' emissive and rim read without it).
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
7. What is left of the design, in the order below (each its own commit; see "Step 7"): the overlay out of DOF's
   blur; hit areas larger than what is drawn; billboards in the vertex stage; the handles' look; the screen layer's
   hover; the scale line; thin handles as screen-space lines.
8. Left over from step 7: hit areas for the clipping gizmos; the lines' smallest pick width.
   - **As built**: the clipping gizmos' geometry builds its handles twice (`gizmoGeometry`, `buildHandles`): as drawn,
     and as a hit area with the shafts, rings and edges four times as thick, the knobs as drawn and the cone's slide
     arrow starting past its move knob (a thick arrow from the apex would take the knob); the binding passes it as
     `GizmoLayout.hitArea`. Test: `clipping/gizmo/_test/geometry.test.ts` (each kind's hit area of the same groups,
     the cone's apex its knob's; fails with the arrow from the apex).
   - The lines' smallest pick width needed no change: a line is at least a pixel of the target it is drawn in
     (`max(1.0, linewidth)`), the pick's included (`uPixelRatio` the pick's ratio there), about 4 CSS pixels at the
     default pick scale, where a mesh shaft as thin falls between the pick's pixels. Test: the headless gizmo picks a
     line handle along its shaft without a hit area.
9. The clipping gizmos in two parts: the clip object's shape with the scene's scale, the handles that move it at a
   fixed size on screen, as the structure-move gizmo (see "Step 9").
   - **As built** (the slots since a list of parts: `gizmo-controls.md` §5): `GizmoLayout.fixedHandles` (a mesh, its hit area, a screen size; its groups the handles') drawn
     beside the layout's handles, and `GizmoLayout.handlePoints` (knobs as points: round, `pointSize` CSS pixels at
     any zoom, a group each); the controller keeps every part standing in for the handles in one list (`standIns`:
     the hit area, the lines, the points, the fixed handles and their hit area), their loci the handles'. The clip
     geometry (`buildHandles`) builds the shape's lines and points and the move handles at once, their groups shared:
     the shape where the clip is, every part of it a set size on screen (edges and rims as lines, `ringLines` in
     `primitives.ts`; a box's corner and face knobs and a cylinder's cap as points; a sphere's ring as lines facing the
     camera, `Lines.Params.billboard`); the move knob, slide arrow and turn ring about the origin in units of
     `MoveScreenSize` (90), the transform handles' proportions (a plane's slide is its knob, its arrow to the side it
     cuts away; the turn ring around the slide's axis). Hit areas: the move handles' (thick shafts and rings, the
     knob's center its own); the shape needs none (lines and points are a pick pixel at least, a world-sized area
     would not match them). The billboard is one function per backend (`billboardModel`, `chunks/billboard.glsl.ts`
     and `.wgsl.ts`) for the mesh and lines vertex stages; lines facing the camera in perspective are widened by
     d / √(d² − r²) (d the eye's distance from the origin, r a vertex's): a ring of a sphere's radius lies on the
     sphere's outline, not inside it. Knobs as points, not one instanced screen-sized mesh: a
     mesh's groups are not an instance's, so each knob would be a representation of its own (seven for a box), where
     points carry a group each. Tests: `clipping/gizmo/_test/geometry.test.ts` (every handle in a part, the move
     part's hit area, the move part the same at any reach, a box's corners where the box is, the cone's apex its
     knob's), the controller (the fixed part, the points: a set size, their picks the handles'), the headless gizmo (the fixed part the same size at two zooms, the
     outline where the world puts it, picked along it, larger when nearer), visual cases
     `overlay/clip-shape/zoom-1` and `zoom-0.5` (edges, corner points and a facing ring the same size at both zooms) and
     `overlay/clip-shape/sphere-outline` (the ring on a sphere's outline).

## Step 9

**Why.** A clip gizmo is sized by what it cuts (`reach`, the radius of its bounds): zoomed out on a large structure its
move knob, slide arrow and turn ring are a few pixels; zoomed in, they fill the view and hide the cut. Its outline
(a box's edges, a cylinder's rims, a sphere's ring, a cone's base rim) is the clip object itself: it has to stay where
the clip is, at the scene's scale. Its move, slide and turn handles only say where it is and how to move it: they read
best at a fixed size, as the structure-move gizmo's (90 CSS pixels, `screenSize`).

**Split.**

- **The shape part, at the scene's scale**: what is the clip object, and the handles that size it, drawn where they are
  in the world: the edges and rims as lines (`handleLines`: a set width on screen, step 7 item 8), the knobs on faces,
  corners and caps as one knob mesh drawn instanced at each knob's place with a screen size (`screenSize` scales about
  each instance's origin: the knobs a set size, their places the clip's). The body (the walls) as now.
- **The move part, at a fixed size on screen**: the move knob, the slide arrow and the turn ring, about the clip's
  origin (the cone's apex, a plane's center): the transform handles' parts (`transformHandles` with only `move` and
  `rotate`, or a clip's own layout of them) at 90 CSS pixels.

**How.** One gizmo of two parts with sizes of their own, not two gizmos: one item, one drag at a time, one frame. A
layout's handles mesh is screen-sized as a whole (`GizmoLayout.screenSize`), so the layout gets a second, screen-sized
handles part beside the world-sized one (e.g. `GizmoLayout.fixedHandles`: a mesh and its groups' ids in the same
`groups`, with its own `screenSize` and hit area), drawn and picked as the handles (`ShapeBuilder.lociShape`: its loci
the handles', so hover, highlight, fade and drags are as now). The knob instances need instanced shapes in a gizmo
part (transforms per knob in the shape: `Shape.create`'s `transforms`, the group per instance by `instanceGranularity`
off).

- Gestures: unchanged. A sizing handle (`scale`, `faces`, the cone's rim) is grabbed where it is in the world (its pivot
  the clip's); the move handles' grab point is the picked one on the screen-sized mesh (the overlay camera unprojects
  it, step 1).
- The sphere cut: its outline ring faces the camera (`billboard`) and is the shape part; its move knob the move part.
- Hit areas: the shape part's lines are a pick pixel at least (step 8); the knobs and the move part as the transform
  handles' (thick shafts and rings, the knob's center its own).
- Prism's clipping spec (`viewer/_spec/clipping.md` §4: "Sized with the scene's bounds") changes with it.

**Tests.** The geometry (unit): the shape part in world units scaling with the clip, the move part the same at any
`reach`. The controller: a gizmo with both parts (both drawn at their sizes, a pick on either the handles' loci, the
fade over both). Headless: a box cut's move knob the same size on screen at two zooms while its edges scale, an edge
picked anywhere along it. Visual cases: a box and a cylinder cut over a structure at two zooms.

## Step 7

Each item: what, why, how, tests. The order is by what a user meets first: a blurred or ungrabbable gizmo before a
new widget.

1. **The overlay out of DOF's blur.** DOF runs on the antialiased image (`DrawPass`: overlay, then antialiasing, then
   DOF), so with DOF on, gizmos and the corner axes are blurred. Swap the two: DOF on the post-processed image, then
   the overlay and the screen layer, then antialiasing (the usual order: an edge DOF softens needs no smoothing, a
   gizmo's does). Answers the first open question. `IlluminationPass` has the same order (overlay on its composite,
   then antialiasing, then DOF) and changes with it.
   - **As built**: `DrawPass._renderOutput` runs DOF, then (DOF on) the overlay into DOF's target, then antialiasing
     or the copy; with DOF off the overlay stays in `_renderHelpers` (the screen layer alone still straight to the
     drawing buffer). `getColorTarget` follows: antialiasing's target first. Tests: the headless `webgpu-overlay.test.ts`
     (an overlay ball's edge 7 pixels wide with DOF, the scene's 19; drawn before DOF it fails), visual cases
     `overlay/passes/dof*`, `postprocessing/dof/canvas`, `postprocessing/dof/antialiasing`.
   - Found on the way: DOF on a canvas never resized drew the background only (its `uTexSize` the null texture's
     size from creation; WebGL2's uniform blocks and WebGPU upload a cell only when it is updated), and a frame held
     for its props prepared before the overlay committed, so on WebGPU it skipped DOF into its own target (the
     prepare now runs after the overlay's commit, `Canvas3D.commitScene`).
   - Tests: a visual case with DOF on and a gizmo out of focus (sharp); the DOF cases unchanged but for the order's
     few pixels (checked on the images, not by the suite's equality).
2. **Hit areas larger than what is drawn.** The pick is drawn at a quarter of the canvas's pixels and searched one
   pick pixel around (`pickPadding`), so a shaft a few CSS pixels wide falls between its pixels (the headless gizmo
   test has to aim at the arrow's head).
   - `RenderableState.pickOnly` (the mirror of `colorOnly`: drawn in the pick variants only) and
     `Representation.State.pickOnly`. A pick-only object is out of the color passes, not out of the bounds (the
     overlay's planes have to hold it to pick it).
   - `GizmoLayout.hitArea?: Mesh`: its groups the handles' (same ids), added as a pick-only representation beside the
     handles, under the same frame, screen size and fade (its transparency layers the handles').
     `transformHandles` gives one: shafts and rings at four times the tube, the arrows' heads and knobs as drawn.
   - Lines: a smallest width in the pick variants (`uPickMinWidth`, in pick pixels), for thin lines anywhere (a
     measurement's).
   - Tests: the pick of a shaft at its middle (headless), a proxy that draws nothing (visual: the image without it),
     the bounds without it.
   - **As built**: `RenderableState.pickOnly` (optional: literal states elsewhere stay valid), `DrawState.PickOnly`;
     the pick passes are flagged (`PassDraws.pick`, the pick and its depth) and `Renderer`'s `isPassDrawn` keeps
     pick-only objects to them (all three paths: table rows, recordings, single objects). `State.pickOnly`, applied by
     shape and generic representations (as `colorOnly`). A pick-only object stays in the bounds, the overlay's planes
     included: out of them it could not be picked. `ShapeBuilder.lociShape`: a hit area's picks are the handles'
     loci (`Shape` identity decides hover, highlight and clicks). `transformHandles` returns `hitArea` (shafts and
     rings four times as thick, the knob half as large again and the shafts starting past it, so the center stays
     the knob's). Tests: the controller (a pick-only representation, the handles' loci, hover through it, the fade),
     the headless gizmo (a shaft picked at its middle; fails without the hit area), the headless overlay (picked with
     its position, not drawn in images or the canvas).
   - Left: the lines' smallest pick width (`uPickMinWidth`), with thin handles as lines (item 8); hit areas for the
     clipping gizmos (`clipping/gizmo/geometry.ts` builds its own meshes, sized with the scene).
3. **Billboards in the vertex stage.** `GizmoLayout.facesCamera` turns a gizmo after a frame (`didDraw`, one frame
   late, a state update each frame the camera moves). A mesh define `dBillboard` (a hidden `billboard` mesh param)
   turns the object about its instance's origin to face the camera in the vertex stage, normals with it, as text
   does: every camera of a frame right (stereo, the pick, XR), no lag, no state updates.
   - Culling: about the origin, so its bounding sphere grows to the origin's distance plus its radius (as screen size
     does with `DrawFact.ScreenSized`).
   - The controller sets it for `facesCamera` layouts and drops `#faceCamera`'s frame turning (keeps the fade's
     update on camera moves).
   - Tests: the sphere cut's ring round from two views in one frame (stereo), the pick where it is drawn.
   - **As built**: `Mesh.Params.billboard` (hidden), the uniform `uBillboard` (first built as a define, `dBillboard`:
     a define doubles the programs of the objects that set it, compiled when a gizmo first shows, for a few vector
     operations a vertex; a uniform branch costs nothing to compile, as `uScreenSize`). The vertex stage turns the model matrix
     about the instance's origin (the smallest turn taking the mesh's +Y to the camera: `uCameraPosition`, or
     `-uCameraDir` orthographic), before the view: positions, normals, `vModelPosition` and the pick follow (GLSL in
     `assign_position`, WGSL `billboardModel`). `DrawFact.ScreenSized` is `DrawFact.CameraBounds` (`hasCameraBounds`:
     screen-sized or billboarded, not culled by its bounds); the overlay's planes take a billboard per instance as a
     sphere about its origin (its reach from there). The controller sets `billboard` for `facesCamera` layouts and
     draws them at their frame; after a camera move it only fades (no state updates for the turn). Its drag maths
     still turn the frame on the CPU (`facingCamera`, from world up): the same place for the sphere cut, whose parts
     are about its origin; a facing layout with a pivot off the origin would differ by the twist. Tests: the
     controller (the prop, the frame kept as the camera moves), the headless gizmo (a ring around +Y round from three
     sides, picked on its rim; fails without the billboard), visual cases `overlay/billboard/side` and `/stereo`.
4. **The handles' look.** Lit by the scene's lights today, a handle on the dark side of a structure reads as a hole.
   Unlit with a rim: a mesh param `rim` (hidden; 0 off): the color brightened by `rim * (1 - |n·v|)²`, in the material
   (GLSL and WGSL), so it reads as 3D under `ignoreLight`. The controller draws handles with `ignoreLight` and a rim,
   the body flat. Answers "handles lit or unlit": unlit.
   - Tests: a visual case of the transform handles on every profile; a rim of 0 identical to today's images.
   - **As built**: `Mesh.Params.rim`, the material uniform `uRim` (no define: no new programs), mixed toward white
     after the light (`mesh.frag`, the WGSL mesh fragment), plain meshes only (texture meshes share the shader without
     the uniform). `GizmoControllerOptions.rim` (0.5): handles `ignoreLight` with the rim and their glow; the body
     flat as before. Tests: visual cases `overlay/handles/rim-0` and `rim-0.5` (the rim on the parts' outlines: 3913
     pixels differ, up to 102 levels).
   - Defines of the layers so far: `dOverlay` only (the occlusion's depth texture and its read, which scene objects
     should not carry). Occluded alpha and saturation, screen size, billboard and rim are uniforms: a gizmo compiles
     the programs of one more define value, the overlay's.
5. **The screen layer's hover.** The corner axes dimmed (their `alpha`) while the pointer is elsewhere, at full alpha
   while it is over them, as navigation widgets do: Canvas3D's hover sets the axes' alpha through `CameraAxes`
   (`setHovered`), no rebuild (`uAlpha`).
   - Tests: the alpha after a hover on and off (headless).
   - **As built**: Canvas3D subscribes to its own hover (`interactionHelper.events.hover`): `CameraAxes.setHovered`
     sets the axes and planes' `uAlpha` (1, or their `alpha`; the labels keep `labelOpacity`), kept through a rebuild
     while hovered. A change updates the screen scenes (`updateChanged`): at alpha 1 an object leaves the transparent
     pass, and the scene plans its opaque pass only if its summary (`hasOpaque`) says so (without it the hovered axes
     vanished). Test: the visual case `screen/axes/hovered` (a canvas capture, hovered through `canvas3d.interaction`;
     checked against `screen/axes/canvas` by eye: the suite compares backends, which would agree on a vanished
     widget). Headless has no readback of the canvas, so no headless test.
6. **The scale line** (`camera.helper.scale`, on by default; the helper's TODO): a line on the screen layer,
   upright, bottom-middle by default (the axes bottom-left), a label above it. Its length the largest of 1, 2, 5 × 10ⁿ
   Å that fits a set width (CSS pixels) at the camera's target depth (`perPixel`): exact in orthographic, at the
   focus in perspective (the label says neither: what a scale line means). After a camera move: a new length moves
   the line's instance scale; a new label rebuilds the text (rare).
   - In images as the axes are (`cameraHelper.scale`).
   - Tests: the length and label at two zooms (unit), a visual case orthographic.
   - **As built** (`canvas3d/helper/scale-line.ts`, `ScaleLine`; `Helper.scaleLine`, the image pass's own):
     `CameraHelperParams.scale` (location, offsets, `maxWidth`, color, line width, label size 20 and bold; on by default; states without it off). Lines
     (width in pixels) of unit length stretched along x by their transform, ticks up at both ends; the label a text
     object (Å below 1000 Å, nm below 10000 Å, µm). Both made in `setProps`, never in a frame (a frame draws only
     what was compiled before it: a label made in the frame was "program not finalized"); `update` rewrites them in
     place (the transforms, the label's text into its buffers with its counts and bounds, as shape representations
     do), before Canvas3D's passes each render (a few operations when unchanged) and in `ImagePass` at its camera.
     Screenshots show it as the canvas does (`viewport-screenshot`, the axes keep their own option). `CameraAxes`
     takes only `axes` (`Pick<CameraHelperProps, 'axes'>`). Tests: `canvas3d/_test/scale-line.test.ts` (the value,
     the label, the world size of a pixel), visual cases `screen/scale-line/image` and `/canvas`.
7. **The view cube**: built, then removed: the corner axes do its work (a click aligns the camera).
8. **Thin handles as screen-space lines.** `transformHandles` with `lines`: shafts and rings as lines of a set width in
   CSS pixels (`lineSizeAttenuation` off), the heads and knobs as meshes, a hit area (item 2) as proxy. Needs
   `GizmoLayout` to take a lines part (`handleLines`, its groups the handles'). Last: the meshes at a screen size
   (step 3) already keep handles thin at any zoom.
   - **As built**: lines take a screen size too (`Lines.Params.screenSize`, the object uniform `uScreenSize`: both
     ends scaled about the instance's origin as the mesh's vertices, GLSL and WGSL), so a screen-sized gizmo's lines
     meet its meshes. `transformHandles` with `lines` returns `handleLines` (the shafts to where the heads start, the
     rings as 64 segments; `arrowHead`, `arrowShaftEnd` in `primitives.ts`); its hit area as without. The controller
     draws `handleLines` at `GizmoControllerOptions.lineWidth` (2.5 CSS pixels), its loci the handles'
     (`ShapeBuilder.lociShape` now decides marking as well as picking: a highlight of a handle marks its line, and the
     fade's layers are on the handles' shape for every part). Opt-in: structure-move keeps its meshes. Tests: the
     handles with lines (unit: which groups are lines, which mesh), the controller (a line's pick and mark are the
     handle's, the fade), visual cases `overlay/handles/lines` and `lines-screen-size/zoom-1`, `zoom-0.5`.

## Decisions

- Occlusion by transparent objects: by the front transparent depth, at the picking threshold (step 2, as built). Why:
  what can be clicked hides what is behind it, and a faint surface does not.
- Gizmos of two controllers overlapping: any order. Rare; a priority waits for a case in practice.
- Several widgets on the screen layer: distributed over the anchors, not stacked in one (the defaults: axes bottom-left,
  scale line bottom-middle). Nothing moves a widget out of another's way. A host app's own widgets
  are out of scope here.
- XR: the screen layer's `isEnabled` is tested, a session end to end is not: enough for now.

## Fewer schema additions

**Built.** What the gizmos added to the schemas, reviewed by who sets each outside tests; the three below removed, the
as-built notes above that name them (step 2's `occludedSaturation`, step 7's rim and `facesCamera`) are history.

**Kept**, used by both gizmos: the mesh's `uScreenSize` (every handle at a size on screen) and `uBillboard` (the
transform handles' facing circles, their hit areas, the trackball's fill); `dOverlay` and `uOccludedAlpha` (the
controller's body and handles); the globals `uCssHeight`, `uOccluderNear`, `uOccluderFar`. The mesh's `uDepthOffset` is
the clip walls' (`clip-volumes.ts`), not a gizmo's. The lines' pick width (`uPickMinWidth`) is gone already (a line is
at least a pick pixel wide).

**Removed:**

1. **`uRim`** (`Mesh.Params.rim`, `GizmoControllerOptions.rim`): a material uniform and a fragment branch in every mesh
   program for one look of the gizmo handles, which the controller no longer asks for by default. Out of the schema,
   both mesh shaders (`mesh.frag`, `mesh.wgsl`), the controller's options and its tests, the visual cases.
2. **`uOccludedSaturation`** (`BaseGeometry.Params.occludedSaturation`): nothing sets it, so every program in the
   overlay carries a material uniform and a mix that is always 1. Out of the schema, both occlusion chunks
   (`apply-occlusion.glsl`, `material.wgsl`), the overlay test (`webgpu-overlay.test.ts`) and visual case
   (`overlay.ts`).
3. **`GizmoLayout.facesCamera`**: left over from the clipping gizmo's round sphere cut, whose only producer went when
   sphere cuts became ellipsoids; one controller test and one headless test set it. Out of the layout, `#frameOf` (the
   frame as is), the tests. With it the lines' perspective outline under `uBillboard` (`lines.vert`, `lines.wgsl`: a
   ring widened to a sphere's outline, "a sphere cut's ring stays round"), wrong for a flat facing circle.

Tests: the quick lane (the controller, the transform handles, the geometry params), the headless overlay and gizmo
tests (the ring facing the camera a part's `billboard` now); the visual cases `overlay/**` (WebGL2 and WebGPU, 41 each:
`overlay/handles` one case, unlit and glowing; the `desaturated` occlusion gone).

## The axis guides near the view axis

**Built.** A drag's axis guide (a ring's, a slide's, a plane's two; `guideLines`) went `reach` (twice the eye's
distance) either way: in perspective the half toward the viewer passed the eye and streaked across the view as the axis
turned toward it. Now that half ends at half the origin's depth (`GuideFrame.depth`): on screen tan(angle to the view)
from the origin, a short stub end-on. An axis within the slides' hidden cone (`Fade.slide`, |cos| ≥ 0.98) is left out:
decided at the press (the axis is the press frame's), nothing fades. Tests: `guideLines` (an axis 20° off the view, one
5° off, one in the screen).

## Open

None.

## Next: the facing parts as points

**Built** (the `ring` style: `pointRingWidth`, `uPointRingWidth`; the depth offset:
`pointDepthOffset`, `uPointDepthOffset`, clamped into the depth range; tests: `points.test`, the headless
`webgpu-points.test`, visual `geometry/points/ring-0.3`, `-0.6`, `depth-offset`; the handles: `GizmoPart.point`,
`transform.test`, the headless `webgpu-gizmo.test`'s pick order unchanged; step 4: visual `overlay/handles/parts`,
`overlay/facing/side`, `/stereo` for the billboard's cases; the line handles' and the sphere outline's cases gone). The transform handles' parts facing
the camera as points: a point is turned to the camera by nature and, without size attenuation, a set size in pixels,
so they need neither a billboard nor a screen size. A `ring` point style for the circles, the `circle` style for the
discs, and a depth offset for points to keep the pick order.

**Why.** The center circle, the roll circle, the trackball's fill and the hit areas of the center, the roll circle and
the trackball face the camera: meshes now, turned in the vertex stage (the mesh's `uBillboard`) at a set size on screen
(`uScreenSize`). As lines (`transformHandles({ lines: true })`, only tests use it) the circles need the lines'
`uBillboard` and `uScreenSize`, which nothing else uses. Points are both by construction, one vertex a part.

**What.**

1. **A `ring` point style** (`Points.StyleTypes`, `dPointStyle`; `points.frag`, `points.wgsl`, the fragment's discard
   as `circle` does): the quad's circle less an inner one, its width a uniform (`pointRingWidth`, a fraction of the
   radius, 1 the full disk: no variant per width). Picked where drawn (the discard in the pick variant too), so a wider
   pick-only ring is a hit band. Opaque, with the occlusion (`apply_occlusion`) as every style. The define's values grow
   (`dPointStyle`, the signatures): no other program changes.
2. **A depth offset for points** (`pointDepthOffset`, CSS pixels toward the camera, negative away; a uniform), in depth
   only as the mesh's `depthOffset` (the point stays where it is on screen): a world length of that many pixels at the
   point's depth, so it scales with the gizmo's size on screen. As built: kept in front of the camera and clamped into
   the depth range (the camera's planes fit the bounds, which know nothing of the offset: past a plane, a point stays
   at it, still in front of or behind what is drawn). What keeps the pick order the meshes have by their
   placement: the center's hit disk in front of every part (1.25 times the handles' size toward the camera), the
   trackball's behind every part (1.25 times away). `pointFront` (twice the point's own radius) stays for knobs.
3. **The transform handles** (`transformHandles`): the center circle and the roll circle `ring` points at the frame's
   origin, the trackball's fill a `circle` point (its alpha, on hover, as now); their hit areas pick-only points: the
   center's a `circle` in front, the roll circle's a wider `ring`, the trackball's a `circle` behind; the scale knobs
   `circle` points, as the clipping gizmo's knobs. Sizes in CSS pixels: a part's radius times the layout's
   `screenSize` (`GizmoPart.pointSize`, a part's own where the controller's one size does not fit). The arrow heads,
   the plane squares and the axis rings stay meshes (a direction, a plane, an axis: they turn with the frame). As
   built: `GizmoPart.point` (`radius`, `style`, `ringWidth`, `depth`, `front`; the radius and depth in the frame's
   units times the screen size); a circle's ring point is its band (the radius plus the width, the band's width its
   fraction), its hit band four times as wide; `front` only for the knobs (it moves a point by twice its radius, not
   clamped: a large ring past the near plane).
4. **Removed**: the mesh's `uBillboard` (no part faces the camera as a mesh anymore) and the lines' `uScreenSize` and
   `uBillboard` with the `lines: true` opt-in (its facing circles were their only use; the shafts and axis rings stay
   meshes). The mesh's `uScreenSize` stays for every mesh part.

**Tests.** Quick: the params and the define (`points.test`), the handles' parts (`transform.test`: the facing parts
points of `ring` and `circle`, the knobs `circle`). Headless WebGPU: a `ring` point picks on its ring, not in its
hole nor outside; a point with a depth offset picks over (or under) a mesh at its depth and is drawn where it was; the
facing parts round and the same size on screen from any side and zoom, the center picked over an edge-on ring, the
trackball under every part (`webgpu-gizmo.test`, the existing pick-order cases). Visual: the points cases for `ring`
(two widths), the overlay handles case again, both backends.

**Decided.** A ring's width a fraction of its radius, not pixels: simple, one param; the controller computes a set
stroke from it where it wants one (`4 / size` for the lines' 4 px).

