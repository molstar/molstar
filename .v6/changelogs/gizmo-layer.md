## Graphics: an overlay layer for gizmos

- [Breaking] The handle helper is removed (`HandleHelper`, `isHandleLoci`, `HandleGroup`, `Canvas3DParams.handle`): its
  handle was never positioned; gizmos are overlay representations made with the gizmo library.
- [Breaking] A cylinder clip object is elliptic: its radii are the scale's x and z halves (the height y), as the
  sphere's three axes; z was ignored, so a cylinder given as (diameter, height, 1) is now flat. GLSL, WGSL and
  `Clip.getSignedDistance`.
- [Breaking] `CameraHelper` is `CameraAxes`, the corner axes in the screen layer, without a camera or scene of its own;
  `Helper.camera` is `Helper.axes` and `Helper.screen`. Its params, loci and axis ids are unchanged.

- `Representation.State.layer` (`'scene'` by default, `'overlay'`): an overlay representation is drawn and picked above
  the scene, with a depth and near and far planes of its own (not cut by the scene's), after post-processing and
  illumination's composite and before antialiasing, without fog. It is out of the scene's bounds (camera reset and
  focus, the scene's near and far planes) and wins the pick wherever it is drawn. A `VisualSet` passes it with the state of
  `add`, `create` and `addShape`; a `setState` that changes it moves the representation.
- Per object, how the scene's depth applies to the overlay: `occludedAlpha` (geometry param, hidden; 1 on top, 0 depth
  tested, between x-ray) multiplies the alpha of the parts behind the scene. Parts too transparent to click leave the
  pick to the scene. GLSL and WGSL, every primitive geometry.
- `screenSize` (mesh param, hidden): CSS pixels per unit of the mesh, the same size on screen at any zoom, in every
  pass (the pick at its scale, stereo, images); such objects are not culled by their bounds.
- `ImageProps.overlay` draws the overlay into images; off by default (screenshots without gizmos).
- A screen layer (`Representation.State.layer` `'screen'`): drawn and picked above the scene and the overlay, anchored
  to the viewport by `State.screen` (the center or a corner, an offset, turning with the camera or upright), its unit a
  CSS pixel; off in XR. The corner axes are drawn in it, at their opacity, and at full opacity while the pointer is
  over them.
- The axes and the scale line take an edge's middle as their place (`top-middle`, `bottom-middle`,
  `middle-left`, `middle-right`) beside the corners; the scale line is centered there.
- A scale line (`camera.helper.scale`, on by default): a line at the bottom middle of the screen layer as long as the largest
  of 1, 2, 5 × 10ⁿ Å that fits its width at the camera target, labeled (Å, nm, µm; outlined, readable over what has its
  color); in screenshots as on the canvas.
  `ImageProps.cameraHelper` takes `scale` beside `axes`.
- Transform handles (`transformHandles`, Prism's structure-move) as Blender's move and rotate gizmo combined: the
  center a small circle facing the camera (move in the screen's plane, or a slide: `centerAxis`), picked over every
  other part; arrow heads past the rings (their shafts drawn while moving; `arrows`); small plane squares; each ring's
  front half (`rings`), turning by the angle the cursor goes around the center on the screen, by steps (`RingRule`
  `around`, Blender's; `tangent` and `across` where the ring is edge-on); a grey circle around them rolling the item about the view axis (`roll`, a new `Handle` kind: the
  angle the cursor goes around the item's origin, by steps); the area inside the rings turning it as the trackball
  turns the scene (`free`), a faint fill on hover.
- Gizmos (`GizmoController`) are drawn and picked in the overlay layer: a handle inside a structure can be grabbed.
  A gizmo is a list of parts (`GizmoLayout.parts`: meshes, lines or points, each a handle, a hit area or the body, at
  its own screen size and facing the camera or not) and of representations of the app's own (`GizmoLayout.extras`:
  placed, faded and removed with it, a handle where their loci are the handles'). `GizmoLayout.occludedAlpha` per part (by default the body x-ray, the handles on top) and `GizmoLayout.screenSize`
  (the same size on screen at any zoom). A slide or plane handle turned edge-on to the view fades out and is not picked.
  During a drag a turn hides every other handle and a move the handles that do not move (the others dim); guides show
  what it does (a slide's or a ring's axis across the view, a plane's two, a ring's radius at the grab and now, a dashed
  line to the cursor and a double arrow there, a cross at the center of a trackball's turn), and a moved handle a grey copy where it started. A part may be shown only in a state
  (`GizmoPart.when`: its group hovered, a kind dragged), at an alpha of its own (`alpha`), or only its half toward
  the viewer (`front`: a ring's, whole while it turns).
- Hit areas larger than what is drawn: `Representation.State.pickOnly` (shape representations; `RenderableState.pickOnly`)
  is drawn in the pick passes only, and `ShapeBuilder.lociShape` reports its picks as another shape's loci.
  `hit` parts (`GizmoPart`), and `transformHandles` gives one: a thin shaft is grabbed anywhere along it.
- Gizmo parts as lines (a set width on screen, `GizmoControllerOptions.lineWidth`; their lengths the world's). A
  representation built with `ShapeBuilder.lociShape` is marked (highlight, overpaint, transparency) by that shape's
  loci, as it is picked.
- Gizmos: the transform handles' parts facing the camera are points (`GizmoPart.point`: a radius, a style, a ring's
  width, a depth): the center and roll circles `ring` points, the trackball's fill and the hit discs `circle` points
  (the center's in front of every part, the trackball's behind), the scale knobs `circle` points.
- `pointDepthOffset` (points param, hidden): CSS pixels toward the camera (negative: away), in depth only (as large
  and where it was on screen), clamped into the depth range; GLSL and WGSL.
- `ring` point style: a ring, its width a fraction of the radius (`pointRingWidth`, 1 the full disk), drawn and
  picked on the ring only; GLSL and WGSL.
- `pointFront` (points param, hidden): each point twice its radius toward the camera along the view ray (as large
  and where it was on screen), in front of what meets it, drawn and picked; GLSL and WGSL. Gizmo knobs take it: a corner
  knob wins over the edges ending in it. Gizmo lines are thicker (4 px) and the transform handles' circles and the guides
  grey, not white (lost on a white background or structure).
- Depth of field runs before antialiasing (the usual order), and the overlay and the screen layer are drawn after it:
  gizmos and the corner axes are never blurred.
- Prism: the clipping walls x-ray through the structure and their handles are on top, with hit areas larger than
  what is drawn (a thin ring or edge is grabbed anywhere along it); the structure-move handles are the same size on
  screen at any zoom.
- Prism: a sphere cut is drawn and edited as the ellipsoid it clips (its three sizes and turn, as Mol\*'s shaders
  have it): its walls and edges, its gizmo (moved and turned as a box's; sized along an axis by its tips, across a
  section by that section, whole by a knob in each octant), its form (up, twist, three diameters), a box's sizes and
  twist kept when it becomes one and back. A box and an ellipsoid get a twist in their form and in `clipping.update`
  (`twist`: the turn about their up); a new direction keeps it.
- Prism: a cylinder cut is elliptic (its two diameters and turn): its walls and rims, its gizmo (a knob on each side
  halfway up and on each cap moving that side), its form (axis, twist, two diameters, height).
- Prism: a clipping gizmo's move and turn handles are the transform handles (Blender's move and rotate gizmo, as
  structure-move's) by cut: a plane's center sliding along its normal with one arrow toward the side cut away and its
  rings about the two axes in it, a cylinder and a cone without a ring about their own axis; the turn ball goes.
- Prism: a clipping gizmo is in two parts: its shape where the clip is (edges and rims as lines, knobs as round points,
  a set size on screen), its move, slide and turn handles at a fixed size on screen, as the move handles (parts of
  their own screen size, points parts).
- A particle list moved and turned as a whole: `Particle.transform` (the coordinates moved, the orientations turned,
  given where it turns a list without them) and the `TransformParticleList` decorator (`ParticlesRegistry`), as
  `TransformStructureConformation` for a structure.
- Prism: `particles.transform` and `particles.transform-reset` (the arguments of `structure.transform`, a command each);
  structure-move's handles and Move section take particle lists as well as structures.
- `Representation.State.previewTransform`: a transient transform after `transform`, for showing a representation moved
  while something is under way (a drag before its release), not a placement: applied on the GPU (a uniform; a merged
  object moves as one, set on what a representation draws; the culling, CPU and GPU, the sphere raster and the drawn
  bounds follow it), no instance, grid or culling layout touched, so setting it every frame costs the same for a
  mesoscale model as for a small one. GLSL and WGSL;
  `Visual.setPreviewTransform`. Prism's structure-move handles preview with it.
- `moveInstanceGrid` (core, `math/geometry/instance-grid`): an instance grid moved as a whole (its spheres and
  transforms; its cells and instance order shared, the grid it was moved from as `layout`, the move as `placement`);
  the GPU culling of merged objects keeps its split of a grid moved so, moving its cells instead.
- `Particle.transform` records the move (`ParticleList.movedFrom`), `Particle.getMove` gives the move between two lists
  of one line of moves; `moveTransformData` (graphics, `transform-data`) moves instances as a whole. A particle list
  moved (the release of the move handles) moves its targets' instances and grids instead of building them again, and
  keeps what it cached not by position (the target groups, the hierarchy colors' index).
- `Renderer.occlusionTest`: the test of `setOcclusionTest`; `Renderer.renderOverlay` and an `occluder` for
  `renderPick`: the scene's opaque depth an overlay is drawn and picked over.
