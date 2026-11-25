# Plan: Gizmo rotation controls, and clip visuals over solid interiors

## Status

Proposed. Causes below are from reading the code; none is reproduced yet (step 1 of each part does that).

Four items:

1. Z-fighting of solid-interior caps with the clip objects' visuals.
2. The rotation rings of the transform handles (structure-move) feel awkward; no free rotation.
3. The clip gizmos' turn handle: not like the camera's trackball, and its ring does not say what it does.
4. The clip gizmos' slide arrow ("normal vector") fades out as it turns toward the viewer: an option, off by default.

## What other tools do

From memory of Blender, Unity, Maya, 3ds Max and three.js's `TransformControls` (check the details by trying them):

- **Three axis rings.** Unity and three.js: the drag along the ring's tangent at the press, at a constant rate; the
  cursor need not stay on the ring. Blender: by default the angle around the gizmo's center on the screen.
- **Free rotation is the area inside the rings**, bounded by a thin circle facing the camera (the silhouette of a
  ball). Drag there and the item turns as a trackball does. Blender, Unity, Maya, 3ds Max and three.js all do this.
- **An outer ring facing the camera** turns about the view axis (roll). Same tools.
- **A plane's direction** (ParaView's plane widget): the normal arrow's tip is dragged to aim the plane. It needs the
  arrow that item 4 hides.

## 1. Z-fighting

**As built.** `Mesh.Params.depthOffset` (hidden; the object uniform `uDepthOffset`, Å toward the camera; only
`gl_Position.z` changes, so the pixels stay), the clip volumes shape's walls at 0.002 of the radius of what they cut. In
`apply_occlusion` a fragment is occluded when behind by more than a thousandth of its distance. Tests: the headless
`webgpu-clip-walls.test.ts` (a plane's walls against the same walls well in front of the cap, scene and overlay: 395 and
341 pixels differ without the fixes), visual cases `overlay/clip-walls/**`. Only the plane showed it in those cases; the
faceted walls (sphere, cylinder, cone) looked even at this offset.

**As planned:**

**Cause.** A cap's depth is the clip surface itself: `fragmentDepth = mix(calcDepth(viewPosition), gl_FragCoord.z,
0.0001)` (`mesh.frag`), and the walls (`clipWallsMesh`) are drawn on that surface. There are two meeting points:

- the clip volumes shape (`ClipVolumesShape3D`), in the scene: the depth test decides per pixel;
- the gizmo's body, in the overlay: `apply_occlusion` compares view z (`fragmentViewZ < occluderViewZ`), so x-ray alpha
  flips per pixel.

A sphere, cylinder or cone is also not an exact match: the walls are faceted, the cut is analytic, so they cross by up
to a facet's sagitta. A bias has to exceed that, or the walls have to be finer.

**Fix.**

- Overlay: a tolerance in the occlusion test, a fragment occluded only when behind by more than a few of a pixel's world
  size at its depth (GLSL and WGSL). It needs no new state.
- Scene: a hidden mesh param, the object uniform `uDepthOffset` (Å toward the camera, in the vertex stage), set on the
  walls by the clip volumes shape, scaled with `reach`. A uniform, as `uScreenSize`: no new programs. Not polygon offset:
  there is none in the neutral API, and WebGPU's `depthBias` is pipeline state.
- Facets: measure the sagitta at the current tessellation. Raise the tessellation or the offset, whichever is cheaper.

**Steps.** Visual cases first, failing: each shape's walls over a solid-interior structure, drawn and as a gizmo's body,
on every profile. Then the two fixes.

## 2. The rotation rings

**As built.** `RingRule` (`gestures.ts`: `grab`, `tangent`, `across`), default `tangent`, to compare in the app;
`rotateAbout` takes a rule (the tests pass each). `tangent` and `across` go at the trackball's rate (`turnRate`),
normalized on the screen, each falling back to the other where its way vanishes (`MinOnScreen`), then to
`angleOnScreen`. Free rotation: `transformHandles` `free` (with `rotate` by default), a `turn` group last (the others'
ids stay), two billboard parts: the circle at the rings' radius (`RingRadius`; lines with `lines`) and a hit part of
a thick circle and a disk `FreeAreaDepth` (1.25 `size`) behind the origin. No fill on hover (as the clip ball: the
circle lights). The roll ring is left until missed. Open: pick the rule, remove the others and `PlaneFacing`.

**As planned:**

**Why it feels awkward.** `rotateAbout` follows the grabbed point (`angleOfGrab`), so the feel depends on the press and
on the gizmo's place:

- the rate depends on how far the cursor is from the center, and gets very fast near it;
- the cursor is projected along the camera's direction, not the ray through the pixel: in perspective a gizmo off the
  view's center sees its ring from another angle than the projection assumes;
- the rule switches at `PlaneFacing` (projection to tangent) as the ring turns edge-on.

**Change.** Tangent at the press, constant rate: the grabbed point's tangent on the screen, taken at the press and kept
for the drag; the angle is the drag along it times a rate. No projection, no distance from the center. Near edge-on, where
the tangent on the screen is short, it is normalized; where it vanishes, the perpendicular of the axis on the screen.

- The press still decides the direction, and has to: the near and far sides of a ring turn opposite ways on the screen for
  the grabbed part to go with the cursor. Nothing else depends on it.
- The rate: the trackball's (`getRotateFactor`), through `DragFrame`, so a gizmo and the camera turn alike.
- Given up: the grabbed point no longer stays under the cursor.

**Three rules to try in the app**, behind a constant, the losers then removed (`PlaneFacing` with the grabbed point's):

- the grabbed point followed (today's);
- the tangent at the press (above);
- one drag direction a ring: the perpendicular of its axis on the screen, whichever part is grabbed, so the press
  decides nothing. Where the axis points at the viewer that perpendicular vanishes: there, the tangent at the press.

**Free rotation.** `Handle` `turn` becomes a trackball (part 3's gesture). `transformHandles` gets:

- a thin circle facing the camera, at the rings' radius;
- the area inside it, picked as free rotation: a disk in the billboard's frame, pushed back along the view by more than
  the radius, so every other handle wins the pick (`GizmoLayout` part with a billboard of its own: `facesCamera` is per
  layout now);
- a faint fill while hovered;
- later, if missed: the outer ring for roll.

The disk takes camera drags over the gizmo's area. Structure-move is a mode entered on purpose, so that is fine there.

## 3. The clip gizmos' turn handle

**As built.** `gestures.turn` is one step, about the axis across it; the controller applies a `turn` to the last value
and the step (`d.next`), every other handle to the press and the total. The rate is `trackballTurnRate`
(`drag-frame.ts`, twice `rotateSpeed` over the canvas's CSS height), through `DragFrame.turnRate`; a test drives the real
`TrackballControls` and compares. The ball: two great circles in the move handles, an outline facing the camera in
a part of its own (`GizmoPart.billboard`; the frame stays the cut's). Both grabs were tried in the app:
the lines stay, the disk inside the outline went. The primitive remains for part 2's free rotation (a facing hit part).

**As planned:**

**Why it differs from the camera.** `turn` composes two rotations from the drag's total: about the screen's up by `dx`,
about its right by `dy`. Any path to the same point gives the same pose, at a fixed 0.5° a pixel. The trackball
(`rotateCamera`) turns per move about the axis across that move, accumulating, at `rotateSpeed`.

**Gesture.** The same: per move, about the axis across the move in the screen's plane, at the trackball's rate, without
damping.

- Gestures take the drag's total from the pose at the press. For a `turn` the controller applies each step to the
  previous result instead (`d.next`, not `d.from`), keeping `from` for cancel and commit. `poseCut` already goes through
  a normalized quaternion.
- Tests: a closed drag path leaves a twist (as the camera), the rate equals the trackball's, cancel restores `from`.

**Visual.** A ring around the axis says "turns about this axis", which it does not. Proposed: a small ball drawn as
wireframe: the circle facing the camera, plus two great circles turning with the cut, so it reads as a ball being turned.
The drag is on those lines (thick hit area); the area between them stays the camera's. Why not the disk of part 2: a cut
often sits at the view's center, and a 90 px disk there takes the camera's drags while the section is open.

Try both in the app, behind a constant: the lines, and the disk of part 2.

- Turning about the axis itself: a plane has none to turn; a box's or cylinder's twist comes with the trackball's roll
  by circling, as the camera's.

## 4. The slide arrow's fade as an option

**As built.** `GizmoControllerOptions.fadeSlides` (the library's default on, as structure-move wants), `setFadeSlides`
while shown; the clip gizmos' off, from `clipping.settings` `fadeNormal` (`ClipGizmos`: the session's, not saved: the
scene's settings are the component manager's options, a plugin type that has no place for it). Unfaded and end-on
(|cos| past `Fade.slide`), a press on the slide is the move handle's (`#grabbedGroup`), where the gizmo has one: no
change to the picks or hit areas. A plane's knob is its slide, so it has nothing to lose.

**As planned:**

The arrow is always drawn, by default also when it points at the viewer. Today the controller fades a slide as its axis
turns toward the view (`Fade.slide`, |cos| 0.9 to 0.98), hidden and not picked at the end.

`GizmoControllerOptions.fadeSlides` (default off; the clip gizmos' from `clipping.settings` `fadeNormal`, the scene's,
saved with it). The planes' fade (`Fade.plane`) stays. Structure-move's three arrows keep their fade: with three, one
along the view is what makes the knob behind it reachable.

- Unfaded, seen end-on: the arrow's head covers the knob at the center, in the hit area too (it starts past the knob,
  toward the viewer). The knob has to stay grabbable: its hit area drawn over the arrow's (a pick priority, or the
  arrow's hit area not covering the knob's disk end-on). A drag on the head then slides at `MaxSlideGain`.
- Tests: the controller both ways (faded, hidden and unpicked; unfaded, drawn and picked end-on), the knob picked end-on
  (headless); the setting's command. `viewer/_spec/clipping.md` §4 changes with it.

## 5. A gizmo as a list of parts

**As built.** `GizmoLayout` is `parts` (`GizmoPart`: a mesh, lines or points; `role` `handle`, `hit` or `body`; its
own `screenSize` and `billboard` over the layout's), the groups, the frame, and `extras`. It replaces the slots
(`body`, `handles`, `hitArea`, `handleLines`, `handlePoints`, `fixedHandles` with its hit area and `facing`): one
path in `#make`, so a new kind of part (part 2's facing disk) is data, not a case. Every handle and hit part stands in
for one shape of no geometry (`ShapeBuilder.lociShape`), whose labels are the tooltips'; a shape representation's
`getAllLoci` is that shape's too. `GizmoExtra` (built with the canvas's graphics and given that shape) is the way in
for a representation of the app's own: placed, faded and removed with the gizmo, a handle where its loci are the
handles' (`#hover` checks the shape), drawn only otherwise. `transformHandles` gives its meshes and `parts`.

## 6. Blender's transform gizmo, for structure-move and the clip gizmos

**What Blender does** (its move and rotate gizmos combined), at rest:

- the center: a small circle facing the camera, no knob: moves in the screen's plane;
- an arrow head on each axis, no shaft, just past the outer circle: slides along it;
- each axis ring's front half only (the half toward the viewer), at a radius inside the outer circle;
- the outer circle, white, facing the camera: turns about the view axis (roll);
- the area inside the rings: the trackball, not drawn; hovered, a faint fill lights and a tooltip says so.

During a drag, what is not dragged is hidden, and guides show what the drag does:

- a ring: only it, whole; a line along its axis across the view; a radius at the start angle (the swept sector);
  a dashed line from the center to the cursor;
- the outer circle (roll): only it; radii at the start and the current angle; the dashed line to the cursor;
- the trackball: nothing;
- an arrow: a line along its axis across the view, a dot and a grey copy of the head where the drag started; the other
  arrows stay (with shafts), the rings and circles go;
- the center: a grey copy of the circle where the drag started; the arrows stay;
- a plane square (Blender's move gizmo has them; ours adds them to the combination): the lines along its two axes
  across the view, crossing at the start, and a grey copy of the square there; the arrows stay.

**Ours.** The transform handles become that layout (`transformHandles`), used by structure-move (move and rotate, the
small plane squares) and by the clip gizmos in place of their move knob, slide arrow and turn ball; a clip keeps its
own sizing handles (the box's edges, corners and face knobs, the rims, the cone's rim, the sphere's ring), so it is
Blender's combination with our scale. The parts (§5) carry it; what is missing is in the controller.

1. **The controller.**
   - Parts shown by state (`GizmoPart.when`: `hover` of its group, `drag` of a kind): the trackball's fill (a drawn disk
     at the rings' radius, faint), the arrows' shafts while a move is dragged. `GizmoPart.alpha` for the fill.
   - The front half of a ring (`GizmoPart.front`): its representation's clip, a plane through the gizmo's origin across
     the view, set as the camera moves (`#fadeOnCameraMove`; a clip-only update is uniforms). The vertex stage clips
     after the screen size and the billboard turn (WebGL `assign-position`, WebGPU `mesh.wgsl`, `lines.wgsl`), so the
     plane in the world is right. The dragged ring: its clip off.
   - Hidden during a drag, by the handle dragged (`Fade.dimmed` becomes hidden, not picked): a turn (`rotate`, `roll`,
     `turn`) hides every other handle, a move (`move`, `slide`, `plane`) the turning ones. The clip gizmo's sizing
     handles: hidden during either.
   - Guides: one lines representation of the controller's, in the overlay, not picked, its geometry made each move
     from the drag (the handle, the frame at the press, the value now, the cursor): axis lines across the view, the
     radii, the dashed line to the cursor (dashes as segments). Start ghosts: the dragged group's parts again at the
     frame at the press, grey, the other groups hidden (transparency), made at the press.
   - A drag's cursor (the double arrow, the four-way): the canvas element's CSS cursor, if the plugin allows it;
     otherwise the dashed line says it. To check.

   **As built** (all but the cursor). `when`, `alpha`, `front` on `GizmoPart`; `#fade` hides by kind (`Turning`:
   `rotate`, `turn`; `Moving`: `move`, `slide`, `plane`; a sizing drag dims as before) and calls `#front` and
   `#reveal`. The front clip is always one plane (clip on and off would change `dClipObjects`, another shader): whole,
   it is pushed `WholeOffset` behind; updated only where it changed. Guides (`guideLines`, exported, tested alone): in
   the world, through the origin at the press, `reach` twice the eye's distance; the ring's radius to the hover's
   grab point; dashes `DashPixels` at the origin's depth. Ghosts: every handle part shown at rest (not `when`), the
   other groups transparent, at the layout's frame at the press. The roll's guides come with `roll`.
2. **The gestures.** `roll` (a new `Handle` kind): the angle the cursor goes around the gizmo's origin on the screen
   (`angleOnScreen`'s rule, the axis toward the viewer). The rings by the rule chosen in §2 (Blender's is the angle
   around the center on the screen, as `grab` facing the viewer; try it beside `tangent`).
3. **The handles** (`transformHandles`, in units of the arrows' length): the center circle (radius 0.1, its hit a disk
   pushed toward the camera past every other part, so it wins over the rings' near sides seen edge-on: the knob's
   problem today), the arrow heads (0.85 to 1, hit areas thick), the rings (0.65, `front`), the outer circle (0.8,
   `roll`), the trackball's hit disk behind everything and its fill (`when: hover`), the plane squares (small, near
   the center between their axes). Options: which axes slide and turn, the center's handle (`move`, or a slide: a
   plane cut's), plane squares, roll, the trackball. Axis colors as Blender's (X red, Y green, Z blue); white
   circles.
   **As built** (2 to 4). `roll` (`gestures.roll`) by steps as `turn` (the controller's step frame: `grab` the last
   cursor from the origin), so it goes on past half a turn; its guides the radius at the press and the dashes.
   `transformHandles`' groups: center, slides, planes, rings, scale, roll, trackball (the earlier ids kept); its parts:
   the heads, squares and knobs (`handles`) and their hit area, the shafts (`when` a move), the rings and their hit
   area (`front`), the circles and the facing hit area (center disk at `CenterHitDepth` 1.25 toward the camera, past a
   scale knob's 1.16; the trackball's at `FreeAreaDepth` 1.25 behind), the fill (`when` the trackball hovered). The
   headless tests pick the center over the edge-on ring and the x ring by its front half (seen along the diagonal of x
   and z the z ring's near half lies on the x ring's far half: the x ring alone there). Structure-move takes the
   defaults as they are. Open: the cursor, the fill's discoverability.

   **The rule chosen** (from the app): `around`, Blender's, the ring following the cursor around the center, by steps
   (`Stepped`: `turn`, `roll`, `rotate`; their `frame.grab` the last cursor from the origin), so on past half a turn;
   `tangent` where the ring is edge-on (`MinOnScreen`), then `across`. `grab` went (its grabbed point does not hold
   by steps). Its guides as Blender's: the radius at the grab and now (the start turned by the frame now over the frame
   at the press), the dashes to the cursor and a double arrow across them there (`ArrowPixels`); a trackball's turn a
   cross at the center. Each move's guides are built into the last ones' `Lines` (`guideLines`' `into`): a shape
   representation keeps its first geometry's buffers, so a new geometry each move left the guides drawn at the press.
4. **Structure-move**: move and rotate with plane squares, no scale.
5. **The clip gizmos**: the move handles from `transformHandles` at `MoveScreenSize`, their groups first, the shape's
   after. By cut:
   - plane: the center slides along the normal (moving it within itself changes nothing), one arrow toward the side cut
     away, no squares, the rings about its two in-plane axes, roll and trackball;
   - sphere: an ellipsoid (§7): everything, as the box;
   - box: everything;
   - cylinder, cone: no ring about their own axis (it changes nothing);
   The ball of §3 (two great circles and the outline) goes: the trackball area replaces it.

   **As built.** `gizmoGeometry`'s `moveParts` (`transformHandles(1, moveOptions(object))`, each part at
   `MoveScreenSize`), its groups first, `buildShape` the lines and points after; the ellipsoid all of them (§7). A cone's
   axis is its opening (−Z): rings about x and y. The geometry test checks the kinds by cut, every handle in a part and
   picked through a hit part, the move groups first.
6. **Tests and docs.** Controller: parts by state, the front clip set as the camera turns, hiding by kind, the guides
   made and removed. Headless: the center picked over an edge-on ring, the back half of a ring not picked, the
   trackball's fill on hover. Visual: the gizmo at rest and in each drag (overlay suite). `structure-move.md`,
   `clipping.md` §4, the changelog entry.

**Open:** the cursor; whether the trackball fill alone (no outline) is found without Blender's
habit (the tooltip says it). White handles on a white background were very hard to use (the maintainers', from the app):
grey instead (the circles, the center, the guides: `0x8c8c8c`) and every line thicker (`lineWidth` 4, guides 3,
the mesh tubes 0.024 `size`); a dark outline was tried and dropped (it did not read well).

## 7. The sphere cut is an ellipsoid

The shaders clip a sphere as an ellipsoid: its half axes the scale's halves, turned by its rotation (`localSD` type 2:
`length(t / size) - 1` times the smallest axis; `clipSphereClass` bounds it so). Prism draws and edits it as a round
sphere of `scale[0]`: the walls (`clipWallsMesh`: `addSphere` of `scale[0] / 2`), the edges (`clipEdgesMesh`: three
great circles of that radius), the gizmo (one rim facing the camera, `facesCamera`, no turning), the form (one
diameter) and `reshape` (the three axes `across`). A scale not uniform (the API's `scale`, MolViewSpec's) cuts
another volume than the one shown. Ours: an ellipsoid everywhere, a sphere its round case.

1. **Walls and edges** (`clip-volumes.ts`): the walls a unit sphere scaled by the half axes (`addPrimitive` with the
   scaling: its normals by the inverse transpose, to check); the edges the three principal sections, ellipses in the
   frame's planes.
2. **The gizmo** (`gizmo/geometry.ts`): turned as the box (§6's rings, roll, trackball; not `facesCamera`); its sizing
   handles, the box's on its surface:
   - a knob at each axis's two tips (points): sizes along that axis (`faces`, the opposite tip stays, as the box's face
     knobs);
   - the three principal sections (lines), the edges joining the tips: each sizes the two axes of its plane (`scale`
     with those `axes`, as the box's edges across an axis);
   - a knob at each octant's middle on the surface (the point along (±a, ±b, ±c) / √3, a box corner's place): sizes it
     whole, keeping its proportions.
   The rim facing the camera goes (an ellipsoid's outline is not one of its sections and changes with the view).
3. **The form** (`ui.tsx`): three diameters (`size`, as the box's) and its turn (`up` and `twist`, as the box's);
   round when placed (`placeCut` from a bounding sphere, unchanged). The box and the ellipsoid get `twist` (degrees):
   the turn about `up`, from the shortest turn of the local +Y onto `up` (`rotationTo`) to the rotation, measured on
   the local +X; set, the rotation is that shortest turn then `twist` about `up`. A form changing `up` alone keeps the
   twist (read, then set with the new `up`: `turnedTo`'s shortest turn from the old `up` would change it). Not for the round shapes (cylinder, cone, plane: a turn about their axis changes nothing).
4. **`reshape`**: a box and an ellipsoid keep each other's three sizes; a cylinder's diameter and height give its axes
   (diameter, height, diameter); a plane or cone gives a round one, as now.
5. **Names**: the shape stays `sphere` (the clip type, MolViewSpec's, saved states); its label and the API's
   descriptions say ellipsoid where it matters (`clipping.add`'s `scale`: three diameters).
6. **Tests and docs**: the walls' and edges' extents per axis, a tip knob sizing one axis, a section its two, an
   octant knob all three in proportion, `reshape` box ↔ sphere, a headless pick of the tips, `twist` read back as set
   (and kept by a change of `up`); `clipping.md`, the changelog entry.

**As built:** all of it but the turn: the ellipsoid moves, slides and turns with the box's handles of today (the
knob, the arrow, the ball); §6 gives it the rings. The walls' normals: `twoFaced` makes its own, a triangle each.

**The cylinder too** (elliptic in the shaders now: its radii the scale's x and z halves, `localSD` and the segment
intervals scaled as the ellipsoid's): its walls and rims, a knob on each side halfway up (`faces`), the rims sizing
both radii, its form (axis, twist, two diameters, height), `reshape` keeping its three sizes. The cone stays round.

## Order

1. Z-fighting: on its own, a bug.
2. Item 4: small.
3. The trackball gesture (item 3's), then its visual.
4. The parts (§5), then the rings and free rotation (item 2), on the same gesture.
5. Blender's gizmo (§6): the controller, the gestures, the handles, then structure-move, then the clip gizmos.
6. The ellipsoid (§7): the walls, edges, form and `reshape` before the clip gizmos of §6 (the sphere's handles are made
   once, as an ellipsoid's).
