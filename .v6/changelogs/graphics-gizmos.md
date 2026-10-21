## Graphics: 3D gizmos

- `@molstar/graphics/util/gizmo`: handles for direct manipulation in the 3D view, drawn outside the state tree. A
  `GizmoController` on a `Canvas3D` shows an app's items (a `GizmoBinding`: their layouts, what a drag makes of them,
  live previews, one commit, cancel), takes the trackball's drags while a handle is hovered, and cancels on Esc, on a
  change under the drag or a commit not taken. Gestures keep the part held under the cursor (a move, a slide along a
  foreshortened axis, a scale by distance from a pivot, a turn about an axis by the angle around it); primitives (knob,
  tube, arrow, ring, facing the camera) and the standard transform handles (axis arrows, plane squares, axis rings). Prism's
  clipping gizmos are built on it.
- Gizmo controllers on one canvas share the trackball's drags (saved once, put back by the last); a commit that fails
  or is refused cancels its own drag; no drag starts while a commit runs; a ring turns by the point grabbed on it.
- `Canvas3D.holdFrames(until, maxMs)`: frames wait (the last one shown) until a change made in steps is whole.
- A clip object's drawing honours its `transform` (`clipObjectFrame`).
