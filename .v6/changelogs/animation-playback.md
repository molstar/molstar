## Animation playback: per-frame cost

- Snapshot transitions apply each frame once: a frame used to be applied again on every animation tick, and during the
  first frame the camera transition restarted on every tick (1cbs: 181 → 30 state updates over 3 s) (#2022)
- Playing or scrubbing a snapshot transition no longer logs the timing of every frame (exosome: 10,515 log entries in 10 s
  → 0); `State.setSnapshot` takes the update options (#2022)
- Snapshot transitions show a frame for its own duration: each frame showed during the slot of the one before it, so the
  second frame never showed and an MVS hold showed the frame after it (also when scrubbing) (#2022)
- [Breaking] The snapshot manager's `changed` no longer emits for the time of a playing transition:
  `events.animationTimeChanged` does, on every tick (the time also moves while a frame holds), so the snapshot list and
  the descriptions no longer re-render at the frame rate (#2022)
- Applying a snapshot with a long transition checks its frames' transformers and theme names once, not on every apply
  (an exosome story step: 18.8 → 0.02 ms) (#2022)
- A frame of a snapshot transition no longer builds every evaluated node's param definition to find its dependencies
  (exosome story: 4.6 → 3.7 ms per frame) (#2022)
- Playing a trajectory updates the model in place and keeps the custom properties that frames of the same topology
  share; every frame disposed their assets (11.5k atoms, the default preset: 13.2 → 9.5 ms and 8.9 → 4.4 MB allocated per frame) (#2022)
- With a focus, a trajectory frame is one state update: the focus representation committed a second one per frame,
  and disposing the plugin right after a frame threw in the focus interactions (#2022)
- A trajectory frame recreates no state cell: the structure's properties, the static components and the overpaint,
  transparency, substance, emissive, clipping and wiggle layers update in place, and the components share the new frame's
  root structure (root-level properties computed once per frame, not once per component) (#2022)
- `Structure.remapModel(m, parent?)` takes the parent already remapped (#2022)
- A representation update that changes the representation's own params but neither the structure nor the themes' params
  keeps the color and size themes (an MVS frame changing opacity or emissive: exosome story 3.7 → 3.0 ms per frame) (#2022)
- A playing trajectory commits a state update only on the ticks that show a new frame; the ticks in between went through
  the task manager and toggled the busy state (#2022)
- A trajectory frame without a selection no longer clears the selection marks of every representation (#2022)
- The state interpolation animation updates a node only when its interpolated params change (it updated every node
  that differs between the snapshots on every tick) and matches the snapshots' nodes by ref in linear time (#2022)
- The overpaint, transparency, substance, emissive, clipping and wiggle bundle layers below a representation are not
  applied again when the representation updates and their layers and structure are unchanged (#2022)
- A snapshot transition ends on its last frame also when that frame is shorter than an animation tick (#2022)
- A visual checks its render object for another graphics context's resources once per context, not on every update
  (exosome story: 3.0 → 2.1 ms per frame) (#2022)
