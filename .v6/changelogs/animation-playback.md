## Animation playback: per-frame cost

- Snapshot transitions apply each frame once: a frame used to be applied again on every animation tick, and during the
  first frame the camera transition restarted on every tick (1cbs: 181 → 30 state updates over 3 s)
- Playing or scrubbing a snapshot transition no longer logs the timing of every frame (exosome: 10,515 log entries in 10 s
  → 0); `State.setSnapshot` takes the update options
- Snapshot transitions show a frame for its own duration: each frame showed during the slot of the one before it, so the
  second frame never showed and an MVS hold showed the frame after it (also when scrubbing)
- [Breaking] The snapshot manager's `changed` no longer emits for the time of a playing transition:
  `events.animationTimeChanged` does, on every tick (the time also moves while a frame holds), so the snapshot list and
  the descriptions no longer re-render at the frame rate
- Applying a snapshot with a long transition checks its frames' transformers and theme names once, not on every apply
  (an exosome story step: 18.8 → 0.02 ms)
- A frame of a snapshot transition no longer builds every evaluated node's param definition to find its dependencies
  (exosome story: 4.6 → 3.7 ms per frame)
- Playing a trajectory updates the model in place and keeps the custom properties its frames share; every frame
  disposed their assets (11.5k atoms, the default preset: 13.2 → 9.5 ms and 8.9 → 4.4 MB allocated per frame)
- With a focus, a trajectory frame is one state update: the focus representation committed a second one per frame,
  and disposing the plugin right after a frame threw in the focus interactions
- A trajectory frame recreates no state cell: the structure's properties, the static components and the overpaint,
  transparency, substance, emissive, clipping and wiggle layers update in place, and the components share the new frame's
  root structure (root-level properties computed once per frame, not once per component)
- `Structure.remapModel(m, parent?)` takes the parent already remapped
- A representation update that changes neither the structure nor the themes' params keeps the color and size themes
  (an MVS frame changing opacity or emissive: exosome story 3.7 → 3.0 ms per frame)
