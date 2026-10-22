## Animation playback: per-frame cost

- Snapshot transitions apply each frame once: a frame used to be applied again on every animation tick, and during the
  first frame the camera transition restarted on every tick (1cbs: 181 → 30 state updates over 3 s)
- Playing or scrubbing a snapshot transition no longer logs the timing of every frame (exosome: 10,515 log entries in 10 s
  → 0); `State.setSnapshot` takes the update options
- Snapshot transitions show a frame for its own duration: each frame showed during the slot of the one before it, so the
  second frame never showed and an MVS hold showed the frame after it (also when scrubbing)
