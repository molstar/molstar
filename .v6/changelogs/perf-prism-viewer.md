## Prism: faster undo checkpoints, capped history

- Undo checkpoints share the transforms a step did not change and measure each transform once: at 1401 cells the size
  pass 16 → 0.2 ms per undoable step, 100 checkpoints 14.8 → 1.6 MB
- Past the hard size limit, the history stops snapshotting the scene for 30 s instead of on every agent call; the time
  budget counts the size pass
