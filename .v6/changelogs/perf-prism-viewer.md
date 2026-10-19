## Prism: faster undo checkpoints, capped history

- Undo checkpoints share the transforms a step did not change and measure each transform once: at 1401 cells the size
  pass 16 → 0.2 ms per undoable step, 100 checkpoints 14.8 → 1.6 MB
- Past the hard size limit, the history stops snapshotting the scene for 30 s instead of on every agent call; the time
  budget counts the size pass
- Undo takes back the latest 10 steps (`options.ui.undoSteps`, 1–100); older steps stay in History, too old to undo,
  their checkpoints let go (`HistoryOptions.checkpoints`)
- Showing a snapshot (`snapshot.apply`, `.next`, `.previous`, an MVS story's steps) is no longer an undo step and takes
  no checkpoint: the steps before it stay undoable, and Undo brings back the snapshot that was current with the scene
- The sequence keeps its document over a trajectory's frames: a frame step no longer rebuilds the Seq\* document,
  empties the band or re-registers the 3D link (20 steps of a 10-model ensemble: 40 rebuilds → 0)
