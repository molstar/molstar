## Prism: faster undo checkpoints, capped history

- Undo checkpoints share the transforms a step did not change and measure each transform once: at 1401 cells the size
  pass 16 → 0.2 ms per undoable step, 100 checkpoints 14.8 → 1.6 MB
- Past the hard size limit, the history stops snapshotting the scene for 30 s instead of on every agent call; the time
  budget counts the size pass
- Undo takes back the latest 10 steps (`options.ui.undoSteps`, 1–100); older steps stay in History, too old to undo,
  their checkpoints let go (`HistoryOptions.checkpoints`)
- Showing a snapshot (`snapshot.apply`, `.next`, `.previous`, an MVS story's steps) is no longer an undo step and takes
  no checkpoint (`undoable: false`): the steps before it stay undoable, Undo is refused while one is applied, and Undo brings
  back the snapshot that was current with the scene; the assistant asks before it shows one
- The sequence keeps its document over a trajectory's frames: a frame step no longer rebuilds the Seq\* document,
  empties the band or re-registers the 3D link (20 steps of a 10-model ensemble: 40 rebuilds → 0)
- Plugin: `StructureInfo` walks the tree for the max model index only when a model or a structure is created or
  removed, not for every created object; a removed model lowers the max at once, and a model or structure replaced in
  place (another frame) gets it at once (100 structures: first overpaint over
  all 101 → 39 ms, transparency 79 → 16 ms)
- Plugin: the measurement manager's sync, on every state change, searches the tree for the measurement group once per
  tree instead of six walks (none when there is no group); it still reports an update on every sync, as before
- Style layers (the Layers list) read all five kinds of a representation in one walk: 300 representations, a read with
  nothing changed 1.15 → 0.24 ms
- Superpose's picked atoms are counted once per selection, not again on every state change (an all-atom selection of
  186k atoms: ~4 ms a change while the section shows)
- `trajectory.frame` is `record: false`: a scrub through a trajectory writes no record entry per frame (as
  `snapshot.seek`)
- The sequence's chain index lists water chains without indexing their atoms (nothing links them): 186k atoms with 12k
  waters, 1.3 → 0.85 ms per new structure
- The sequence's link takes a 3D select over several structures once, not once per structure
- Fix: the sequence of the followed structure reads its CIF block again: SIFTS numberings onto UniProt (mapping,
  alignment), a predicted model's styling and the auth labels of unobserved residues were silently missing
- Fix: the marks' paint lets go of a removed representation (its structure and theme) at once, not at the next paint
- Fix: a large sequence hover or selection over ions of one Mol\* unit lights the ions it touches, not every chain
  between them
- Prism: a story's own navigation (a description's `!apply-snapshot` / `!next-snapshot` links, a clickable primitive's
  `snapshot_key`, the keyboard's next and previous) goes through the snapshot operations: earlier steps stay undoable;
  Undo to a scene shown with no snapshot marks none current
- Plugin: `managers.snapshot.setNavigator({ applyKey, applyNext })` routes the manager's `applyKey` / `applyNext` through
  a host's own path (unset: as before), and `clearCurrent()` marks no snapshot current
