## Seq* (next): editing

- `@molstar/seqx-viewer`: `SpecModel` applies patches as commits with the language's checks of a change, keeps recent
  roots for a host's rebase of a stale patch through `merge3`, treats previews as patches, and tells the session the
  moves on every change of head (undo and redo included). The history goes through Seq\*'s history port; the adapters
  `uiHistory`, `hostHistory` and `HistoryParticipants` are at `@molstar/seqx-viewer/ui-history`, and a mounted viewer
  gets its own stack. A new `editing` behavior provides commands, each declaring its change kind; edits can nest in a
  host's run (`within`). Language 0.3.0 adds the `snapshots` level and `track.hidden`.
