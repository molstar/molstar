## Seq* (next): bookmarks from the selection, the history toolbar, links before load

- New commands: `selection.bookmark` (Cmd/Ctrl+D with the view focused) adds a bookmark for the selection, one undo
  step; `refs.label` labels a bookmark; `highlight.clear` clears one highlight channel or all; `track.add` adds a track.
- Back and Forward as a toolbar a host places (`SeqViewerUI.toolbar`), disabled when there is no step, their keys in
  the tooltip; the default viewer no longer bundles a toolbar.
- A link followed before a spec is loaded waits and runs when the spec arrives; another jump or a spec that no longer
  names it ends the wait.
