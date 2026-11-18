## Seq* (next): several ranges in one selection, followers that lead

- Several ranges in one selection: Cmd/Ctrl+drag adds a range, Alt/Option+drag removes one, Shift+click, Shift+drag
  and Shift+arrows extend only the range the anchor is in (`selection.extendScope`: `range`; `selection` for the old
  whole-selection extend); Cmd/Ctrl+Shift+arrows start a new range at the cursor; `selection.add` and
  `selection.subtract` commands take a range or an address. Esc during a drag restores the selection.
- A follower section leads: a pan, zoom or pinch on it moves its leader (`follow.followerMoves`: `lead`, the default;
  `ignore` and `detach` kept), so a protein and its gene move together from either section.
- The hover no longer blinks at far zoom (a cell narrower than a device pixel keeps an outline).
