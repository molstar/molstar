## Seq* (next): accessibility

- The exported SVG is a `role="figure"` with a title and a description (sections, systems, track counts), a named group a
  track and a named item a range, variant, secondary element, row or pair (capped at 200 a track and 2,000 a figure,
  "n more not named"); its boxes never take the pointer; ids unique per export by default (`idPrefix` to override).
- The viewer announces the cursor, the hover (after it rests) and the selection in a live region, selection first, quiet
  during a drag; a repeated message is announced again.
- Forced colors apply to any theme source; `rows.order-next` (Shift+O) cycles the rows' order from the keyboard.
