## Seq* (next): section layouts and the text layout

- Sections are laid out by a layout: `tracks` or `text` (lines of whole groups of ten wrapped to the width, numbers per
  line and per group by a layer's numbering, thin bands under the letters collapsed where empty, scrolling inside a box
  capped at `text.maxHeight`, `pieceBreak: 'line'`, a `frame` layer outlining regions across tracks).
- A text section is navigated through its layout: Cmd/Ctrl+wheel, pinch, - and = step the letter size about the
  pointer's line, 0 returns to the preset; the wheel and Space+drag scroll lines; ↑↓ move the cursor a line, Shift+↑↓
  extend, Shift+Alt+↑↓ scroll; a drag selects in reading order and scrolls the box when held past its edge.
- `selection.copy-letters` (Cmd/Ctrl+C) copies the selection's letters as shown and reports how many; a selection made
  elsewhere scrolls a text to it (`reveal.linked`); the viewer config's `cursor.reveal` is now `reveal.cursor`.
- A `window` behavior (drags, slider steps); `renderSvg`'s `window` option; a text view's export holds every line.
