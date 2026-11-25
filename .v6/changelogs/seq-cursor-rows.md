## Seq* (next): the keyboard cursor reaches an alignment's rows

- ↑ ↓ step rows (the pinned band, a folded clade one stop; from the section's other tracks onto its rows), Page Up and
  Page Down a page, Cmd/Ctrl+↑ ↓ the first and last row, ← onto a row's header and → back, Shift+arrows extend a block
  and Cmd/Ctrl+Shift+arrows start one, Space toggles a cell or (on the header) folds the row's clade, Enter opens a
  folded clade or focuses the row; a row never loaded is stepped over by its index, a hidden row's cursor stands on the
  nearest shown. The cursor is drawn on its row, announced and carded by the row's label. On a grid, ↑ ↓ step its y,
  which scrolls to show it. Page Up and Page Down are now `cursor.page-up`/`cursor.page-down` (`rows.page-*` keep no
  default key).
