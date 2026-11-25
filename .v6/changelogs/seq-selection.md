## Seq* (next): selection, highlights and the keyboard cursor

- Click selects (again: clears), Shift+click extends from the anchor, Cmd/Ctrl+click toggles, a drag on the data
  selects a range (whole items when it begins on one), Esc or a click on empty space clears; the arrow keys move a
  keyboard cursor (Shift extends, Space toggles, Enter focuses), shown once the keyboard moved it. The hover, the
  selection and named highlight channels are drawn and carried onto linked sections (a protein on its gene) in the same
  look, lighter. The selection is one band per section by default (`column`; `tracks` draws a box per track), a soft
  tint with a thin outline; the hover stands just outside its mark. SVG exports draw the selection and highlights.
  New Storybook page: Features/Selection.
