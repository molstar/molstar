## Seq* (next): fixes from use

- The vertical overview is back in the Globins alignments, with stories of the rows rail and the overview.
- A filling rows band no longer grows when its frame is resized.
- A label that fits neither inside nor beside its mark is cut with an ellipsis, keeping at least `labels.minChars`
  letters (4); the hover card keeps it whole.
- Ring letters sit on their cells; a moved label in a lower lane is placed instead of dropped; an empty track label
  draws nothing.
- Zooming 10,000 rows builds the new level across frames, the rows' pending placeholder drawn meanwhile, within
  `rows.buildBudget` cells a frame.
- A null inline `strand` reads as none; the cursor card stays on a re-hit; a folded clade's summary row is named to
  assistive technology.
