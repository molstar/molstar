## Seq* (next): multiple sequence alignments

- Aligned FASTA, Stockholm and Newick readers; the `alignment` and `tree` kinds.
- A `rows` layer: a tree beside the rows (widening the label column), clades folded into consensus summary rows with a
  partial selection mark; pinned and hidden rows; order by tree, data or a field; row zoom on `]` and `[`, scrolling
  with the wheel, Shift+click extending in display order; first and last residue numbers; held or lazy rows (up to
  10,000,000); image tiles from a shared reduction.
- A sequence `logo` layer; the rows in the text figure; SVG export carries image tiles as PNG `<image>` elements.
- Globins examples and stories.
- Hovering a tree node or a folded clade's summary row lights the clade: its branches and one band over its rows.
- A text figure's ruler draws one label a position (a group's end number wins) and the line numbers sit on the ruler's
  label baseline.
- A row's first and last residue numbers follow the view (the residues in view, never a gap).
- Rows layers show data columns beside each row's header (`columns`: `field`, `side` left or right, `as` text, number,
  bar or chip, `width`, `title`), from the alignment's own fields (Stockholm `#=GS` tags, an inline table's extra
  columns), the measures `length`, `gaps` and `identity`, or a records file joined by row key; a column's title orders
  the rows by it (again: the other way); on narrow views the columns yield, right ones first. The Globins family shows
  species and identity.
