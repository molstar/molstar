## Seq* (next): rings and blocks

- A section can be laid out as rings (`layout: 'rings'`, `rotation`, `labels`): its tracks as concentric rings round its
  whole axis, labels outside with leaders, the name and length in the hole, a linear detail's window drawn as a sector
  (quiet when it covers the whole ring) with grips on a 12 px rim.
- A linear section follows a ring's window and a ring follows a linear section, across a cyclic origin; a cyclic
  viewport wraps in every layout. The ring's window is navigated, zoomed about the pointer's angle, dragged, selected
  across the origin by the shorter arc, and focused on an item's arc.
- Sections take `block` children (`of`, `range`, `label`; `arrange: row | stack`, `copies: merge | each`) for a complex's
  chains and ligands, in tracks and on rings; blocks on several systems select and focus per block.
- A one-position block (a ligand, an ion) is a chip in every layout: a fixed slot its label's width, never a
  position's share, so it reads at the fit view and is never clipped; in a row of blocks where its block stands, in a
  stack the chips share one row, on a ring a fixed-angle arc, in text flowing like words; its letters band draws its
  label in a pill. Chips take at most half the width, their labels elided past it.
- A row of blocks pans and zooms through its own geometry: no jump as a chip or a gap crosses the edge, and a zoom keeps
  the point under the pointer; a block label that does not fit is elided.
