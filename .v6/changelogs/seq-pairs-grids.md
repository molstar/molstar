## Seq* (next): pairs and grids

- A `pairs` data kind (two ends each, kinds, strands), drawn as arcs standing on the sequence, ring chords, or dots in
  a grid; crossing pairs split above and below; colours from themes (by kind by default); pairs follow the layer's
  data filter.
- A `grid` data kind and a `grid` section layout with a second axis (`of_y`, `follow_y`, `aspect`; cells square where
  the room allows): dot plots and PAE as rects or cached image tiles, with min, mean or max reduction and triangles.
- A crosshair hover lights both positions and the arcs ending there; rectangle selection (click, Shift, Cmd/Ctrl,
  drag, Alt+drag subtracts); rects subtract exactly, also in queries (`(and :of [a (not :of b)])`); `(row :key k)`
  addresses an alignment's rows; a zoom or scroll over a grid moves the sections its two axes follow.
- Examples: the tRNA (arcs, split, chords, a dot plot of synthetic base-pair probabilities) and the Complex's
  synthetic PAE and contacts.
