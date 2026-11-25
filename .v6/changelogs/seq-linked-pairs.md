## Seq* (next): linked pairs

- The `linkedPairs` behavior, in the default viewer spec: hovering or selecting a range element (a feature, a CDR)
  lights the pairs with an end in it on every pairs layer (arcs, chords, dots), on its own `linked-pairs` highlight
  channel; remove it from the spec to opt out. Params `on` (`hover`, `selection`) and `layers`.
- The query language: a pairs annotation answers `items`, a pair one item with a row per end (`QueryItem.system`), so
  `(items :overlapping <location>)` gives the pairs with an end in the location, on either system.
- A highlight channel carries element keys (`EmphasisHighlight.elements`); the pairs layer draws the pairs a channel
  lists in the channel's colour, in the canvas and the SVG export. Highlight channels keep the order they were first set.
