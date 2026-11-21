# Fixes: the fifth slice's small findings

Small findings that need no design collect here for the slice's fix points, as the fourth slice's did (the author,
2026-11-21: "lets have a 'fix phase' for slice 5 too where we accumulate these smaller things?"). An FX lane takes the
open items at each fix point (plan, Fix points); a finding that needs a design becomes a drift note instead.

## Open

- **The navigator's position readout crowds its edge** (look; the author, 2026-11-21): on a rows navigator (Globins,
  the alignment's map), "globins/columns:112 · row 47" sits flush against the map's top-left corner and wraps onto a
  second line over the data. Give the readout the chrome's padding (the hover card's inset), a backdrop that keeps it
  legible over data, one line (truncate with an ellipsis before wrapping), and a space after the colon. Check every
  navigator form (bar, strip, rows map, ring) in light and dark.

- **Letters larger than their labels** (decided design; the author, 2026-11-21: "yes, do that"): ux-foundations
  section 2 now caps a letter at the chrome's label size, optically matched (cap heights equal; about 0.9× for the mono
  face), regular weight by default. Build `letters.maxPx` (a `RenderConfig` key, default the label size; a host can
  make it follow the cell) in every letters path (tracks, text, blocks, rows, rings, the logo's letters if any), canvas
  and SVG alike; its row in core's Configuration table; a test that a letter never exceeds the label size at the
  closest zoom. Look across every example with letters, light and dark, at the fit and the closest zoom.

- **A legend hover sticks after a fast pass** (bug; the author, 2026-11-21, on Protein Over Its Coding DNA): moving fast
  across the legend sometimes leaves a legend entry's highlight on the view ("Highlight: legend Clear" stays, exon-2's
  region tinted); only hovering that entry again clears it. Likely a `pointerleave` lost when the legend redraws under
  the pointer (`viewer/src/mount-legend.ts`'s per-item listeners). FX1, item 8.

## Fixed
