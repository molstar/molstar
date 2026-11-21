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

## Fixed
