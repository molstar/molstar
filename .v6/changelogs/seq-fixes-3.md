## Seq* (next): third-slice fixes

- The conservation logo's letters (and its bars below 2 device px) keep their class hue, moved in lightness to 3:1
  against the track, so they read as glyphs in both themes.
- A double-click gives back only the selection its own first click cleared: a double-click on a rail, after a drag, or
  after the selection changed otherwise no longer restores an earlier cleared selection.
- Hovering a residue in a section of blocks outlines that position only on blocks of the same chain (its other copies
  too), no longer on other chains' blocks that happen to hold the position.
- The logo's letters are saturated as glyphs: each keeps its hue, its oklab chroma raised to a floor
  (`--seq-logo-chroma`, 0.12 in light, 0.10 in dark), then its lightness moved to 3:1 against the track; greys stay
  grey.
- The navigators' cursor says the gesture: a window's grip shows a resize cursor along its axis (`ew-resize`, or
  `ns-resize` on the rows rail and a vertical navigator) through its drag, the body `grab` and `grabbing`, outside
  `pointer`.
- Overlapping highlight channels no longer blend into a third color: where channels overlap, the band splits into
  equal stripes in channel order (up to `emphasis.highlightStripes`, 4, the rest sharing the last), in tracks, text,
  rings and grids; a channel can be a `tint`, a `band` or an `outline` (set with `address.highlight`'s `style` or
  `highlight.set`'s value, kept when carried onto linked sections). The Legend part lists the active highlight
  channels with their colors, each removable.
