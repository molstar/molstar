## Seq* (next): third-slice fixes

- The conservation logo's letters (and its bars below 2 device px) keep their class hue, moved in lightness to 3:1
  against the track, so they read as glyphs in both themes.
- A double-click gives back only the selection its own first click cleared: a double-click on a rail, after a drag, or
  after the selection changed otherwise no longer restores an earlier cleared selection.
- Hovering a residue in a section of blocks outlines that position only on blocks of the same chain (its other copies
  too), no longer on other chains' blocks that happen to hold the position.
