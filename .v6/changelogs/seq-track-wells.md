## Seq* (next): track wells

- A track's `padding` (px inside its well: a number top and bottom, or `(:top n :bottom n)`; the render config's
  `trackPadding` by default) insets its data, hit tests and emphasis; a track's label is centred on its well.
- Each track's well is a squircle its data is clipped to (`trackShape: 'square'` keeps boxes); a grid over blocks
  rounds only its plot's outer corners; Canvas2D and SVG draw the same curve.
- Every emphasis is a squircle, round where it really ends and square where cut or below 6 px: the hover and its
  projection, highlights in each style, the selection, the keyboard cursor, and a cell's, a chip's, a helix's and a
  grid cell's hover outline, concentric with what it surrounds.
