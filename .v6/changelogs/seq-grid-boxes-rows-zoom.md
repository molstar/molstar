## Seq* (next): grid selection boxes, the rows' zoom limits

- A selection's boxes on a grid meet the cells' edges (clipped to each cell, the shared cell tinted once).
- A lazy alignment's rows stop zooming out at a device pixel a row while its source serves no row levels, so the plot never goes blank; `rows.zoom` (`minPx`, `maxPx`, `enabled`) sets the rows' zoom limits for a host, `enabled: false` leaving a pan-only rail.
