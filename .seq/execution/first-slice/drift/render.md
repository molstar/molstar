# Drift: render (lane F)

Where lane F's code and `design/layout-and-draw-lists.md` disagree or the design is silent. Recommendations are the
lane's; each is to be folded into the design at the stage 4 gate unless the author decides otherwise.

1. **Reusing buffers and lanes**: `layout(..., prev?)` takes the layer's last layout (option: a state object in the
   context). Recommend `prev`.
2. **Emphasis reaching `draw`**: `draw(lay, pass, out, emphasis?)`; selection and cursor per track (layout), hover per
   layer (the hit decides granularity). Recommend recording it.
3. **Selection tint**: drawn translucent over the data, not under it (option: a `selection-under` pass). Recommend the
   overlay.
4. **Squircles**: `rect` keeps a radius and a round-ends mask (square cut ends need it); the curve is defined once
   (`squirclePath`, quartic superellipse) and every backend draws that path. Recommend the design say so.
5. **Density thresholds in device px**: letters need 14 CSS px at DPR 1 but 7 at DPR 2. Recommend letters by CSS px,
   the other steps by device px. **Author's call** (made configurable meanwhile: `RenderConfig`).
6. **Follower sections on a minus strand** show the leader's hull in genomic order (columns run opposite to the
   protein). Recommend genomic order now, a config item later.
7. **Length of a source-read system**: tier 3 reads the source coverage's length once attached, else
   `seq.render.axis.pending`. Recommend the design allow it.
8. **Inline wrapped rows (EC-9)** share `id`/`row`; readers dedupe by key, start and end. The data spec's "dedupe by id
   or row" should say so.
9. **`renderSvg`'s result**: `{ svg, settled, issues, lists }`; the SVG carries `data-*` attributes for settled and open
   domains.
10. **`query` layers** read the whole annotation they start at; filtering belongs to the sublanguage slice.
