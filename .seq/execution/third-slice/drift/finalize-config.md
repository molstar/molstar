# Drift: configuration the designs name and the code does not have (SR, 2026-11-19)

Found by SR (`reconcile.md`). Each a question for the author at the next gate, with a recommendation.

1. **`rows.pinInPlace` or `rows.pin.inPlace`** (grids-and-rows.md:549; the code flat). Recommend: the design takes the
   code's flat spelling (one item, no group of one).
2. **`reveal.opensClades` and `rows.clades.cutLine`** (grids-and-rows.md:553-554, layouts.md): the code always opens
   the clades holding a revealed row and has no cut-line drag. Recommend: `reveal.opensClades` fixed (no second use
   for "off"); the cut line moved to `TODO.md` as planned, the row dropped.
3. **`navigator.idle: quiet | hidden`** (navigators.md): idle is quiet only. Recommend: fixed quiet, `hidden` to
   `TODO.md` until a host asks (nothing configurable without a second use).
4. **`rows.reduce: theme | letter`** (grids-and-rows.md:549): the reducer is the majority only. Recommend: fixed, the
   tension kept open in the design.
5. **The `figure` preset for `text.*`** (layouts.md): presets are not built. Recommend: the row says planned, with
   presets.
6. **A Configuration table for `packages/seq/examples`** (examples.md promises one; no `_spec` there). Recommend:
   reword examples.md, "documented in the module that owns it", the design's table the one list.
7. **`hover.card` and the `hoverCard` behavior's port**: check ports.md says the port is the behavior's param (the
   author's "Accept").
8. **Sizes** (reading 211, headless 677, viewer 1064 kB against 160 / 300 / 600 kB): reset at finalize, after MV.
