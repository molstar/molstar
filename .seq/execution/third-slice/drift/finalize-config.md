# Drift: configuration the designs name and the code does not have (SR, 2026-11-19)

Found by SR (`reconcile.md`). Finalize (CL, 2026-11-20) resolved items 1 (the design takes `rows.pinInPlace`), 5 (the
`figure` preset marked planned in layouts.md), 6 (examples.md reworded: no `_spec` for `seqx-examples`), 7 (ports.md
says the hover card's port is the `hoverCard` behavior's param) and left 8 (sizes) to stage 7's M3. The designs now
mark the three items below "not built"; whether each stays planned or is dropped as fixed is the author's call.

2. **`reveal.opensClades` and `rows.clades.cutLine`** (grids-and-rows.md, Configuration; layouts.md): the code always
   opens the clades holding a revealed row and has no cut-line drag (`clades.collapseTo` is built). Recommend:
   `reveal.opensClades` fixed (no second use for "off"), its row dropped; the cut line stays planned in `TODO.md`.
3. **`navigator.idle: quiet | hidden`** (navigators.md): idle is quiet only. Recommend: fixed quiet, `hidden` in
   `TODO.md` until a host asks (nothing configurable without a second use).
4. **`rows.reduce: theme | letter`** (grids-and-rows.md): the reducer is the theme's classes only. Recommend: fixed,
   the tension kept open in the design (`letter` stays in `TODO.md`).
