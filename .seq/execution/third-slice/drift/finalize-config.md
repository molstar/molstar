# Drift: configuration the designs name and the code does not have (SR, 2026-11-19)

Found by SR (`reconcile.md`). Items 1–7 resolved: 1, 5–7 at finalize (CL, 2026-11-20); 2–4 decided by the author,
2026-11-20 ("Fix all three": `reveal.opensClades`, `navigator.idle` quiet and `rows.reduce` fixed, their design rows
dropped, the alternatives in `packages/seq/TODO.md`). Open:

8. **Sizes** (reading 211, headless 677, viewer 1064 kB at stage 6 against 160 / 300 / 600 kB): reset after lane RR
   reduces `seq/core/src/render` (the author, 2026-11-20: "reduce render first, then move the size targets").
