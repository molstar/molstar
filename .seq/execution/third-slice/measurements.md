# Stage 7: measurements of the third slice (M3)

Measured on 2026-11-20 on `packages-seq` after every stage 6 lane merged (`ab68a2ba0`, plus plan-only commits).
Two environments:
- **Node** (the heavy lane's settings);
- **Chrome** in the app's visible browser pane at dpr 2, with the canvas draw included but not the GPU's raster.

To rerun:
- `pnpm test:heavy packages/seq/`;
- `pnpm -C packages/seq/core test:perf`;
- `pnpm -C packages/seq/storybook test:perf`;
- `pnpm -C packages/seq/core size`;
- the browser: a script in the pane wraps `requestAnimationFrame`, which `DomHost` looks up at call time, so the
  viewer's ticks are timed. It sends one wheel event a frame (240 pan, 120 zoom, 120–160 row scroll) and observes
  `long-animation-frame` entries.

| Target (performance.md) | Measured (median / p95 / max) | |
| --- | --- | --- |
| Heavy lane: wall-clock bounds, scaling, oracles | 27 files, 109 tests, 24 s | pass |
| Draw and hover loops allocate nothing warm (core) | 135 cases within baseline + 4 B | pass |
| Draw and hover loops allocate nothing warm (Storybook) | **the viewer's warm emphasis frame (p53 on its gene: selection, highlight and hover, projected): 800 B a frame** against 0; the other 4 cases pass | **fail**: lane EA |
| Browser pan and zoom under 8 ms, never over 16: the 10,000 × 300 alignment (Globins › Rows › Scale) | pan 0.1 / 0.2 / 0.2; row scroll 2.1 / 2.9 / 3.1; zoom 1.2 / 2.5 / 6.2 (warm). **Two frames of about 100 ms in every 120-step zoom**: a frame callback outside the viewer's ticks, most likely the rows' image tiles rebuilt at a new level | ticks pass; **the 100 ms frames fail** |
| The same: a ring (Plasmid › Ring) | pan 0.2 / 0.3 / 0.5; scroll 0.3 / 0.7 / 4.2; zoom 1.0 / 1.7 / 2.1; no long frame | pass |
| Text reflow at a changing width (Features › Text, 120 widths) | 1.6 / 3.3 / 4.1; no long frame | pass |
| A story's transitions (Features › Story, six steps) | 647 frames 0.6 / 2.1 / 4.0; no long frame | pass |
| SVG export under 200 ms | four steps exported in one click, 3.3 ms synchronously | pass |
| PNG export (1644 × 224 device px) | core's `encodePng` 4–6 ms warm (11 ms cold), **1.47 MB** (stored, by design: synchronous, `png.ts`); the canvas's own `toBlob` 6 ms, 18 kB | time passes; size noted |
| Bundle: reading a spec under 160 kB | 175.4 kB min (56.4 kB gzip) | **fail** |
| Bundle: headless rendering under 300 kB | 712.0 kB (244.2 kB gzip) | **fail** |
| Bundle: the viewer, default, under 600 kB | 1142.4 kB (386.0 kB gzip) | **fail** |
| Bundle: viewer, minimal, under 600 kB | 853.1 kB (286.6 kB gzip) | **fail** |
| Quick lane under 30 s | 2463 tests (`pnpm test:unit packages/seq/`) | pass |

## The misses, for the targets gate

- **The emphasis frame's 800 B**: a regression since `0a92bdb64` added the case. Lane EA (fast) is finding and
  removing it; the budget is not raised.
- **The 100 ms zoom frames on 10,000 rows**: rare (two in 120 steps) but over 16 ms. This is a fix lane, not a target
  move. The likely remedy is building the new level's tiles across frames, keeping the old level drawn meanwhile
  (performance.md, Levels).
- **Sizes**: every row fails. `seq/core/src/render` alone is 370 kB in headless rendering and in the default viewer.
  Rendering grew with the third slice's layers, layouts, navigators and rows; the targets date from the first slice.
  The author's call: move them to the measured sizes with room, or reduce `render` first (per-layer and per-layout
  entries lazy, as MV made the viewer's registries injectable).
