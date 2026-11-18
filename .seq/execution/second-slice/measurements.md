# Stage 5: measurements of the second slice

Measured 2026-11-18 on `packages-seq` after lanes T and U (`3561f0a02`). Node 22 headless (the heavy lane's settings:
1440 CSS px wide, dpr 2, warm) and Chrome in the app's visible browser pane at dpr 2 (the canvas draw included, the
GPU's raster not). To rerun:
- `pnpm test:heavy packages/seq/` (43 passed; `npx vitest run --mode heavy <files> --disableConsoleIntercept` prints);
- `pnpm -C packages/seq/core test:perf` (allocation: 54 passed);
- `pnpm -C packages/seq/core size` (bundles);
- the browser: a script in the pane wrapping `requestAnimationFrame` to time each frame's callbacks, a wheel event a
  frame (240 pan frames, 120 zoom frames) on each story.

| Target (performance.md) | Measured (median / p95 / max) | |
| --- | --- | --- |
| An edit: patch, resolve, redraw under 16 ms; undo the same | edit 0.44 / 0.90 / 2.22 ms; undo 0.00 / 0.29 / 0.92 ms (p53 Tracks, 50 each, Node) | pass |
| Selection drag and cursor step under 4 ms | drag-select 0.04–0.15 ms; additive drag 0.09; cursor 0.02–0.03; follower pan 0.17 | pass |
| Decode MessagePack: columns 200 MB/s, node-heavy 20 MB/s | columns without a copy (aligned views: 24 MB in about 2.5 ms); 20,000 small nodes 60 MB/s | pass |
| Draw and hover loops allocate nothing warm | 0 bytes, every scenario (the allocation tier) | pass |
| Layout bounded a frame | within the tier's bounds | pass |
| Hover under 2 ms | about 0.001 ms everywhere | pass |
| Pan and zoom frame under 8 ms, never over 16 (Node) | p53 views 0.02–0.20 ms; 3M letters 0.03; 3M cells 0.12; 3M bars 100k pan 2.18 / 2.40 / 3.44, zoom 1.61 / 2.18 / 2.38; 3M bars 200k zoom 2.68 / 3.64 / 3.97 | pass |
| **Browser**, pan and zoom frame under 8 ms (Chrome, dpr 2) | p53 Tracks at letters: pan 0.8 / 1.1 / 1.3, zoom 0.8 / 1.0 / 2.0; Protein on its gene: pan from the protein 1.1 / 1.7 / 2.1, from the gene 1.0 / 1.4 / 3.1, zoom 0.9 / 1.4 / 2.3; Gene by region: pan 1.0 / 1.5 / 1.6, zoom 1.3 / 2.3 / 5.5; the frame gap 16.7 ms median (no frame dropped) | pass |
| Validate and resolve 1,000 nodes under 20 ms | 8.07 ms (validate 1.40, resolve 0.23, views 6.45); from JSON 4.26 ms | pass |
| SVG export under 200 ms | p53 Tracks 0.66 / 0.79 / 1.25 ms | pass |
| First paint of a dataset view under 100 ms after the data is in memory | not isolated in the browser: p53 Tracks settles 212 ms after navigation, Storybook's own load included; headless, decode + validate + resolve + views is about 12 ms and a frame about 1 ms | pass by the parts, not measured whole |
| Bundle: reading a spec under 150 kB | **176.4 kB** min (55.4 kB gzip) | **fail** |
| Bundle: headless rendering under 300 kB | 293.7 kB (94.7 kB gzip) | pass, 6 kB left |
| Bundle: the viewer, default, under 600 kB | 504.8 kB (165.9 kB gzip) | pass |
| Quick lane under 30 s, a folder under 2 s | 796 tests, about 1.5 s | pass |

Not measured: the 3M fixtures in the browser (no Seq\* (next) story shows them; the first slice measured them headless
only too), bigWig and the 10,000 × 300 alignment (not in this slice).

## The miss: reading a spec, 176.4 kB against 150

- **Where it comes from**: the first slice's 133.4 kB grew by 43 kB. 26 kB of it is `@molstar/core/data/immutable`,
  which only the editing validator uses (`lang/check-changes.ts`, the ref index's persistent maps): the read path
  imports `@molstar/seqx-core/lang`, whose index re-exports that validator, and the packages declare no
  `"sideEffects": false`, so the bundler keeps the collections' top-level code. Measured with that export removed:
  **150.6 kB** (48.3 kB gzip). The rest (about 17 kB) is this slice's language: the security bounds, refs by key, the
  address and query forms.
- **Options**:
  1. move the editing validator to its own entry (`@molstar/seqx-core/lang/changes`, imported by the viewer's
     `SpecModel`), so reading a spec carries no editing: 150.6 kB; and move the target to 160 kB for the language this
     slice added;
  2. 1 without moving the target: find 0.6 kB (and the next slice's language has no room);
  3. keep the code as is and move the target to 180 kB;
  4. `"sideEffects": false` on the packages instead of a new entry: a repository-wide change for `@molstar/core`.
- **Recommendation**: 1. Reading a spec should not pay for editing; the entry split is small (one export path, the
  viewer's import, a hot `package.json` the orchestrator edits); and 160 kB keeps the read path guarded with room for
  one slice of language.

## For the author (the stage 5 gate)

Confirm or move the targets the slice touched. One needs a call: the read bundle (above). Headless rendering has 6 kB
left under its 300 kB: worth knowing before the third slice's layers and themes.

## Decided (the author, 2026-11-18: "Yes, as recommended")

Option 1: the change validators moved to `@molstar/seqx-core/lang/changes` (the viewer's `SpecModel` imports them
there), the read target is 160 kB. Remeasured: reading a spec 150.6 kB (48.3 kB gzip), pass; headless rendering
267.9 kB (87.8 kB gzip), which also dropped editing, 32 kB under its target; the viewer 504.9 kB. The other targets
stand as measured.
