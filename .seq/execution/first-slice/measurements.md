# Stage 7: measurements of the first slice

Measured 2026-11-18 on the merged tree (`seq/measure`, merged at `ec6808f3b`'s parent). Node 22, headless, 1440 CSS px
wide, dpr 2, warm (1,000 unmeasured steps first), 200 measured steps. To rerun:
- `pnpm test:heavy packages/seq/` (bounds; `--disableConsoleIntercept` prints the numbers);
- `pnpm -C packages/seq/core test:perf` (allocation, about 10 minutes);
- `pnpm -C packages/seq/core size` (bundles).

| Target (performance.md) | Measured (median / p95) | |
| --- | --- | --- |
| Pan frame under 8 ms | p53 Tracks 0.04–0.09 / 0.10 ms; p53 Protein on gene 0.07–0.22 / 0.23; 3M letters 0.05; 3M cells 0.21; **3M bars, 100k positions in view: 9.0 / 9.2** | fail at bars |
| Zoom frame under 8 ms | p53, 3M letters and cells 0.03–0.17 ms; 3M bars 100k 6.5 / 9.0; **3M bars 200k: 11.2 / 14.8** | fail at bars |
| Never over 16 ms | max 15.5 ms (3M bars, 200k zoom) | pass, barely |
| Hover under 2 ms | about 0.001 ms | pass |
| Frame time flat across scale fixtures | letters and cells flat; bars not | fail at bars |
| SVG export (p53 Tracks) under 200 ms | 0.74 / 1.03 ms warm, 0.8 ms cold (Node) | pass |
| Validate and resolve 1,000 nodes under 20 ms | 7.1 ms (validate 1.0, resolve 0.2, views 6); from JSON 3.7 ms | pass |
| Draw and hover loops allocate nothing warm | pan 47–170 KB a frame (p53), 96 KB (3M letters), 422 KB (3M cells), **44 MB (3M bars)**; hover 0.9–1.2 KB | fail |
| Bundle: reading a spec under 150 kB | 133.4 kB min (42.9 kB gzip) | pass |
| Bundle: headless rendering under 300 kB | 228.9 kB min (75.2 kB gzip) | pass |
| Bundle: the viewer, default, under 600 kB | 318.2 kB min (106.7 kB gzip) | pass |
| Quick lane under 30 s, a folder under 2 s | 471 tests, about 1.2 s; the largest folder (viewer) 0.61 s | pass |

Not measured: browser frame times (canvas draw included), first paint in the browser, bigWig, the 10,000 × 300
alignment (not in this slice).

## Where the misses come from

- **3M bars**: the synthetic source's levels are 1, 64 and 4096. At 38 positions a device pixel the values layer
  takes the finer level (1) and places 100,000 positions one by one (`place` in `render/values.ts` 3.2 ms,
  `toSection` 1.7, `axis.spans` 1.0). With a level that fits, the same frame is 0.66 ms. Cost follows the gap between a
  source's levels, not the pixels. Options:
  1. allow a level up to k× coarser than a pixel, k configurable (a values layer's resolution drops by at most k);
  2. denser levels in sources (the store's summaries use `factor` 16; the synthetic fixture's 64 is its own choice);
  3. a per-position fast path in `values` (no `toSection`/`spans` per position on an identity map).
  Recommendation: 3 (it is waste) and 1 with k = 2 as a config item.
- **Allocation**: the layout path allocates (the values `place` closure's `pieces` and `spans` arrays, ranges `layout`,
  `layoutView`, sorts and `Int32Array` copies); `letters.draw` allocates about 12 KB; the hover allocates the frame key
  string (a join of nine parts). The design says "the draw and hover loops are allocation-free"; performance.md's
  "frames allocate nothing once warm" reads stricter (layout included). Options: hold layout to the same bar (reusable
  scratch arrays, a numeric key), or count only draw and hover (the hover key and `letters.draw` still need fixing).
  Recommendation: draw and hover allocation-free now (small); the layout's allocation bounded per frame, not zero.

## For the author (the stage 7 gate)

Confirm or move the targets in performance.md. The two misses above each need a call: the level choice (k) and what
"allocation-free" covers.

## After lane J (2026-11-18, merged `47a27cca7`)

The author's decision ("as recommended"): the targets stand; a values direct path and `values.coarser` (default 2);
draw and hover allocate nothing, a layout a bounded amount. Remeasured on the merged tree, same machine and settings:

| Target | Before | After |
| --- | --- | --- |
| Pan frame under 8 ms (3M bars, 100k) | 9.0 / 9.2 ms | 0.78 / 0.90 ms, max 1.7 |
| Zoom frame under 8 ms (3M bars, 200k) | 11.2 / 14.8 ms, max 15.5 | 1.01 / 1.61 ms, max 1.98 |
| Frame time flat across scale fixtures | bars not | every scenario under 2.3 ms max; the heavy bounds now hold all of them |
| 3M bars at 100k, `coarser: 1` (the direct path alone) | — | pan 2.05 / 2.19 ms, zoom 1.49 / 1.98 ms |
| Draw, hover, selection and cursor, warm | 12 KB (letters draw), 0.9–7.1 KB (hover) | 0 bytes, every scenario (a +88 B control is caught) |
| Layout and draw a pan or zoom frame | 47 KB–44 MB | 25–40 KB (p53, 3M letters and cells), 0.73 MB (3M bars 100k: about 300 B a ranges item in view), each under a stated bound |

Checks: 32 heavy tests, 45 allocation cases, 552 quick tests. The review compared the old and new renderers over 520
seeded frames a spec: no difference in any draw list.

Open: `values.coarser`'s default (gaps 15: with the direct path, 1 costs about 2 ms; recommendation 1). Browser frame
times and first paint are not measured (the second slice, in a visible pane).
