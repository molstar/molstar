# M4 browser measurements (fourth slice, 2026-11-21)

Chrome in the Browser pane, emulated 1440 x 900, **dpr 2** (confirmed in page), page visible, load average 1.3 to 1.4 (`uptime`)
before the runs. Console errors: none.

**Method.** The pane is not reliably in view (rAF stalled when it was hidden; vsync-bound at 60 Hz when shown, which hides cost), so
`requestAnimationFrame` was replaced in the page by a queue the harness flushes right after each synthetic input. A frame is
the input dispatch (pointermove/mousemove, or a wheel: shift-less dx for pan, ctrl+dy for zoom, dy for row scroll) plus the
draw callbacks it queued, timed with `performance.now()`; 20 warm frames, then 300 measured. This is main-thread CPU
(handler, layout, canvas 2D commands): not raster, GPU or compositor time. Fake structure lit chains checked changing
(38 distinct states across the sweep).

| Case | Story | dpr | Frames | Median ms | p95 ms | Max ms | Target | Verdict |
|---|---|---|---|---|---|---|---|---|
| 2OCJ hover, P04637 row, lights 4 chains in the fake | Features/Structure link: Hover both ways | 2 | 300 | 0.7 | 1.2 | 2.3 | hover under 2 ms, frame under 8 | pass |
| same, hover cycling three rows | same | 2 | 300 | 0.9 | 2.6 | 7.5 | same | median pass; p95 over 2 ms, max under 8 |
| 200 copies merged, hover | Examples/Complex/Views: Copies | 2 | 300 | 0.0 | 0.1 | 0.4 | same | pass |
| 200 copies merged, pan | same | 2 | 300 | 0.1 | 0.2 | 0.4 | pan under 8, never over 16 | pass |
| 200 copies merged, zoom | same | 2 | 300 | 0.1 | 0.2 | 0.2 | same | pass |
| 200 copies split (200 blocks), hover | same, after Split | 2 | 300 | 0.2 | 0.3 | 0.6 | same | pass |
| split, pan / zoom / row scroll | same | 2 | 300 each | 0.1 / 0.1 / 0.1 | 0.2 / 1.2 / 0.2 | 0.3 / 5.7 / 0.7 | same | pass |
| split patch, one click to drawn | same | 2 | 1 | 33.6 | | | edit under 16 ms | over (one sample, includes a 500 ms wait's timer slack; re-measure) |
| gene-by-region pan / zoom / hover | Examples/p53/Gene by region | 2 | 300 each | 0.2 / 0.1 / 0.0 | 0.5 / 0.5 / 0.1 | 3.1 / 1.1 / 0.2 | same | pass; NOT a bigWig |
| signal by levels pan / zoom (proxy) | Examples/p53/Chromosome | 2 | 300 each | 0.4 / 0.1 | 0.5 / 0.4 | 0.7 / 4.7 | same | pass; NOT a bigWig |

**Gaps.** (1) No story reads a bigWig: `runtimes/bigwig.ts` is tested in Node only, and the p53 stories have no conservation
signal over the gene region; Gene by region draws transcripts, Chromosome a generated GC-like signal by levels, so the bigWig
pan is unmeasured in the browser (needs a story). (2) 200 copies exist only in Examples/Complex/Views: Copies (CDK2 x 200),
not 2OCJ. (3) Raster and compositor time are not in the numbers; a real-vsync run needs the pane in view.
