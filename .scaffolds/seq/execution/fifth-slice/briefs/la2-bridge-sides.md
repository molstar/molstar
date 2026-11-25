# Lane LA2: a bridge reads each side over assemblies; LA's rest

**Base**: `packages-seq` at the orchestrator's head after LA merges (`git reset --hard packages-seq` first). Branch
`seq/la2`. Model: strong (the access plan, the store's request, ownership). Alone in resolve, the store and
`render/pairs-ribbons-*`.

**Why**: LA's strong review rated it a blocker: over assemblies or blocks a bridge reads once, on side `a`
(`render/pairs-ribbons-read.ts`, `perSide`). In the peach/grape stories, zoom peach in and leave grape at the fit:
pairs that end in grape's window but not in peach's are never read; a grape zoomed in only ever gets peach's coarse
bins. across-sections.md Default 3 ("each side reads its own view") is built only for one system a side.

**Task**:

1. **A side in a request**: carry the side through the access plan, `LayerSource`, `wantOf` and the store's request,
   so a sided source (the `seq` runtime's sided pairs, a BEDPE by either end) reads by `b`'s positions too.
   `AssemblyRegion.unitOf` is keyed by name today, so `chr1` on both sides collides: key it by (side, name). Over
   assemblies each side reads per block of its own system (GA2's per-sequence plan).
2. **Ownership over multi-system sides**: the drawn-once rule (the finer side owns; `R` its window widened to the far
   level's grid) with `R` through `dataRegions` (normalized, merged, projected through `src.back`; bounded on a
   collapsed axis, EC-33), `owns()` not linear per pair, and the near side's re-read never coarser than the far level.
3. **LA's left items**: ribbons inside the SVG bridge group (`data-bridge`, `data-a`, `data-b`; the export skips the
   bridge in the passes); tests for the hatch while a side reads items and for ownership when a side is refused.
4. **Stories**: the peach/grape stories gain "one side zoomed in over assemblies" (peach zoomed in, grape at the fit;
   then the reverse), each with its headless twin asserting every pair in either window is drawn exactly once.

**Tests**: the drawn-once property over assemblies with both sides at every pair of levels (items and bins), windows
overlapping or not, a refused side, a pending side; a region read per side equals the whole read clipped to that
side's window; EC-209 and EC-33 in titles where they apply; both perf tiers (a frame 0 B warm, the hover's card
2,168 B a step today: do not raise it).

**Out of scope** (LA3, after): grids over two sections (`follow_y` a second section), categorical grid axes, the
liftover (EC-212).

**Where**: a worktree on `seq/la2`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files never
staged: list them. Anything in `@molstar/core` or `@molstar/io` stops the lane. A file past 800 lines is a question.
Sizes (headless 930 — at about 929.6 —, default 1,480, minimal 1,060 kB, the bridge chunk 24 kB — at 23.5): a bound
passed stops the lane with a report; do not edit bounds.

**Checks**: the quick lane, `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`, each
`typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf tiers to completion, size,
`pnpm test:medium packages/seq/`, the heavy `bridge-time` file when load is under 6, the `.scaffolds/seq/_checks`. A fix's test
fails on its parent.

**Report**: under 350 words.
