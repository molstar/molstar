# State: the third slice

Work id `third-slice`. Plan: `plan.md`; loop: `loop.md`.

## Done

- Gate 0 decided (2026-11-18); the size targets are guides, not bounds ("don't overdo it on the limits ... but
  obviously try to minimize the amount of code").
- The designs drafted, reviewed for their decisions, self-checked as a set (the four steps), decided by the author and
  agreed into `.seq/design/` (axis, themes, layouts, grids-and-rows, navigators, the datasets in examples) with their
  Storybook Design pages.
- C0 merged (`6d4b96898`): review's six findings fixed, the look's disabled buttons fixed. Stage 0's gate passed:
  quick 822, heavy 44, typechecks, lint, format, sizes (150.6 / 267.9 / 507.8 kB, all under target); no medium tests
  in `packages/seq` yet.
- AX merged (`1e7e87870`): the review's seven findings fixed; the look's solid conflict ink fixed (the toy gene now
  encodes p53, R175H its one variant). Quick 892, typechecks, lint, format green on the merged tree; headless
  rendering 290.4 of 300 kB (a guide).
- SL merged (`2457f3ab9`): the review's eight findings fixed; the look's struck-through labels fixed. Quick 973 on the
  merged tree; headless rendering 313.0 kB against the 300 kB guide (reset at the close).
- TH merged (`1772436b9`) with the `./theme` export: the review's findings fixed but the core module split (a gap); the
  look's two empty stories fixed. Sizes against the guides: reading 170.1, headless 409.4, viewer 668.7 kB.
- Stage 1's gate passed: quick 2,745 (core color, model, mvs, seq), heavy 65, typechecks, lint, format, `pnpm build`,
  `check:workspace`; sizes over the guides (reset at the close; the core split, a gap, saves ~17 kB); no medium tests.
- TX's interface merged (`c6282f553`, with the section-cache fix `9cfc6801f`): its review's seven findings fixed. TX
  then built the text layout, its targets and the p53 text figure (branch `seq/tx-text`), leaving the viewer work:
  that is lane TV. An allocation regression found and bisected to TH (42.6 → 72.4 KB a pan).
- TA merged (`24db27a58`): the regression's cause was `Color.contrast` per label ink (boxed doubles); "p53 on its
  gene, whole protein" back to 49.9 KB against 64; the review found no correctness bug, its test gap fixed.
- TV merged (`be56ae0be`): the text layout and its viewer; the review's findings and the look's two fixed (the export
  cut at the canvas height; the domains band); quick 1,155, the allocation tier 84/84 on the merged tree.
- Found by TV: a range starting in an earlier tile is not drawn when the view lies wholly past the tile edge (the
  synthetic source answers only items starting in a region): a lost item, to fix test-first after RB.

## In progress

- RB (strong): `briefs/rb-blocks-rings.md`: its fix round's Part A done; Part B rebases onto `be56ae0be` and adopts
  TV's mechanisms.

## Next

Stage 2's gate after TV and RB. The early use gate is deferred.

## Decided

- At the design review: plain defaults (`pieces` `whole`, `direction` the writer's); a layer's default theme the
  language's; a collapsed intron zero width on every layout; hiding columns waits; the bookmark command and
  `highlight.clear` in C0; clade folds session state with `clades.keep`, a node click toggles, a reveal opens.
- At the set's self-check: the ruler rail the default linear navigator (a ruler drag pans; an alignment's ruler selects
  above its grab zone), with its grips and edge cases; `layer.filter` with a chip (TH); export choices as optional spec
  params; `selection.copyLetters` copies what is shown; p53's generator on core's `PCG`; figure-meaning settings are
  spec params.

## Findings by kind (gates, use gates, novice pass)

| Gate | Bug | Missing design | Look |
|---|---|---|---|
| stage 2, TV (review and look) | 5 (groups out of phase; end numbers untested; text allocation; a cyclic window drag; the copy message) + 1 found (an item from an earlier tile lost) | 0 | 2 (the export cut short; no domains band) |
| stage 2, RB (review and look) | 6 (an empty ring throws; >24 blocks; blocks on several systems select wrongly; label hits; overlapping label boxes; the start grip across the origin) | 0 | 2 (a whole-ring window loud; stacked rows empty) |
| stage 2, TA (review) | 0 (a test gap) | 0 | 0 |
| stage 2, TX's interface (review) | 5 (rings drew no emphasis; contract parts unwired; place not inverse of at; stale glyph values; a `.seq/` edit) | 0 | 0 |
| stage 1, TH (review and look) | 8 (log binned; clamp unapplied; the MolViewSpec round trip; mid not diverging; value-color overwriting; two legends off their scale; layer.filter called blocked) | 1 (a safe substitute that is not safe) | 3 (two stories empty; the view not on SL's layers; dark letters dim) |
| stage 1, SL (review and look) | 4 (heads on the wrong piece off-tile; hidden items untick'd, a design departure; secondary codes misread; an untested stem hit) | 0 | 1 (labels struck through) |
| stage 1, AX (review and look) | 6 (Shift+C on a lone gene; focus outside a range; pieces outside a range; a deopt hidden by scenario order; the history key's size; a reversed cyclic follower) | 0 | 1 (the protein row solid conflict ink) |
| stage 0 (C0's review and look) | 3 (a test not checking the drawing; a kept focus outliving a load; untested parts) | 0 | 1 (disabled buttons not dimmed) |

## Running

- 2026-11-18: overnight, without the author (loop, "Running while the author is away"): gaps decided as recommended
  into `gaps.md`; the use gates deferred.

## Blockers

None.

## How to verify

`loop.md`, Checks.
