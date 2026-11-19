# State: the third slice

## Paused (2026-11-19): resume here

The loop paused for a hand-over to a new orchestrator session, with no lane running. To resume: read `loop.md`, this
file, `gaps.md`, `drift/`, then:

1. **PG** and **NV**: briefs written (`briefs/pg-pairs-grids.md`, `briefs/nv-navigators.md`); set their base to
   `packages-seq` at the latest commit (RW merged at `98d1af9e2`) and start both (two strong lanes at most).
2. Then stage 4's GATE (build, `check:workspace`, heavy run alone, sizes, gaps), stage 5 onward per the plan.
3. The author's open questions: `gaps.md` 1–5. The use gates stay deferred to finalize.

RW's leftovers for the author (from its fix round): the rows labels now use the items' font (10 px) so Pfam labels fit
(keep, or widen the column instead); `rows.pinInPlace` kept against the design's spelling `rows.pin.inPlace`; the heavy
lane's theme change (36 ms alone) reached 57–80 ms under a load average of 4, so run the heavy lane alone. The look of
RW's fixes was checked in light only (Scale's tiles in canvas and SVG, labels whole, wedges before labels).

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
- TE merged (`78f20c5fd`): the item lost past a tile edge was the test sources' contract, not the product's (the reader
  keeps an item by its anchor); fixtures fixed, a tile-edge test, the contract in the core `_spec`. Rings to check
  after RB.
- RB merged (`4e9a5ded6`): the review's findings, Part B's (a wide window flipping; a test changed to fit the code,
  restored to the design: tracks wrap on a cyclic axis) and the look's three (a loud window, stacked rows empty, a blank
  view after zooming out past the fit) fixed; rechecked in the browser.
- Stage 2's gate on `4e9a5ded6`: quick 1,228, typechecks, lint, format, `pnpm build`, `check:workspace`, the
  allocation tier 88/88; heavy: one max bound (`frame-time`, 150k pieces) failed under load beside the build and passed
  3/3 alone (to G2: best of three). Sizes against the guides: 173.2 / 468.8 / 744.1 kB. The look's two text findings
  (gutter labels overlapping; the domains count strip a heavy rule) to G2; the canvas/SVG spacing matched.
- G2 merged (`e2aac7989`): the three gate findings fixed, looked at in light and dark. Stage 2 closed.

- CH merged (`a5252eabd`, changelog `662d04133`): the clade hover (red once in the real story: the viewer's hover
  record dropped the target; fixed with a test through the viewer), the ruler's one label a position, the line numbers
  on the ruler's baseline (gap 7, decided). Review: findings 2–7 fixed; 1 (clade names) to the author (gap 6). Looked
  at in light and dark (Family, Folded, Figure). Post-merge: unit 145 files / 1395 passed, tsc, typechecks, storybook
  typecheck, lint, format, the `.seq` checks: exit 0.

- PG merged (`e99c48af6`, changelog `seq-pairs-grids.md`): review findings fixed but A7 (to GD); the look's five
  fixed, a regression from the look's own "square" finding reverted (the grid's x is the view's x); a tile seam fixed.
  Post-merge: unit 150 files / 1465, typechecks, storybook, lint, format, perf 100, the `.seq` checks: exit 0.

- RF merged (`1e0782e9e`): the rows' numbers in view, the generated labels' ranges (the 247), the right gutter
  clipped at the band's foot. Post-merge: unit 151 / 1469, typechecks, storybook, lint, format: exit 0.

- NV merged (`96a649d1f`, changelog `seq-navigators.md`): four rounds (review's 13 bugs; stories added; the ideogram
  and plasmid strip blank at dpr 2; placement; the overview's width and rows window; the rows rail, decided; the
  wash). Post-merge: unit 161 / 1578, typechecks, storybook, lint, format, perf 102, the `.seq` checks: exit 0.

- GD merged (`b69771348`, changelog `seq-gap-decisions.md`): gaps 1-6 and the individual colour-list exports (-13.4 kB);
  review findings fixed; looked at in dark. Every ColorLists importer typechecks; MVS and the base sequences pass.

- SO merged (Features apart from Examples; a story-link test). FS's design agreed (the author, 2026-11-19: "Ok as recommended"): `.seq/design/files.md`, the Design page
  "Files and sources", sources, protocol, examples, base-library, themes and performance amended; FS split into FS-A1
  (strong), FS-A2 (fast), FS-B (conditional). Next: FS-A1. The heavy lane capped at 2 workers (the author, 2026-11-19: "Yes, cap at 2"; QT's measurement).

## In progress

- Stages 3 and 4's gate (2026-11-19): medium, `pnpm build`, `check:workspace` exit 0 (no known failure this time);
  heavy: 2 bounds failed in the full run (rows-time's theme change, pairs-grids-time's arcs hover), both pass alone
  twice at load 5 (the parallel run's pattern, a third bound after frame-time: QT measures a serial heavy run); quick
  lane: seq/core 2.05 s and seq/storybook 2.25 s over the 2 s guide, seq 5.3 s in all: lane QT (fast) fixing. Sizes:
  198.5 / 603.5 / 935.9 kB. QT merged (`61f351a30`): core 1.56, viewer 1.28,
  examples 1.56, storybook 2.00 s; 8 tests to medium. **The gate passed** (2026-11-19). Open to the author: the heavy
  lane's pool (QT: `maxWorkers: 2` in heavy mode, 1/3 → 5/6 passing, same wall time; `vitest.config.ts`, repo-wide).

- MD drafted (`drafts/molecules.md`): a `molecule` kind, `sequence.components`, `complex`/`unit`, a `group` map, 0.5.0
  for the fourth slice (lanes ML, MC, MR, MX); 9 gaps. Reviewed (8 must-changes, all taken) and revised; agreed by the author 2026-11-19 ("Other recommendations from
  the design ok"; gap 10 off; gap 11 superseded by the seq data ports design, fourth slice): `.seq/design/molecules.md`
  with its Design page; maps.md, base-library.md, examples.md updated.
- PG and NV reported (branches `seq/pg-pairs-grids`, `seq/nv-navigators`, base `51787bb2b`); checks rerun by the
  orchestrator: PG green (unit 1433); NV green but heavy, whose two failures (rows-time's theme change, themes' 1.1×
  ratio) also fail on `packages-seq` under load and pass alone on NV: flaky bounds, to stage 4's gate. Reviews running
  (independent, opus) done: NV 12 bugs (double history step, wrapped grips, the miniature's body, rows history,
  the view's images and colours, the idle key, ...), the ring's second window geometry, config beyond the table; PG
  9 (y never moves on the PAE, a self grid's rect painted as bands, pairs ignoring themes and filters, antiparallel
  only, triangles' diagonal, the hit path's boxing), the data not marked synthetic. Both in one fix round each
  (2026-11-19), with the look's findings below; NV's worktree was briefly switched by PG's reviewer and restored.
  Design calls written: navigators (slim rail under the tracks; what replaces the rail), grids-and-rows (height
  floor, two-rect subtraction, `aspect: none`, the `type` column). The look's findings:
  - PG: tRNA arcs spring from the band's foot, away from the letters, feet clipped; the RNA legend says T (not U) and
    shows a gap swatch; dot plots not square (`aspect` 1); the generated PAE mostly low with a high stripe (the design:
    low across the interface, high elsewhere); contacts dots barely visible in dark. The crosshair works.
  - NV: no Storybook stories for its five views (the brief asked); the idle rail at the fit view a heavy rule (the
    design: quiet); the SVG export draws a line where the rail is (NV: omitted by default) — check against the base.
  - Sizes: headless 573 (PG) / 566 (NV) kB and the viewer 879 / 890 kB against the 300 / 600 kB guides.

## Next

Decided 2026-11-19: one-position blocks are chips in every layout (layouts.md), contacts optional (the host's);
MD (the molecules design) scheduled in stage 5, built in the fourth slice. Scheduled in `plan.md` (the author, 2026-11-19: "schedule these appropriately in the third slice plan"): CH and GD in
stage 4; SO first in stage 5 (Features apart from Examples); DP's features as arrows by strand. After PG and NV: lane GD (gaps 1–5 decided as recommended, 2026-11-19: the core palette split, `classes-safe`,
the dark cell edge, `view.motion`, the pyramid's tolerance), then stage 4's GATE, then stage 5. The early use gate is deferred.

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
| stages 3–4's gate | 1 (heavy bounds flaking in the parallel run) | 0 | 0 (quick-lane time: QT) |
| stage 4, NV (review, four looks) | 12 + 6 (blank pictures at dpr 2; parts over the view; the overview's window and width; the wash; the stories missing) | 1 (the rows rail) | 2 (the idle rail loud; the wash) |
| stage 3, PG (review, looks) | 9 + 1 (a regression from the look's own finding) | 0 | 5 (arcs off the letters; T in RNA; non-square; the PAE; dark contacts) |
| stage 4, GD (review) | 5 | 0 | 0 |
| stage 4, RF (author) | 3 (numbers not in view; a label's range; the last row) | 0 | 0 |
| stage 4, CH (look, review) | 5 (the hover target dropped by the viewer; a band over the pinned band; a double tint with `pinInPlace` false; a unary chain's stem; a ruler blanking another numbering's label) + 1 to the author (clade names) | 1 (the line numbers' baseline) | 0 |
| stage 2, TV (review and look) | 5 (groups out of phase; end numbers untested; text allocation; a cyclic window drag; the copy message) + 1 found (an item from an earlier tile lost) | 0 | 2 (the export cut short; no domains band) |
| stage 2, RB (review and look) | 9 (an empty ring throws; >24 blocks; blocks on several systems select wrongly; label hits; overlapping label boxes; the start grip across the origin; a wide window flipping; tracks not wrapping, a test changed to fit; a blank view past the fit) | 0 | 2 (a whole-ring window loud; stacked rows empty) |
| stage 2's gate (look, heavy) | 1 (a max bound flaking under load) | 0 | 2 (gutter labels overlap; a heavy count strip) |
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

- SY ↔ RW bands seam (decided by the orchestrator, within stories.md and grids-and-rows.md): the band record's `bandTop`
  is in display rows (SY's language unit), the layout turns it into px; history keeps the top row by key plus `rowScale`,
  restored clamped (RW's px drift rejected). The ViewModel's bands (RW, with folds) are the source; a shown story frame's
  `bandTop`/`rowScale` merge over them, keeping folds. Merge order: SY, then RW rebased onto it.
- RW review: 15 findings (4 confirmed bugs, the bands units, config beyond the table, untested logic); look: SVG export
  drops the image tiles, labels cut to 3 characters, wedges over summary labels. All in RW's fix round.
- SY merged (6b36ad68d; the ./story export and changelog 2f0e4e417). Post-merge: unit 132 files / 1282 passed, tsc,
  storybook typecheck, lint, format all exit 0. RW rebasing onto 2f0e4e417 in its fix round.
- RW merged (`98d1af9e2`, changelog `08c6880de`): the review's 15 findings and the look's three fixed, the bands seam
  with SY as decided. Post-merge: unit 143 files / 1379 passed, tsc, storybook typecheck, lint, format exit 0. Sizes
  185.2 / 532.5 / 836.5 kB against the 160 / 300 / 600 kB guides.
- LG merged (the logo's letters at 3:1 against the track; looked at in light and dark: dark reads well, light reads but
  muted, pale class colours darkened at low chroma). DX merged after one fix round (the restore guarded by a count of selection changes; the drag reset removed). HB and FS-A1 running.
