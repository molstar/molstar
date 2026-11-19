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
  muted, pale class colours darkened at low chroma). DX merged after one fix round (the restore guarded by a count of selection changes; the drag reset removed). HB merged (looked at in light and dark on Complex › Blocks › Split; the ring of blocks draws no letters, so untested there). LS merged (looked at in light and dark: the logo reads colourful). FS-A1 reported green, under review. CR merged (checked in the browser: a grip `ew-resize`, the body `grab`). IZ (the zoom strip, decided 2026-11-19) running.
- Scheduled 2026-11-19 (the author, "ok as recommended"): HL (overlapping highlights), QS (Features/Query), the grid
  over a complex (in DC, its matrix read in FS-A2). Decided: the fourth slice's structure integration is a fake 3D
  viewer; the Mol\* integration waits for the new mol viewer (molecules.md decision 15, plan).
- FS-A1 merged after one fix round (prototype keys, values at their rows, row issues reported through a new optional
  `OpenedSource.issues`, the own system lazy, the seq runtime's sniff and bound, MessagePack loaded lazily). Checks in
  the merge: tsc -b 0; quick lane seq+core+io 3201 passed; mvs and sequences 4698 passed; check:workspace 0; lint 0;
  format 0. Sizes (split measure): reading a spec 201.0, headless 606.3, viewer default 973.8 kB. The author on sizes,
  use gates and spellings recorded (plan stage 7, USE GATE, drift/themes-spellings.md).
- IZ merged after one fix round (cyclic arcs, export settings, small heights, touch in the gap, own-axis sections, host
  heights, contexts other than bands, hysteresis); a conflict with CR in the viewer's `navigator.ts` resolved keeping
  both (the part's axis and the context check). Quick lane 1719 passed. FS-A2 in its fix round; HL in its fix round;
  RC running; PR starting.
- HL merged after its review round (legend redraws, the grid striped, outlines whole, pixel snapping, breaks once, API
  tidied); the grid rows' stripes still allocate about 96 B a frame (bounded, in `packages/seq/TODO.md`; was 536 B).
- FS-A2 merged after its review round (legend addresses on a base, one parser for base and filter, wraps clamped, the
  matrix parsed once, the variants fixture UniProt's own records at 392 KB as Biome formats it, 281 KB compact: the
  `biome.json` ignore for `packages/seq/examples/fixtures/**/*.json` left for the author's call). A spec conflict with
  HL resolved keeping both. Quick lane 1827 passed. QS can start.
- PR merged (`drag-session.ts`: window listeners for a drag's life, a buttonless move, blur, a hidden document).
  Quick lane 1847 passed. Running: RC, LB; QS next.
- The examples' JSON fixtures stored compact (the author, 2026-11-19: "Yes, add the ignore"; `01063c2cd`). RC, LB and QS
  in fix rounds after their reviews.
- QS merged after its review round (looked at: Residue 175 → `/data/p53 175`); the Files stories' head indents compact
  JSON (the fixtures are compact now). Quick lane 1866 passed. RC and LB in fix rounds.
- LB and RC merged after their review rounds (conflicts with LB and FS-A2 resolved keeping both). Quick 1924, medium 8,
  perf 120 passed. Identity is PID1-like (either row a residue; gap-gap and insert cells left out), the author may
  choose another. Sizes after RC: reading a spec 203.9, headless 626.8, viewer default 995.5 kB. DP and DT running; DG
  and DC start (strong slots).
- Decided (the author, 2026-11-19): pairs, ranges and maps on an alignment's columns wait for 0.5.0 (followups).
- DP merged after its review round (sites asserted under every rotation, IUPAC letters, honest summaries; its four core
  gaps are CY's). Running: DG, DC, DT's fix round, and OP, AY, SW (the author's exception); queued: NB, DA, CY.
- DT merged after its review round (pairs over a source's letters, the rows legend by alphabet, the tertiary column,
  generated data formatted); record.ts conflicts with DP resolved (both recorders, one dispatch table). Quick 1989,
  medium 9. CY gains the grid-from-a-source gap. Running: DG, DC, OP, AY, SW; queued: NB (next strong slot), DA, CY.
- SW merged (a pan about the pointer, a ring chip's label once, the minimal viewer measured at 959 kB: what keeps it
  large in plan stage 7). Quick 1990.
- AY merged after its review round (an SVG `role="figure"` with named groups and items, boxes that never take the
  pointer, ids hashed per export, announcements by priority and quiet in a drag, forced colors for any theme source;
  rows' and pairs' items untested: followups). OP merged after its review round (the core and the viewer spec conflict
  kept both; mono fonts resolved to Liberation before resvg, `image-rendering="optimizeSpeed"`, links outside the root
  refused, exit 5 internal, atomic writes, a pixel cap; the font files unverified against upstream: followups); its hot
  files applied (the resvg catalog entry, the inventory, the root references, the lockfile). Quick 2043 (with
  seqx-render), lint, format, workspace green. DG in review; DC running; NB waits for DG's slot.
- DG merged after its review round (attribution shared by the family and columns views with UniProt's changes stated,
  the recorder's rules tested offline, `#=GF SQ` and the seed's hash fixed, a summary row keeps its `+n`: core's
  rows layer elides the first label instead; the Query stories' alignment examples over four inline rows). Conflicts
  with DP and DT (`record.ts`, `draw-checks.ts`) resolved keeping both. Quick 2070, medium 9. Its core gaps (6–8) are
  CY's. DC in review; NB starts in DG's slot.
- Decided (the author, 2026-11-19, "Yes, as proposed"): a track's vertical `padding` inside its well and squircle
  wells, data clipped (ux-foundations §3, §5; the Design page). TW queued after CY. DC in its fix round (merging
  `packages-seq` first); NB running; queued strong: CY (DC's slot), then TW, DA (after NB).
- Decided (the author, 2026-11-19, "Yes, into TW"): a highlight takes the shape of what it marks; the hover outline
  around a chip made concentric (ux-foundations §5). TW's brief gains it.
- Amended (the author, 2026-11-19, "Yes, round all"): every emphasis is a squircle as the selection already is (the
  earlier "a rect over positions" dropped); TW's item 3.
- DC merged after its review round (emphasis drawn once a strip, `gridSummary` counting each cell once and keeping the
  triangle, a grid over stacked or ring blocks refused, the rect over copies by the hit's block, the over-blocks drag
  kept inside the select behavior; one mmCIF reader for the three recorders; ligand arcs by ligand). Public beyond
  `gridSummary`: `Hit.yBlock` (a split copy told from its twin), `GridLay.off` and the track layout's `filter.dy`
  (fields of existing layout types), accepted. Quick 2131, medium 9. CY starts in DC's slot.
- NB merged after its review round (a cap on a source's numbering with an issue past it, held labels kept through
  eviction, a pending numbering resolves as waiting and the address model wants its tiles, blank and mismatched
  numberings reported, labels matched by spelling then case; `shownNumbering` removed, the reader internal). No view
  shows it yet: DA's Antibody does. Quick 2166, medium 9. Designs (features, files, molecules) say the kind is built.
  TW starts in NB's slot; DA after.
- Decided (the author, 2026-11-19, "Yes, accept it"): CY's lengths through resolve (`sourceLengths`,
  `ResolveOptions.lengths`, `Resolver.resolve(spec, lengths?)`, `StaticModel.lengths`, `SpecModel.setLengths` replacing
  `ViewModel.setLengthSource`/`refreshLengths`, `overlayQueryModel`'s `RowsReader`), subject to its review.
- TW merged after its review round (a grid over blocks rounds only its plot's outer corners and takes its padding,
  outlines drawn after the well's clip, the secondary, bands, rows, cursor and grid cell hovers squircles, labels
  centred on the well). Quick 2206, medium 9, perf 120. CY in its fix round (will meet TW's layout hunks at merge);
  DA starts in TW's slot.
- CY merged after its review round (one `openModel` for the open-resolve-open protocol, lengths keyed by source so an
  edited uri never keeps the old length, `renderSvg` assigning its model only when the opening wins, row keys pending
  and requested, the label skip, the text guard by name; `setLengths`, `lengthsEqual`, `StaticModel.lengths`,
  `RowsReader` gone; `mapsEqual` added to `@molstar/core/util/map`). No conflicts with TW. Quick 2553 (DP's `it.fails`
  now passes), medium 9, perf 127, workspace green. DA running; the strong queue is empty after it.
- DA merged after its review round (the against view said to be the edit list's alignment, SAbDab's licence field and
  the entry that picks its segments recorded, the design's rules marked unchecked, `kabatFromImgt` refusing another
  framework occupancy, the record script running from built `lib/` for every dataset; core: a headless export opens a
  ruler's numbering source). examples.md corrected: IMGT from SAbDab, Kabat derived from IMGT, ANARCI later. Quick
  2594, medium 9, workspace green. Sizes: reading 211.2, headless 676.7, viewer default 1063.6, minimal 999.6 kB.
  Stage 5's gate passed; stage 6: the NOVICE PASS next.
- Stage 5's gate reopened: the quick lane's examples 2.04 s and storybook 2.45 s are over the 2 s guide (core 1.91,
  viewer 1.86, seq 6.37 s in all); QT2 (fast) fixing. The NOVICE PASS (fast) starts beside it.
- QT2 merged (its diff moves tests only: slow renders and sweeps to medium with quick representatives; no waste found):
  core 1.98, viewer 1.96, examples 1.08, storybook 1.96 s, seq 5.4 s (at load 9, the novice pass running); quick
  2568, medium 105. Core, viewer and storybook sit near the guide. **Stage 5's gate passed.** The novice pass runs.
- The NOVICE PASS reported (`novice-findings.md`): 21 findings, bug 5, missing design 8, look 8 (2 good). Placed: NP
  (fast) the bugs, the stories and the looks; the hover card and a keys hint asked of the author.
- Decided (the author, 2026-11-19, "Yes, lane HC now"): the hover card and the keys panel this slice (ux-foundations
  §9, viewer-models, features; the Design pages). HC (strong) starts; NP runs beside it.
- Decided (the author, 2026-11-19: "implement the hover card via a port too … so apps can easily provide their own
  styles"): the hover card's content and place computed headless, drawn through a `hoverCard` port (ports.md); the
  keys' list as data (`describeKeys`). HC's brief updated before it started.
- The author's finding (2026-11-19, Globins › Family): a row column's expanded label drawn under the next column; long
  labels: grow to the view's edge, ellipsis past it, the full text in the hover card (recommended). Given to NP.
  HC started (strong).
- NP merged after its review round (Page Down with pinned rows, `followSvg` on revisions only, the tip canvas 1×1 until
  a tip, the tip walk skipped without a rows layer, the tip clamped to its band, the tree target stopping at the
  labels). Not built: the keyboard cursor's rows (it is column-only; grids-and-rows' cursor designs `ySystem`,
  `yPosition`): the Design page and the design say so; asked of the author. Quick 2585, medium 105, perf 127.
- Decided (the author, 2026-11-19): HC's port stays the `hoverCard` behavior's param (answered "Accept"; the review's
  top-level field offered, unanswered); the default adapter through the default spec, an ownership flag. A ruler in
  its own numbering has no label (ux-foundations §9; the examples' "Position" dropped in HC's round). RK (the cursor's
  rows) after HC ("Yes, after HC").
- The author (2026-11-19): "let's pause the loop after this task finishes": HC's fix round is merged and recorded, then
  the loop pauses; RK (briefed) waits for the author's restart.
- HC merged after its review round (the port supplied by the default spec and owned, no stale card, the keys panel's
  keys and focus, the API trimmed, units and "–" in the card, rulers in their own numbering unlabelled). Quick 2606,
  medium 105, perf 127, workspace green; the quick lane at load 4.2 over the guide (core 2.33, viewer 2.53, storybook
  2.46 s): re-measure quiet. Looked at: the card on p53's domain. **The loop paused** (the author); `resume.md` holds
  the state for a new orchestrator.
- **The loop resumed** (2026-11-19, `/loop`). The quick lane re-measured near quiet (load 3.1–4.8, two runs): core
  2.40 / 2.69, viewer 2.47 / 2.22, storybook 2.63 / 2.08, examples 1.54 s: over the 2 s guide. QT3 (fast,
  `briefs/qt3-quick-lane.md`) starts; RK waits for it (its brief: nothing beside it in `packages/seq`).
- Decided (the author, 2026-11-19: "we should also add a "genome-level" example(s)"; "let's do chromosome scale this
  slice; and the whole genomes in appropriate slice"): **CS**, p53's chromosome (chromosome 17 by region, generated,
  TP53 in place; examples.md, its Design page), strong, beside RK after QT3 (`briefs/cs-chromosome.md`); the Genomes
  dataset is the fifth slice (version-1.md, plan). The data's choices (generated, not recorded) taken as recommended:
  in `gaps.md`.
- QT3 merged (`6ad0195c1`; test-only, reviewed by the orchestrator): the large figure (14 × 400) and the per-command
  sweep to medium with quick representatives (3 × 260; one command a kind), the navigator's warm-up 8 frames. Quick
  2253 in `packages/seq` (222 files), medium 107. At load 5.5: core 2.08, viewer 2.07, storybook 1.93 s: re-measure
  quiet at stage 7's gate. RK and CS (strong) start beside each other.
- Decided (the author, 2026-11-19): SF and SR (fast) start beside RK and CS ("Small fixes", "Spec reconcile prep"); the
  minimal viewer "Injectable, opt-in" (viewer-composition.md, Creating a viewer; lane MV, strong, after RK:
  `briefs/mv-minimal-viewer.md`); the open gaps "Yes, as recommended", found already decided and built by GD (gaps.md's
  Open heading was stale, corrected).
- SF merged (`2afa74936`; reviewed by the orchestrator, 48 lines): the size script's `--files` (inputs summed over the
  chunks loaded first) and `--only` (an unknown entry refused), its quick test; `@molstar/ui`'s disabled toolbar button
  at opacity 0.4 (1.7–1.9:1, computed; not looked at in a browser), Seq\*'s own rule removed; changelog
  `.v6/changelogs/ui-toolbar-disabled.md`. On the merged tree: quick 3190 (ui and seq), tsc, lint, format green.
- SR reported (read-only): `reconcile.md` (Configuration tables 17, features.md 41 rows, status.md 7, edge-cases 6
  with the third slice's outcomes missing, the themes drift's 3 edits) for finalize; its design questions in
  `drift/finalize-config.md` for the next gate.
- The author (2026-11-19, three screenshots of rondo.app): a genome ring like Rondo's for the fifth slice, "perhaps
  even as an extension to show this is "easy""; GR (fast, read-only research) started.
- GR reported (`research-genome-ring.md`): a Circos-like ring is close with what is built (blocks of several systems
  on a ring are tested, which GR missed); the Rondo look needs a bundled-pairs layer, a pairs level by bin pair, `ext:`
  layer kinds, a navigator over systems and a ring-zoom channel. Six questions for the author.
- Decided (the author, 2026-11-19: "Rondo-like, extension"): the fifth slice's genome ring is Rondo-like (bundled,
  count-weighted) and an extension package; `ext:` layer kinds before it (examples.md, version-1.md, the Design page,
  `research-genome-ring.md`). Its other five questions wait for the fifth slice's plan.
- RK's review (bugs 5, design and API 4, minimal code 6, test gaps): RK's fix round running; ↓ from the pinned band
  keeps the band's scroll (decided as recommended). CS's review (frames at the fit growing with each gene visited,
  about 460 ms after five; the ruler; density answers): CS's fix round running in examples and storybook.
- Decided (the author, 2026-11-19): a ruler's `units: none | bases` (axis.md); template instances hold their node only
  while in view (dynamic-nodes.md); both Design pages. CF (strong) briefed for them and inline data's tiles; it starts
  when RK or CS merges. EO (fast, read-only) started.
- EO reported (`ec-outcomes-draft.md`): of 56 rows, kept 45, changed 2 (EC-92, EC-138), not met 9; three
  contradictions in `drift/ec-contradictions.md`. AT (fast) running; ET (fast) starts for the ECs built but untested.
- AT merged (`f88979d0c`, test-only): the exported SVG names no pairs arcs (it.fails; fix added to CF's brief); a
  folded clade's summary row unnamed (drift/ec-contradictions.md 4). ET (fast) started.
- **Lanes at once** (the author, 2026-11-19: "let's keep 1heavy+2 light or heavy active after current stuff
  finishes"): three lanes, one strong and two of either kind, in place of loop.md's two strong and two fast.
- RK merged (`e76bbf40f`) after its fix round; its API accepted by the author ("Accept all"); looked at (Rows › Folded
  light, Family dark). ET merged (`50956fce0`): six ECs tested; EC-135 not wired (to CF). CS merged (`ce6566143`)
  after its fix round (frames at the fit 0.40 ms after browsing, from about 460 ms); looked at (fit and letters, light
  and dark). On the merged tree: quick 2308 (+1 expected fail), medium 108, tsc, storybook typecheck, lint, format
  green. grids-and-rows.md and its Design page say the cursor's rows are built; changelogs `seq-cursor-rows.md`,
  `seq-chromosome.md`.
- CS's PCG finding: core's `PCG` (`@molstar/core/data/util/hash-functions.ts:588`) is a 32-bit-state variant with a
  64-bit output step, biased (CS: G 33%, C 20% over 1 Mb); fixing it changes Mol*'s distinct colors and every seeded
  Seq* fixture: the author's call.
- CF, MV (strong) and SC (fast) started from `6883ee3c7`. Decided (the author, 2026-11-19: "Fix PCG in core"): a real
  PCG32 in `@molstar/core`, Mol*'s distinct colors and every seeded Seq* fixture regenerated once (base-library.md);
  lane PG (strong) after CF.
- SC merged (reviewed by the orchestrator; an inline numeric strand added): the hover card and announcer say a range's
  strand. Quick 2310 (+1 expected fail) on its branch.
- The author (2026-11-19, Storybook): Examples › Antibody › Index › "Antibody Stories" broken: the index's map was
  exported from a stories file (and Complex's); fixed by the orchestrator (`24967738f`) with a quick guard over every
  stories file. "The numbering on the rail is broken": the ruler's labels in a numbering with insertions and gaps
  (52A and 82C at the fit, every insertion ticked, the step not held across 100A–C, IMGT's 60 and 110 missing, the
  rail's window on the tick row): sent to CF (it owns `ruler.ts` this round), with a recommended rule for the author.
- SD reported the fourth slice's structure syncing draft (`drafts/structure-syncing.md`), seven gaps for the author.
  MF (fast) starts: a Files story of one file over several sequences (the author's question; built, not shown).
- The author (2026-11-19, using the Complex's PAE): the plot resizes as one axis zooms, a second double-click
  re-centres, the wheel falls through to the page; the chain boundaries a gap plus an off-centre hairline. Decided
  ("Fixed box", "Zoom a step about it", "Never scroll the page"; the boundaries as the blocks' gap, recommended):
  grids-and-rows.md and its Design page; lane PZ (strong) after CF and MV. The use gate's count: bug 2 (the rail's
  labels, the index stories), missing design 2 (the plot's navigation, its boundaries).
- MF merged: Features › Files › Many Sequences (one TSV, a target column, two sequences). Decided (the author,
  2026-11-19: "Behavior, on by default"): `linkedPairs` lights the pairs with an end in a hovered or selected range,
  a behavior on the public API in the default spec (viewer-models.md, the Design page); lane LP (strong) after MV.
- HP (strong, a fourth-slice draft: track header parts) started. Decided (the author, 2026-11-19: "the full labels
  should only display when you actually hover on the species column"): a row column's cut value whole only under the
  pointer (grids-and-rows.md, the core `_spec`); lane RX (fast) at the next free slot. Use findings so far: bug 2,
  missing design 4, look 1.
- Decided (the author, 2026-11-19: the clade visuals "quite ugly", the line widths uneven; "Light wedge"): branches
  snapped to device pixels, the wedge a light outlined triangle to the deepest leaf (grids-and-rows.md, the Design
  page); folded into RX (now strong: the rows layer's look). Use findings: bug 2, missing design 4, look 3.
- The author (2026-11-19): "2 columns hidden" misleading (reads as alignment columns; does not say which): the chip
  names them, its card lists them, on their side of the title band (decided as recommended; grids-and-rows.md, the
  Design page); folded into RX. Use findings: bug 2, missing design 5, look 3.
- The author (2026-11-19): "also the positioning of the label is wrong": read as the hidden columns' chip hanging over the data's top-right corner; RX places it in its side's title band (to confirm with the author).
- The author (2026-11-19, Rows › Scale): the rows rail dragged to the whole band cannot be grabbed again: the code idles on 'whole' alone against navigators.md (idle when whole and zoom-in at its limit); the core _spec corrected; folded into RX. Use findings: bug 3, missing design 5, look 3.
- Confirmed (the author, 2026-11-19: "Yes, leave it empty"): the label column empty below rows.labelMinPx (a dark slab drawn there on Rows › Scale); folded into RX. Use findings: bug 4, missing design 5, look 3.
