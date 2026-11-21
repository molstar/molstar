# State: the fourth slice

## Ready to start

The third slice is closed (2026-11-20; its record in `.seq/records/third-slice/`): RR merged, the size bounds set
(195 / 780 / 1250 / 880 kB), its fold-back applied, its worktrees cleaned (CE's kept, parked). Start with
`/loop run .seq/execution/fourth-slice/loop.md`.

**Base**: `packages-seq`. Every brief's base is set to its head when the lane starts.

## Decided

- Gate 0 (the author, 2026-11-19): "Accept all" (plan, Gate 0).
- Stage 0's designs agreed and promoted: structure syncing (a library, not a port), track header parts, seq data ports,
  extension openings (DR, SP, XD; plan, Stage 0).
- Lanes at once: two (the author, 2026-11-20: "also we can run 2 lanes in parallel again now"; one for a while
  before it); a lane's review and fix round count as that lane. Before it: two or three, each strong or fast by its work (the author, 2026-11-20: "2-3 lanes, all either heavy
  or light depending on context"); limited by files: three tracks (molecules, formats and data, views; plan, Order and
  parallelism).
- DP kept as a strong lane (nothing else builds the seq data ports); the order of HP1 and ST1 is the orchestrator's
  (the author, 2026-11-20: "order doesn't really matter, decisions as recommended").
- `@molstar/core`'s standard components (ML1): "Accept" (2026-11-20). BW's range reads: the surface in
  `drafts/bw-range-reads.md` "Accept", 64 KiB pages, "fail on change" (2026-11-20). The first fetch of a service's
  fixtures is still asked when SV reaches it.
- A unit's positions count from 0, as every system's ("From 0, as built", 2026-11-20); molecules.md reworded.
- Small findings collect in `fixes.md` for the fix points; the author's findings from checking the third slice go
  there as they come (the author, 2026-11-20).
- A label that does not fit its box: an ellipsis inside before it is dropped, at least 4 characters by default,
  configurable (`labels.minChars`; the author, 2026-11-20; `fixes.md`).
- Module size (§module-size, `.seq/_checks/file-size.ts`): lane SZ splits the four files over 1500 lines after FX1,
  while no other lane is in them (the author, 2026-11-20: "Yes, add all three").
- MX's drifts (the author, 2026-11-20: "ok for the drifts"): a section's `via` counts every path of the section, its
  projections included (maps.md; MX fixes it with the projection bug); a component query on a source-backed sequence
  stays pending until the letters arrive (as built; TODO.md).
- AC's and ST1's drifts (the author, 2026-11-20: "Yes, as recommended"): an alignment's columns named
  `<alignment>/columns`; a map onto a unit a source describes resolves as an overlay edge, pending until the
  description arrives. Both to FX3 (`fixes.md`).
- Stage 1's gate passed (the author, 2026-11-20: "Pass, FX3 fixes times"): checks green on the merged tree, medium,
  heavy on a quiet machine (122, QE's `:index` bound in it), the tiers, build, workspace, sizes 221.3 / 824.8 /
  1289.9 / 922.6 kB; the quick-lane folders over 2 s (core, viewer, storybook) are FX3's first item. Submenus from a
  command's select param ("Commands own choices"); T-DM1's CCD ids and the pairs layer's numeric columns ("Both to
  FX3").
- **Public API changes allowed without stopping** (the author, 2026-11-20: "you can change public api of packages/seq,
  core/repr, packages/ui if needed. this has no other consumers beyond the work here yet"): a lane may change the
  public surface of `packages/seq/*`, `@molstar/core`'s `repr`, and `@molstar/ui` where its work needs it, each change
  in its changelog entry and listed in `gaps.md` for the author's review; `pnpm build` and `check:workspace` at its
  merge when core or ui changed. The rest of `@molstar/core` and `@molstar/io` still stop the lane.
- Size targets are bounds once reset after RR (the author, 2026-11-20: "reduce render first, then move the size
  targets"). Moved for the query language: reading a spec 230, headless 820, minimal 920 kB (the author,
  2026-11-20), the default viewer 1,290 kB ("Move to 1290 kB", 2026-11-20); for molecules drawn, headless 830 and
  minimal 930 kB ("Move to 830 / 930 kB", 2026-11-20); for the structure link and the rest of the slice, the default
  viewer 1,330 kB ("Move to 1330 kB", 2026-11-20).

- **The overnight gaps** (the author, 2026-11-21: "for gaps go as recommended. download what you need for fixtures
  etc"): every `gaps.md` recommendation taken. The minimal viewer's bound 940 kB; a hover through the structure link at
  most 2 kB a lit unit (performance.md); `header.slots`' rationale amended (one slot in the measured column); HD's drift
  to the next slice's view-model work (TODO), the issue until then; the viewer folder's budget unchanged, fewer jsdom
  files in a fix lane; F1 and F22 (the hover card on the structure stories' viewers and compact legends) to FX5, the
  other novice questions to the fifth slice's gate; F5 kept. Fetches allowed: SV's first service fetch, T-DM1's CCD
  ids, real PAE data, each with its licence recorded. The fifth slice's drafts still wait for the author.
- **The author's look, 2026-11-21** (to FX5): the extension's lollipops show the card but no hover emphasis; the PAE in
  AlphaFold DB's colours, on real data where a permissive one exists.

- **T-DM1's placeholders kept** (the author, 2026-11-21: "keep UNL placeholders"): no CCD component exists for DM1
  or the MCC linker (FX5's lookup); `UNL` stays, the finding and its source in `adc.ts`.
- **The use gates passed** (the author, 2026-11-21: "i have looked at things here and there in the storybook, so for
  now user gates passed -- we will need to do followups anyway as small things appear"): stages 2 and 5's; findings
  from here on are follow-ups (`fixes.md`).
- **The loop runs to the close** (the author, 2026-11-21: "you can run slice 4 loop all the way to the close. I will
  review the foldback.md file before we apply it"): FX5, SV and PC merged, stage 3's gate, M4, FINALIZE; it stops once
  `foldback.md` is written, for the author's review, then applies it, writes the record and closes.

## In progress

- 2026-11-20: stage 1's gate passed. Next: FX3 and ST2, two lanes.
- 2026-11-21: stage 2's gate passed (overnight, as recommended): ST1–ST3, HD, HP1–HP3 merged and reviewed; quick,
  medium, heavy (123, load 3.7) and both perf tiers green; sizes 222.4 / 826.5 / 1289.3 / 925.3 kB under their bounds;
  gaps in `gaps.md`; matrix cells checked by each lane's review. Its use gate deferred to the author. The fifth slice's
  stage 0 done (six drafts). Next: FX4 (the fix point after stage 2's gate) and EX, two lanes; SV waits for the author
  (its first service fetch is asked).
- 2026-11-21: FX4, EX, LE, CN merged; stage 4's gate passed (overnight, as recommended): quick, medium, heavy (123 at
  load 2.6; the chromosome view's pan-frame ratio missed its 1.9× once at load 4.9, 1.97×, and passed alone and in the
  full run after), both perf tiers; sizes 223.1 / 806.0 / 1296.4 / 929.5 kB. Next: the novice pass and M4's
  measurements without the browser's fetches; SV and stage 3's gate wait for the author; stop before FINALIZE.
- 2026-11-21: the novice pass (29 findings) and NP merged (`34a944a17`). The overnight work is done: what is left
  waits for the author — SV's first service fetch and stage 3's gate, M4's targets (the minimal viewer at 929.9 of
  930 kB), the use gates, `gaps.md`, the fifth slice's drafts; then FINALIZE.
- 2026-11-21 (morning): the author took the gaps as recommended and allowed the fetches. Next: FX5 (the author's look,
  F1/F22, T-DM1) and SV, two lanes; then stage 3's gate and M4; the fifth slice's drafts wait for the author.

## Next

1. **Overnight (the author, 2026-11-20: "you can continue with slice 4 too, just don't close it until i review it in
   the morning")**: the loop runs "Running while the author is away" (loop.md): decisions a gate would ask are made as
   recommended and listed in `gaps.md`; a public API change past a design in `packages/seq`, core's `repr` or `@molstar/ui` goes ahead (Decided, the
   author), listed in `gaps.md`; elsewhere in `@molstar/core` or `@molstar/io` the lane stops; the use
   gates are deferred, not passed. Two lanes at once in all, shared with the fifth slice's stage 0
   (`../fifth-slice/state.md`, "Overnight"): one fifth-slice draft (or its review) and one fourth-slice lane at a time,
   while both have work.
2. The fourth slice's order: ST2 and FX3 merged; then ST3 → HD (after FX3's map onto a described unit) → HP3; SV (after
   ST3); stage 2's and stage 3's gates; stage 4 (EX → CN, LE); the novice pass; M4. **Stop before FINALIZE**: no
   fold-back, record or close until the author has reviewed (the morning).
3. In the morning, with the author: the fifth slice's drafts agreed; the fourth slice's review, then its finalize.

## Findings by lane

One row a merged lane (loop, step 5, record): review findings and look findings by kind.

| Lane | Merged | Bug | Missing design | Look | Notes |
| --- | --- | --- | --- | --- | --- |
| AR | 2026-11-20 | 3 (review: duplicate names past the bounds; a written archive failing the ratio; case-folded names) | 0 | 0 | six low (a Configuration table, host spec bounds, version, exports, the fuzz's asserts, the story's reads); the look clean |
| NG | 2026-11-20 | 3 (review: public data access called missing in the _spec; a partial copy of layerSource; tile-edge ranges by label only) | 0 | 0 | three low (the import scan's dynamic imports, the variant heads' width, the missing values on core's side) |
| ML1 | 2026-11-20 | 1 (a section's range unbounded on units and molecules) | 1 (a unit's positions: from 0, decided) | 0 | the core table duplicated `AminoAcidNames` (moved to `util/sequence-constants`); five low |
| FX1 | 2026-11-20 | 4 (review: the levels fix's wiring untested; a whole level still built when none was drawn; an empty band while building, against EC-169; a second build budget) | 0 | 2 (look: the "2 hidden" chip over the overview; a narrow frame squeezing the alignment, decided for FX2) | 10 items fixed, 4 partly, the rest to FX2 with the decided designs |
| SZ (1) | 2026-11-20 | 0 | 0 | 0 | review: no behavior change; an import cycle grown to 7 modules, the ranges modules unlisted, three nits |
| DP | 2026-11-20 | 6 (review: malformed descriptions crashed; values on the wrong unit; row-order content; tile nodes losing identity; narrowed runtimes compared by object; a removed node's selection unreported) | 0 | 0 | public names trimmed to the design's |
| BW | 2026-11-20 | 6 (review: range reads unused by any host; memory exhaustion on malformed files; If-Range breaking servers; an abort shared across readers; maxPages not a memory bound; offsets past 2^53) | 0 | 0 | core helper reuse left (new public names) |
| ML2 | 2026-11-20 | 2 (review: a reference to a unit without a system threw; a ranged unit over a source drew no letters) | 2 (the atom form's meaning, decided: monadic, then `elements`/`groups`; equivalents hooks) | 0 | internals exported; model.ts 1410 lines |
| MC | 2026-11-20 | 5 (review: a crafted operator expression exhausted memory; quadratic lookups and a rethrown error; branch links on the complex; PO4 classed by what was modeled; the fuzz's own crash) | 0 | 0 | the io allowlist narrowed; the re-count made independent; the lazy chunk's imports reported (90 kB) |
| FX2 | 2026-11-20 | 3 (review: the rail not back when a narrow frame hid the overview; a large group's card cut; a possible resize loop) | 2 (no tree toggle; model issues unwired) | 0 | the ruler's allocation; the look clean |
| MK | 2026-11-20 | 2 (review: provenance on the wrong run; overlaps on `to` unchecked) | 2 (the ranking rule and an unused via, decided) | 0 | a weak property test; model.ts 1483 → 1370 |
| SZ (2) | 2026-11-20 | 0 | 0 | 0 | review: no behavior change; a test anchored on shims, viewer-motion's viewports, two nits (to the next fix lane) |
| QE | 2026-11-20 | 3 (review: an empty column printed as everything; a pending property dropped in node filters, EC-169; column lengths unbounded before zipping) | 0 | 0 | `:index` scanned rows × residues (4.2 s → indexed); monadic laws untested; internals exported; PYL; the default viewer 2.7 kB past its bound, moved to 1,290 kB (the author) |
| ST1 | 2026-11-20 | 6 (review: the entry unresolvable by name; colorRuns not compiled as the track; past the tile bound read complete-empty, EC-169; an overlay unit's maps ignored; the gate lighting unobserved-only units; toUnits quadratic) | 1 (resolve cannot hold a map onto a source's unit: drift) | 0 | test gaps (units past 0, EC-4, atoms); re-exports internal; Seq*'s own surface accepted by the orchestrator; nothing drawn |
| AC | 2026-11-20 | 0 | 2 (the `alignment` map kind owed, a `segments` stopgap; one fixed name for every alignment's columns: drift) | 0 | weak tests (a path compared with itself, a test passing on the parent, the hover's section unchecked, no gap end drawn); EC-48 claimed whole; a duplicated kinds list; the look deferred to MX; core/render 2.28 s at load 8.9 (the gate measures it quiet) |
| HP1 | 2026-11-20 | 4 (review: slots hit while not drawn; rows tracks' states over the tree; navigator pictures given headers; marks at negative x) | 3 (a label hit for `Emphasis.header`; the "⋯" order and backdrop; marks in a narrow column: the author amended the design) | 0 | `HeaderFacts` renamed; `nodeAt` a private copy of repr/tree; test gaps; no new part visible without a host (HP3's stories); minimal viewer 915.8 of 920 kB |
| MR | 2026-11-20 | 4 (review: EC-95's label over the start number; the allocation tier order-dependent, a test helper's megamorphic copy; an atoms-form block lit as residues; a view after an edit never opening its end plans) | 1 (`block.form` froze the unit list: a block child overrides one unit, the author) | 3 (look: "Fc A1"; the code mode never shown; the story's data against its text) | sizes past 820 / 920 twice trimmed (lazy molecule letters), the bounds moved to 830 / 930 (the author); three brief items built late |
| CU | 2026-11-20 | 2 (review: the real example as on the parent, every contact on the first heme; a private copy of the unit step) | 1 (an unnamed pair across split copies: never on a guessed copy, the author) | 1 (look: the split hemoglobin's contacts all on the first heme) | weak issue tests; no tile-edge test; `ResolvedBlock.copy` internal; CE's parked branch deleted |
| MX | 2026-11-20 | 4 (review: chip pills cut below `labels.minChars`, a core cause; the 1FIN twin counting arcs only; module state in molecule letters; rival vias a silent missing link) | 3 (a section's `via` scope and a component query pending: drifts decided; the pairs layer reads text columns only: TODO, the author's) | 2 (look: "AT…", "2…" chips; transplants' contacts not said to come from other entries) | the projection bug (a blocks section unlinked) fixed in core; half the brief built late; T-DM1's CCD ids placeholders; storybook folder 3.0 s |
| HP2 | 2026-11-20 | 7 (my check: a legends race, imports outliving teardown; review: settle not waiting for the ports; clicks lost before the late chunk; `seq:filtered` on one layer; DOM parts drawn under chips; the hidden chip a mark at 112 px, core) | 3 (the default viewer past 1290: `@molstar/ui` fixed, the author; lazy chunks allowed, the author; the chip in the minimal viewer, the author) | 1 (look: "1 hidden" a bare mark) | red twice at my check; internals exported; submenus from items pending for the author; default viewer 1289.9 of 1290 |
| ST2 | 2026-11-20 | 2 (review: a double-click focus never reached the structure; ST1's paths empty over an empty overlay) | 2 (unnamed jumps recorded as focus: kinds added; the link's hover allocation target: for the author, gaps.md) | 0 | the default viewer 0.4 kB past 1,290 (moved to 1,330, the author); minor copies; nothing drawn (ST3's stories) |
| FX3 | 2026-11-21 | 4 (review: the columns alias not serving old sources; the pending state unwired; three quick representatives dropped; resolve/model.ts at 1433) | 0 | 0 | 9 items fixed, 1 partly (the viewer folder 2.78 s, gaps.md); T-DM1's CCD ids open (a fetch, gaps.md); the look's items and SZ's four to the next fix point |
| ST3 | 2026-11-21 | 5 (review: the camera frame misplaced; a late non-null 3D hover taking the pointer's; focus-and-back asserting nothing of where back lands; color not toggling; unplaced selectors widening to whole units) | 0 | 1 (the color command's icon drawn as a letter) | `StructureSide.onHover` takes `at` (gaps.md); the host command's change kind to the TODO; within 5 Å faked by residues; the folder 2.2–2.5 s (FX3's partial) |
| HD | 2026-11-21 | 3 (review: the edits counter's listener dropped on mount; a dropped section silent; an unchecked cast in `withAdapter`) | 1 (no section can show a source's complex: drift/host-complex-sections.md) | 1 (chain D's track lingering after Drop chain: not reproduced, blanks without input) | `withAdapter` generic (gaps.md); a ligand's author number from `pdb_seq_num` (1FIN's ATP 299); 2OCJ shows `b_factor`, pLDDT by a synthetic table |
| HP3 | 2026-11-21 | 5 (review: three hints contradicting the code; keyboard entry into DOM parts losing focus, HP2's; the DOM-parts, hidden-list and ring twins not testing what they named) | 1 (one slot in the measured 80–160 px column against the design's "the eye and the color fit in 112 px": gaps.md) | 2 (host badges a row too high; a chip beside a mark in a narrow column: core now all or nothing) | the stories' screenshots stale in the pane at times: checked by the canvas's pixels and the DOM |
| FX4 | 2026-11-21 | 1 (review: the short CDR labels lost the scheme from the card) | 1 (a rows gap mark: TODO) | 0 | 4 items fixed (SZ's four, the IMGT ruler, the CDR labels, the ruler switch); every item card names its track (gaps.md); the viewer folder within noise; the last column cut open (a browser look) |
| EX | 2026-11-21 | 7 (review: the minimal viewer grown 4.6 kB; versions mixing 0.x and 1.x; a failed extend unmapping the viewer's behaviors; a newer type's params refused; late parts dropped silently; the navigator picture unsaid; a property test that could not fail) | 0 | 1 (the demo scheme colouring by index: TAD2 in proline-rich's colour) | the minimal viewer +3.7 kB to 929.3 of 930 (gaps.md); public surface for acceptance (gaps.md); the size script measures the `_measure` entry once a lazy chunk exists |
| LE | 2026-11-21 | 3 (review: the compact form's names unread by keyboard, pre-existing; the 0.5.0 history missing `legend`; the shared-legend rule untested) | 0 | 1 (the story's keys shown twice) | compact steps untested (the binned params unconfirmed); the minimal 929.5, headless 829.1 |
| CN | 2026-11-21 | 4 (review: the minimal viewer at its bound with LE; no public position↔x for a host navigator; exportSvg's two sources of pictures; a host navigator's missing picture unsaid) | 0 | 1 (the host bar's drag: not reproduced in jsdom, pans after the fix round) | headless 829.8 → 806.0 (navigators lazy in renderSvg); the minimal 929.5 (docs trimmed from language.ts) |
| NP | 2026-11-21 | 4 (review: the channel suffix before the merge; the Mol* converter's title; the DOM backdrop untested and without the ellipsis; the specs) | 0 | 0 | 5 bugs and 9 stories from the novice pass; F5 as designed, F7 not reproduced; the minimal viewer 929.9 of 930 |
| SV | 2026-11-21 | 5 (review: a dropped file's rows and length trusted; no size check before reading; maps keyed on the author's chain, not `label_asym_id`; my look: the structure section unaligned with P04637; the back hover asserted missing where the spec's per-chain sequences serve it) | 2 (SIFTS lengths unchecked; the drop target by name only) | 1 (the hint naming the wrong icon) | fixtures 28 kB with licences; the structure section's fit alignment to the TODO; no core change |

## How to verify

The loop's Checks on `packages-seq`; the plan's ticks against `git log --merges`.
