# Gaps: for the author's review (the fourth slice, overnight)

Decisions made while the author is away (loop.md, "Running while the author is away"), and public API changes made
under the author's standing permission (`state.md`, Decided). Each a yes or no question with what was done.

## Public API changes (allowed: `packages/seq`, core's `repr`, `@molstar/ui`)

- **ST2**: `JumpKind` gains `zoom` and `move`; `view.jump(move, record?, kind = 'move', location?)`;
  `NavigateContext.location`; core `toUnits`/`fromUnits` over an empty overlay read the model alone. Keep? Recommended
  yes (the review: unnamed jumps were recorded as `focus`).
- **ST3**: `StructureSide.onHover`'s callback takes the event's time (`at?: number`); `linkStructure` drops a side hover
  older than the pointer's last own hover (the review: a late non-null 3D hover took the pointer's). Keep? Recommended
  yes (EC-141 covered for hovers, not only leaves).

- **HD**: `withAdapter` is generic over the registry's entry type (a viewer's registry passes through); the examples'
  chain tables gain `ChainResidue.b`, `ChainTable.model`, `HostModel`. Keep? Recommended yes.- **EX**: the extension surface. core's repr: `KindSpec.extensionParam`, `KindSpec.fallbackParam` (any `x:y` key kept
  in every repr language). Seq* core: `SeqExtension`, `defineExtension` (it checks and composes), `SeqContextSpec.extensions`,
  `SeqContext.kinds`/`schemes`/`extensions`, a `scheme` registry entry, `LayerType.shape`, `LayerKinds`, `LayerShape`,
  a `kinds` argument on `resolveView(s)`, `ResolvedLayer.shape`, `SeqSpecOptions.extensions`, issues
  `seq.lang.layer.fallback`, `seq.view.layer.extension`, `seq.view.extension.version`. Viewer: `SeqViewerExtension`,
  `defineViewerExtension`, `SeqViewerSpec.extensions`, `extend()` (async, a lazy chunk, resolves to its undo),
  `BehaviorHandle.commands`. Language 0.5.0: `spec.extensions`, `layer.fallback`, `layer.type` a word, `<id>:<name>`
  for `layout`, `heights`, `content`. Drifts: `extend` async with no public `context.extend`; params checked by the
  descriptor's `check`; no type-table version in the draw-list key (a new resolved view is a new renderer); late
  header entries and settings refused. Keep? Recommended yes.

## Decisions as recommended

- **ST2**: the link sends every jump that names a Location (`focus`, `reveal`, and back/forward over them); `fit`,
  `zoom`, `move` never go out. Keep `reveal` going out? Recommended yes (a reveal is a navigation to a place the
  structure can show).
- **FX4**: every item's hover card and announcement now name its track (a Track row: "Domain A, p 3–8, Domains"), so
  two tracks' items of one label are told apart (Kabat's and IMGT's CDR-H1); ux-foundations amended. Keep it for every
  card, or only where two tracks share a label? Recommended: every card (one rule, and a row under the cap).

## For the author (no default taken)

- **T-DM1's CCD ids** (FX3, item 6, open): the repo holds no CCD data for DM1 or the SMCC/MCC linker. Fetch their
  chemical component pages from the PDB once (two small files, recorded with their licence), or give the ids? A
  fetch is asked before it is made.
- **The viewer's quick-lane folder still over 2 s** (FX3, item 1): 2.78 s, its time imports (28%) and jsdom setup
  (20%) across 65 files. Fewer jsdom files, a lighter DOM environment, or a viewer-folder budget of 3 s?
  Recommended: fewer jsdom files first (the next fix lane), the budget unchanged.

- **The link's hover allocation**: a 2OCJ hover allocates 11.3 kB through the link (9.2 kB in `toUnits`, 1.7 kB per
  lit unit; was 31.8 kB), against performance.md's 0 B for warm hover loops; the target and the hover's Location are
  new each hover. Set a target for the link's path (performance.md names none), or ask for pooled results? Recommended:
  a bound per lit unit (2 kB) in performance.md, pooling later if a large assembly shows it.

- **HD's drift**: no section can show a source's complex (`drift/host-complex-sections.md`): keep static sections with
  an issue, sections over a source's complex, or the source bringing `values` nodes? Recommended: sections over a
  source's complex in the next slice's view-model work, the issue until then.

- **HP3: one slot in the measured label column.** `actionColumns` gives `floor(labelWidth/20) − 3` slots, and the
  measured 80–160 px column usually leaves one (Hide; Color only in the menu), while track-header-parts' default says
  "the eye and the color fit beside a label in 112 px". Keep one slot and amend the rationale, or reserve two slots
  (the column wider, or the slot rule changed)? Recommended: amend the rationale (the menu holds every entry; a wider
  column costs every view).
- **EX: the minimal viewer grew 3.7 kB** (925.6 → 929.3 of 930; the brief: "the minimal viewer is not grown"): 2.4 kB
  are the reader's fallback parts (the fallback chain, the issue once per extension, version reading, the navigator
  issue, language and catalogs), 1.2 kB composition entry points; `checkExtensions` and the tables moved behind
  `defineExtension`, late `extend` lazy. 0.7 kB is left under the bound, and CN and LE come next. Accept and move the
  minimal bound (e.g. 940 kB), or trim design items (the info-level version notice, the navigator issue, the builder's
  `extensions`)? Recommended: accept, move the bound to 940 kB (a reader without the extension must say what it lacks).
