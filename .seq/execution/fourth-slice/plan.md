# Plan: the fourth slice

**Reconciled 2026-11-20, after the third slice's last lanes merged** (its use gate's look, M3 and finalize remain).
Stage 0 is done and Gate 0 is decided (below); the build starts when the third slice closes. "What the third slice
changes in how this one runs" is provisional until the third slice's `foldback.md` is decided, and is rewritten from it
then.

Work id `fourth-slice`. Version 1's third milestone, **integration** ([version-1](../../design/version-1.md)): structure
maps and the structure showcase (on a fake 3D viewer, decided), services, bigWig and range reads, the resolved view's test consumer (Nightingale
dropped), archives. With it, what the third slice and the author carried here: the molecules design built
([molecules](../../design/molecules.md): units, copies, the 0.5.0 bump), copy ends on units (CE's parked branch and its review), the designs
agreed in stage 0 and built here (structure syncing as a library, track header parts, seq data ports, the extension
openings: `ext:` kinds, custom navigators, legends' surface).

## Where this is heading

| Slice | Milestone | In one line |
|---|---|---|
| first (closed) | the end-to-end slice | the tree core, the language, resolve, draw lists, Canvas2D and SVG, the viewer's models, data by region, measured |
| second (closed) | 1, core editing and addressing | patches, history, previews, addresses, refs and links, selection, security bounds, MessagePack, navigation |
| third (closing: use gate's look, M3, finalize) | 2, figures | the axis, text, ring and block layouts, alignments, dot plots, arcs and logos, themes and legends, stories, navigators, PNG and the CLI, the datasets |
| **fourth (this plan)** | 3, integration | molecules and units, structure maps and syncing (a library) on a fake 3D viewer, track header parts, seq data ports, services, bigWig with range reads, the resolved view's test consumer, archives, extension kinds, custom navigators, legends |
| fifth | 4, Genomes | layers across sections, synteny, contact maps, MAF, two assemblies on one axis, the genome formats, the Rondo-like ring as an extension package |

**No Mol\* adapters, no live 3D view** (decided, the author, 2026-11-19: "let's only do a demo with fake 3D viewer ...
There will be a new mol viewer work after we are done with the sequence viewer which will replace prism, so i don't
want to waste work"). Milestone 3's "Mol\* showcase" is the fake structure viewer in Storybook on the real library (translation and
`linkStructure`), maps and messages; Mol\* `Model` adapters wait for the new mol viewer.

**No MolViewSpec types in Seq\*** (decided at the draft, 2026-11-19: "this should be provided by the application layer
anyway"): the library speaks Seq\*'s locations; a host's adapter translates to MVS and back. The fake viewer is a
host, and its adapter does the translation. **No structure port** (decided, 2026-11-20: "this is just a wrapper around
events the viewer should expose normally ... not call it a port"): translation and wiring as a library.

## What the third slice changes in how this one runs

From the third slice's state and its fold-back (decided by the author, 2026-11-20: merged; the rules in the
principles, "Folded back from the third slice"); the bullets below are checked against them.
- **Lanes: two or three at once, each strong or fast by its work** (the author, 2026-11-20: "2-3 lanes, all either
  heavy or light depending on context. that way we can do things in parallel better"; it replaces "1+1-2"): a lane
  is strong when its work is design, a model or ambiguous, fast when it is specified (FX, stories, tests, read-only).
  What runs together is limited by files, not by kind: never two lanes on the language, never two on core's store or
  resolve, never two touching one file, never two heavy test runs at once. A lane's review and fix round count as that
  lane; they run before a lane that builds on it starts.
- **Designs agreed before their lanes**: the third slice's stages 1–4 ran unattended on agreed designs; its later
  items waited on drafts. Here every design was drafted, self-checked and decided in stage 0, before a build lane.
- **The author's findings arrive mid-slice**: a finding that needs a design goes to a stage 0 style draft; one that needs
  none goes to [fixes.md](fixes.md) and is taken by the next FX lane (below). Why: the third slice's twenty finding
  lanes (CH to DB; bug 6, missing design 7, look 6 by the end) each cost a brief, a review and a merge.
- **A lane whose fix needs an unbuilt design is parked, unmerged, with its review kept** (CE), not merged with known
  high findings.
- **A lane's public API is accepted by the author before merge** (MV, DB): a brief names the public surface it expects;
  anything past it stops the lane.
- **Read-only fast lanes prepare the author's questions and the close** (SR, EO, DA, GR): cheap, kept.
- **Use gates do not block the loop**: the author's use comes when they are back, its findings fixed after; the novice
  pass runs before the close.
- **The third slice's fold-back, applied** (principles, "Folded back from the third slice"):
  - *Defaults designed with a picture* (§defaults-designed): each stage 0 design that adds something drawn (the
    molecules' sections, the header parts, the fake 3D viewer) has a Defaults section, reviewed at stage 0's gate.
  - *Drift notes are gaps* (§drift-notes-are-gaps): every gate's message lists `drift/`'s open items with the gaps.
  - *Edge cases tested by name* (§ec-tested-by-name): a lane's tests name the catalog items they test;
    `node .seq/_checks/ec-ids.ts` runs with the gate's checks, `--worklist` at the close.
  - *Configuration tables checked* (§config-tables-checked): the core's and the viewer's `config-table` tests run in the
    quick lane; a lane adding a config item adds its row.
  - *Features move with the merges* (§features-move-with-merges); *a changed test expectation explained*
    (§test-expectation-explained); *a matrix cell built when a test names its gesture* (§matrix-cell-by-test).
  - *Use is continuous* (§use-continuous): the author's use counted by kind as it comes; the novice pass before the
    close plays a host's integration task (sync a selection with the fake 3D view, wire a seq data port), not a
    reader's.
  - *The heavy lane alone, at two workers* (§heavy-lane-alone); *a hot loop's representative in every package that
    composes it* (§hot-loop-representatives).
  - **The findings table is kept per lane** (B6): the merge step adds the lane's row (bug, missing design, look; review,
    look, use) to `state.md`'s table before the merge note, for every lane, not only the early stages. Why: the third
    slice's table stopped at stage 4 and its last forty lanes could not be counted at the close.
  - **Each lane keeps its own quick-lane time** (B7): a lane's checks include its packages' quick-lane time against the
    2 s guide; a lane over it moves its own slow tests to medium (a quick representative kept) before its review. Why:
    three lanes (QT, QT2, QT3) moved what earlier lanes had added.
  - **Reviewers in their own worktrees** (B8): every reviewer and every lane works in a worktree of its own, never in
    another lane's. Why: a reviewer switched another lane's worktree in the third slice.
  - **The lane limit is one line** (B8): `loop.md` holds the author's current limit as one line, edited when the
    author changes it (now: two or three lanes, each strong or fast by its work); `resume.md`'s shape is the template for any pause.
    Why: the author changed the limit five times in the third slice.

## Done when

- **Molecules**: `molecule`, `complex`, `unit`, `modified`, `compounds` in language 0.5.0 (one bump, unreleased); units
  and atoms systems, copy groups derived, the atoms map at a path's ends, the modified index, the query forms; the
  mmCIF runtime (`runtimes/mmcif`, opt-in) with its conformance suite on 1EHZ, 1FIN, 4HHB and a glycoprotein; a
  section of a complex drawn with `form`, the letters' mark and code, hover; DT, DC and DG on the new kinds (changed
  once), the Fc glycan and ADC views ([molecules](../../design/molecules.md), Consequences).
- **Copy ends on units**: a pair's end on split copies by unit, every finding of CE's review a test (hemoglobin split
  keeps its contacts; copies counted as the resolver's).
- **Structure maps**: the `structure` map kind kept (decided), `segments` built, a map's `priority` and a section's
  `via` to choose between two maps joining one pair of systems ([structure-syncing](../../design/structure-syncing.md),
  Choosing between two maps; `StaticModel.path(from, to, via?)`); EC-22 to EC-27 on every encoding with a ranked pair,
  a tie and a `via`.
- **Structure syncing as a library** ([structure-syncing](../../design/structure-syncing.md): not a port, no behavior):
  `@molstar/seqx-core/structure` (`structureFacts`, `unitRows`, `toUnits`, `fromUnits`, `colorRuns`) and
  `@molstar/seqx-viewer/link-structure` (`linkStructure`), with the viewer's four additions (`HoverModel.setLocation`,
  `AddressModel.revealLocation`, `view.lastJump`, `viewer.retiredSignal`); hover, selection and focus both ways, "color
  by this track", a structure's chains mapped onto sequences, on the fake structure viewer (its own chain table, its
  host's adapter translating to MolViewSpec); p53 with 2OCJ and the Complex with 1FIN; the stories with headless twins
  on a fake clock (one with latency). No MolViewSpec types in Seq\*.
- **Track header parts** ([track-header-parts](../../design/track-header-parts.md)): entries in fixed slots drawn on
  the canvas and a host's DOM parts over the same geometry (`headerPartsDom()`), the hidden-tracks chip, the menu
  behind its port, the secondary click a behavior param, the cursor on headers; `Features/Track headers`.
- **Seq data ports** ([seq-data-ports](../../design/seq-data-ports.md)): Seq\*'s own reading is the default adapter
  of the one source protocol; `withAdapter` replaces one per format; `source.ref`, `values.field`, `changes` and
  `reattach` (a host's updates reopen a source in place, keys kept); `StructureDescription` and `structureNodes` in
  `seqx-core`'s `io/structure` (shared with MC's converter); the fake host's adapter and a `Host data` story. No
  Mol\* adapter, no package of its own.
- **Data**: bigWig read by region and zoom level through range reads of assets (p53's conservation signal over the
  gene's region); services (UniProt, PDBe SIFTS, AlphaFold DB) composed in Storybook from recorded fixtures with a live
  switch; archives read and written over the assets port with EC-91's bounds; the predictions-onto-a-structure and
  addressing showcases ([examples](../../design/examples.md), Showcases).
- **The resolved view's second consumer**: a test-only plain-SVG consumer draws a subset from the resolved view and
  the store alone (EC-183). Nightingale is dropped (decided, the author, 2026-11-19): nothing of it is a dependency, a
  story or a lane.
- **Extension openings** ([extension-openings](../../design/extension-openings.md)): part 1, a `SeqExtension`
  descriptor and namespaced layer types (`<id>:<name>`; plain params on an extension's own type, `<id>:<name>` keys on
  others' nodes), registered atomically with undo, the fallback drawn and said when absent (EC-235); part 2 (cut first),
  a navigator's picture as a registry entry with the window and `NavigatorModel` public; part 3 (cut second), legends'
  library surface (the `legend` layer param, `mountLegend`'s options, a compact form).
- Measured: a hover message for 2OCJ, `maxUnits`-scale hover, a bigWig region decoded, a header pointer move over 200
  tracks, the bundles (the mmCIF and bigWig runtimes as their own rows, the minimal viewer not grown).

## Gate 0 (the author, before any build lane)

Decided (the author, 2026-11-20: "Accept all"): every question as recommended; the lane default one strong and one
or two fast at once (not one at a time).

1. **When it starts.** The build lanes after the third slice closes (finalize, fold-back, record); DR and the stage 0
   designs before, since they write only this folder. DR copies the two drafts into `fourth-slice/drafts/` and revises
   them there (the third slice's folder is closed after its record). Recommend **yes**.
2. **Scope under the lane limit.** About 30 lanes, the strong ones one after another, is long. Keep milestone 3 and
   the carried items, and if usage runs short cut in this order: custom navigators (CN), legends' surface (LE), the
   resolved view's test consumer (NG), archives' writer (AR's write half). Recommend **this order**.
3. **`ext:` layer kinds here** (version-1: "`ext:` layer kinds ... land before" Genomes' ring; GR's question 3): an
   extension descriptor and `ext:` layer kinds with their fallback, as lane EX near the end, so the fifth slice opens on
   its datasets. Recommend **yes, here**. **Answered** by XD's design (2026-11-20): here, as lane EX on
   [extension-openings](../../design/extension-openings.md) part 1.
4. **Range reads in `@molstar/core`** (a public API change, [base-library](../../design/base-library.md), ports.md: a
   page cache by asset and page, ref-counted, `Blob.slice` for files; `maxBytes` on `file://` honoured, a TODO item):
   BW drafts the exact surface as its first commit and stops for your acceptance. Recommend **yes**.
5. **bigWig's home and data.** An opt-in entry `runtimes/bigwig` of `seqx-core` (as mmCIF's, molecules gap 12), its
   fixture generated by a test-only writer on core's `PCG` (both zoom levels, edges, an empty region), p53's signal
   over the gene's region generated the same way; a real UCSC slice only if you want one (licence to check). Recommend
   **opt-in entry, generated**.
6. **Services' home.** In the Storybook as application-layer runtimes (features.md: "the storybook package first"),
   SIFTS read there; base-library's `@molstar/io` SIFTS row stays planned until a second host. Tests read recorded
   fixtures only; the live switch is the Storybook's, through the asset policy's allowed hosts. Recommend **yes**.
7. **Nightingale**: decided (the author, 2026-11-19): dropped; the check is our own test consumer (NG).
8. **Copy ends**: rewritten on units (lane CU) with `seq/ce-copy-ends`'s tests and its review's six findings as the
   cases, the branch deleted after, rather than rebasing it. Recommend **rewrite**.
9. **Seq data ports, how far**: the design (SP) and a fake host's adapter in Storybook (reading the fake viewer's chain
   table into Seq\*'s nodes); no Mol\* adapter. Whether MC reads through the same converter is SP's decision.
   Recommend **design and the fake adapter only**. **Answered** by SP's design (2026-11-20): the design and the fake
   adapter only; MC reads through the shared `StructureDescription` converter.
10. **The third slice's not-met edge cases** (`ec-outcomes-draft.md`): EC-115, 116, 135 and 195 (built, untested) to
    lane UT (**removed 2026-11-20: all four are tested**: EC-115, 116 and 195 by ET, merged `50956fce0`, EC-135 by the rows
    layer's tests); EC-99 (ligand chips flowing in text) to MR; EC-59 (two grids on one block pair) and EC-209 (ribbons
    between sections) to the fifth slice (contact maps, layers across sections); EC-63 as the third slice's drift note
    decides (slices kept, glyphs to `TODO.md`); EC-103 (track groups and summaries) see 11. Recommend **as listed**.
11. **Version 1's planned rows outside milestone 3**: track groups and summaries (EC-103), search and results, the
    context menu beyond what HP2 brings, patches on streamed things, keep this region, typed configuration and presets,
    the query forms `search`, `onto`, `expand`, `group`, EC-29's indel carets (a `pairwise` map). None here; listed at
    this slice's finalize for you to place (a version-1 close-out, or the fifth slice). Recommend **none here**.
12. **A user's own file in the Storybook** (carried from the third slice with services): drop a file onto a view for
    the formats built (`seq`, FASTA, Stockholm, the records, mmCIF through MC), in SV; GenBank deferred. Recommend
    **yes, without GenBank**.
13. **Use gate**: a short one after ST3 and HP3 (the first integration drawn), not blocking; the novice pass before
    finalize with a task script for the structure stories. Recommend **yes**.
14. **Datasets**: 2OCJ recorded for p53 (examples.md), an IgG1 Fc glycan entry chosen at recording (CC0), the ADC
    construct generated; recorded in MX. Recommend **yes**.

Not questions (decided): the fake 3D viewer; no MVS types in Seq\*; the `structure` map kind kept; one language bump,
0.5.0, opened by ML1 and added to by later lanes in turn (one lane at a time: no conflict), closed at finalize; the
secondary click a `headers` behavior param; canvas slots and DOM parts both; `links.hover` off by default, on in the
examples.

## Testing strategy

- **Every lane**: as the third slice (quick-lane tests beside the code; a story composing more than one consumer has a
  headless twin on a fake clock; drawing checked against the data at the fit view and across a tile edge; a fix lands
  with a test shown failing on its parent; the two consumers its shared pieces serve, live in one test).
- **Conformance suites, offline**: the mmCIF runtime's (every entity, chain, branch link, `struct_conn` row and author
  number read, against counts recorded from the file; an entity without canonical letters; an altloc); bigWig through
  the shared source suite with its synthetic twin; each service runtime against its recorded fixture; archives' bounds
  (EC-91) by fuzz. No test reaches the network.
- **Units and maps**: property tests (copies derived equal for any unit order; atoms collapsed then expanded cover what
  they started from; a path's result unchanged by a molecule's atoms form; EC-22 to EC-27 equal on implicit, runs and
  expanded encodings).
- **Structure syncing, two consumers**: Seq\* and the fake viewer in one test: both translation round trips; every lit
  cell in the fake viewer in the image of Seq\*'s hover and back; the fake viewer's colors equal `colorRuns` at each
  observed residue's preimage; echoes never sent back; a stale hover ignored.
- **Performance representatives**: one 2OCJ hover message within a bound and none for an unchanged image; a
  `maxUnits`-scale hover gated to chains (heavy); 200 merged copies and a glycan chip at frame time (EC-67); a bigWig
  region decoded under 50 ms (heavy, `itWithin`); a header pointer move over 200 tracks at 0 bytes (the opt-in tier).
- **Sizes**: rows for `runtimes/mmcif` and `runtimes/bigwig` (under 30 kB, loaded only when composed), the minimal
  viewer not grown by the `headers` behavior (opt-in in a minimal composition) nor the structure library (its own entries).
- **At every gate**: quick, medium and heavy for the touched packages (heavy alone), the allocation tier, sizes, the
  Storybook typecheck, lint, format, the `.seq` checks; `pnpm build` then `pnpm check:workspace` when `@molstar/core`
  or `@molstar/io` changed (ML's standard components, BW's range reads, MC's allowlist).

## Stages and lanes

Each lane: a brief in `briefs/` (named below; every lane of stages 1–5 has one), its `_spec` and Configuration table
kept, a v6 changelog entry (in its report), an independent review and one fix round, the look in a visible pane in
light and dark before the merge. Briefs in `briefs/` (all written 2026-11-20): `sp-seq-data-ports`, `xd-extension-drafts` (stage 0, done); `fx-fix-lane`; `ml1-language`, `dp-data-ports`, `ml2-model`, `mc-mmcif`, `mk-structure-maps`, `ac-alignment-targets`, `mr-molecules-render`, `cu-copy-ends-units`, `mx-examples`; `hp1-header-line`, `st1-structure-core`, `hp2-header-behavior`, `st2-structure-link`, `st3-fake-viewer`, `hd-host-data`, `hp3-header-stories`; `bw-bigwig`, `sv-services`, `ar-archives`, `ng-test-consumer`; `ex-ext-kinds`, `cn-custom-navigators`, `le-legends`; `novice-pass`. Strong: the strongest model, high effort (design, model, ambiguous). Fast: the fast
model (well specified, stories, tests, read-only). Each brief cites the design sections it builds; a lane whose design
disagrees with its code writes a drift note and stops.

### Stage 0: designs (done)

- [x] DR (done 2026-11-19, `3e1cc3320`; promoted 2026-11-19 `b85b59e6d`, structure syncing again as a library
      2026-11-20 `64dd776ca`; strong; `../third-slice/briefs/dr-drafts-revision.md`) the two drafts revised by the
      author's decisions and agreed: [structure-syncing](../../design/structure-syncing.md) (a library, no port, no MVS
      types; the `structure` kind kept; `priority` and `via` for two maps of one pair) and
      [track-header-parts](../../design/track-header-parts.md) (canvas slots and DOM parts, the secondary click a param).
- [x] SP (done 2026-11-20, gaps decided `f8561bb26`, promoted `2fb469bc5`; strong; `briefs/sp-seq-data-ports.md`):
      [seq-data-ports](../../design/seq-data-ports.md) (`FormatRuntime` kept, one adapter a format with `withAdapter`,
      no package, `source.ref`, `changes` and `reattach`, `values.field`; MC's converter shared).
- [x] XD (done 2026-11-20, `f8561bb26`, promoted `2fb469bc5`; strong; `briefs/xd-extension-drafts.md`):
      [extension-openings](../../design/extension-openings.md) (gaps 1, 3–8 yes; 9 no, no extension loading in the CLI;
      2 changed: no `custom.<id>`, plain params on an extension's own type, `<id>:<name>` keys on others' nodes).
- [x] UT **removed** (2026-11-20): EC-115, 116 and 195 are tested by ET (`50956fce0`), EC-135 by the rows layer's tests
      (`rows-layer.test.ts`, `rows-model.test.ts`). The third slice's drift left open at its close, if any, goes to
      [fixes.md](fixes.md) (FX1) or the stage it touches.
- [x] GATE (2026-11-20): the three designs agreed; Gate 0 decided (`cb740554d`). Nothing of stage 0 remains but the
      third slice's close (finalize, fold-back, record), which opens stage 1.

### Fix points (fast; [fixes.md](fixes.md); `briefs/fx-fix-lane.md`)

- [x] FX1 (merged 2026-11-20) before stage 1 (beside ML1 when its files are not the language's): the open items of fixes.md (the ring
      letters, the three DS looks, LX's dropped label) and the third slice's use gate's findings.
- [ ] FX2, FX3, FX4 after stages 1, 2 and 3's gates; FX5 before M4 (the use gates' findings). A fix point with nothing
      open is skipped. An FX lane takes only what needs no design.

### Module split (strong; §module-size; `briefs/sz-module-split.md`)

- [ ] SZ (part one merged 2026-11-20: tracks.ts 2871 → 865, rows.ts 2520 → 413, ranges.ts 1521 → 464; part two, `viewer.ts`, after DP) after FX1 and before HP1, while no molecules lane is in the files it splits (else after stage 1's gate):
      the four files over 1500 lines (`render/tracks.ts`, `render/rows.ts`, `render/ranges.ts`, `viewer/viewer.ts`)
      split by what changes together, and `models/view.ts` (1466) looked at; behavior unchanged, the baseline lowered
      by the orchestrator (`file-size.ts --update`) at the merge. The author, 2026-11-20: "Yes, add all three".

### The query language's elements (strong; `briefs/qe-query-elements.md`; designs: locations-and-queries.md, properties.md)

- [ ] QE after MC (its 3SGJ fixture) and ML2: `elements` and `groups`, a generic `filter` over typed properties with
      positional operands, modifiers and combinators (`merge`/`intersect`/`except`), the `query-form` and
      `query-property` registry kinds, short and long aliases, selections as records or columns; `atom`, `entity` and
      `component` removed (the author, 2026-11-20: "we should not really have (atom), we should have (elements ....)").
      Molecules track, beside MK when their files differ.

### Stage 1: molecules, maps and copies (strong lanes in order; [molecules](../../design/molecules.md), The fourth slice's lanes)

- [x] ML1 (merged 2026-11-20) (strong, first) the language: the kinds, params, checks and builder (molecules.md, Decision and Consequences,
      `seqx-core` `lang/`); the `seq` format's tables; `source.ref` and `values.field` (seq-data-ports, Decided, gaps 2
      and 4); the 0.5.0 bump opened and owned (`briefs/ml1-language.md`). The core's standard component per letter
      (`@molstar/core`, a base-library row, decided).
- [x] DP (merged 2026-11-20) (strong, after ML1) the seq data ports built (seq-data-ports, Decision and Consequences): `OpenedSource.changes?`
      and `RegionRequest.field?`, the store's `reattach` (node content replaced by key), `withAdapter` (one adapter a
      format, conflicts listed), `io/structure` (`StructureDescription`, `structureNodes`), the viewer's store model
      reopening in place; the shared adapter conformance suite (`briefs/dp-data-ports.md`). No package of its own.
- [x] ML2 (merged 2026-11-20) (strong, after DP: both touch resolve and the store) the model: unit and atoms systems, the atoms map at a
      path's ends, held pairs with per-row units, copy groups, the modified index, `maxUnits`, the query forms
      (`briefs/ml2-model.md`).
- [x] MC (merged 2026-11-20) (strong, after ML2 and DP) the mmCIF runtime (`runtimes/mmcif`, opt-in): CIF categories to
      `StructureDescription` to `structureNodes`, `values` by `field` (B-factor, pLDDT), the allowlist, the conformance
      suite split as seq-data-ports says, the fixtures re-recorded (`briefs/mc-mmcif.md`). Parallel with ML2 is gone:
      two strong lanes at once break the lane limit.
- [ ] MK (strong, after MC) structure maps: `segments`, the `structure` kind with per-residue status, the map's
      `priority` and the section's `via` with `StaticModel.path(from, to, via?)` (structure-syncing, Choosing between
      two maps; both params join the bump), a numbering's `offset` and `segments` methods (NB's leftovers); EC-22 to
      EC-27 on every encoding (`briefs/mk-structure-maps.md`).
- [ ] AC (strong, after MK; followups: "Yes, slice 4 (0.5.0)") pairs, ranges and maps targeting an alignment's columns:
      the tRNA's `SS_cons` arcs over the alignment, a chain mapped onto its closest seed row
      (`briefs/ac-alignment-targets.md`).
- [ ] MR (strong, after AC) drawing and interaction: a section of a complex, `form` through the chips, the letters' mark
      and code, the pairs layer's per-row ends, hover, `links.hover`, `block.form`, ligand chips flowing in text
      (EC-99), a block's label inline on its first text line (EC-95, LB's leftover), the Design page
      (`briefs/mr-molecules-render.md`).
- [ ] CU (strong, after MR) copy ends on units: CE's parked branch `seq/ce-copy-ends` (`189097018`) and its review's six
      findings as the cases (Gate 0, 8); grids-and-rows.md's "end on the first copy" and its caret replaced
      (`briefs/cu-copy-ends-units.md`). Parked until units exist, so the branch is not rebased: rewritten, then deleted.
- [ ] MX (strong, after CU) the examples changed once: DT, DC, DG on the new kinds; the Fc glycan and ADC views; 2OCJ
      recorded for p53 (its units and the fake viewer's chain table from one mmCIF, an `alignment` map of the same
      pair); DC's leftovers (`@molstar/io` as a dev dependency of seqx-examples, a ligand chip coloured on its own);
      stories, twins, drawing against the data (`briefs/mx-examples.md`).
- [ ] GATE: merged, reviewed, checks green, sizes (`runtimes/mmcif` its own row), gaps, matrix cells (molecules.md's).

### Stage 2: structure syncing and track headers (after stage 1; the agreed designs' lanes)

- [ ] HP1 (strong) core's header line, slots, hits and chips; `TrackLayout.header` replaces `TrackLayout.filter`
      (track-header-parts, Consequences; `briefs/hp1-header-line.md`).
- [ ] ST1 (strong, after HP1; MK done) the translation library, `@molstar/seqx-core/structure` (a new entry, not a
      port): `structureFacts`, `unitRows`, `toUnits` and `fromUnits` (the size gate, `seq.structure.*` issues),
      `colorRuns`; core's `SeqPorts.structure` stub removed; the round-trip test against DP's `structureNodes`
      (`briefs/st1-structure-core.md`).
- [ ] HP2 (strong, after ST1) the viewer: the `headerActions` registry kind, built-in entries, the `headers` behavior
      (reveal, clicks, the menu and its port, long press, the secondary click param), `headerPartsDom()` over
      `viewer.headers`, the cursor on headers, words (`briefs/hp2-header-behavior.md`).
- [ ] ST2 (strong, after HP2; smaller: no behavior, no port, nothing in the composition) the viewer's four additions
      (`HoverModel.setLocation`, `AddressModel.revealLocation`, `view.lastJump`, `viewer.retiredSignal`; accepted as
      the design lists them), then `@molstar/seqx-viewer/link-structure`: `linkStructure` with echoes, coalescing,
      clearing on either end, focus with the jump's duration, `via` (`briefs/st2-structure-link.md`).
- [ ] ST3 (fast, after ST2; HP2 done) the fake structure viewer (`_stories/support/`, importing nothing of Seq\*, its
      own chain table, MolViewSpec validated by `mvs-builder`) with its host's adapter (translation on `structureFacts`
      and `unitRows`, its `side` for `linkStructure`); the stories: hover, selection, focus, chain maps (SIFTS against an
      alignment by `priority` and `via`), "color by this track" as a host's header entry; twins; the Design page
      (`briefs/st3-fake-viewer.md`).
- [ ] HD (fast, after ST3 and DP) the fake host's data adapter beside the fake viewer (it serves 2OCJ and a pLDDT track
      from its chain table through `withAdapter`) and the `Host data` story with its twin (no asset read; equal nodes
      to the default adapter); the Extensibility page's adapter list (`briefs/hd-host-data.md`).
- [ ] HP3 (fast, after ST3) `Features/Track headers` and its twins, the Design page (`briefs/hp3-header-stories.md`).
- [ ] USE GATE (short, not blocking): the author tries the structure and header stories; findings counted by kind.
- [ ] GATE: merged, reviewed, checks green, sizes (the minimal viewer), gaps, matrix cells (both designs').

### Stage 3: data and references (after stage 0; independent of stages 1–2 but for SV's SIFTS, after MK)

- [x] BW (merged 2026-11-20) (strong) range reads in `@molstar/core` (the surface accepted first, Gate 0, 4), the bigWig runtime
      (`runtimes/bigwig`, opt-in) with its levels and synthetic twin, p53's conservation signal over the gene's region
      (`briefs/bw-bigwig.md`). Takes a strong slot after stage 2's gate.
- [ ] SV (fast, after MK and ST3) services in Storybook: UniProt, PDBe SIFTS, AlphaFold DB as runtimes from recorded
      fixtures with a live switch; the predictions-onto-a-structure and addressing showcases on the fake viewer; a
      user's own file (Gate 0, 12) (`briefs/sv-services.md`).
- [x] AR (merged 2026-11-20) (fast, any time after stage 0) archives over the assets port: read and write zip with `@molstar/core/util/zip`,
      entries resolved as assets, EC-91's bounds; every dataset saved and opened (`briefs/ar-archives.md`).
- [x] NG (merged 2026-11-20) (fast, any time after stage 0) the resolved view's second consumer (the author, 2026-11-19: Nightingale
      dropped, "our own test consumer only"): a test-only plain-SVG consumer of letters, ranges, values and variants
      from the resolved view and the store alone; the field test (EC-183) (`briefs/ng-test-consumer.md`).
- [ ] GATE: merged, reviewed, checks green, `pnpm build` and `check:workspace` (core changed), sizes, gaps.

### Stage 4: extension openings (after stage 0; after HP2's registry kind)

- [ ] EX (strong) part 1: `SeqExtension` and namespaced layer types (`<id>:<name>`, plain params, namespaced keys on
      others' nodes), registered atomically with undo, the fallback drawn and said when absent (EC-235), the
      `@molstar/seqx-core/extension` entry; `Features/Extensibility` with an example extension (`briefs/ex-ext-kinds.md`).
- [ ] CN (strong, after EX; cut first) part 2: a navigator's picture as a `navigator-picture` entry; `NavigatorModel`
      and the window public; the four built-ins as entries (`briefs/cn-custom-navigators.md`).
- [ ] LE (fast, after EX; cut second; beside CN) part 3: legends' library surface, the `legend` layer param (a language
      change through the bump's owner, after CN has none open), `mountLegend`'s options, the compact form
      (`briefs/le-legends.md`).
- [ ] GATE: merged, reviewed, checks green, sizes, gaps.

### Stage 5: use, measure, finalize

- [ ] NOVICE PASS (fast; `briefs/novice-pass.md`, the third slice's shape) a task script per new story (structure
      syncing both ways, headers, host data, services, archives).
- [ ] USE GATE (not blocking; at finalize a look at what changed): findings fixed or placed, counted by kind.
- [ ] M4 the heavy lane and the opt-in tiers; the browser at dpr 2 (a 2OCJ hover round trip, a complex of 200 copies,
      a bigWig pan); the bundles; the targets confirmed or moved by the author (below: the minimal viewer's size).
- [ ] FINALIZE: reconcile (designs, `_spec`s, status, features, edge cases' outcomes, `TODO.md`, drift emptied); the
      fold-back; the record; the fifth slice's plan (Genomes, on extension kinds).

### The third slice's carried items, placed

| Item | Placed |
|---|---|
| CE, parked (`seq/ce-copy-ends`, `189097018`), its review's six findings | CU (after MR); the findings are CU's cases, its brief lists them |
| DR, the drafts revision | done (stage 0) |
| Structure integration; track header parts; seq data ports; custom navigators; legends' surface | ST, HP, DP and HD, CN, LE (the third slice's followups, "decided for the fourth slice") |
| Pairs, ranges and maps on an alignment's columns | AC |
| NB's leftovers: `offset` and `segments` methods | MK; the rest (a text section's numbers in a source's numbering, labels keyed other than by position, a records runtime converting labels) stays in `TODO.md` |
| DC's leftovers (`@molstar/io` dev dependency, a ligand chip coloured on its own) | MX |
| LB's leftover: a block's label inline on its first text line (EC-95) | MR |
| EC-99 (ligand chips in text) | MR |
| EC-59 and EC-209 | the fifth slice |
| Query forms `search`, `onto`, `expand`, `group`; EC-29's indel carets; track groups and summaries (EC-103) | none here (Gate 0, 11) |
| Small looks and bugs from the third slice (disabled toolbar look, HC's, SC's, CY's, DA's) | [fixes.md](fixes.md) |
| `render/_test/frame-time.test.ts`'s flaky max bound | fixes.md (test fix) |
| QT2's `itWithin` option, DA's ANARCI oracle, OP's font and `image-rendering` checks, AY's coarser naming, CS's strand and p53 gene link, the toolbar's `ToolbarItem.title` | the author's call or `TODO.md`; not lanes here |
| The third slice's size failure | below |

**Open target question, for the author at the third slice's M3 gate (not decided here):** the minimal viewer is about
850 kB (830 kB after MV, 2026-11-19) against the 600 kB guide, and the default viewer about 1 MB; the third slice's
gate moves or confirms the targets. This slice adds the structure library and the `headers` behavior as opt-in entries
and `runtimes/mmcif` and `runtimes/bigwig` as their own rows; its sizes gate reads "not grown" against whatever the M3
gate decides. If M3 leaves the target open, M4 asks again.

## Order and parallelism

Two or three lanes at once, each strong or fast by its work (the author, 2026-11-20; loop, Lanes at once). Three
tracks, each in order, run beside each other; a lane waits only for its own track and the named dependencies:

- **Molecules** (the language and the store, one at a time): ML1 → DP → ML2 → MK → AC → MR → CU → MX → stage 1's gate.
- **Formats and data** (core's assets, runtimes, test packages): MC after ML1's node shapes (beside ML2); BW (stops at
  its `@molstar/core` surface for the author); AR; NG. Any of them beside any molecules lane.
- **Views** (render and the viewer): FX1 first (its files are render's and the stories'); SZ after FX1 and before any
  header lane (it splits the files they touch), beside the molecules track when stage 1 is not in those files; then
  HP1 → HP2 → HP3; ST1 after MK (map priority) → ST2 → ST3; HD after ST3; SV after ST3 and MK.

Then EX (after the language bump's lanes) → CN, LE beside CN; FX2–FX5 at their fix points beside whatever runs; the
novice pass after stage 4's gate; M4 → finalize.

Never two lanes on the language at once: ML1, MK (`priority`, `via`), AC, EX (the extension params) and LE (`legend`)
go through the 0.5.0 bump in turn, each added by the lane that owns it. Never two lanes on core's store or resolve
(DP, ML2, MK). Never two lanes in one file: SZ runs alone in the files it splits. The author's findings queue as
fixes.md lines for the next fix point.

## Hot files (the orchestrator's only)

- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/`;
- `packages/seq/*/package.json` and `tsconfig*.json` (the new `./structure`, `./io/structure`, `./extension` and
  `./runtimes/mmcif`, `./runtimes/bigwig` entries of seqx-core, the `./link-structure` entry of seqx-viewer, NG's dev
  dependency, `sideEffects` listings), `packages/core/package.json`, `packages/io/package.json`;
- `packages/seq/*/_spec/README.md` sections shared by lanes (Language, Context, Configuration): each lane edits only
  its own sections; the bump's owner merges the Language section;
- `.seq/`.
No new package is planned (seq-data-ports: no package; structure syncing and extensions are entries).

## Known risks

- **Units change what a system is**: paths, copies, hover, selection, follow and export all read it; ML2 lands it with
  a twin per consumer, and CU and MR rest on it.
- **The 0.5.0 bump changes every fixture once** (MX): a recording again for DT, DC, DG and p53; the records
  (`recordDrawLists`) change once, named in MX's report.
- **A public API change in `@molstar/core`** (range reads, the standard components): stopped for the author before the
  lane builds past the surface.
- **Designing for the fake viewer**: the library is plain data and minimal (translation and wiring, no port), so the new mol viewer can use it; the
  fake viewer imports nothing of Seq\* and holds its own chain table, so the mapping crosses two representations.
- **Length**: about 30 lanes, the 16 strong ones one after another; Gate 0's cut order decides what goes first.
- **The third slice's size failure**: see "The third slice's carried items, placed"; the sizes gates wait on the M3 gate's answer.
- **Two lanes on the store** (DP, ML2): ordered, not parallel, because both change what a source's nodes are keyed by.
- **The third slice's tail**: DB, DS, PZ and M3 may move the files stage 1 touches (layout heights, labels); the first
  briefs take their base after the third slice's close.

## Out of scope

Mol\* adapters and a live 3D view (the new mol viewer); Genomes and its formats (the fifth slice: layers across
sections, EC-59, EC-209, the genome ring); the version-1 rows of Gate 0's question 11 unless the author places one
here; GenBank; settings layers, a generated settings form; agent tools; a React adapter; writers beyond SVG, PNG and
archives; the compute port; a glycan's tree drawing (SNFG), a PRD group's chip, a conjugate drawn at its site
(molecules.md's standing tensions). Event composition across parts stays the application's (behaviors offered, a host
composes them).
