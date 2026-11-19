# Plan: the fourth slice

**A draft (2026-11-19), written while the third slice runs its last lanes** (DB, DS, PZ, its use gate's look, M3,
finalize). Gate 0 below is for the author; nothing here is decided until they answer it. The section "What the third
slice changes" is provisional until the third slice's `foldback.md` is decided, and is rewritten from it then.

Work id `fourth-slice`. Version 1's third milestone, **integration** ([version-1](../../design/version-1.md)): structure
maps and the structure showcase (on a fake 3D viewer, decided), services, bigWig and range reads, the resolved view's test consumer (Nightingale
dropped), archives. With it, what the third slice and the author carried here: the molecules design built
([molecules](../../design/molecules.md): units, copies, the 0.5.0 bump), copy ends on units (CE's parked branch), the
two drafts revised (DR) and built (structure syncing, track header parts), the seq data ports design, and `ext:` layer
kinds if Gate 0 places them here.

## Where this is heading

| Slice | Milestone | In one line |
|---|---|---|
| first (closed) | the end-to-end slice | the tree core, the language, resolve, draw lists, Canvas2D and SVG, the viewer's models, data by region, measured |
| second (closed) | 1, core editing and addressing | patches, history, previews, addresses, refs and links, selection, security bounds, MessagePack, navigation |
| third (closing) | 2, figures | the axis, text, ring and block layouts, alignments, dot plots, arcs and logos, themes and legends, stories, navigators, PNG and the CLI, the datasets |
| **fourth (this plan)** | 3, integration | molecules and units, structure maps and syncing on a fake 3D viewer, track header parts, seq data ports, services, bigWig with range reads, the resolved view's test consumer, archives, `ext:` kinds |
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

Provisional, from the third slice's `state.md`; replaced by its decided fold-back.
- **One lane at a time by default** (the author, 2026-11-19: "slow down to just 1 lane", usage): one lane runs; its
  review and fix round run before the next starts. When the author allows more ("you can also run another 1-2 lanes"),
  the pairs under Order and parallelism run together, never two heavy test runs at once. Why: usage limits paused the
  third slice once, mid-lane.
- **Designs agreed before their lanes**: the third slice's stages 1–4 ran unattended on agreed designs; its later
  items waited on drafts. Here every design is drafted, self-checked and decided in stage 0, before a build lane.
- **The author's findings arrive mid-slice and become lanes** (the third slice: about twenty lanes, CH to DB, from use;
  counted bug 6, missing design 7, look 6 by the end). Each queues as a fast lane with its own brief between the
  strong lanes; the plan is edited, not bypassed.
- **A lane whose fix needs an unbuilt design is parked, unmerged, with its review kept** (CE), not merged with known
  high findings.
- **A lane's public API is accepted by the author before merge** (MV, DB): a brief names the public surface it expects;
  anything past it stops the lane.
- **Read-only fast lanes prepare the author's questions and the close** (SR, EO, DA, GR): cheap, kept.
- **Use gates do not block the loop**: the author's use comes when they are back, its findings fixed after; the novice
  pass runs before the close.

## Done when

- **Molecules**: `molecule`, `complex`, `unit`, `modified`, `compounds` in language 0.5.0 (one bump, unreleased); units
  and atoms systems, copy groups derived, the atoms map at a path's ends, the modified index, the query forms; the
  mmCIF runtime (`runtimes/mmcif`, opt-in) with its conformance suite on 1EHZ, 1FIN, 4HHB and a glycoprotein; a
  section of a complex drawn with `form`, the letters' mark and code, hover; DT, DC and DG on the new kinds (changed
  once), the Fc glycan and ADC views ([molecules](../../design/molecules.md), Consequences).
- **Copy ends on units**: a pair's end on split copies by unit, every finding of CE's review a test (hemoglobin split
  keeps its contacts; copies counted as the resolver's).
- **Structure maps**: the `structure` map kind kept (decided), `segments` built, a way to choose between two maps
  joining one pair of systems (DR's gap, as decided); EC-22 to EC-27 on every encoding.
- **Structure syncing** (the revised draft, agreed): hover, selection and focus both ways, "color by this track", a
  structure's chains mapped onto sequences, on the fake structure viewer; p53 with 2OCJ and the Complex with 1FIN; the
  stories with headless twins on a fake clock (one with latency).
- **Track header parts** (the revised draft, agreed): entries in fixed slots drawn on the canvas and a host's DOM parts
  over the same geometry, the hidden-tracks chip, the menu behind its port, the secondary click a behavior param, the
  cursor on headers; `Features/Track headers`.
- **Seq data ports**: designed (a fake host's adapter reading Seq\*'s nodes from its in-memory representation, and
  annotations given back as the host's), the package layering decided; built as far as Gate 0 says.
- **Data**: bigWig read by region and zoom level through range reads of assets (p53's conservation signal over the
  gene's region); services (UniProt, PDBe SIFTS, AlphaFold DB) composed in Storybook from recorded fixtures with a live
  switch; archives read and written over the assets port with EC-91's bounds; the predictions-onto-a-structure and
  addressing showcases ([examples](../../design/examples.md), Showcases).
- **The resolved view's second consumer**: a test-only plain-SVG consumer draws a subset from the resolved view and
  the store alone (EC-183); Nightingale dropped (decided, the author, 2026-11-19).
- **`ext:` layer kinds** (if Gate 0 says here): an extension's layer named in a spec, drawn, its fallback when absent
  (EC-235), registered atomically through a descriptor.
- Measured: a hover message for 2OCJ, `maxUnits`-scale hover, a bigWig region decoded, a header pointer move over 200
  tracks, the bundles (the mmCIF and bigWig runtimes as their own rows, the minimal viewer not grown).

## Gate 0 (the author, before any build lane)

Decided (the author, 2026-11-20: "Accept all"): every question as recommended; the lane default one strong and one
or two fast at once (not one at a time).

1. **When it starts.** The build lanes after the third slice closes (finalize, fold-back, record); DR and the stage 0
   designs before, since they write only this folder. DR copies the two drafts into `fourth-slice/drafts/` and revises
   them there (the third slice's folder is closed after its record). Recommend **yes**.
2. **Scope under the one-lane limit.** About 25 lanes at one at a time is long. Keep milestone 3 and the carried items,
   and if usage runs short cut in this order: custom navigators (CN), legends' surface (LE), the resolved view's test
   consumer (NG), archives' writer (AR's write half). Recommend **this order**.
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
    lane UT; EC-99 (ligand chips flowing in text) to MR; EC-59 (two grids on one block pair) and EC-209 (ribbons
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

Each lane: a brief in `briefs/` (named below), its `_spec` and Configuration table kept, a v6 changelog entry (in its
report), an independent review and one fix round, the look in a visible pane in light and dark before the merge.
Strong: the strongest model, high effort (design, model, ambiguous). Fast: the fast model (well specified, stories,
tests, read-only).

### Stage 0: designs and carried tests (before the build; the author decides each design)

- [x] DR (done 2026-11-20: structure syncing (as a library) and track header parts promoted; strong; `briefs/dr-drafts-revision.md`, carried, retargeted to this folder) the two drafts revised by their
      "Decided by the author" sections: structure syncing without MVS types (the port on Seq\*'s locations, the host's
      adapter translating; the `structure` kind kept; the MVS writer moves out of Seq\* and out of the seq data ports),
      and a design for choosing between two maps of one pair of systems (EC-22 to EC-27) as one new gap; track header
      parts with canvas slots and DOM parts over the geometry, the secondary click a param. Self-checked; agreed into
      `design/` with their Design pages.
- [x] SP (strong, design only, after DR; `briefs/sp-seq-data-ports.md`) seq data ports: Seq\*'s nodes read from a
      host's in-memory representation (sequences, entities, components, units) and annotations given back as the
      host's, as registry entries by representation kind (ports.md's rule); the layering against `@molstar/io` and
      `@molstar/model` (a package of its own or not); whether MC's converter is shared; a fake host's adapter.
      Reads: molecules.md (gap 11, Tensions), the revised structure syncing. **Done 2026-11-20**: agreed as
      [seq-data-ports](../../design/seq-data-ports.md) (gaps decided: `FormatRuntime` kept, one adapter a format with
      `withAdapter`, no package, `source.ref`, `changes` and `reattach`, `values.field`; MC's converter shared).
- [x] XD (strong, design only; `briefs/xd-extension-drafts.md`) the extension openings, one draft: the extension
      descriptor and `ext:` layer kinds with their fallback (extensibility.md, EC-235; GR's research); a navigator's
      picture as a registry entry with the window public (`render/window.ts`, `NavigatorModel`); legends' library
      surface (placement, hiding, a compact form; editing stays the application's). Each part cut by Gate 0's answer 2.
      **Done 2026-11-20**: agreed as [extension-openings](../../design/extension-openings.md) (gaps 1, 3–8 yes; 9 no,
      no extension loading in the CLI; 2 changed: no `custom.<id>`, plain params on an extension's own type and
      `<id>:<name>` keys on others' nodes).
- [ ] UT (fast, beside the designs; `briefs/ut-untested-ecs.md`) tests for the built, untested edge cases: EC-115 (a
      stack's zoom about a point), EC-116 (a row's text for all its columns or none), EC-135, EC-195; and the third
      slice's drift left open at its close, if any.
- [ ] GATE: the designs agreed (each with its self-check and its matrix cells), UT merged.

### Stage 1: molecules, maps and copies (strong, in order; [molecules](../../design/molecules.md), The fourth slice's lanes)

- [ ] ML1 (strong, first) the language: the kinds, params, checks and builder; the `seq` format's tables; the 0.5.0
      bump opened (`briefs/ml1-language.md`). The core's standard component per letter (a base-library row, decided).
- [ ] ML2 (strong, after ML1) the model: unit and atoms systems, the atoms map at a path's ends, held pairs with per-row
      units, copy groups, the modified index, `maxUnits`, the query forms (`briefs/ml2-model.md`).
- [ ] MC (strong, after ML1 and SP) the mmCIF converter and runtime as SP decides: the allowlist, the tables read, the
      conformance suite, the fixtures re-recorded (`briefs/mc-mmcif.md`).
- [ ] MK (strong, after ML2) structure maps: `segments`, the `structure` kind with per-residue status, the map choice
      as DR's gap decides, a numbering's `offset` and `segments` methods (NB's leftovers); EC-22 to EC-27 on every
      encoding (`briefs/mk-structure-maps.md`).
- [ ] AC (strong, after MK; followups: "Yes, slice 4 (0.5.0)") pairs, ranges and maps targeting an alignment's columns:
      the tRNA's `SS_cons` arcs over the alignment, a chain mapped onto its closest seed row (`briefs/ac-alignment-targets.md`).
- [ ] MR (strong, after ML2 and MC) drawing and interaction: a section of a complex, `form` through the chips, the
      letters' mark and code, the pairs layer's per-row ends, hover, `links.hover`, `block.form`, ligand chips flowing
      in text (EC-99), the Design page (`briefs/mr-molecules-render.md`).
- [ ] CU (strong, after MR) copy ends on units: CE's branch and review as the cases (Gate 0, 8); grids-and-rows.md's
      "end on the first copy" replaced (`briefs/cu-copy-ends-units.md`).
- [ ] MX (strong, after CU) the examples changed once: DT, DC, DG on the new kinds; the Fc glycan and ADC views; 2OCJ
      recorded for p53 (its units and the fake viewer's chain table from one mmCIF); stories, twins, drawing against
      the data (`briefs/mx-examples.md`).
- [ ] GATE: merged, reviewed, checks green, sizes (`runtimes/mmcif` its own row), gaps, matrix cells (molecules.md's).

### Stage 2: structure syncing and track headers (after stage 1; the revised drafts' lanes)

- [ ] HP1 (strong) core's header line, slots, hits and chips; `TrackLayout.header` (`briefs/hp1-header-line.md`).
- [ ] ST1 (strong, after MK) the translation library, `@molstar/seqx-core/structure` (a new entry): `structureFacts`,
      `unitRows`, `toUnits` and `fromUnits` (the size gate, `seq.structure.*` issues), `colorRuns`; core's
      `SeqPorts.structure` stub removed (`briefs/st1-structure-core.md`).
- [ ] HP2 (strong, after HP1) the viewer: the `headerActions` registry kind, built-in entries, the `headers` behavior
      (reveal, clicks, the menu and its port, long press, the secondary click param), DOM parts over `viewer.headers`,
      the cursor on headers, words (`briefs/hp2-header-behavior.md`).
- [ ] ST2 (strong, after ST1; smaller: no behavior, no port, nothing in the composition) the viewer's API gaps the
      design lists, as the author decides them (`HoverModel.setLocation`, `revealLocation`, `view.lastJump`, a retire
      signal; each public surface accepted before merge), then `@molstar/seqx-viewer/link-structure`: `linkStructure`
      with echoes, coalescing, clearing on either end, focus with the jump's duration (`briefs/st2-structure-link.md`).
- [ ] ST3 (fast, after ST2 and HP2) the fake structure viewer (`_stories/support/`, importing nothing of Seq\*) with its
      host's adapter (MVS translation on `structureFacts` and `unitRows`, its `side` for `linkStructure`); the stories: hover, selection, focus, chain maps, "color by this track" as a
      host's header entry; twins; the Design page (`briefs/st3-fake-viewer.md`).
- [ ] HP3 (fast, after ST3) `Features/Track headers` and its twins, the Design page (`briefs/hp3-header-stories.md`).
- [ ] USE GATE (short, not blocking): the author tries the structure and header stories; findings counted by kind.
- [ ] GATE: merged, reviewed, checks green, sizes (the minimal viewer), gaps, matrix cells (both drafts').

### Stage 3: data and references (after stage 0; independent of stages 1–2 but for SV's SIFTS, after MK)

- [ ] BW (strong) range reads in `@molstar/core` (the surface accepted first, Gate 0, 4), the bigWig runtime with its
      levels and synthetic twin, p53's conservation signal over the gene's region (`briefs/bw-bigwig.md`).
- [ ] SV (fast, after MK and ST3) services in Storybook: UniProt, PDBe SIFTS, AlphaFold DB as runtimes from recorded
      fixtures with a live switch; the predictions-onto-a-structure and addressing showcases on the fake viewer; a
      user's own file (Gate 0, 12) (`briefs/sv-services.md`).
- [ ] AR (fast) archives over the assets port: read and write zip with `@molstar/core/util/zip`, entries resolved as
      assets, EC-91's bounds; every dataset saved and opened (`briefs/ar-archives.md`).
- [ ] NG (fast) the resolved view's second consumer (the author, 2026-11-19: Nightingale dropped, "our own test
      consumer only"): a test-only plain-SVG consumer of letters, ranges, values and variants from the resolved view and
      the store alone; the field test (EC-183) (`briefs/ng-test-consumer.md`).
- [ ] GATE: merged, reviewed, checks green, `pnpm build` and `check:workspace` (core changed), sizes, gaps.

### Stage 4: extension points (after XD's design and HP2's registry kind)

- [ ] EX (strong) the extension descriptor and `ext:` layer kinds: named in a spec, registered atomically with undo,
      the fallback drawn and said when absent (EC-235); an example extension layer in `Features/Extensibility`
      (`briefs/ex-ext-kinds.md`).
- [ ] CN (strong, after EX; cut first) a navigator's picture as a registry entry; the window, gestures and keys public
      (`briefs/cn-custom-navigators.md`).
- [ ] LE (fast, after EX; cut second) legends' library surface as XD decides (`briefs/le-legends.md`).
- [ ] GATE: merged, reviewed, checks green, sizes, gaps.

### Stage 5: use, measure, finalize

- [ ] NOVICE PASS (fast) a task script per new story (structure syncing both ways, headers, services, archives).
- [ ] USE GATE (not blocking; at finalize a look at what changed): findings fixed or placed, counted by kind.
- [ ] M4 the heavy lane and the opt-in tiers; the browser at dpr 2 (a 2OCJ hover round trip, a complex of 200 copies,
      a bigWig pan); the bundles; the targets confirmed or moved by the author.
- [ ] FINALIZE: reconcile (designs, `_spec`s, status, features, edge cases' outcomes, `TODO.md`, drift emptied); the
      fold-back; the record; the fifth slice's plan (Genomes, on `ext:` kinds).

## Order and parallelism

One lane at a time: DR → SP → XD (UT between) → stage 0's gate → ML1 → ML2 → MC → MK → AC → MR → CU → MX → gate →
HP1 → ST1 → HP2 → ST2 → ST3 → HP3 → use gate → gate → BW → SV → AR → NG → gate → EX → CN → LE → gate → novice pass →
M4 → finalize. The author's findings queue as fast lanes between them.

If the author allows two (one strong, one fast, as "1+1"):
- UT beside DR, SP or XD;
- MC beside ML2 (molecules.md: it needs only the node shapes);
- BW beside any stage 1–2 strong lane (its files are core's assets and a new runtime);
- HP1 beside ML2 or MK (the header line against the language and resolve);
- AR, NG and SV (after ST3) beside any strong lane.
Never two lanes on the language at once (ML1, AC and MK's params go through the bump in turn).

## Hot files (the orchestrator's only)

- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/`;
- `packages/seq/*/package.json` and `tsconfig*.json` (MC's and BW's opt-in entries, NG's dependency),
  `packages/core/package.json`, `packages/io/package.json`;
- a new package's inventory entry (if SP decides one);
- `.seq/`.

## Known risks

- **Units change what a system is**: paths, copies, hover, selection, follow and export all read it; ML2 lands it with
  a twin per consumer, and CU and MR rest on it.
- **The 0.5.0 bump changes every fixture once** (MX): a recording again for DT, DC, DG and p53; the records
  (`recordDrawLists`) change once, named in MX's report.
- **A public API change in `@molstar/core`** (range reads, the standard components): stopped for the author before the
  lane builds past the surface.
- **Designing for the fake viewer**: the port is plain data and minimal, so the new mol viewer can implement it; the
  fake viewer imports nothing of Seq\* and holds its own chain table, so the mapping crosses two representations.
- **Length**: about 25 lanes at one a time; Gate 0's cut order decides what goes first.
- **The third slice's tail**: DB, DS, PZ and M3 may move the files stage 1 touches (layout heights, labels); the first
  briefs take their base after the third slice's close.

## Out of scope

Mol\* adapters and a live 3D view (the new mol viewer); Genomes and its formats (the fifth slice: layers across
sections, EC-59, EC-209, the genome ring); the version-1 rows of Gate 0's question 11 unless the author places one
here; GenBank; settings layers, a generated settings form; agent tools; a React adapter; writers beyond SVG, PNG and
archives; the compute port; a glycan's tree drawing (SNFG), a PRD group's chip, a conjugate drawn at its site
(molecules.md's standing tensions). Event composition across parts stays the application's (behaviors offered, a host
composes them).
