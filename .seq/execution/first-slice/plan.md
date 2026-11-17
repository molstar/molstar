# Plan: the first slice

Work id `first-slice`. Designs: `.seq/design/version-1.md` (Order of work) and the designs each stage names. The
smallest end-to-end slice that tests the riskiest decisions first; after it, the milestones in version-1.

## Done when

- The two views, **p53 tracks** and **p53's protein on its gene**, run end to end: built with the builder, validated,
  resolved, laid out, drawn on a Canvas2D with projected hover, exported to SVG, in Storybook and in tests.
- A source read by region (synthetic, then FASTA with its index) carries tables and dynamic nodes, shown by a track
  template.
- Frame time, allocation, test time and bundle size are measured, and the performance targets confirmed or moved by
  the author.
- Designs, `_spec`s, `status.md` and the features catalog say what was built; the record is written.

## Gate 0: before any lane

- [x] GATE 0 (2026-11-17, "approve both"; range reads wait for bigWig, the slice reads its FASTA whole): the author confirms the provisional choices the slice builds on (each with its recommendation):
  - meso\*'s tree core: fork into `repr` now, migrate meso\* later ([language-description](../../design/language-description.md));
  - the slice's kinds and param names (same);
  - the three result types and labels 1-based by default ([locations-and-queries](../../design/locations-and-queries.md));
  - tiles of 1,024 bins and the table columns ([protocol](../../design/protocol.md));
  - the base-library changes the slice needs: `@molstar/core/repr` (agreed); range reads in `AssetManager` (stage 6,
    FASTA by region; or the slice reads a small FASTA whole and range reads wait for bigWig).

## Stages and lanes

Lanes run in worktrees, each on its own branch, merged by the orchestrator (§commit-often, §lanes-report-not-write).
At most two heavy lanes and two light ones at once (§parallelism). Model per lane (§model-per-task): **strong** for
design-heavy or ambiguous work, **fast** for well-specified implementation.

### Stage 1: the tree core in `@molstar/core/repr` (lane A, strong)

- [x] A1 `repr/tree`: node with `custom`; values with typed columns (identity by type and bytes, cached hash); Merkle
      hashes; issues with severity; namespaced kind pattern; addresses with the form-resolver table (meso\*'s child
      selection as the default)
- [x] A2 `repr/schema`: languages as data with `column` and `table` params; validation with every issue; defaults;
      normalization on read (columns to typed arrays)
- [x] A3 `repr/builder`: the typed builder's engine
- [x] A4 `repr/codec`: the codec interface and its suite; `codecs/json`, `codecs/sexpr` (columns elided in print)
- [x] A5 `repr/_spec` (rules and reasons, what is built); `package.json` exports and `_test` nulls; a v6 changelog entry
      (`.v6/changelogs/core-repr.md`)
- [x] GATE (2026-11-17): quick lane green within budget (71 tests, about 0.2 s); review (strong) done, its blockers
      fixed in one round; merged `fbc27cde5`

### Stage 2: the language, locations, the protocol (three lanes in parallel, after stage 1)

- [x] B (fast) `seqx-core/lang`: `SeqLanguage` with the slice's kinds; the builder (`createSeqSpec`); the issue catalog
      with its test; the p53 example built, printed, encoded both ways, validated, misuses giving their codes
- [x] C (strong) `seqx-core/query`: Location (parts and rects), element keys, the three result types, the first forms
      as resolvers (`range`, `at`, `all`, `item`, `node`, `extent`, `items`), numberings `index` and `ordinal`, labels
      EC-1 to EC-8
- [x] D (strong) `seqx-core/data`: the protocol's types, the `seq` format (round trip), `FormatRuntime` for `seq`, the
      store (keys, states, versions, `KeyedState` per entry, requests ordered and bounded, cancellation, byte-bounded
      eviction) with a fake clock
- [x] GATE (2026-11-17): B, C, D merged (`77d03dfc1`, `5fa7c43d4`, `6e885495a`), each reviewed (strong) and fixed in one
      round; the merged tree 233 tests in about 0.3 s, tsc (tests too), lint, format, scaffold checks green; designs
      updated from the drift; the orchestrator's decisions await the author (gaps.md)

### Stage 2b: range sets into core (lane C2, fast, beside stage 3)

- [x] C2 `@molstar/core/data/int`: half-open range sets (normalize, union, intersect, subtract, overlap by binary
      search) moved from `seqx-core/query/location.ts` with its property tests; Location uses them; core's `_spec`,
      exports and a changelog entry (the author approved the move, 2026-11-17)

### Stage 3: resolve (lane E, strong, after B and C)

- [x] E1 the static model: systems, units, numberings; `offset` and `codon` maps in their encodings (implicit, runs);
      paths and projection with the copy rule and closure; sources bound to runtimes; annotations' access plans
- [x] E2 the resolved view: views, sections, tracks, layers, themes as values, fallback chains, each layer's data and
      static path
- [x] E3 maps' edge cases on every encoding (EC-9 to EC-32 that apply to offset and codon) with equal answers
- [x] GATE (2026-11-17): lane E merged after one review round; the merged tree 500 tests (repr, data/int, packages/seq)
      in about 0.4 s, tsc with tests, lint, format green

### Stage 4: render and models (two lanes in parallel, after E's types)

- [x] F (strong) `seqx-core/render`: the layer-type contract; `letters`, `ranges`, `values`; the tracks layout; draw
      lists (units, passes, primitives in reused buffers); hit indexes; stable lanes; the token table; the text metrics
      port with the Node table; the recording and SVG backends; `renderSvg` (asynchronous, waits to settle)
- [x] G (strong) `seqx-viewer` models: `SpecModel`, `ViewModel` (a scroll group, one follow through the codon map),
      `HoverModel` (ownership, projection), `StoreModel`; the `hover` and `navigate` behaviors; the gesture-binding table
- [x] GATE (2026-11-17): review of F and G (strong), each fixed in one round and merged; the merged tree 581 tests
      (about 0.4 s), tsc with tests, lint, format, scaffold checks green; drift folded into the designs (two items marked
      to confirm: the selection tint's order, letters' density in CSS px); nothing drawn in Storybook until stage 5

### Stage 5: the viewer and the views (lane H, fast, after F and G)

- [x] H1 the Canvas2D backend with the overlay canvas for emphasis passes
- [x] H2 `createSeqViewer`, `DefaultSeqViewerSpec`, the three ports (assets with the policy wrapper's defaults, text
      metrics, structure unset), `createSeqContext`
- [x] H3 `packages/seq/examples` (orchestrator adds it to the workspace, inventory and tsconfig): fixtures and the two
      views' builder modules, with their tests (build, validate, resolve, lay out, SVG); two draw-list records
- [x] H4 stories: `Seq* (next)/Examples/p53/{Tracks, Protein on its gene}`, canvas and SVG side by side, the printed
      spec beside them; a visual pass in light and dark
- [x] GATE (2026-11-18): integration: lane H merged after one review round; both p53 views in Storybook (light and
      dark, canvas and SVG side by side) and in tests; two bugs the visual pass found fixed (layers' buffers dropped
      contents when they grew; the mounted viewer ignored the UI theme); 678 quick tests, medium lane, typechecks, lint,
      format, `pnpm build` green; `check:workspace` red on a pre-existing `@molstar/ui` export, not this slice's

### Stage 6: by region (lane I, strong, after stage 5)

- [ ] I1 (4a) a synthetic asynchronous `seq` source with a fake clock: coverage, levels, tiles, pending vs empty,
      errors, cancellation, eviction, summaries (value pyramids); planning through the static codon map; then FASTA
      with its index (read whole: Gate 0)
- [ ] I2 (4b) dynamic nodes: a transcript with its item coordinate system and codon map; the overlay; planning through
      a dynamic map; node lifetime against eviction
- [ ] I3 (4c) the transcript's domains as a nested source, shown by a track template
- [ ] I4 a story: p53's gene region by region, transcripts arriving
- [ ] GATE: review (strong)

### Stage 7: measure (orchestrator, light)

- [ ] M1 frame time (pan, zoom, hover) at the scale fixtures the slice supports; allocation of draw and hover (opt-in
      tier); the quick lane's time per folder and in total; bundle sizes (reading a spec; headless SVG; the viewer)
- [ ] GATE: the author confirms or moves the performance targets ([performance](../../design/performance.md))

### Finalize

- [ ] FINALIZE: reconcile: designs edited where the code taught otherwise (drift notes resolved), `_spec`s with what is
      built, `status.md`, the features rows' status, the edge cases' outcomes for what was built, `TODO.md`; write the
      record in `.seq/records/first-slice/`; end `state.md` with its closing line

## Hot files (the orchestrator's only)

`pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json` (the allowlist rule for
`packages/seq` lands with its first allowed import), `storybook/packages.ts`, `vitest.config.ts`, `.v6/changelogs/`,
and every file in `.seq/`.

## Testing strategy ([testing](../../design/testing.md))

- **In the edit loop**: the touched package's quick lane (`pnpm test:unit packages/core/src/repr`,
  `pnpm test:unit packages/seq/`), each folder under 2 s.
- **At each gate**: the repository's quick lane for the touched packages, `npx tsc -b` for them, `pnpm lint`,
  `pnpm format:check`, the scaffold checks (`.seq/AGENTS.md`, Checks); at the integration gate the medium lane too.
- **Properties** at small fixed seeds: codecs, normalization, addresses, map encodings agreeing.
- **Records**: two draw-list records (one per view); the rest by properties.
- **Heavy and opt-in tiers** only in stage 7.

## Order and parallelism

```
Gate 0 → A ─┬─ B ─┐
            ├─ C ─┼─ E ─┬─ F ─┐
            └─ D ─┘     └─ G ─┴─ H → I → M → Finalize
```
