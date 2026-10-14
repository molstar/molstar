# Plan: the sequence viewer (Seq\*) into Mol\* v6, then Prism over it

Supersedes the order of work in [../prism-ui/handoff.md](../prism-ui/handoff.md) (its "How the user works" and
"Practical notes" still hold). Background:
[../prism-ui/analysis/harness-and-seqstar.md](../prism-ui/analysis/harness-and-seqstar.md).

Reconciled 2026-10-14 with Seq\*'s state after its review gate in Stories Untold (branch `seq-view-spec-proto`; the
plan, findings and fix rounds in that repository's `designs/seqstar/review/`, `2026-10-14-*.md`; the handoff in
`designs/seqstar/handoff.md`; the progress log `designs/progress/052-seqstar-move-prep.md`).

## 1. The order

1. **Port Seq\*** from `stories-untold/lib/seqstar` (§2), from the commit Stories' handoff names as the port base (the
   pushed head after the closing pass). The move is one way: from then on this repository is Seq\*'s source of truth,
   and Stories Untold's copy is frozen (Stories keeps reading it until it switches to the published packages).
2. **Review, merge** Seq\* as a PR of its own (about 225k lines with tests), after Mol\* UI v-next phase 1 (Prism and
   the AI packages, `v6-ui-next`) has merged, and before the Prism viewer is written. The code review already happened
   in Stories (eight area reviews, five gap reviews, three fix lanes each reviewed again, every must-fix and
   should-fix fixed; §2.6): this PR's review is of the port itself (layout, imports, tooling, the upstreamed io and zip
   changes), not a second code review.
3. **Prism over Seq\*** (§3), on a branch of its own.
4. **Then the agent in Prism**: the handoff's former step 2 (command metadata, the selection param, the core viewer
   commands, the checkpointer, `AgentModel` over `AgentSession`, History over fibers), decided with the user.

## 2. The port

### 2.1 What arrives

Seq\* at the port base (2026-10-14; lines with tests): `core` (77k, no dependencies), `viewer` (123k: the canvas
engine, SVG export and a CLI, DOM and React UIs, rxjs), `nightingale` (6k, the Nightingale elements, small systems
only), `molstar-poc` (5k, Mol\* v5), `scripts` (1k: the spec, header and export checks, `ported-lines`), `storybook`
(15k: the stories and the `Seq*/Docs` pages, React), `_spec/` (the spec; tables generated into `_spec/reference/`),
`todo.md`, `AGENTS.md`. Already in Mol\*'s shape (done in Stories, keep it):

- **Package names**: `@molstar/seq-core`, `seq-viewer`, `seq-nightingale`, `seq-storybook` (private), `seq-scripts`,
  `seq-molstar-poc`. The `seqstar:` system namespaces, the `seqstar-tiles/1` protocol, the `seqspec:` citations and the
  `--seq-*` tokens are data or contracts: unchanged.
- **Headers**: Mol\*'s header with the Anyfold line on every source file; `scripts/spec.ts` fails a file without it
  (`HEADER_SKIP` lists the vendored and generated ones).
- **Layout**: tests in `_test/` beside each package's `src`, stories under `storybook/src/_stories/`; every package's
  `exports` hides `./_test/*` and `./_stories/*` and its `files` leaves them out.
- **Exports**: named export lists per entry (no `export *`); `scripts/exports.ts` fails when a public signature names
  an unexported type, for every package (the `--contract` walk), and `--unused` lists names no consumer uses (a check,
  not the list's source). Tests reach internals by source path or `/testing`.
- **Names**: snake_case for every name a SeqViewSpec document holds, camelCase for the viewer's (spec, behaviors,
  config, commands, events), kebab-case for DOM events, CSS and CLI names.

### 2.2 Steps

Decided with David on 2026-10-14, after an analysis of both trees (the decisions are in the rows).

| Step                   | What |
| ---------------------- | ---- |
| Layout                 | `packages/sequences/{core,viewer,nightingale,molstar-poc,scripts,storybook}`, with `_spec/`, `todo.md` and `AGENTS.md` beside them, copied verbatim from the port base's tracked files (`git archive`; not the untracked `viewer/__debug__/`). Seq\*'s `AGENTS.md` adapted: keep its code rules (§2.4), drop Stories' shared-tree commit rules and its pointers into Stories |
| Imports                | relative imports `.ts`/`.tsx` → `.js` (static, dynamic, re-exports, `vi.mock` paths), except in `scripts/`: the private tooling package is never built, Node's type stripping runs it, and it keeps `.ts` imports |
| Packages               | exports in the `molstar-src` / `types` / `import` shape with a `lib` build; CLI `bin`s point at `lib/*.js`; the inventory, the root tsconfig's references, `check:workspace`. `check:workspace` skips `_stories` folders as it skips `_test`, so tests and stories keep their relative imports into the other Seq\* packages' tests and fixtures (no `./testing` export for now) |
| tsconfig               | per package: a composite build config (`src` only) and a `noEmit` one for `_test`, `bench` and `_stories`. Overrides of the base: `useDefineForClassFields: true` (Seq\* was built and guarded under ES2022 define semantics; Mol\*'s base keeps v5's `false`), `lib: ["ES2023"]` with DOM only for `/canvas`, `/dom`, `/react` and their tests, `noUncheckedIndexedAccess` as Seq\* has it (false; true for core's `query`). Switching the whole repository to define semantics is a PR of its own (§4) |
| React                  | `react` / `react-dom` `>=18` as optional peers (`peerDependenciesMeta`, as prism-ui: only `/react` imports React); built and tested against the catalog's React 18 and its types (only basic hooks are used) |
| Formatting and lint    | Mol\*'s Biome and Prettier settings and rules, applied in one format-only commit. Not formatted: the generated `_spec/reference/*.md` (test snapshots), the Storybook's MDX pages, JSON and binary fixtures (exact bytes; `case9-pixels.json`), and the vendored code kept as is. The layer rules stay GritQL plugins (`layers.grit`, `core/builder-imports.grit`) with their paths updated, plus Seq\*'s `noRestrictedImports` / `noRestrictedGlobals` overrides moved into the root `biome.json` |
| msgpack                | Seq\*'s copy (`core/src/model/msgpack`, `VENDORED.md`) is Mol\*'s plus `STORIES-UNTOLD(M1–M5)`: M2 extension types, M3 JSON's rules for what JSON cannot hold, M4 characters outside the BMP, M5 untrusted lengths. M2–M5 go into `@molstar/io`'s msgpack (`packages/io/src/common/msgpack`) and M4 into `@molstar/core/util/utf8`, unconditionally: BinaryCIF output loses the nil `mask` key (M3; fine, David), ASCII and BMP text encodes as before (M4). M5's checks are one comparison per header; they stay unless a benchmark on a large `.bcif` decode shows a cost (then a trusted flag). Seq\*'s `index.ts` (the typed-array extension) stays in `seq-core` |
| zip                    | Seq\*'s copy (`core/src/format/zip`) is Mol\*'s UZIP port (identical to v6's but for formatting and imports) plus `STORIES-UNTOLD(Z1)`. Z1 goes into `packages/core/src/util/zip` (a security fix; stricter for every caller: data-source, the plugin's file actions, MVS formats); the copy and its shims are deleted and `seq-core` imports `@molstar/core` and `@molstar/io` (io already depends on core) |
| `molstar-poc`          | keeps its name, a private package for now; its v5 imports (25 paths, 18 files) moved to the v6 packages through `.v6/plans/migration-map.json` and `migration-symbols.json`; its plugin built by v6's explicit composition with the basic setup and MolViewSpec. Its limits and the `molstar:*` providers (`pae`, `plddt`, `sifts_coverage`, `secondary_structure`) stay as they are |
| Tests                  | the `test` tier joins the repository's one Vitest run (`pnpm test:unit`) as a `sequences` project with `isolate: false` and an `isolated` project for files that mock modules (as `vitest.shared.ts`); the `slow` tag declared at the root and filtered out (`!slow`); browser folders (`_test/**/{canvas,dom,react}`) and `*.heap.ts` never matched. `test:slow`, `test:browser`, `test:heap`, `test:browser-heap`, `test:heap-deep` stay per-package scripts (Vitest 5, `@vitest/browser-playwright` 5); the slow tier builds `lib` first (its CLI tests spawn the built CLIs). Test names in the fast tier equal the port base's (a baseline list, §2.7) |
| Storybook              | `packages/sequences/storybook` registered in `storybook/packages.ts` with no group (its titles already start with `Seq*/`); its React stories mounted by a decorator of that entry (the Storybook is `web-components-vite`), their `Meta`/`StoryObj` types local; the MDX pages on the docs addon; `seqstarFixtureFs`'s path regex updated. Done when `build-storybook` passes and the smoke tests pass; David checks it visually |
| Dependencies           | to the catalog: the Nightingale elements by major (as the AI SDK); `@vitest/browser-playwright` pinned exactly; Seq\*'s `@types/node` 24 and React 19 pins dropped for the catalog's |
| CI                     | Node 24 in `.github/workflows` (Seq\*'s tooling and scripts use it); CI and `pnpm test` run the fast tier only |
| Docs                   | the root `AGENTS.md` table, `.v6/plans/checklist.md`, the changelogs (msgpack M2–M5, utf8 M4, zip Z1) |

### 2.3 Tests and CI (Seq\*'s policy, carried over)

- `test`: node only, fast (about 4–5 s per package on an M4, the viewer about 13 s), deterministic: no heap guards, no
  Chromium, no wall-clock bounds (a bounded test runs without the clock here, through `itWithin`). The only tier CI and
  `pnpm test` run.
- `test:slow` (the `slow` tag: fuzz at full size, the SeqQL oracles, million-row cases, CLIs as processes, PNG through
  Chrome, timing bounds, the heaviest smoke stories), `test:browser`, `test:heap`, `test:browser-heap`,
  `test:heap-deep`: opt-in, by hand during heavy work on Seq\*; heap guards only on a quiet machine. Never in CI (the
  heap tiers take about half an hour).
- Heap guards follow Seq\*'s template: every `*.heap.ts` has a positive control through `leak()` / `controlCaught`
  run through the file's own warm-up and `settleFrame` (the spec check fails a guard without one; browser controls plant
  garbage sized to the run). Keep that check.
- The port does not re-run the slow, browser or heap tiers (David, 2026-10-14): every tier passed in Stories' closing
  pass at the port base. The port is checked by the fast tier (names equal to the baseline), the type check, lint and
  `check:workspace`; the heap guards' numbers are re-measured on Mol\*'s Node the next time someone works on a hot path.

### 2.4 What to keep from Seq\*'s `AGENTS.md`

The spec is `_spec/` and says what, never how (tables generated; `seqspec:<file>#<slug>` citations checked by
`scripts/spec.ts`); the code is the truth; the performance contract (no allocation per frame, pointer event or run; no
boxed doubles across helper calls; every hot path guarded); interaction is viewer-spec behaviors; the replay rule
(`interaction.md#replays-and-carries`: a state's replay or a graph carry is never anyone's move); severities come from
the registries (`svs.*`, `seqql.*`, `seqv.*`: error, warning, info); `todo.md` holds what is not built, each with its
trigger; size budgets are recorded, not stopped for.

### 2.5 Done in Stories, not to redo

The composer and `streams.ts` splits; the guard template on every heap file (none pending); the `info` severity
(core and UI); the input bounds (zip, MessagePack, JSON and SeqQL depth, BGZF, BigWig/BigBed, 2bit); the Clustal
reader; the File formats stories (`Seq*/File formats`, recorded fixtures with a Live toggle for UniProt, PDBe and
AFDB).

### 2.6 Known open items at the port base

Seq\*'s `todo.md` is the list. Worth knowing for the PR: the conformance kit does not yet catch a hit naming an item
where nothing is drawn (LC-17's remainder); a few composer branches are uncovered by the proof corpus (LC-41); a guard
whose control reports a leak without planting one would still pass the spec check (LC-23's remainder); `y` with its own
window for dot plots (PAE's interface blocks) is the next dot-plot step.

### 2.7 Delegation

Sonnet does the port; Opus (this session) copies, formats, commits, writes the adapted `AGENTS.md` and decides; a fresh
Opus agent reviews. Rules for every agent: no commits (the coordinator reads each diff and commits it, one commit per
step); one writer per package folder; checks scoped to what the step touches (the package's fast tier through
`pnpm test:unit packages/sequences/<pkg>`, its type check, `biome lint <paths>`), never a full run while another agent
runs one, never the slow, browser or heap tiers; stop and report when the work differs from this plan.

| #   | Step                                                                                                         | Who                   | After      |
| --- | ------------------------------------------------------------------------------------------------------------ | --------------------- | ---------- |
| A0  | Baseline: the fast tier's test names per package at the port base (`vitest list`, sorted), into the scratchpad | Sonnet, low           | —          |
| A1  | Verbatim copy (`git archive 9b2742c2 lib/seqstar`), committed as is                                         | Opus                  | —          |
| A2  | Imports codemod (a script, run once), manifests, tsconfigs, inventory, references, `check:workspace`'s `_stories` skip, `pnpm install` | Sonnet, medium | A1         |
| A3  | Format-only commit: the formatting exclusions, then `pnpm format`                                           | Opus                  | A2         |
| A4  | Vitest: the root config's `sequences` projects, the tag, the per-package tier scripts; names against A0      | Sonnet, medium        | A3         |
| B1  | msgpack and zip into `@molstar/io` / `@molstar/core`: tests, the M5 decode benchmark, changelog entries      | Sonnet, medium        | A1 (parallel with A2–A4; touches only `packages/io`, `packages/core`) |
| B2  | core: imports from `@molstar/io` / `@molstar/core` (copies deleted), type check, lint, fast tier             | Sonnet, medium        | A4, B1     |
| B3  | viewer: type check (React 18 types), lint, fast tier                                                         | Sonnet, medium        | A4 (parallel with B2, B4, B5) |
| B4  | nightingale and scripts: type check, lint, fast tier; scripts' paths (`SEQSTAR`, the repository root, `HEADER_SKIP`) | Sonnet, low     | A4         |
| B5  | molstar-poc on v6: the migration map, the plugin by explicit composition with MolViewSpec, fast tier         | Sonnet, high          | A4         |
| C1  | Storybook: registration, the React decorator, the stories' types, MDX, fixture shim; `build-storybook` and the smoke tests | Sonnet, high | B2–B5 |
| C2  | Biome: GritQL paths, Seq\*'s overrides into the root `biome.json`; CI on Node 24                             | Sonnet, low           | B2–B5      |
| C3  | Docs: `packages/sequences/AGENTS.md` adapted, the root table, checklist, changelogs                           | Opus                  | C1, C2     |
| R   | Migration-correctness review, read-only: what the port changed against the port base (imports, manifests, configs, tests kept and passing, molstar-poc on v6, the io and zip changes), not Seq\*'s code again (reviewed in Stories); findings fixed by the step's agent, resumed | fresh Opus agent, high | C3 |

Then the full checks once, in this order: `pnpm lint`, `pnpm test:unit`, `npx tsc -b`, `pnpm check:workspace`,
`pnpm format:check`; then the PR.

## 3. Prism over Seq\*

The first slice is designed in [prism-sequence.md](prism-sequence.md) (2026-10-14, in the mocked viewer mode). prism-ui
is the runtime for generic Mol\* UIs; the Prism viewer is an opinionated composition on it and owns the sequence view
directly (`packages/prism/viewer/src/sequence/`, no extension mechanism, no new package): its model owns the Seq\*
contexts (by id, so several later), follows the canvas ink and manages the themes; its control, in the top dock
without a card, subscribes the two-way link to the 3D view while it is rendered; compact text and tracks; hover and
selection; Prism's own tooltip and legend.

## 4. Follow-ups

An item leaves this list in the commit that does it. Each says what triggers it. Seq\*'s own deferred features stay in
`packages/sequences/todo.md`.

| Item                                                                                                                                           | Trigger                                      |
| ---------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------- |
| Seq\*'s React UI (`viewer/react`): kept for React hosts, or dropped once Prism has its own                                                     | Prism's sequence UI built                    |
| The subscription guard (`prism/ui/src/_test/subscriptions.test.ts`) covers prism-ui only; the other Prism packages want it too                 | a leak found there, or a shared check script |
| `useDefineForClassFields: true` for the whole repository (Seq\* has it per package): an audit for re-declared fields (`declare`), then the full tests and apps | a class-field bug, or the next tooling pass |
| **Seq\*'s input as rxjs streams, reconciled with Mol\*'s** (David, 2026-10-14): Seq\*'s DOM input on Mol\*'s `InputObserver`, one binding vocabulary with Mol\*'s, public interaction streams on the viewer like Canvas3D's `interaction` (hover, click, drag, context menu with the hit, place and modifiers), configurable; behaviors the built-in reactions on the same streams. A Seq\* core change: by Opus, design first (`packages/sequences/todo.md`) | a PR of its own, after the one bringing Seq\* and the Prism prototype |
| **Tracks in Seq\*'s text (wrapped) view as lanes** (David, 2026-10-14): ranges and features (secondary structure, domains, sites) as underline lanes with tidy heights, continuous data (pLDDT, hydrophobicity) as a thin colored strip under each line, letters as now; not dot plots, pairs, trees or multi-row alignments. A Seq\* renderer change under its performance contract (allocation-free lane drawing, heap guards, tests, spec): by Opus, design and implementation, not Sonnet. Until then Prism's Text view shows letters only | a PR of its own, after the one bringing Seq\* and the Prism prototype (which only reviews the current renderer changes) |
| **Performance pass on Seq\*'s new host options** (`wrappedAlign`, `wrappedMaxWidth`, `wrappedColors`, `unobservedGray`, `focusViewport`'s `reveal`; a9648833e): precompute the unobserved mask per unit instead of scanning `observed` per position per frame, heap guards (the template, with positive controls) for `unobservedGray` and `wrappedColors`, then run the letters and layout heap guards, `test:heap`, `test:browser-heap` and `test:slow` on a quiet machine | the same later PR as the text-view lanes, after the one bringing Seq\* and the Prism prototype (David, 2026-10-14) |
| Seq\*'s heap tiers in a scheduled CI job (opt-in today)                                                                                        | a hot-path regression missed between checkpoints |
| The model catalog's prices are secondhand (ported from Stories, with sources and dates)                                                        | credits shown to people                      |
| The Inspector's mock Trace link always lands on the latest turn                                                                                | a session feeding the Inspector              |
| Whether browsers encode AAC and Opus at Gemini's 24 kHz (a runtime check falls back to WAV)                                                    | narration in the story mode                  |
