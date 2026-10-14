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

| Step                   | What                                                                                                                                                                                                                                                                                                                                                                        |
| ---------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Layout                 | `packages/sequences/{core,viewer,nightingale,molstar-poc,scripts,storybook}`, with `_spec/`, `todo.md` and `AGENTS.md` beside them. Seq\*'s `AGENTS.md` adapted: keep its code rules (below), drop Stories' shared-tree commit rules and its pointers into Stories                                                                                                         |
| The repository's rules | relative imports `.ts` → `.js` (Seq\* writes `.ts` because Node's type stripping runs its scripts; here the scripts run as the repository's do), exports in the `molstar-src` / `types` / `import` shape with a `lib` build, the repository's Biome formatting, the inventory, tsconfig references and `check:workspace`, `noRestrictedImports` per package (from Seq\*'s `layers.grit` and `core/_test/layers.test.ts`) |
| msgpack                | Seq\*'s copy (`core/src/model/msgpack`, `VENDORED.md`) is Mol\*'s plus local changes `STORIES-UNTOLD(M1–M5)`: M2 extension types (encode and decode), M3 JSON's rules for what JSON cannot hold, M4 characters outside the BMP, M5 untrusted lengths (an array or map claiming more than the bytes left fails, arrays grow by `push`, no `__proto__` prototype). They go into `@molstar/io`'s msgpack (`packages/io/src/common/msgpack`; additions and fixes, in the changelog); `seq-core` imports `@molstar/io`, its one dependency (the spec's "no dependencies" line says so) |
| zip                    | Seq\*'s copy (`core/src/format/zip`, `VENDORED.md`) is Mol\*'s UZIP port plus `STORIES-UNTOLD(Z1)`: the end record searched only in the last 64 KB, the central directory read first with optional `UnzipLimits` (entry count, total size), sizes and CRC from it (Finder's data-descriptor archives), each entry checked against its size and CRC-32, `_inflate` bounded. Z1 goes into `packages/core/src/util/zip` (a security fix, in the changelog); `seq-core`'s `.svsx` reader imports it. If `seq-core` should not depend on `@molstar/core`, keep the copy and note it |
| `molstar-poc`          | keeps its name, a private package for now; its v5 imports (`molstar/lib/...`) moved to the v6 packages (`model`, `io`, `query-language`, `plugin`, `mvs/runtime`). Its `NIGHTINGALE_MAX_POSITIONS`-style limits and the `molstar:*` providers (`pae`, `plddt`, `sifts_coverage`, `secondary_structure`) stay as they are                                                         |
| Tests                  | Seq\*'s tiers (§2.3): its `test` tier joins the repository's one Vitest run (`pnpm test:unit`, per checklist.md); `test:slow`, `test:browser` (Chromium through Playwright), `test:heap`, `test:browser-heap` and `test:heap-deep` become scripts of their own, never in `pnpm test` or default CI. Keep the test counts equal before and after the port (by name) |
| Storybook              | `packages/sequences/storybook` registered in `storybook/packages.ts` under a `Seq*` group; its React stories mounted by a decorator of that entry (the Storybook is `web-components-vite`); its paths into other packages' `_test` become exports (`./testing`). The smoke test (`storybook/_test/smoke*.test.ts`) opens every story headless and stays in the `test` tier (its four heaviest stories are `slow`) |
| Dependencies           | to the catalog: the Nightingale elements by major (as the AI SDK), Vitest's browser provider and Storybook's React pieces pinned exactly                                                                                                                                                                                                                                  |
| Docs                   | the root `AGENTS.md` table, `.v6/plans/checklist.md`, the changelog (msgpack M2–M5, zip Z1)                                                                                                                                                                                                                                                                                |

Mechanical steps go to Sonnet subagents at low or medium effort, each checked against the tests its change touches.

### 2.3 Tests and CI (Seq\*'s policy, carried over)

- `test`: node only, fast (about 4–5 s per package on an M4, the viewer about 13 s), deterministic: no heap guards, no
  Chromium, no wall-clock bounds (a bounded test runs without the clock here, through `itWithin`). The only tier CI runs
  by default.
- `test:slow` (the `slow` tag: fuzz at full size, the SeqQL oracles, million-row cases, CLIs as processes, PNG through
  Chrome, timing bounds, the heaviest smoke stories), `test:browser`, `test:heap`, `test:browser-heap`,
  `test:heap-deep`: opt-in, by hand, at checkpoints or in a manual or scheduled job; heap guards only on a quiet machine.
- Heap guards follow Seq\*'s template: every `*.heap.ts` has a positive control through `leak()` / `controlCaught`
  run through the file's own warm-up and `settleFrame` (the spec check fails a guard without one; browser controls plant
  garbage sized to the run). Keep that check.
- Stories pinned Node 24.19.0 for the heap tiers (the guards' numbers differ between V8 13.6 and 14.1); Mol\* runs on
  the repository's Node: re-run the heap tiers there once after the port and fix what moves, as the closing pass did.

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

## 3. Prism over Seq\*

- **The context**: Seq\*'s view context on Prism's context, at least in the viewer mode (the story mode may differ), per
  [contexts.md](../prism-ui/design/contexts.md).
- **The UI**: written in prism-ui's own framework, not Seq\*'s React UI, in `prism/ui/src/components` for now. So
  prism-ui depends on `@molstar/seq-viewer` (Seq\* is not one of the opt-in modules); the viewer mode's package is the
  alternative if every Prism app should not carry it. The default UIs' state plumbing is not worth polishing before
  this (David, 2026-10-13: the immediate renderer hosts the UI).
- **The mocks**: the sequence strip replaced by the real viewer. Start with a plain text view or a single track with a
  toggle. No viewer frame: its theme follows the canvas background.

## 4. Follow-ups

An item leaves this list in the commit that does it. Each says what triggers it. Seq\*'s own deferred features stay in
`packages/sequences/todo.md`.

| Item                                                                                                                                           | Trigger                                      |
| ---------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------- |
| Seq\*'s story support vendors Prism's tokens (`storybook/src/_stories/ui/prism`): use prism-ui's own                                           | the port's review                            |
| Seq\*'s React UI (`viewer/react`): kept for React hosts, or dropped once Prism has its own                                                     | Prism's sequence UI built                    |
| The subscription guard (`prism/ui/src/_test/subscriptions.test.ts`) covers prism-ui only; the other Prism packages want it too                 | a leak found there, or a shared check script |
| Seq\*'s heap tiers in a scheduled CI job (opt-in today)                                                                                        | a hot-path regression missed between checkpoints |
| The model catalog's prices are secondhand (ported from Stories, with sources and dates)                                                        | credits shown to people                      |
| The Inspector's mock Trace link always lands on the latest turn                                                                                | a session feeding the Inspector              |
| Whether browsers encode AAC and Opus at Gemini's 24 kHz (a runtime check falls back to WAV)                                                    | narration in the story mode                  |
