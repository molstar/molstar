# Seq\* (packages/sequences)

The entry point for work in the Seq\* packages. Seq\* is SeqViewSpec, a format for sequence views built on one
coordinate model (positions from a single residue to a whole genome, and the maps between them: a protein onto its
gene's codons across introns, an alignment's rows onto its columns, a design onto its wild type, a chain onto its
structure; data placed on any of these shows on the others), with its data model, SeqQL, its query language, and the
builder. The Seq\* viewer renders it, and so does Nightingale. It imports nothing of Mol\*'s apps or UI; Mol\*'s
packages it may use as any package here does (the restriction to a few was the separate repository's, with its vendored
Mol\*, lifted 2026-10-14), except in the core's leaves, `query/` and `links/`. The Mol\* adapter, which reads the
structure packages, is Mol\*'s: it lives in the Prism viewer (`@molstar/prism-viewer/sequence/molstar`), and the
Storybook's Mol\* demo imports it from there. This file, the repository's [AGENTS.md](../../AGENTS.md) and the spec are
all it needs.

- **The spec** is `_spec/` (start with [_spec/README.md](_spec/README.md): what Seq\* is, the packages, the principles,
  the glossary). It says what, never how: contracts and meaning. Tables (nodes, types, layer types, themes, issue codes,
  SeqQL's symbols, commands, behaviors, config) are generated from the code into `_spec/reference/`; never edit those by
  hand.
- **The code is the truth.** The spec holds the contracts and what they mean; where the two disagree, fix whichever is
  wrong (often the spec), in the same commit.
- **Packages:** `core` (platform-free; of Mol\* only its MessagePack, zip, synchronous task context, alignment and value
  palettes) and `viewer` (with the interactions extension), plus the adapter `nightingale` and the private `storybook`
  (the table in `_spec/README.md`). Each folder of a package is an entry: inside a folder import relatively; across
  folders name the entry (`@molstar/seq-core/model`, the package importing itself), never a path into another folder.
  `core/_test/layers.test.ts` holds the core's folder graph. `scripts/` (`@molstar/seq-scripts`) holds Seq*'s own checks
  and measures, and each package's `_test/fixtures/` its recorded data.
- **Headers:** every source file (`.ts`, `.tsx`, `.css`, in `src`, `_test`, `_stories` and scripts) starts with Mol\*'s
  header, as Mol\* v6's Prism does (the spec check fails on a file without it); vendored, copied and generated files
  keep their origin's header (their `VENDORED.md` says so; `HEADER_SKIP` in `scripts/spec.ts`):
  ```
  /**
   * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
   *
   * Developed in collaboration with Anyfold <https://anyfold.ai/>
   *
   * @author David Sehnal <david.sehnal@gmail.com>
   */
  ```
- **What is not built** is `todo.md`: work in progress, next, deferred features (each with its trigger), follow-ups and
  known costs. An item leaves it in the commit that builds it; a deferral comes in with its trigger.
- **Phase ids** in comments (`LY1`, `QL2b`, `CS2`, …) name the development phases; `_spec/dev-phases.md` says what each
  was.

## Tests: run the tier the change needs

The repository's test lanes (quick, medium, heavy: the root [AGENTS.md](../../AGENTS.md#tests-three-lanes), Tests) hold
Seq\*'s Node tests; this section adds Seq\*'s own tiers. The quick lane runs in the repository's one Vitest run
(`vitest.config.ts`, the `sequences` projects): about 3 s for core, 10 s for the viewer (M4, 2026-10-14). The other
tiers are each package's own scripts (`vitest.shared.ts` and the package's `vitest.*.config.ts`):

| Command                                                                  | What                                                                                                                                                                                                                                                                                                                                     | Run it                                                                                                                                         |
| ------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| `pnpm test:unit packages/sequences/<pkg>`                                | the quick lane: node unit tests, fast and deterministic, no wall-clock timing, no Chrome, no heap guards                                                                                                                                                                                                                                 | on every edit loop, in the packages you touched; a folder or file narrows it (`pnpm test:unit packages/sequences/core/_test/format`)           |
| `pnpm -C packages/sequences/<pkg> test:perf` (core, viewer, nightingale) | the quick allocation tier: every hot path's bytes a step, read exactly (`v8.GCProfiler`), against the package's `_test/perf-baseline.json` (fails over it by 4 bytes a step; `--update` records an intended move); warm-up by count under lowered tier-up thresholds, so seconds (the viewer's about 1 min), any load, files in parallel | first in the edit loop when you touch a hot path (draws, input, streams, query kernels), beside the quick lane; pass the touched `*.quick.ts`  |
| `pnpm -C packages/sequences/<pkg> test:heavy` (core, viewer, storybook)  | the heavy lane, node tests tagged `heavy`: wall-clock bounds and ratios (`itWithin` and the scaling tests), the storms, fuzz and oracle runs at full size, a million rows and large fixtures, Chrome from Node (PNG export), the CLIs in a process (built first: `tsc -b`), the four largest smoke stories                               | when you change what such a test covers, and at checkpoints                                                                                    |
| `test:browser`                                                           | Chromium: pixels, UI, surface (the viewer's `_test/**/{canvas,dom,react}`; Nightingale's `mount`, `_test/dom`)                                                                                                                                                                                                                           | when you change drawing, the DOM or the UIs; pass the affected files                                                                           |
| `test:heap`, `test:browser-heap`                                         | the deep guards: heap growth over n and 4n steps under V8's default tuning, warmed by count until optimized (`settleFrame`), each with its positive control                                                                                                                                                                              | at checkpoints and before a release, and to confirm a quick-tier finding; pass the touched guard files (the viewer's whole tier: about 25 min) |
| `test:heap-deep` (viewer)                                                | the deep guards as before the 2026-10-12 trim: full steps (catch about 22 bytes a step, the trimmed run about 43) and the dropped cases                                                                                                                                                                                                  | before a release, or when chasing a small leak the quick tier cannot place; never in an edit loop                                              |
| `test:full`                                                              | all of the above                                                                                                                                                                                                                                                                                                                         | at reviews and checkpoints only, one package at a time                                                                                         |

**CI policy** (David, 2026-10-14; the lanes, 2026-10-20): every push runs the quick lane (and the medium lane, which has
no Seq\* tests); the heavy lane runs nightly and by hand (`tests-heavy.yml`, built first, with Chrome): a wall-clock
failure there is re-run alone before it is diagnosed. `test:perf`, `test:browser`, `test:heap`, `test:browser-heap` and
`test:heap-deep` are opt-in, by hand during heavy work on Seq\*, never in CI (the heap tiers take about half an hour).
The root config never matches the browser folders, `*.heap.ts`, `*.quick.ts` or the benches. Two quick-lane tests run a
CLI and skip without a build (`lib`). The Mol\* adapter's link guards (Prism's,
`packages/prism/viewer/src/sequence/molstar/_test`) run on this harness: `pnpm -C packages/prism/viewer test:perf` and
`test:heap`.

- **The `heavy` tag** (vitest's test tags, `vitest.shared.ts`, the root's lanes): a test slow by nature carries
  `{ tags: ["heavy"] }` (a `describe` can carry it for all its tests) and keeps a quick representative untagged (the
  fuzz at 100 mutations, the oracles at a fifth, a storm of 50, 20,000 rows for a million, most smoke stories). A
  wall-clock bound goes in `itWithin` (`core/_test/model/timing.ts`): the test runs in the quick lane without the bound
  and again in the heavy lane with it. A loop of thousands of checks records its disagreements (`mismatches`,
  `core/_test/model/mismatches.ts`) and checks once: an `expect` costs microseconds.
- **Files share a worker's modules** (`isolate: false`): a file that mocks or resets modules (`vi.mock`, `vi.doMock`,
  `vi.resetModules`) is found by reading it and runs isolated, in the `isolated` project; module state a file changes (a
  limit, a config default) it restores in `afterEach`, or the next file in the worker sees it.

- **Every hot path has a quick case; the deep guards run at checkpoints** (David, 2026-10-14). A new hot path gets its
  quick case first (`quickCase`, `core/_test/model/quick.ts`: the step, then its positive control `leaking(step)` on the
  same state, built in; `control: false` only with a `why`). Both tiers warm up by count, so neither depends on the
  load: run them beside other work. If a deep guard fails, re-run it alone once before you diagnose it.
- **The quick tier's flags** (`PERF_FLAGS`, `vitest.shared.ts`): a 1 MB young generation (the measurement accounts every
  collection: Σ(heap before a collection − the previous after) + the final delta, the code and trusted spaces apart),
  `--single-threaded` (V8 compiles on the main thread: tier-up falls on the same step under any load) and lower tier-up
  thresholds (`--invocation-count-for-turbofan=300 --minimum-invocations-after-ic-update=50`). Why the thresholds: under
  V8's defaults a step stays in its middle tier (Maglev) for 15,000 to 60,000 steps, allocating there (the wrapped
  lanes' draw 130 bytes a frame, a hover 9 to 26); with the flags most steps are optimized in a few thousand and the
  steady state is the defaults'. The warm-up ends only on two quiet windows in which V8 compiled nothing (a steady read
  in a lower tier, or a lull between tiers, never ends it), so a case needs no tuned warm-up; a case that allocates by
  design never reads quiet and sets a short `maxWarmSteps`. A case over its baseline by 4 bytes a step fails (a tenth of
  a baseline above 40: a run of kilobytes reads a few percent apart). Exact to the byte: `leak` reads 88.0, a 32-byte
  object every 64 steps 0.5. Call `hitTestAt` (the point in the record), not `hitTest(section, x, y)`, from a case's
  step: a fractional coordinate passed from test code to a call V8 does not inline is boxed, 16 bytes, the test's (three
  findings of the roll-out were this). The quick flags can change what V8 optimizes: the dot plot's drag allocated 18
  bytes a move under the defaults and none under them (2026-10-15), so a hot path's deep guard still runs at
  checkpoints.
- In a per-frame path, read enum-like object members (`SkinToken.x`, `PathVerb.x`) once into module constants: under the
  test runner's module wrapper each read allocated (about 41 bytes a frame), a false "allocates only under vitest"
  result.
- A guard's cases each load the viewer afresh (`vi.resetModules`), for independence. Cases do not "share V8's feedback"
  so that one allocates only after another (the old reading, 2026-10-14): under V8's default tuning a step climbs to
  TurboFan over 15,000 to 60,000 steps, allocating in its middle tier until then, longer after another case's draws (the
  wrapped lanes' hover after their draw: 16 bytes a step for 13,000 steps, then nothing, in either order), and the timed
  warm-up (`FRAME_WARM_MS`, 3 s) ended inside that climb whenever the machine was loaded. The deep guards now warm by
  count until optimized: `FRAME_WARM` steps, then `settleFrame` (exact windows of `FRAME_N` reading under 2 bytes a step
  and compiling nothing over 6,000 steps, at most `FRAME_WARM_CAP`), or `frameScaling(step, n)` for both and the
  measurement, `cheapScaling` for a step of tens of microseconds (60,000 steps first: its climb has quiet spans of
  thousands of steps); `viewer/_test/heap.ts`. A step that allocates once optimized never settles, runs the cap and is
  caught. A positive control allocates through `leak` (garbage, the same bytes every step), never onto a kept list: a
  list's reallocation inside the shorter measured run cancels the longer's growth, and the control goes blind.
- **The guard template:** every `*.heap.ts` runs its positive control as its guards run their step:
  `expect(controlCaught(() => <the file's own warm-up, settle and measurement> of leaking(step))).toBe(true)`
  (`core/_test/model/heap.ts`, re-exported by `viewer/_test/heap.ts`; `canvas/heap-browser.ts` in the browser).
  `controlCaught` counts a run the measurement refuses because a collection ran during it as caught (in the browser that
  is how a leak often shows). The spec check fails on a guard file without it unless `scripts/guards-pending.txt` lists
  it; converting a file removes its line. Every Node `x.heap.ts` has its quick twin `x.quick.ts` beside it, the two
  sharing their cases (a `x-cases.ts` module or the guard's fixture: no copies), and every `*.quick.ts` calls
  `quickCase` (a `control: false` with its `why`): the spec check fails otherwise. The browser's guards (`canvas/`,
  `dom/`, `react/`) stay deep-only: Chrome has no `v8.GCProfiler`.
- Pixel hashes are exact and the machine's: `canvas/case9-pixels.test.ts` keeps a record per user agent
  (`case9-pixels.json`) and fails on one with none, naming it; record a setup's hashes with `VITE_SEQ_PIXELS=record`
  (and `VITE_SEQ_PIXELS_MACHINE`) only after checking its render, and again after an intended look change.
- Run the guards on the repository's Node (`mise.toml`, 24): another Node's V8 tiers up and inlines differently (a guard
  can pass on one and fail on the other). The guards' numbers were taken on Node 24.19 in Seq\*'s former repository.
- Before every commit: `pnpm lint` and `pnpm format:check` (the repository's), the package's `typecheck` script (its
  build and its noEmit configs), the export guard (`node packages/sequences/scripts/exports.ts`, with `--contract`) and
  the spec check (`node packages/sequences/scripts/spec.ts check`: every `seqspec:` citation names a heading, the spec's
  identifiers exist in the code, headers, issue codes, guard templates). The scripts run on Node's type stripping and
  keep `.ts` imports.

## Size, specs and stories

- **Size:** phases carry size budgets in production lines
  (`node packages/sequences/scripts/ported-lines.ts --base <commit>`). Record an overrun and what took the room; don't
  stop for size (David, 2026-10-12: "the budget overflow is ok"). Do stop and report when the work turns out different
  from the plan.
- **Spec:** the spec says what, never how. A change to a contract updates `_spec/` in the same commit; a change to a
  table updates its reference (`pnpm test:unit packages/sequences/<pkg> -u`: the reference tests write
  `_spec/reference/`, which no formatter touches). Cite `seqspec:<file>#<slug>` (the slug is GitHub's anchor of the
  heading, e.g. `seqspec:model#closure`) once per file or per contract, where code implements it; never section numbers,
  never a generated reference, never for how something is built. Nothing in Seq\* points outside `packages/sequences`
  for design or history. The code decides what the docs say.
- **Stories and docs:** every major feature a user sees gets a Storybook story, or extends a case story, checked in the
  toolbar's light, dark and blue themes. Evals and internal tooling don't need one. Feature work also updates the
  `Seq*/Docs` pages it touches (Features, Concepts, Extending, Streaming, Color schemes, Query language). A new story
  brings its entry in `storybook/src/_stories/support/seqstar-smoke.ts` (the documents it opens, the patches its buttons
  send, the warnings it shows on purpose): `storybook/_test/smoke.test.ts` and `smoke-more.test.ts` (half each,
  `smoke-stories.ts`) open each on a headless viewer in the quick lane (the four largest, `HEAVY_STORIES`, in
  `test:heavy`) and fail on an error, an unlisted warning, an empty visible section, or a story export with no entry. A
  story's builder lives in `storybook/src/_stories/support/` (not in the `.stories.tsx`), so the test can import it
  without React.
- **Where the Storybook content lives:** `storybook` (`@molstar/seq-storybook`, private), all under `src/_stories/`
  (exports and `files` hide `_stories` and `_test` in every package): the stories, `support`, `docs` (the `Seq*/Docs`
  pages), `ui` (its own small components on prism-ui's own `--msp-*` tokens, `@molstar/prism-ui/styles/base` adopted by
  `support/seqstar-styles.ts`), `_test` (the smoke test and the stories' search and patch tests). It is registered in
  the repository's Storybook (`storybook/packages.ts`, no group: its titles start with `Seq*/`), whose framework is
  web-components: `src/preview.ts` (registered as the entry's `preview`) is the decorator that mounts each React story
  with `react-dom/client` and the page tokens, the stories' `Meta` and `StoryObj` are the local ones
  (`support/seqstar-story.ts`), and `readPromise` (`support/seqstar-read.ts`) is React 19's `use` on the catalog's
  React 18. `mise run storybook` serves it on :6006.
- **The stories' look is Mol\* Prism's** (David, 2026-10-13): every story is a `PrismPage` (title, caption, cards:
  `seqstar-prism.tsx`), themed by the Storybook toolbar (one viewer per story, as Mol\*'s Storybook); the viewer draws
  on a `PrismCanvas` (`ViewerCard`'s default), pinned to a theme only where comparing themes is the point (the color
  schemes pages). Toggles are `Choices`; spacing, type and colors are the `--msp-*` tokens, no borders (tone and the
  inline block's hairline).

## Design

- **Domain-independent core:** biology meaning lives in tables, presets and app-level code, never in core types.
  Complement is a named letter map, and an exon is a range.
- **Performance is a contract:** no allocation on hot paths (per frame, per pointer event, per run), and no boxed
  doubles across helper calls. Every hot path has a quick case (its positive control built in); the deep guards run at
  checkpoints. The public input streams (`ctx.interaction`) may allocate a payload per event, and build nothing while a
  stream has no subscriber.
- **snake_case is the format's** (David, 2026-10-13): every name a SeqViewSpec document holds — fields, enumerated
  values, registry entry names it refers to (themes, color schemes, letter maps, aligners and their modes) and the
  format ids — is snake_case, as in MolViewSpec. The viewer's own names stay camelCase, as code: the viewer spec,
  behaviors and their params, config, commands' params, events. Web and CLI conventions stay kebab-case (DOM events, CSS
  classes, data attributes, `--seq-*` tokens, CLI names).
- **Interaction is viewer-spec behaviors** (David, 2026-10-13): what happens on hover, click, double click, Shift-click,
  drag, autoscroll, a header's hover, and how a large set is handled, is a behavior of the viewer spec (`SeqViewerSpec`,
  the viewer's runtime behavior, as Mol\*'s PluginSpec behaviors are) with params and reasonable defaults; a host
  overrides, replaces or removes it there. Input delivery (the DOM engine, the dispatcher) only delivers events; nothing
  reacts outside a behavior. Linked views (the Mol\* adapter) compose their reactions the same way, in their own
  options.
- **Small-agent design:** id references, flat schemas, instructive errors that say what to write instead.
