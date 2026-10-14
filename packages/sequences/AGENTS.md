# Seq\* (lib/seqstar)

The entry point for work in the Seq\* packages. Seq\* is SeqViewSpec, a format for sequence views built on
one coordinate model (positions from a single residue to a whole genome, and the maps between them: a protein onto its
gene's codons across introns, an alignment's rows onto its columns, a design onto its wild type, a chain onto its
structure; data placed on any of these shows on the others), with its data model, SeqQL, its query language, and the
builder. The Seq\* viewer renders it, and so does Nightingale. It is self-contained: it imports nothing of the app
around it, and this file and the spec are all it needs.

- **The spec** is `_spec/` (start with [_spec/README.md](_spec/README.md): what Seq\* is, the packages, the
  principles, the glossary). It says what, never how: contracts and meaning. Tables (nodes, types, layer types,
  themes, issue codes, SeqQL's symbols, commands, behaviors, config) are generated from the code into
  `_spec/reference/`; never edit those by hand.
- **The code is the truth.** The spec holds the contracts and what they mean; where the two disagree, fix whichever
  is wrong (often the spec), in the same commit.
- **Packages:** `core` (platform-free, no dependencies) and `viewer` (with the interactions extension), plus the
  adapters `nightingale` and `molstar-poc` and the private `storybook` (the table in `_spec/README.md`). Each folder
  of a package is an entry: inside a folder import relatively; across folders name the entry
  (`@molstar/seq-core/model`, the package importing itself), never a path into another folder.
  `core/_test/layers.test.ts` holds the core's folder graph. `scripts/` (`@molstar/seq-scripts`) holds Seq*'s own checks
  and measures, and each package's `_test/fixtures/` its recorded data: nothing outside `lib/seqstar` is needed.
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
- **What is not built** is `todo.md`: work in progress, next, deferred features (each with its trigger), follow-ups
  and known costs. An item leaves it in the commit that builds it; a deferral comes in with its trigger.
- **Phase ids** in comments (`LY1`, `QL2b`, `CS2`, …) name the development phases; `_spec/dev-phases.md` says
  what each was.

## Tests: run the tier the change needs

Each package's `test` is the fast tier: about 5 s a package on the M4 (2026-10-14), the viewer's 1,700 tests about 13 s.
The slow tiers are separate scripts:

| Script | What | Run it |
|---|---|---|
| `test` | node unit tests (`vitest run`): fast and deterministic, no wall-clock timing, no Chrome, no heap guards; add `--changed` or `related <files>` for only the affected tests | on every edit loop, in the packages you touched; a folder narrows it (`vitest run _test/format`) |
| `test:slow` (core, viewer, storybook) | node tests tagged `slow`: wall-clock bounds and ratios (`itWithin` and the scaling tests), the storms, fuzz and oracle runs at full size, a million rows and large fixtures, Chrome from Node (PNG export), the CLIs in a process, the four largest smoke stories | when you change what such a test covers, and at checkpoints; `mise run seqstar:test:slow` runs every package |
| `test:browser` | Chromium: pixels, UI, surface (the viewer's `_test/**/{canvas,dom,react}`; Nightingale's `mount`, `_test/dom`) | when you change drawing, the DOM or the UIs; pass the affected files |
| `test:heap`, `test:browser-heap` | heap guards | when you change a hot path (draws, input, streams, query kernels); pass the touched guard files (viewer heap: about 4 min) |
| `test:heap-deep` (viewer) | the heap guards as before the 2026-10-12 trim: full steps (catch about 22 bytes a step, the trimmed run about 43) and the dropped cases; about 10 min | when chasing a small leak, or before a release; never in an edit loop |
| `test:full` | all of the above | at reviews and checkpoints only; `mise run seqstar:test:full` runs every package one at a time |

**CI policy** (David, 2026-10-14): CI runs each package's `test` only (Node, fast, deterministic: no heap guards, no
Chromium, no wall-clock timing). `test:slow`, `test:browser`, `test:heap`, `test:browser-heap` and `test:heap-deep`
are opt-in: by hand, at checkpoints, or in a separate manual or scheduled job, never in the default workflow. There is
no CI workflow here yet; Mol\*'s gets this rule with the move (designs/seqstar/handoff.md).

- **The `slow` tag** (vitest's test tags, `vitest.shared.ts`): a test slow by nature carries `{ tags: ["slow"] }` (a
  `describe` can carry it for all its tests) and keeps a quick representative untagged (the fuzz at 100 mutations, the oracles at
  a fifth, a storm of 50, 20,000 rows for a million, most smoke stories). A wall-clock bound goes in `itWithin`
  (`core/_test/model/timing.ts`): the test runs in the fast tier without the bound and again in the slow tier with
  it. A loop of thousands of checks records its disagreements (`mismatches`, `core/_test/model/mismatches.ts`) and
  checks once: an `expect` costs microseconds.
- **Files share a worker's modules** (`isolate: false`): a file that mocks or resets modules (`vi.mock`, `vi.doMock`,
  `vi.resetModules`) is found by reading it and runs isolated, in the `isolated` project; module state a file
  changes (a limit, a config default) it restores in `afterEach`, or the next file in the worker sees it.

- Heap guards are load-sensitive: run them only when no other suites are running. If one fails, re-run it alone before you diagnose it.
- In a per-frame path, read enum-like object members (`SkinToken.x`, `PathVerb.x`) once into module constants: under the test runner's module wrapper each read allocated (about 41 bytes a frame), a false "allocates only under vitest" result.
- A new draw guard gets a file of its own, or a fresh module per case (`vi.resetModules`): cases in one file share V8's feedback, so one can allocate only after another has run. Keep the shared warm-up (`FRAME_WARM`, `FRAME_WARM_MS`); shortening it made guards flaky. A guard that measures many frames apiece (the conformance kit's, a layer per zoom) warms each with `warmFrame` (`viewer/_test/heap.ts`: until the frame stops allocating, at most `FRAME_WARM_MS`): a draw called once a frame climbs to TurboFan over thousands to tens of thousands of frames, allocating in its middle tier until then, and a count alone measured that climb. A guard with a timed warm-up settles after it (`settleFrame`, `FRAME_SETTLE_MS`). A positive control allocates through `leak` (garbage, the same bytes every step), never onto a kept list: a list's reallocation inside the shorter measured run cancels the longer's growth, and the control goes blind.
- **The guard template:** every `*.heap.ts` runs its positive control as its guards run their step: `expect(controlCaught(() => <the file's own warm-up, settle and measurement> of leaking(step))).toBe(true)` (`core/_test/model/heap.ts`, re-exported by `viewer/_test/heap.ts`; `canvas/heap-browser.ts` in the browser). `controlCaught` counts a run the measurement refuses because a collection ran during it as caught (in the browser that is how a leak often shows). The spec check fails on a guard file without it unless `scripts/guards-pending.txt` lists it; converting a file removes its line.
- Pixel hashes are exact and the machine's: `canvas/case9-pixels.test.ts` keeps a record per user agent (`case9-pixels.json`) and fails on one with none, naming it; record a setup's hashes with `VITE_SEQ_PIXELS=record` (and `VITE_SEQ_PIXELS_MACHINE`) only after checking its render, and again after an intended look change.
- Run the guards on the pinned Node (`mise.toml`, 24): where it is not installed, `mise run` falls back to the PATH's (Node 25 on the M4 with 16 GB, 2026-10-13), whose V8 tiers up and inlines differently (a guard can pass on one and fail on the other).
- `mise run check`, `mise run typecheck`, the export guard (`mise run seqstar:exports`, `scripts/exports.ts`) and the spec check (`mise run seqstar:spec`, `scripts/spec.ts`: every `seqspec:` citation names a heading, the spec's identifiers exist in the code) must pass before every commit.

## Working in the shared tree

Several agents often work in this tree at once.
- Commit with a temporary git index built from HEAD plus your own hunks. Never `git commit -- <paths>` on a file another agent has staged: it sweeps their hunk in.
- After that commit, `git reset -q -- <your paths>` (exact paths): the shared index still holds their old versions, which `git status` then shows as staged reversals of your commit.
- Give each temporary index a fresh, literal name and build it (`git read-tree HEAD`) in the same command as the commit: HEAD moves between commands, and a reused index from an older HEAD silently reverts other lanes' commits.
- Never move the branch: no `git reset` of HEAD, `update-ref`, `rebase`, `commit --amend` or force. If a commit of yours is wrong, fix it forward with a new commit and say so in your report; if it touched another lane's work, stop and report instead of repairing history.
- Before you `git add` a whole file into the temporary index, check `git diff -- <file>` shows only your hunks. If another agent has hunks in it (or the file is behind HEAD), stage only yours: `git diff -- <file>` filtered to your hunks, then `git apply --cached` (with context lines: never `-U0`/`--unidiff-zero`, which can place a hunk at the wrong spot). Whole-file adds have swept other lanes' unfinished hunks into a commit and reverted others' committed text.
- Verify in a clean worktree of HEAD plus your hunks, with the relative `node_modules` links copied.
- No more than two agents write in the same package folder at a time (the core's folders, `viewer/src/interactions` and the viewer's entries count separately).
- Scratch files go in the session scratchpad, in a folder named for the phase. Never delete files with `rm`: leave scratch files in place, and give each temporary index a fresh, literal name.

## Size, specs and stories

- **Size:** phases carry size budgets in production lines (`mise run seqstar:ported-lines`, `scripts/ported-lines.ts`). Record an overrun and what took the room; don't stop for size (David, 2026-10-12: "the budget overflow is ok"). Do stop and report when the work turns out different from the plan.
- **Spec:** the spec says what, never how. A change to a contract updates `_spec/` in the same commit; a change to a
  table updates its reference (`vitest run -u` in that package: the reference tests write `_spec/reference/`). Cite
  `seqspec:<file>#<slug>` (the slug is GitHub's anchor of the heading, e.g. `seqspec:model#closure`) once per file or
  per contract, where code implements it; never section numbers, never a generated reference, never for how
  something is built. Nothing in Seq\* points outside `lib/seqstar` for design or history. The code decides what
  the docs say.
- **Stories and docs:** every major feature a user sees gets a Storybook story, or extends a case story, checked in the toolbar's light, dark and blue themes. Evals and internal tooling don't need one. Feature work also updates the `Seq*/Docs` pages it touches (Features, Concepts, Extending, Streaming, Color schemes, Query language).
  A new story brings its entry in `lib/seqstar/storybook/src/_stories/support/seqstar-smoke.ts` (the documents it opens, the patches its buttons send, the warnings it shows on purpose): `lib/seqstar/storybook/_test/smoke.test.ts` and `smoke-more.test.ts` (half each, `smoke-stories.ts`) open each on a headless viewer in the `test` tier (the four largest, `SLOW_STORIES`, in `test:slow`) and fail on an error, an unlisted warning, an empty visible section, or a story export with no entry. A story's builder lives in `lib/seqstar/storybook/src/_stories/support/` (not in the `.stories.tsx`), so the test can import it without React.
- **Where the Storybook content lives:** `lib/seqstar/storybook` (`@molstar/seq-storybook`, private), all under `src/_stories/` (Mol\* v6's layout; exports and `files` hide `_stories` and `_test` in every package): the stories, `support`, `docs` (the `Seq*/Docs` pages), `ui` (its own small components on Mol\* Prism's `--msp-*` tokens, a local copy in `ui/prism` until the move: it imports no `@untold/*` or `stories-*`, and Biome enforces that), `_test` (the smoke test and the stories' search and patch tests). `apps/storybook` adds the globs (`@molstar/seq-storybook/vite`) and imports `@molstar/seq-storybook/tokens.css`. The one story that mocks Stories' own workspace screens, `Seq*/Cases/Workspace drafts`, stays in `apps/storybook/src/seqstar/` (the registry lists it as skipped).
- **The stories' look is Mol\* Prism's** (David, 2026-10-13): every story is a `PrismPage` (title, caption, cards: `seqstar-prism.tsx`), themed by the Storybook toolbar (one viewer per story, as Mol\*'s Storybook); the viewer draws on a `PrismCanvas` (`ViewerCard`'s default), pinned to a theme only where comparing themes is the point (the color schemes pages). Toggles are `Choices`; spacing, type and colors are the `--msp-*` tokens, no borders (tone and the inline block's hairline).

## Design

- **Domain-independent core:** biology meaning lives in tables, presets and app-level code, never in core types. Complement is a named letter map, and an exon is a range.
- **Performance is a contract:** no allocation on hot paths (per frame, per pointer event, per run), and no boxed doubles across helper calls. Every hot path has a heap guard with a positive control.
- **snake_case is the format's** (David, 2026-10-13): every name a SeqViewSpec document holds — fields, enumerated
  values, registry entry names it refers to (themes, color schemes, letter maps, aligners and their modes) and the
  format ids — is snake_case, as in MolViewSpec. The viewer's own names stay camelCase, as code: the viewer spec,
  behaviors and their params, config, commands' params, events. Web and CLI conventions stay kebab-case (DOM events,
  CSS classes, data attributes, `--seq-*` tokens, CLI names).
- **Interaction is viewer-spec behaviors** (David, 2026-10-13): what happens on hover, click, double click,
  Shift-click, drag, autoscroll, a header's hover, and how a large set is handled, is a behavior of the viewer spec
  (`SeqViewerSpec`, the viewer's runtime behavior, as Mol\*'s PluginSpec behaviors are) with params and reasonable
  defaults; a host overrides, replaces or removes it there. Input delivery (the DOM engine, the dispatcher) only
  delivers events; nothing reacts outside a behavior. Linked views (the Mol\* adapter) compose their reactions the
  same way, in their own options.
- **Small-agent design:** id references, flat schemas, instructive errors that say what to write instead.
