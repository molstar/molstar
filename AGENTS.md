# Mol\*

Mol\* is a library and set of apps for viewing and analyzing molecular structures and volumes. This branch line is
**v6**: a pnpm workspace of ESM packages, explicit plugin composition, per-app builds, and a new UI (Prism) with an
optional AI assistant. Start with the [v6 summary](.v6/designs/summary.md), the
[workspace guide](.v6/plans/workspace-usage.md) and the [checklist](.v6/plans/checklist.md) (what is done and what is
not).

## Where things are

| Path                                         | What                                                                                                                                                                                                                     |
| -------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `packages/{core,io,model,graphics}`          | library layers, imports one way: `plugin → graphics → model → io → core`; state [_spec](packages/core/src/state/_spec/README.md), model [_spec](packages/model/_spec/README.md), io [_spec](packages/io/_spec/README.md) |
| `packages/plugin/{core,ui,headless}`         | the plugin runtime, its classic React UI, Node capture                                                                                                                                                                   |
| `packages/mvs/{builder,runtime}`             | MolViewSpec: standalone builder, Mol\* runtime                                                                                                                                                                           |
| `packages/query-language`                    | MolScript and the PyMOL, VMD and Jmol transpilers                                                                                                                                                                        |
| `packages/prism/*`                           | Prism, the new UI ([AGENTS.md](packages/prism/AGENTS.md), [_spec](packages/prism/viewer/_spec/README.md))                                                                                                                |
| `packages/ai/*`                              | the agent harness and what runs on it ([AGENTS.md](packages/ai/AGENTS.md), [_spec](packages/ai/_spec/README.md))                                                                                                         |
| `packages/sequences/*`                       | Seq\*, the sequence viewer and SeqViewSpec ([AGENTS.md](packages/sequences/AGENTS.md), [_spec](packages/sequences/_spec/README.md))                                                                                      |
| `packages/extensions/`, `apps/`, `examples/` | each owns its dependencies and build                                                                                                                                                                                     |
| `servers/`, `cli/`, `distributions/`         | servers, command-line packages, the `molstar` distribution                                                                                                                                                               |
| `smoke/`                                     | packed-consumer and browser checks                                                                                                                                                                                       |
| `apps/docs`                                  | the documentation site (Astro Starlight): the v5 and v6 docs in `src/content/docs/`                                                                                                                                      |
| `storybook/`                                 | Storybook for the packages registered in `storybook/packages.ts` (stories in `_stories/` folders)                                                                                                                        |
| `.v6/`                                       | v6-wide designs and plans (`designs/`, `plans/`); a package's own decisions are in its `_spec/`                                                                                                                          |
| `scripts/workspace/`                         | the package inventory (`inventory.json`) and the workspace checks                                                                                                                                                        |

## Commands

Node 22.13+ and the pnpm version in `package.json`; `mise install` sets up both (`mise.toml`, which also wraps the
scripts below as tasks: `mise tasks`).

```sh
pnpm install
pnpm build                 # libraries, apps, distribution
pnpm dev:apps -- viewer    # one app in watch mode
pnpm lint                  # Biome
pnpm format                # Biome (code, JSON) and Prettier (Markdown, YAML, SCSS); format:check in CI
pnpm test:unit             # the quick test lane (Vitest, ESM; one run for every package, `vitest.config.ts`)
pnpm test:unit packages/ai # a folder or file filter; `pnpm test:watch` for watch mode
pnpm test                  # lint and the quick lane
pnpm test:medium           # the medium lane; `pnpm test:native -- medium` with the native GL backend
pnpm test:heavy            # the heavy lane (after `npx tsc -b`: its CLIs run built)
pnpm test:live             # live agent runs on real models (opt in: keys in .env.local; costs money)
pnpm check:workspace       # package boundaries, exports, imports, the inventory
pnpm deploy:local          # build, then the GitHub Pages folder in deploy/data (apps, prism/, ui-storybook/)
npx tsc -b                 # every project, by references
```

## Rules for every change

- **Packages own their dependencies.** Import only what a package's `package.json` declares, through its exports; never
  reach into another package's `src` by path. A new package goes into `scripts/workspace/inventory.json`, the root
  `tsconfig.json` references; `pnpm check:workspace` checks them. Tests need no mapping: Vitest resolves workspace
  packages to their sources through the `molstar-src` export condition.
- **Reuse before writing helpers.** Before inlining a generic routine (hashing, sorting, array and typed-array helpers,
  sets, maps, memoization, vector and matrix math), look in `@molstar/core`: `util/`, `data/util/` (hash functions,
  sort, array, grouping), `data/int/` (`OrderedSet`, `SortedArray`, `Interval`), `math/` (`linear-algebra/`,
  `geometry/`) and `algorithms/` (`sequence/`: alignment). If nothing fits and the code is generic, suggest adding it
  there with a test instead of keeping a private copy; ask before changing a public API. Inlining in a hot loop on
  purpose is fine: say why in a comment.
- **ESM with NodeNext**: relative imports end in `.js`. Source is resolved through the `molstar-src` export condition in
  development and tests.
- **Formatting and lint** are the repository's (`biome.json`, `.prettierrc.json`): run `pnpm format` and `pnpm lint`
  before committing. Biome's per-package `noRestrictedImports` rules are architecture, not style: do not loosen them to
  make an import pass.
- **Headers**: every source file starts with the Mol\* copyright header and its authors. Every file a change modifies
  gets its header updated in the same change: an `@author` line for the change's author (the git user, e.g.
  `@author David Sehnal <david.sehnal@gmail.com>`), added after the existing ones unless it is already there, and the
  copyright year range extended to the current year (`2019` → `2019-2026`, `2019-2024` → `2019-2026`).
- **Tests** are in one of three lanes ([Tests](#tests-three-lanes)) and live in `_test/` folders next to the code, run
  by Vitest (`import { describe, expect, it, vi } from 'vitest'`, no globals; DOM tests opt in with
  `// @vitest-environment jsdom`); they never reach the network. A new `_test/` folder in a published package needs
  `"./<dir>/_test/*": null` in its `package.json` exports, or its `./*` wildcard publishes the tests and
  `check:workspace` fails in CI; `pnpm lint` checks it from the sources (`scripts/workspace/test-exports.mjs`, `--fix`
  adds the entries).
- **User-facing changes** go into the `[Unreleased]` section of `CHANGELOG.md`. v6 collects them in
  `.v6/plans/changelog.md`, but a branch never edits that file: it adds its own entry file to `.v6/changelogs/` (see
  [its AGENTS.md](.v6/changelogs/AGENTS.md)), and the entries get folded into the draft from time to time. A v6 break of
  a public API is recorded in `.v6/plans/breaking-v6-changes.md`.
- **Keep the docs true**: a change that moves what a package's `_spec/` (or a `.v6` design or plan) says updates it in
  the same commit; code comments cite the spec section they follow (`harness.md §6`, `api.md §2`). Specs are for agents:
  terse, a rule and its reason, no history; open items go to the package's TODO.
- **Question the specs**: a spec records decisions, some made in passing by agents. Working in a package, re-read its
  foundational choices and raise any that the code, measurements or the maintainers' intent contradict, or that have no
  reason recorded; ask, don't silently follow or rewrite.

## Tests: three lanes

A test's lane is its Vitest tag; the root `vitest.config.ts` picks a lane by `--mode` and imports only the files that
can hold its tests.

| Lane                                                     | What goes there                                                                                                                                                                                                    | Command            | CI                                                                         |
| -------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------ | -------------------------------------------------------------------------- |
| quick (untagged)                                         | deterministic, no clock, small fixtures: well under a second a test on CI, seconds a package                                                                                                                       | `pnpm test:unit`   | every push, in 4 shards (`Quick tests`)                                    |
| medium (`{ tags: ['medium'] }`, `// @module-tag medium`) | a second or more a test, deterministic: large fixtures (1cbs over two cells, the exosome animation), the rest of a sweep (axe over every tab, every shader, every layer), a whole app session, the native GL suite | `pnpm test:medium` | every push, its own job (`Medium tests (native GL)`, with `gl` under Xvfb) |
| heavy (`{ tags: ['heavy'] }`, `// @module-tag heavy`)    | wall-clock bounds and ratios, scaling by size, full-size fuzz and oracles, a million rows, CLIs in a process (built first), Chrome                                                                                 | `pnpm test:heavy`  | nightly, by hand, and on a PR labelled `heavy tests` (`tests-heavy.yml`)   |

- **A test slow by nature moves, and keeps a quick representative**: one layer of a sweep, the smaller fixtures, a short
  fuzz, a smaller scale. A wall-clock bound goes in `itWithin` (`packages/sequences/core/_test/model/timing.ts`): the
  test runs quick without the bound and heavy with it. Coverage is never dropped: medium runs on every PR. Why: the PR
  checks finish in minutes, and the edit loop in seconds.
- **A test slow by waste is fixed, not moved**: a fixture computed once per file, not per `describe`; fake timers where
  real ones only wait; no rebuild per test.
- **Tag in the test's own file**, inline or by `@module-tag` on a line of its own: the config reads the files to decide
  what a lane imports (a filtered-out file is still imported, and imports are about half the run). Helpers that tag are
  heavy only (`itWithin`, Seq\*'s smoke stories): the heavy lane imports every file. A test file imported by another
  runs its medium tests in its own file only. Tags inherit (a heavy test in a medium `describe` is heavy), and the lanes
  are exclusive. A medium test may build on the state an earlier medium test left, never on a quick one's.
- `--mode pr` runs quick and medium (what a PR checks; `pnpm test:native` runs it with the GL backend), `--mode all`
  every lane; `pnpm check:publish` runs quick and medium (the GL suite skips without `gl`). Run the quick lane in the
  edit loop, the medium lane when you touch what it covers, the heavy lane when you touch what it bounds and at
  checkpoints.
- **Not lanes**: the packages' opt-in tiers, by hand at checkpoints and never in CI: `test:perf` (allocation, the
  `*.quick.ts` cases: not the quick lane), `test:scaling`, `test:browser` and the heap guards (`test:heap`,
  `*.heap.ts`). See the [Seq\*](packages/sequences/AGENTS.md) and [Prism](packages/prism/AGENTS.md) guides.
- **Files share their worker's module registry** (`isolate: false`: importing the libraries again per file was about a
  third of the run's CPU). So a file restores the module state it changes (a registry, a config default);
  `vitest.setup.ts` restores spies, stubbed globals and env and real timers after every file. No `export *` inside an
  import cycle (`pnpm lint`, `scripts/workspace/export-star-cycles.mjs`): under Vite's module runner the names it
  re-exports depend on which file imported the cycle first. A file that mocks or resets modules (`vi.mock`,
  `vi.resetModules`) runs isolated, found by reading it; one that needs a fresh registry otherwise goes in
  `vitest.config.ts`'s `Isolated`, with why.
