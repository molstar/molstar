# Mol\*

Mol\* is a library and set of apps for viewing and analyzing molecular structures and volumes. This branch line is
**v6**: a pnpm workspace of ESM packages, explicit plugin composition, per-app builds, and a new UI (Prism) with an
optional AI assistant. Start with the [v6 summary](.v6/designs/summary.md), the
[workspace guide](.v6/plans/workspace-usage.md) and the [checklist](.v6/plans/checklist.md) (what is done and what is
not).

## Where things are

| Path                                 | What                                                                                                                                |
| ------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------- |
| `packages/{core,io,model,graphics}`  | library layers, imports one way: `plugin → graphics → model → io → core`; state [_spec](packages/core/src/state/_spec/README.md)    |
| `packages/plugin/{core,ui,headless}` | the plugin runtime, its classic React UI, Node capture                                                                              |
| `packages/mvs/{builder,runtime}`     | MolViewSpec: standalone builder, Mol\* runtime                                                                                      |
| `packages/query-language`            | MolScript and the PyMOL, VMD and Jmol transpilers                                                                                   |
| `packages/prism/*`                   | Prism, the new UI ([AGENTS.md](packages/prism/AGENTS.md), [_spec](packages/prism/viewer/_spec/README.md))                           |
| `packages/ai/*`                      | the agent harness and what runs on it ([AGENTS.md](packages/ai/AGENTS.md), [_spec](packages/ai/_spec/README.md))                    |
| `packages/sequences/*`               | Seq\*, the sequence viewer and SeqViewSpec ([AGENTS.md](packages/sequences/AGENTS.md), [_spec](packages/sequences/_spec/README.md)) |
| `extensions/`, `apps/`, `examples/`  | each owns its dependencies and build                                                                                                |
| `servers/`, `cli/`, `distributions/` | servers, command-line packages, the `molstar` distribution                                                                          |
| `smoke/`                             | packed-consumer and browser checks                                                                                                  |
| `apps/docs`                          | the documentation site (Astro Starlight): the v5 and v6 docs in `src/content/docs/`                                                 |
| `storybook/`                         | Storybook for the packages registered in `storybook/packages.ts` (stories in `_stories/` folders)                                   |
| `.v6/`                               | v6-wide designs and plans (`designs/`, `plans/`); a package's own decisions are in its `_spec/`                                     |
| `scripts/workspace/`                 | the package inventory (`inventory.json`) and the workspace checks                                                                   |

## Commands

Node 22.13+ and the pnpm version in `package.json`; `mise install` sets up both (`mise.toml`, which also wraps the
scripts below as tasks: `mise tasks`).

```sh
pnpm install
pnpm build                 # libraries, apps, distribution
pnpm dev:apps -- viewer    # one app in watch mode
pnpm lint                  # Biome
pnpm format                # Biome (code, JSON) and Prettier (Markdown, YAML, SCSS); format:check in CI
pnpm test:unit             # unit tests (Vitest, ESM; one run for every package, `vitest.config.ts`)
pnpm test:unit packages/ai # a folder or file filter; `pnpm test:watch` for watch mode
pnpm test                  # lint and the unit tests
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
- **ESM with NodeNext**: relative imports end in `.js`. Source is resolved through the `molstar-src` export condition in
  development and tests.
- **Formatting and lint** are the repository's (`biome.json`, `.prettierrc.json`): run `pnpm format` and `pnpm lint`
  before committing. Biome's per-package `noRestrictedImports` rules are architecture, not style: do not loosen them to
  make an import pass.
- **Headers**: every source file starts with the Mol\* copyright header and its authors; add yourself to the headers of
  files you change substantially.
- **Tests** live in `_test/` folders next to the code, run by Vitest
  (`import { describe, expect, it, vi } from 'vitest'`, no globals; DOM tests opt in with
  `// @vitest-environment jsdom`); they never reach the network.
- **User-facing changes** go into the `[Unreleased]` section of `CHANGELOG.md` (v6 collects them in
  `.v6/plans/changelog.md`); a v6 break of a public API is recorded in `.v6/plans/breaking-v6-changes.md`.
- **Keep the docs true**: a change that moves what a package's `_spec/` (or a `.v6` design or plan) says updates it in
  the same commit; code comments cite the spec section they follow (`harness.md §6`, `api.md §2`). Specs are for agents:
  terse, a rule and its reason, no history; open items go to the package's TODO.
- **Question the specs**: a spec records decisions, some made in passing by agents. Working in a package, re-read its
  foundational choices and raise any that the code, measurements or the maintainers' intent contradict, or that have no
  reason recorded; ask, don't silently follow or rewrite.
