## Workspace and tooling

- Test lanes (`AGENTS.md`, Tests): every test is quick (untagged), medium (`{ tags: ['medium'] }`) or heavy
  (`{ tags: ['heavy'] }`). `pnpm test:unit` and `pnpm test` run the quick lane, `pnpm test:medium` and
  `pnpm test:heavy` the others (`vitest run --mode <lane>`; `--mode pr` quick and medium, `--mode all` every lane);
  `pnpm test:native [-- <lane>]` runs a lane with the native GL backend (default `pr`). The native GL suite, the
  largest fixtures and the rest of each sweep are medium, with a quick representative each. Seq\*'s `slow` tag and
  `test:slow` scripts are now `heavy` and `test:heavy`.
- CI: the quick lane in 4 shards on Node 24 (`Quick tests (n/4)`) beside one medium job (`Medium tests (native GL)`), in place of
  the single `Native tests` job; the heavy lane in its own workflow (`tests-heavy.yml`: nightly, by hand, and on a PR
  labelled `heavy tests`).
- Mol\*'s tests share their worker's module registry, as Seq\*'s did (a third less CPU). `vitest.setup.ts` restores
  spies, globals, env and timers after each file; files that mock modules, or need a registry of their own, run
  isolated. `pnpm lint` fails on an `export *` inside an import cycle (`scripts/workspace/export-star-cycles.mjs`): under
  Vite's module runner its names depend on import order. The 21 such re-exports (`model/structure`, `core/state`, the
  polymer iterators) now name their exports. `to_mmCIF` from `@molstar/model/model/structure` was undefined after
  `formats/structure/mmcif` had been imported first.
- `mesh-builder.test.ts` compares any NaN as one bit pattern (Node 24 on x64 writes another).
