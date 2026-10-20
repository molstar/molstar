## Workspace and tooling

- Test lanes (`AGENTS.md`, Tests): every test is quick (untagged), medium (`{ tags: ['medium'] }`) or heavy
  (`{ tags: ['heavy'] }`). `pnpm test:unit` and `pnpm test` run the quick lane, `pnpm test:medium` and
  `pnpm test:heavy` the others (`vitest run --mode <lane>`; `--mode pr` quick and medium, `--mode all` every lane);
  `pnpm test:native [-- <lane>]` runs a lane with the native GL backend (default `pr`). The native GL suite, the
  largest fixtures and the rest of each sweep are medium, with a quick representative each. Seq\*'s `slow` tag and
  `test:slow` scripts are now `heavy` and `test:heavy`.
- CI: the quick lane in 4 shards (`Quick tests (n/4)`) beside one medium job (`Medium tests (native GL)`), in place of
  the single `Native tests` job; the heavy lane in its own workflow (`tests-heavy.yml`: nightly, by hand, and on a PR
  labelled `heavy tests`).
