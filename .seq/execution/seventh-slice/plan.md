# Plan: the seventh slice (draft)

**Draft, 2026-11-23**: not agreed and not started. It starts after the sixth slice moves Prism onto the new Seq\*, so
nothing in the repository depends on the old one.

Work id `seventh-slice`. **One Seq\***: the old implementation removed and the new one takes its names.

- The author, 2026-11-23: "slice 7 will then be removal of packages/sequences and renaming of the @molstar/seqx to just
  /seq."

## Stages

- [ ] **Stage 0: what the old Seq\* holds that nothing else does.** A sweep of `packages/sequences` (its packages,
      stories, docs, specs, fixtures, test-harness users, the Nightingale adapter, CLI and smoke entries, CI jobs,
      `storybook/packages.ts`, the inventory, the docs site's v6 pages): each item either already has its new form,
      moves (a fixture, a doc page, a story worth keeping), or is dropped, the author agreeing the list.
- [ ] **Stage 1: remove `packages/sequences`**: the packages, their inventory rows, root references, storybook
      registration, CI and lint rules (`biome.json`'s `packages/sequences` overrides), the AGENTS.md and README rows;
      `pnpm check:workspace`, the build and every lane green.
- [ ] **Stage 2: rename `@molstar/seqx-*` to `@molstar/seq-*`** (core, viewer, examples, storybook, `ext-rondo`, the
      `cli/seqx-render` CLI and any `seqx` binary or story title): package names, imports, exports, folder names where
      they say `seqx`, the `seq` storybook title "Seq\* (next)" to "Seq\*", biome rules, the inventory, docs and
      changelog entries; the language's own name and file format unchanged unless the author says otherwise. A
      repository-wide check that no `seqx` is left outside the history.
- [ ] **Stage 3: the docs and the changelog**: the v6 docs site's Seq\* pages, `.v6/plans/breaking-v6-changes.md`
      (the old packages gone, the new names), `.v6/changelogs/`.
- [ ] **Record**.

## Questions for the author (stage 0)

1. **Package names**: `@molstar/seq-core`, `@molstar/seq-viewer`, … reuse the old packages' names: say so in the
   breaking changes as a new API under old names, or pick fresh ones.
2. **The CLI**: `seqx-render` becomes `seq-render`?
3. **The storybook**: the old Seq\* stories go; any kept as references for comparison until then?
