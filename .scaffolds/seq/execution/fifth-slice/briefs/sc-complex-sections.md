# Lane SC: sections over a source's complex

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/sc`. Model:
strong. Beside GC2 (render only): SC owns `resolve/` and `view/`'s section resolution; GC2 owns `render/`. Touch a file
the other lane owns only to add an export, and list it.

**Task**: HD's drift, decided 2026-11-21 (`packages/seq/TODO.md`, "Sections over a source's complex"): a section `of`
a **dynamic complex** (a complex node an opened source puts in the overlay, `resolve/overlay.ts`) resolves when the
source opens and **reflows on its `changes`**, so a host's own complex shows in a view. Today a section reads only the
static model, and a dropped unit leaves an empty section and its issue.

- the static-model rule's exception for opened sources, said once in core's `_spec` (Resolve, Resolved view, the
  overlay section) with its reason; seq-data-ports' matrix row "a section of a complex after an update: reflows" built
  and its test named;
- a section of a complex not yet held is pending (as a genome assembly before its chunk lands, GC1's
  `StaticModel.loading` pattern), never an error; a unit the update drops leaves its place by the same reflow, and a
  complex the source stops holding says its issue once;
- the reflow is one rebuild a version (no rebuild a hover or a frame); the warm frame stays at 0 bytes (both perf tiers
  are on this lane's path: `pnpm -s -C packages/seq/{core,storybook} test:perf`, B6);
- a story: a host source whose complex gains and loses a chain on a button (the Complex dataset's stand-in is fine),
  light and dark, with its headless twin on a fake clock;
- the TODO line closed (list it as a hot-file need, never edit `TODO.md`).

Read: `packages/seq/core/_spec/README.md` (Resolve, Resolved view, the overlay), `.scaffolds/seq/design/seq-data-ports.md`
(its matrix), `.scaffolds/seq/design/host-data.md` if present, the viewer's spec model and settle (GC1 made them resolve again
when a chunk lands: the same seam).

**Where**: a worktree on `seq/sc`. Headers updated (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files
never staged (`package.json`s, lockfile, `biome.json`, `.scaffolds/seq/`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`): list them. Public API changes allowed in `packages/seq`; anything in `@molstar/core` outside
`repr` or in `@molstar/io` stops the lane and comes back as a question.

**Checks**: as `fx1-checks-and-carried.md` lists them, both perf tiers, `pnpm -s -C packages/seq/core size` (reading
has 0.8 kB left: keep the reflow off the reading path; stop and report if a bound would move), the four `.scaffolds/seq/_checks`
and `features-planned.ts --plan .scaffolds/seq/execution/fifth-slice/plan.md`. A fix's test fails on its parent. Do not claim a
failure is the base's without running it on the base.

**Report**: under 300 words: commits; checks with numbers; the story ids for the orchestrator's look; the public
surface added; sizes; hot-file needs; the `.v6/changelogs/` text.
