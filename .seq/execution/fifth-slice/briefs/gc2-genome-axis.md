# Lane GC2: the genome axis and layouts

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/gc2`. Model:
strong. Beside SC: GC2 owns `render/` and the view's block parts (`blockAxes`, groups, the places index); SC owns
section resolution over the overlay. Touch a file the other lane owns only to add an export, and list it.

**Task**: build [genome-coordinates](../../../design/genome-coordinates.md)'s axis and layouts part (Consequences,
"render" and "view" under `seqx-core`, marked GC2), reading what GC1 declared:

- view: `blockAxes` with a block's direction and scale (blocks on one system merge only on one direction), derived
  blocks keyed by name, overrides by `of`; groups (`section.of` a list: a wider gap and a group label);
  `maxAssemblyBlocks`' closing run; equal block lists follow, unequal refused (`seq.render.follow.blocks`' message);
  the places index (`placesOf(system, position)`, built once per resolved section, for LX, RR and RN);
- render: block strips with their own direction and scale (each block's ruler, `axis.directionMark`, `×N`); the gap
  share on rows and on rings (at most a quarter of the width, 12 px or 1° each; `rings.ts` amended) and runs of small
  blocks in both layouts; group gaps and labels; `section.center`; the `block` and `system` built-in fields and the
  block-ordered domain for themes (EC-218); the SVG's `data-system`, `data-direction`, `data-scale`, `data-assembly`;
  `renderSvg`'s viewport by `{ system, position }`;
- the design's Defaults 1 and 3 (a row, a ring) as drawn there; the Configuration rows marked GC2 built, the rest left
  "Not built" with their lane;
- GC1's declared params this lane reads lose their `seq.resolve.param.unread` info (`block.direction`, `block.scale`,
  `section.center`, whatever else names GC2);
- tests: hit-testing the inverse of placement over mixed directions and scales, runs included (property); offsets
  round-trip through `{ system, position }`; drawing against the data at the fit, across a block gap, inside a run, on
  a reversed block; a ring of 3,000 blocks keeps three quarters of its circle for positions; EC-102, EC-206 (the draw
  part), EC-210, EC-218 named in titles; a whole-genome layout of 25 blocks × 3 layers and of 3,400 blocks with runs
  under `itWithin`; a pan across a block gap at 0 bytes warm (perf tier, B6);
- stories: a Design page "Genome coordinates" in the Seq* storybook: one assembly in a row, a ring, two assemblies on
  one ring (groups), B in A's order with an inverted block reversed, 3,000 generated contigs in runs; light and dark,
  canvas and SVG, each with its headless twin. Toy assemblies generated in the story (GX brings the real dataset
  later); **no downloads**.

Read the design whole (Decision: A section of an assembly, A direction per block, Scaled segments and the zoomed ring,
Positions … when blocks reverse, What layers across blocks read from here; Defaults; Configuration; Tensions: "Reads at
the whole genome" closes on your measure of the 3,000-contig fixture), `axis.md`, `layouts.md`, core's `_spec` (Render:
Blocks, Rings, Export, Configuration; Resolved view), GC1's commits on `packages-seq` for what it declared.

**Where**: a worktree on `seq/gc2`. Headers updated (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files
never staged (`package.json`s, lockfile, `biome.json`, `.seq/`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`): list them. Public API changes allowed in `packages/seq`, core's `repr` and `@molstar/ui`;
anything else in `@molstar/core` or in `@molstar/io` stops the lane and comes back as a question.

**Checks**: as `fx1-checks-and-carried.md` lists them, both perf tiers (`pnpm -s -C packages/seq/{core,storybook}
test:perf`, B6), `pnpm -s -C packages/seq/core size` (reading 0.8 kB left, minimal 4.6 kB: the render growth lands in
the viewers; stop and report if a bound would move), `pnpm test:heavy packages/seq/` for the `itWithin` bounds, the
four `.seq/_checks` and `features-planned.ts --plan .seq/execution/fifth-slice/plan.md`. A fix's test fails on its
parent. Do not claim a failure is the base's without running it on the base.

**Report**: under 300 words: commits; checks with numbers; the story ids for the orchestrator's look (light and dark,
canvas and SVG); the public surface added; the 3,000-contig measure; sizes; what GC3, LX, RR, RN must read; hot-file
needs; the `.v6/changelogs/` text.
