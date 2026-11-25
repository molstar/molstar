# Lane GC1: the genome language and resolve

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/gc1`. Model:
strong. The language changes **once** this slice, in this lane (the minor bump to the next language version; the
version file and its history are this lane's).

**Task**: build [genome-coordinates](../../../design/genome-coordinates.md)'s language and resolve part (its Decision,
Configuration, Consequences "`seqx-core` (GC1)"), and **declare the other designs' params** with their docs, checks and
builder entries, no reader or drawing: [genome-formats](../../../design/genome-formats.md)'s (`source.features`,
`samples`, `map.source`, `line: 'gene'`, the BED family's `fields`), [across-sections](../../../design/across-sections.md)'s
(`bridge`, `pairs.shape`, a map as pairs data, `follow_y` a second section, a rows system on a grid axis,
`member_order`, `member_order_y`) and [region-reads](../../../design/region-reads.md)'s (`pairs.source`, `pairs.field`,
`grid.source`, `grid.field`, `rows` and `values` then optional, the `matrix` kind and layer type; the author,
2026-11-21: "Yes, in GC1"). Each declared param's reader is the later lane's (SC, GC2, GC3, BE, GG, VC, RG, LA); a
declared param a reader lacks says so by the existing issue, never silently.

- the `assembly` kind; `section.of` an assembly or a list; `section.include`; `block.direction`, `block.scale`;
  `navigator.switch`; the locus and block forms (`(block :name …)` quoted names); the builder's `genome()`,
  `assembly()`, `colinearBlocks()`; the names rule (aliases, `match: chr`, `MT` against `chrM` refused); derived
  sequences from a source's coverage; a sequence without a length an error (`seq.resolve.assembly.lengths`); a map
  between assemblies as one map;
- the Configuration tables copied into core's `_spec` (the config-keys test); the language docs and catalogs; the
  version history entry; every new issue in the catalog;
- edge cases named in test titles (`EC-<n>` the design assigns to GC1: read its "Edge cases" section);
- property tests the plan's Testing strategy names for coordinates that are resolve's (the names rule never maps an
  alias twice; a block's direction reversed twice is the identity at the model level);
- the minimal viewer is at **939.4 of 940 kB**: the language's growth there is the risk — keep docs and catalogs out
  of the reading path (as CN did), measure with `pnpm -s -C packages/seq/core size` at each step, and stop and report
  if a bound would move (sizes are the author's).

**Where**: a worktree on `seq/gc1`. Headers updated (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot
files never staged (`package.json`s, lockfile, `biome.json`, `.scaffolds/seq/`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`): list them. Read `packages/seq/core/_spec/README.md` (Language, Resolve), the language's
existing 0.5.0 bump (ML1's work) for the pattern, the repo `AGENTS.md`, `packages/seq/AGENTS.md`.

**Checks**: as `fx1-checks-and-carried.md` lists them (both perf tiers if resolve's hot path changes), plus
`node .scaffolds/seq/_checks/defaults-designed.ts`. Nothing is drawn, so no look.

**Report**: under 300 words: commits; checks; the public surface added (names, against the designs'); the version and
its history line; the sizes; what each later lane must read; hot-file needs; the `.v6/changelogs/` text.
