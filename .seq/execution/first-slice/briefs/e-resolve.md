# Lane E: resolve (the static model and the resolved view)

**Goal**: in `packages/seq/core/src/resolve/` and `packages/seq/core/src/view/`: tier 1 of the compiled forms
([composition.md](../../../design/composition.md)): the static model and the resolved view, for the first slice's kinds.

- **The static model** (`resolve/`): every coordinate system the spec declares (sequences, with length, alphabet,
  cyclic; the built-in numberings `index` and `ordinal`); the slice's maps resolved to their encodings
  ([maps.md](../../../design/maps.md)): `offset` implicit, `codon` as runs (exons, strand, step 3, a direction on every run);
  one map interface (image of a range forward and back as runs, the boundary rules EC-28, EC-29, element status); paths
  composed without expanding (bounded length, at most one codon step, shortest unambiguous, ambiguity reported); projection
  of a Location along a path, with closure (a nucleotide lifts to its residue's codon, EC-30) and the copy rule (EC-25,
  for units when present); sources bound: each `source` node to a `SourceEntry` (derive it from the language's source
  params: one definition) and to a runtime by `format` from a registry given by the caller (the `seq` runtime for inline
  data); each annotation's access plan (inline columns as a `seq` source of one tile; a source by region). It implements
  `QueryModel` (`seqx-core/query`), including `fields` and `chain` (no structures yet: `chain` answers none with an issue),
  so the query forms run on a real spec. Pure, keyed by address, cached by the spec's hashes.
- **The resolved view** (`view/`): views, sections (axis: a system and its range), tracks, layers with their params,
  themes passed through as values (the theme language comes later), each layer's fallback chain, and each layer's data
  as an address or query with the static path from its data's system to its section's system (none when the same).
  Follow links between sections resolved to their paths (scale 1, 3, 1/3). No data, no target, no pixels.
- **Edge cases** on every encoding of the same map (implicit offset, runs, expanded for a small region) with equal
  answers: EC-9 to EC-12 (cyclic, where they apply to offset), EC-15 to EC-21, EC-25 to EC-30 (codon, paths,
  boundaries, closure). Property tests at small seeds: forward then back is the identity where the map is total; a path
  equals the composition of its steps.
- **Tests on the p53 example** (`lang/_test/p53.ts` builds it): its static model, its resolved view; a query over it
  (`/data/p53/(range :from 94 :to 292)`) through `QueryModel`; projection of a protein range onto the gene through the
  codon map (both strands).

**Read first**: `.seq/design/composition.md`, `maps.md`, `building-blocks.md` (the nouns), `locations-and-queries.md`,
`protocol.md` (sources, the store's keys; the store is `seqx-core/data`), `rendering.md` (the resolved view's contract),
`viewer-models.md` (follow), the edge cases named above in `.seq/sources/seqstar/edge-cases.md` (old evidence paths
are guidance only), the `_spec` of `packages/seq/core` (Language, Queries, Data) and of `packages/core/src/repr`,
`.seq/design/testing.md`.

**Model**: strong.

**Base**: `packages-seq` at the commit named in the launch instructions. Branch `seq/e-resolve`; commit on it only.

**Touches**: `packages/seq/core/src/resolve/**`, `packages/seq/core/src/view/**` (new), their `_test/`; the core
`_spec`'s new sections "Resolve" and "Resolved view" only. Read-only: `lang/`, `query/`, `data/` (ask in the report
for changes there). Not in scope: the overlay and dynamic nodes (stage 6), layout and draw lists (lane F), the viewer's
models (lane G), structures and computed maps.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`,
`vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/`, `packages/seq/*/package.json` (ask in the report for
the entry exports and `_test` nulls you need), anything in `.seq/` and in `packages/core/`. Ask in the report for
what you need there.

**Headers**: the Mol\* header, `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more
info.`, with `@author David Sehnal <david.sehnal@gmail.com>`. No mention of `packages/sequences` or `.seq/` in code,
comments or tests.

**Checks**, each its own command, green only on exit 0 with a non-zero test count:
- `pnpm test:unit packages/seq/core` (your folder under 2 s, a test under 100 ms)
- `npx tsc -b packages/seq/core`
- `pnpm lint`
- `pnpm format:check`

**Stop when** the checks are green, or a check has failed twice after a fix, or a permission is refused (report the
action and why; do not work around it).

**Report**: the commits; each check with its exit code and test count; the quick lane's time for your folder; the
exports and `_test` nulls you need in `packages/seq/core/package.json`; any place the designs were unclear or wrong (as
drift: what, where, options, your recommendation); followups.
