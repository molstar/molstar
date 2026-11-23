# Lane C: locations, elements and queries

**Goal**: in `packages/seq/core/src/query/`: the Location type (per-system sorted merged ranges in `Int32Array`s, rects
pairing two parts; immutable, comparable, hashable), element keys, the three result types (Nodes, Elements, Location)
and their explicit conversions; the first forms as `repr` form resolvers (`range`, `at`, `all`, `item`, `node`,
`extent`, `items`, and `and`/`or`/`not` over one type) against a minimal interface of systems, numberings and
annotations that a test fixture implements (the static model comes later, in lane E); numberings `index` and
`ordinal` (1-based default) and labels with the rules EC-1 to EC-8; staged results (`complete`, `pending`).

**Read first**: `.scaffolds/seq/design/locations-and-queries.md`, `.scaffolds/seq/design/addressing.md`, the edge cases EC-1 to EC-14 in
`.scaffolds/seq/sources/seqstar/edge-cases.md` (old evidence paths there are only guidance), the `_spec` of
`packages/core/src/repr` (the form-resolver table, `Named`, segments, bounds), `.scaffolds/seq/design/testing.md`.

**Model**: strong.

**Base**: `packages-seq` at the commit named in the launch instructions. Branch `seq/c-locations`; commit on it only.

**Touches**: `packages/seq/core/src/query/**` (new), its `_test/`. Not in scope: the language's kinds (lane B: use your
own small test language on `repr`), the protocol and store (lane D), the static model (lane E).

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`,
`vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/`, `packages/seq/*/package.json` (ask in the report for
the entry exports and `_test` nulls you need), anything in `.scaffolds/seq/` and in `packages/core/`. Ask in the report for
what you need there.

**Headers**: the Mol\* header, `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more
info.`, with `@author David Sehnal <david.sehnal@gmail.com>`. No mention of `packages/sequences` or `.scaffolds/seq/` in code,
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
