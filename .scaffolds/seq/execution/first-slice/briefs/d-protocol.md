# Lane D: the source protocol and the data store

**Goal**: in `packages/seq/core/src/data/`: the protocol's types (`FormatRuntime`, `OpenedSource`, `Coverage`,
`RegionRequest`, `RegionAnswer` with the `letters`, `ranges` and `values` table shapes on `repr`'s typed columns,
dynamic nodes as a type only); the `seq` format (a coverage and request/answer pairs, inline or a JSON file) with its
round trip and a `FormatRuntime` for it; the data store as a `ReactiveModel`: string keys (`<source>|<system>|<level>|<tile>`
and `<source>#<key>`), entry states (pending, loaded, empty, failed), versions, `KeyedState` per entry, requests
deduped, ordered (view, margin, prefetch) and bounded in flight, cancellation with a grace, retries with backoff for
`failed` and none for `not-served`, byte-bounded LRU eviction never evicting what is in view, settle; all tested with a
fake clock and a synthetic source (no network, no real timers).

**Read first**: `.scaffolds/seq/design/protocol.md`, `.scaffolds/seq/design/data-store.md`, `.scaffolds/seq/design/dynamic-nodes.md` (types only now),
the edge cases EC-169 to EC-180 in `.scaffolds/seq/sources/seqstar/edge-cases.md`, `packages/core/src/reactive/_spec/README.md`
(`State`, `KeyedState`, `ReactiveModel`, `batch`), the `_spec` of `packages/core/src/repr` (columns), `.scaffolds/seq/design/testing.md`.

**Model**: strong.

**Base**: `packages-seq` at the commit named in the launch instructions. Branch `seq/d-protocol`; commit on it only.

**Touches**: `packages/seq/core/src/data/**` (new), its `_test/`. Not in scope: request planning across maps and the
overlay (lane E and stage 6), FASTA and bigWig runtimes (stage 6), the language's kinds (lane B).

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
