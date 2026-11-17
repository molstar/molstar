# Lane C2: half-open range sets into `@molstar/core/data/int`

**Goal**: move the range-set operations of `packages/seq/core/src/query/location.ts` (sorted, merged, half-open
ranges in `Int32Array`s: normalize, union, intersect, subtract, overlap by binary search, the int32 refusal) into a new
`@molstar/core/data/int/half-open-ranges.ts` (the name is yours to propose if a better one fits core's `data/int`
family: `Interval`, `SortedRanges`, `OrderedSet`), with their property tests (moved, at the same small seeds); then
make Location use them, with no change in its behavior or tests. The author approved the move (2026-11-17).

**Read first**: `.seq/design/base-library.md` (Promote), `packages/core/src/data/int/` (its conventions: namespaces,
inclusive vs half-open, `_test`), `packages/seq/core/src/query/location.ts` and its tests, root `AGENTS.md` (reuse,
headers, `_test` null exports).

**Model**: fast (well specified).

**Base**: `packages-seq` at the commit named in the launch instructions. Branch `seq/c2-range-sets`; commit on it only.

**Touches**: `packages/core/src/data/int/**` (the new module and its test), `packages/seq/core/src/query/location.ts`
and `compare.ts` if needed, `packages/core/package.json` (an export only if `data/int` modules need one; check how
`sorted-array` is exported). Not in scope: any other change to Location's API.

**Checks**: as below, plus `pnpm test:unit packages/core/src/data` and `npx tsc -b packages/core`.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`,
`vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/`, `packages/seq/*/package.json` (ask in the report for
the entry exports and `_test` nulls you need), anything in `.seq/`, and in `packages/core/` anything outside `src/data/int/` and its `package.json` export. Ask in the report for
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
