# Lane R: `@molstar/core/data/immutable`

**Goal**: build the persistent collections exactly as `.seq/design/immutable.md` agrees: its API section is the agreed
surface. `HashMap`, `HashSet` (CHAMP, flat below a measured threshold, demoted on delete), `OrderedHashMap`,
`OrderedHashSet`, `Vector` (32-way with a tail), transients (they read like the persistent collections;
`asMutable()` on a transient returns itself; `asImmutable()` returns a new persistent collection and the transient
handle is dead), cursors, per-collection `hash`/`equals`/`valueEquals`/`valueHash` (SameValueZero keys by default),
`defaultHash` with the finalizer, equality and `hashCode` as the design states.

**Also**:
- a **type-check-only spike**: `packages/core/src/state/tree/immutable.ts` and `transient.ts` type-checked against the
  API in a scratch file under `_test/` (no code of the state tree changes); report what does not fit;
- **benchmarks against immutable.js** (already a dev dependency) in the heavy lane: get, set, delete, iterate, build at
  10, 1,000 and 100,000 entries, keys including 22-character strings; ratios; a miss is reported, not hidden;
- the **allocation cases** in the opt-in tier (as `packages/seq/core`'s `test:perf` does: ask for a `test:perf`
  script in core's `package.json` if none fits);
- the **flat threshold measured** (report the numbers that chose it);
- the bundle size of the entry, recorded in the report.

Not in scope: `repr` (lane K builds the child index; the persistent ref index is a later gated step), moving the state
tree or model off immutable.js, the RRB vector.

**Read first**: `.seq/design/immutable.md`, `base-library.md`; the root `AGENTS.md` (tests, lanes, headers, exports);
`packages/core/src/data/` (how `int`, `generic` are laid out: a folder plus an entry file, `data/index.ts`);
`packages/core/src/data/util/hash-functions.ts`; `packages/core/src/state/tree/{immutable,transient}.ts`,
`state/_spec/tree.md`; immutable.js 5's source where the design names it (`node_modules/immutable`).

**Model**: strong (opus).

**Base**: `packages-seq` at `adf4a3018`. Branch `seq/r-immutable`; commit on it only, at each small green step.

**Touches**: `packages/core/src/data/immutable/**`, `packages/core/src/data/immutable.ts`, `packages/core/src/data/index.ts`
(the `Immutable` namespace), a `_spec` for the module (`packages/core/src/data/_spec/immutable.md` or the folder's
convention).

**Hot files, never staged**:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/` (write the changelog entry's text in your report);
- every `package.json` and `tsconfig*.json` (ask in the report: exports, `_test` nulls, dependencies);
- anything in `.seq/`.

**Tests**: model-based property tests against JS `Map`/`Set`/`Array` (with −0 and NaN keys), persistence, the
transient contract (including the colliding-keys delete regression), an all-zero hash, equality and hashCode across
edit orders and configs; quick lane short and seeded, heavy lane long; each folder under 2 s in the quick lane.

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/core/src/data`
- `pnpm test:heavy packages/core/src/data` (the ratios)
- the allocation tier you add
- `npx tsc -b packages/core`
- `pnpm lint`, `pnpm format:check`

**Headers**: a new file has the Mol\* header `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for
more info.` with `@author David Sehnal <david.sehnal@gmail.com>`; a modified file keeps its header, its year range
extended to 2026 and that author line added if missing. No mention of `packages/sequences` or `.seq/` in code, comments,
tests or `_spec`s; code cites its own `_spec` sections.

**Fixes land with a test shown failing on the fix's parent**: the report names the command.

**Stop when** the checks are green, or a check has failed twice after a fix, or the work needs a public API past the
agreed surface, or a target cannot be met without changing a design (report the measurement and the options), or a
permission is refused (report the action and why; do not work around it).

**Report**: the commits; each check with its exit code and test count; the hot-file changes you need; the changelog
entry's text; config items with their defaults (and which you were unsure of); drift (what, where, options, your
recommendation); followups.
