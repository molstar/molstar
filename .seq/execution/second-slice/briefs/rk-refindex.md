# Lane R+K: the ref index across revisions

**Goal**: measure first, then build only what the numbers ask for (plan.md, R+K).
1. **Measure** (a heavy-lane test, kept): a spec of 10,000 bookmark refs (a `refs` container of refs over the p53
   spec, or the conformance language at that size), edited by 100 moves and 100 one-key sets, each applied to the
   previous result. Report the per-apply time split: building the old root's ref index vs the rest.
2. **Derive the next revision's index** (lane K's followup): `applyPatch` returns a root whose `RefIndex` is derived
   from the old one by what the ops wrote, made, moved and removed, cached on the new root the way the old root's is,
   so a chain of applies never rebuilds from scratch. On plain `Map`s first.
3. **Only if** step 2 still leaves the derive itself above budget (16 ms an edit, patch + resolve + redraw; say what
   share the patch may take): put the index on `@molstar/core/data/immutable`'s `HashMap` behind `RefIndex` and
   measure again. Else record "not needed" with the numbers.

The `RefIndex` seam is internal (patches.md, The ref index): no public API change. Any change to what `applyPatch`
returns or takes stops the lane.

**Read first**: `.seq/design/patches.md` (Applying, The ref index), `immutable.md` (In `repr`),
`performance.md`; `packages/core/src/repr/patch/{refs,apply,work}.ts` and the `repr` `_spec`; `packages/core/src/data/immutable`'s
`_spec` if step 3 runs.

**Model**: fast (sonnet).

**Base**: `packages-seq` at `595c6a26d`. Branch `seq/rk-refindex`; commit on it only, at each small green step.

**Touches**: `packages/core/src/repr/patch/{refs,apply}.ts`, their tests, the `repr` `_spec` (The ref index). Lane K2
adds `patch/merge.ts` at the same time: do not touch it or `repr/patch.ts`.

**Two consumers**: the derived index read by link integrity in `applyPatch` and by `diff`'s ref handling (or the
stored-patch dry run), in one test: a derived index equals a fresh build on the same root (by a structural compare,
seeded over the patch fuzz's multi-op patches).

**Performance representative**: one move on the 10,000-ref spec with the old root's index derived (not built): a quick
case, and a heavy bound in `itWithin` (propose it; lane K measured about 0.7 ms with the index cached).

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/core/src/repr packages/seq/`
- `pnpm test:heavy packages/core/src/repr`
- `npx tsc -b packages/core packages/seq/core`
- `pnpm lint`, `pnpm format:check`

**Headers**: a new file has the Mol\* header `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for
more info.` with `@author David Sehnal <david.sehnal@gmail.com>`; a modified file keeps its header, its year range
extended to 2026 and that author line added if missing. No mention of `packages/sequences` or `.seq/` in code, comments,
tests or `_spec`s; code cites its own `_spec` sections.

**Fixes land with a test shown failing on the fix's parent**: the report names the command.

**Hot files, never staged**:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/` (write the changelog entry's text in your report);
- every `package.json` and `tsconfig*.json` (ask in the report: exports, `_test` nulls, dependencies);
- anything in `.seq/`.

**Stop when** the checks are green, or a check has failed twice after a fix, or the work needs a public API past the
agreed surface, or a target cannot be met without changing a design (report the measurement and the options), or a
permission is refused (report the action and why; do not work around it).

**Report**: the commits; each check with its exit code and test count; the hot-file changes you need; the changelog
entry's text; config items with their defaults (and which you were unsure of); drift (what, where, options, your
recommendation); followups.
