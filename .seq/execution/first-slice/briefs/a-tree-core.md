# Lane A: the tree core in `@molstar/core/repr`

**Goal**: a codec-independent tree core in `packages/core/src/repr/`, written from meso\*'s tree core with the agreed
delta, enough for the first slice: `repr/tree`, `repr/schema`, `repr/builder`, `repr/codec`, `repr/codecs/json`,
`repr/codecs/sexpr`.

**Read first** (in this order): `.seq/design/tree-core.md`, `.seq/design/language-description.md` (the delta table, the
description's shape, issue codes), `.seq/design/format.md` (typed columns, identity by bytes, NaN, printing elided),
`.seq/design/addressing.md` and `.seq/design/locations-and-queries.md` (the address grammar, atoms, the form-resolver
table), `.seq/design/testing.md` (budgets, properties at small seeds). The source to write from: meso\*'s
`packages/meso/core/src/{tree,schema,builder}` and `packages/meso/core/_spec/{schema,builder}.md` in the molstar-ludo
repository (its path is given in the launch instructions). Root `AGENTS.md` for the repository's rules.

**Model**: strong, high effort (porting with design changes).

**Base**: `packages-seq` at the commit named in the launch instructions. Worktree and branch `seq/a-tree-core`; commit
on that branch only, at each small green step; never merge.

**Touches**:
- `packages/core/src/repr/**` (new): `tree/` (node with `custom`; values with typed columns identified by element type
  and bytes, hash cached on the node; Merkle hashes; addresses with lists, strings and addresses as values, atoms as
  numbers when they parse as JavaScript numbers; the form-resolver table with meso\*'s child selection as the default;
  issues with `severity`; the namespaced kind pattern), `schema/` (meso\*'s param types plus `column` and `table`;
  validation reporting every issue; defaults; normalization on read turning `column` params into typed arrays),
  `builder/` (the engine), `codec/` (the interface and its conformance suite), `codecs/json/` (NaN as `null` inside
  float columns, meso\*'s `$number` elsewhere), `codecs/sexpr/` (reader and printer; columns elided past a short head;
  the reader refuses an elided column); `_spec/README.md` (rules and reasons of what is built); `_test/` per folder.
- `packages/core/package.json`: the new `_test` folders' `null` exports only (`pnpm lint` checks them).
- Not in scope: meso\*'s `locators`, `context`, `units` (units become a hook in the schema, no import), `patch` and
  `descriptors` (later), any Seq\* package, any change to `@molstar/io`.

**Headers**: the Mol\* header, `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more
info.`, with `@author David Sehnal <david.sehnal@gmail.com>` (meso\*'s authors carried where a file is written from
theirs).

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/`, anything in `.seq/`. Ask in the report
for what you need there (the changelog entry text, a dependency).

**Checks**, each its own command, green only on exit 0 with a non-zero test count:
- `pnpm test:unit packages/core/src/repr` (each folder's tests under 2 s, a test under 100 ms)
- `npx tsc -b packages/core`
- `pnpm lint`
- `pnpm format:check`

**Stop when** the checks are green, or a check has failed twice after a fix, or a permission is refused (report the
action and why; don't work around it).

**Report**: the commits; each check with its exit code and test count; the quick lane's time per folder; what you left
out or changed from meso\*'s code and why; any place the designs were unclear or wrong (as drift: what, where, options,
your recommendation); followups; the changelog entry's text for the orchestrator.
