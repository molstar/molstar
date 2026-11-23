# Lane K3: law 3 with stored patches in the tree

**Goal**: make law 3 (an inverse applied with its `inverseNeeds` gives the old root back, patches.md "The inverse and
the allowances") hold when the tree holds stored patches (story steps) whose typed refs are judged over their dry runs.
- The failing case: `packages/core/src/repr/patch/_test/fuzz.test.ts`, the long derived run (`fuzz(3000, 4, true, false)`),
  seed 4, case 43: the inverse refuses with `repr.patch.result` (it gives another revision than it says). Turn its
  `law3` argument on to see it (`pnpm test:heavy packages/core/src/repr/patch/_test/fuzz.test.ts`); the case's JSON is
  printed. Reduce it to a small hand case first and say what goes wrong.
- Fix it in `patch/apply.ts`, `patch/refs.ts` or `patch/work.ts`, in the spirit of the last fix (`ff2ca5476`: a stored
  patch's typed `set` typed again before a move). Then turn `law3` on for the long run and run it with more seeds
  (say 5 seeds × 3,000) to find the rest; fix what it finds, each with a small hand case.
- If a case shows the design itself cannot hold law 3 for stored patches (a ref whose meaning depends on the dry run
  changing under a patch in a way no inverse op can restore), stop and report it with the case and the options.

Not in scope: `merge.ts`, the ref index's derive (only if the fix needs it), anything outside `packages/core/src/repr`.

**Read first**: `.seq/design/patches.md` (all of it: Applying, link integrity, stored patches, the inverse);
`packages/core/src/repr/_spec/README.md` (Patches); `packages/core/src/repr/patch/{apply,refs,work}.ts` and their tests.

**Model**: strong (opus).

**Base**: `packages-seq` at `ff2ca5476`. Your worktree may be cut from master: first run
`git checkout -B seq/k3-stored-law3 ff2ca5476` with a clean tree, and check `git merge-base --is-ancestor ff2ca5476 HEAD`.
Commit on that branch only, at each small green step.

**Touches**: `packages/core/src/repr/patch/**` and its tests, the `repr` `_spec`.

**Performance representative**: the heavy bounds in `patch/_test/{scale,refindex}.test.ts` stay green.

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/core/src/repr packages/seq/`
- `pnpm test:heavy packages/core/src/repr` (twice)
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
