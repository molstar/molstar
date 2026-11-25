# Lane K2: `merge3` in `@molstar/core/repr`

**Goal**: build `.scaffolds/seq/design/merge.md` exactly. Its "The function" section is the agreed surface: `merge3`,
`MergeOptions`, `MergeResult`, `Conflict`, `ConflictKind`, exported from `@molstar/core/repr/patch` (no new entry).
- `diff` both sides from the base, key each change by the node it touches, apply theirs on ours as one patch under the
  options (link integrity and the validators judged on the merged root against ours);
- independence by node (target, parent, refs written, the sibling a `place` names when theirs made it), dependents
  left out with their conflict;
- the conflict table complete for the six ops; `move` reserved (no branch, no test); orders compared as orders;
- placement: the nearest standing predecessor in theirs, else first; never a failure;
- `partial` = ours plus theirs' independent changes, applied and validated as one patch, else ours;
- without a language an info issue says the merge was not validated;
- `changes` on success net from ours to root;
- rebase as the design's one-liner, with a helper only if the design's text needs none past the agreed surface (else
  a test that does it inline).

Not in scope: `SpecModel` and the viewer (lane N), merging written patches (the reserved `move`), the persistent ref
index (lane R+K).

**Read first**: `.scaffolds/seq/design/merge.md`, `patches.md` (all of it: Applying, link integrity, "The inverse and the
allowances"), `tree-core.md`; `packages/core/src/repr/_spec/README.md` and `packages/core/src/repr/patch/`.

**Model**: strong (opus).

**Base**: `packages-seq` at `595c6a26d`. Branch `seq/k2-merge`; commit on it only, at each small green step.

**Touches**: `packages/core/src/repr/patch/merge.ts` (new), its tests in `patch/_test/`, the export lines in
`repr/patch.ts` and the `repr` `_spec`. Lane R+K works in `patch/refs.ts` and `patch/apply.ts` at the same time:
change them only if `merge3` cannot be built without it, and say so in the report.

**Two consumers**: one test where a merge's result is consumed by `applyPatch` (its `changes` and root) and by
`diff` (rebased patch read back), over `SeqLanguage` (`@molstar/seqx-core`'s language; a test in
`packages/seq/core/src/lang/_test/` is fine) and over the conformance language.

**Performance representative**: a merge of two 10-op edits of the p53 spec: a quick case, and a heavy bound in
`itWithin` (propose it from measurement; well under the 16 ms edit target). And one heavy case: two edits of a
20,000-node tree.

**Tests**: the laws as seeded property tests (merge.md, Laws), every conflict kind both ways round, the R2-B2b and
R2-B2c cases (a ref to a node only theirs inserted when ours inserted another with that id; two follows making a
cycle), placement when the predecessor is gone, orders (one side, same, different), `invalid` with `at` and
`partial` = ours, a rebase of a stale patch with a move whose refs follow.

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/core/src/repr packages/seq/`
- `pnpm test:heavy packages/core/src/repr`
- `npx tsc -b packages/core packages/seq/core packages/seq/viewer packages/seq/examples`
- `pnpm -C packages/seq/core typecheck:test`, the same for `viewer` and `examples`
- `pnpm lint`, `pnpm format:check`

**Headers**: a new file has the Mol\* header `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for
more info.` with `@author David Sehnal <david.sehnal@gmail.com>`; a modified file keeps its header, its year range
extended to 2026 and that author line added if missing. No mention of `packages/sequences` or `.scaffolds/seq/` in code, comments,
tests or `_spec`s; code cites its own `_spec` sections.

**Fixes land with a test shown failing on the fix's parent**: the report names the command.

**Hot files, never staged**:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/` (write the changelog entry's text in your report);
- every `package.json` and `tsconfig*.json` (ask in the report: exports, `_test` nulls, dependencies);
- anything in `.scaffolds/seq/`.

**Stop when** the checks are green, or a check has failed twice after a fix, or the work needs a public API past the
agreed surface, or a target cannot be met without changing a design (report the measurement and the options), or a
permission is refused (report the action and why; do not work around it).

**Report**: the commits; each check with its exit code and test count; the hot-file changes you need; the changelog
entry's text; config items with their defaults (and which you were unsure of); drift (what, where, options, your
recommendation); followups.
