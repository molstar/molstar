# Lane K: revisions, patches and diff in `@molstar/core/repr`

**Goal**: build `.scaffolds/seq/design/patches.md` exactly: its API section is the agreed surface (entry
`@molstar/core/repr/patch`, and in `repr/schema` the `Changes` type, `validateChanges` and `normalizeValue`).
- the six ops with their inverses and `inverseNeeds`; `applyPatch` in the design's eight steps (result then base,
  ops, normalization of what was made, refs following a move, validation by how each object came to be, link
  integrity old root against new, the language's validators, the result check); `Changes` net;
- the core's pass over refs (what a ref names: the node at its id run) and the `RefIndex` seam as written, on plain
  `Map`s (`derive` rebuilds); built only for ops that change an address;
- `diff` with its fallbacks and emission order; `patchLanguage(language)` (memoized, id `<id>.patch`, `at`/`parent`
  plain strings) and `AnyPatchLanguage`; `patchToNode`/`patchFromNode`; stored patches' refs by the stored-patch rule
  over a dry run;
- the child index in `resolveAddress` (a native `Map` per children array, above the flat threshold; immutable.md, In
  `repr`);
- `repr/patch/conformance`: `PatchCases` with options, `checkPatches(apply, cases?)`;
- the companion change in `packages/seq/core/src/lang`: `layer` ids required, `seq.lang.id.required` for every child
  of the spec (not slot bodies), the builder and every example given ids, the language version moved (0.x, no
  migration); every test that builds a layer updated.

Not in scope: `merge3` (lane K2), `SpecModel` and the viewer (lane N), the persistent ref index.

**Read first**: `.scaffolds/seq/design/patches.md` (all of it), `tree-core.md`, `format.md`, `addressing.md`,
`changes.md`, `language-description.md`; `packages/core/src/repr/_spec/README.md` and the code under
`packages/core/src/repr/`; `packages/seq/core/_spec/README.md` (Language) and `packages/seq/core/src/lang/`.

**Model**: strong (opus).

**Base**: `packages-seq` at `adf4a3018`. Branch `seq/k-patches`; commit on it only, at each small green step.

**Touches**: `packages/core/src/repr/**` (new `patch/`, `schema` additions, `tree/address.ts` for the child index) and its
`_spec`; `packages/seq/core/src/lang/**`, the examples' builders (`packages/seq/examples/src/**`), and tests that build
layers across `packages/seq`; `.scaffolds/seq/design/format.md` and `tree-core.md` are the orchestrator's (say in the report
what they should say).

**Hot files, never staged**:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/` (write the changelog entry's text in your report);
- every `package.json` and `tsconfig*.json` (ask in the report: exports, `_test` nulls, dependencies);
- anything in `.scaffolds/seq/`.

**Two consumers**: `checkPatches` run with `PatchCases` and with a case set over `SeqLanguage` (a few cases: a layer
hidden, a track moved with a ref following, a dangling ref refused), so the suite runs over the app's language too.

**Performance representative**: a 10-op patch on the p53 spec, applied and diffed back: a quick case, and its wall-clock
bound (under 1 ms) in `itWithin` (heavy lane).

**Tests**: the conformance cases (every op and inverse at the root, middle, leaves; every refusal by code; the design's
Tests list); seeded property tests of diff's three laws; structural sharing; a refusal leaves the hash; validation
does not re-report an existing error; a kind-changing `replace` validates its kept children.

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

**Stop when** the checks are green, or a check has failed twice after a fix, or the work needs a public API past the
agreed surface, or a target cannot be met without changing a design (report the measurement and the options), or a
permission is refused (report the action and why; do not work around it).

**Report**: the commits; each check with its exit code and test count; the hot-file changes you need; the changelog
entry's text; config items with their defaults (and which you were unsure of); drift (what, where, options, your
recommendation); followups.
