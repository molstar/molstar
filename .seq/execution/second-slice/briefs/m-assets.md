# Lane M: the asset manager's interface and its conformance suite

**Goal**: the amendment agreed at Gate 0 (`.seq/design/security.md`, What the policy enforces today). The author: "we
make sure the interface has all we need, and then app composition still provides that. we should still have basic
conformance tests for the asset manager that any implementation needs to satisfy."
- `AssetManager.ResolveOptions` (`packages/core/src/util/assets.ts`) gains `credentials`
  (`'omit' | 'same-origin' | 'include'`), `redirect` (`'follow' | 'error'`), `signal` and `maxBytes` (enforced while
  the body streams); `fetchData`'s params (`util/data-source.ts`) gain the same. A manager honors them or refuses the
  read. `DefaultAssetManager` honors them; **its defaults stay what they are today** for every Mol\* app (no option
  given, the same fetch as now). A public API addition to `@molstar/core`, no break.
- `checkAssetManager`: the conformance suite, data and a function (as `repr/codec`'s `checkCodec`), on a fake fetch,
  never the network: the four options honored or the read refused; references and release; one download shared by
  concurrent reads; decompression; a failed download not cached. `DefaultAssetManager` runs it in core's tests.
- Seq\*'s asset policy (`packages/seq/core/src/context` or wherever `policyAssets` lives) passes its choices on every read:
  browser `credentials: 'omit'` unless the composition says otherwise (a policy option), `redirect: 'error'`, the
  read's signal, the size bound; its documentation then claims all four.
- the assets `_spec` (`packages/core/src/util/_spec/assets.md`, a new section) and Seq\*'s security notes in its `_spec`.

Not in scope: range reads (they wait for bigWig), changes to how the cache works.

**Read first**: `.seq/design/security.md`, `ports.md`, `base-library.md`; `packages/core/src/util/assets.ts`,
`data-source.ts`, `util/_spec/assets.md`, `util/_test/assets.test.ts`; Seq\*'s asset policy and its tests.

**Model**: fast (sonnet).

**Base**: `packages-seq` at `adf4a3018`. Branch `seq/m-assets`; commit on it only, at each small green step.

**Touches**: `packages/core/src/util/{assets,data-source}.ts` and a new conformance module beside them, their `_spec` and
tests; Seq\*'s asset policy, its `_spec` and tests.

**Hot files, never staged**:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/` (write the changelog entry's text in your report);
- every `package.json` and `tsconfig*.json` (ask in the report: exports, `_test` nulls, dependencies);
- anything in `.seq/`.

**Two consumers**: `checkAssetManager` runs on `DefaultAssetManager` and on a minimal second implementation in the test
(a derived manager an app could write, adding credentials), so the suite is shown to be implementation-independent.

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/core/src/util packages/seq/`
- `npx tsc -b packages/core packages/seq/core packages/seq/viewer packages/seq/examples`
- `pnpm -C packages/seq/core typecheck:test`
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
