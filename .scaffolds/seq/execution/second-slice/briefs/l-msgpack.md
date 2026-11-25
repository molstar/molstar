# Lane L: the MessagePack codec for `repr`, in `@molstar/io`

**Goal**: a `repr` codec (`packages/core/src/repr/codec/codec.ts`'s interface) in `@molstar/io`, built on
`common/msgpack`'s `encodeMsgPack`/`decodeMsgPack` and their extension hooks, with **no change to io's MessagePack
itself** (tree-core.md, base-library.md). Typed columns as extensions `0x11`–`0x18` (and the category column),
decoded without a copy where alignment allows (else one copy, measured). It passes `checkCodec` (the conformance
suite in `repr/codec/conformance.ts`) from its own test. The p53 specs round-trip JSON ↔ MessagePack by hash.

**Performance**: decode throughput as a heavy-lane bound (performance.md: at least 200 MB/s) on a large spec of columns.

Not in scope: changes to `common/msgpack`, Seq\*'s use of the codec beyond the round-trip test.

**Read first**: `.scaffolds/seq/design/tree-core.md` (Codecs), `format.md`, `performance.md`; `packages/core/src/repr/_spec/README.md`
(Codecs); `packages/core/src/repr/codec/`, `codecs/json` (as the model to follow); `packages/io/src/common/msgpack/`.

**Model**: fast (sonnet).

**Base**: `packages-seq` at `adf4a3018`. Branch `seq/l-msgpack`; commit on it only, at each small green step.

**Touches**: a new folder in `packages/io/src/` for the codec (choose by io's layout; say why), its `_spec` note, tests.

**Hot files, never staged**:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/` (write the changelog entry's text in your report);
- every `package.json` and `tsconfig*.json` (ask in the report: exports, `_test` nulls, dependencies);
- anything in `.scaffolds/seq/`.

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/io/src packages/core/src/repr`
- `pnpm test:heavy packages/io/src` (the throughput bound)
- `npx tsc -b packages/core packages/io`
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
