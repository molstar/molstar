# Lane B: the language, first slice

**Goal**: `SeqLanguage` in `packages/seq/core/src/lang/`, described on `@molstar/core/repr`, with the first slice's kinds,
the language builder `createSeqSpec` (handles are addresses), and the issue catalog of `seq.lang.*` codes with its
test; the p53 example (built in a test) printed, encoded in JSON both ways, validated, and each documented misuse giving
its code.

**Read first**: `.seq/design/language-description.md` (the description's shape, the slice's kinds table, the builder by
example, issue codes), `.seq/design/format.md`, `.seq/design/composition.md` (the four levels), the `_spec` of
`packages/core/src/repr` (what the tree core offers: `defineLanguage`, `column` and `table`, `createBuilder` returning
addresses, the codecs), root `AGENTS.md`.

**Model**: fast (well specified).

**Base**: `packages-seq` at the commit named in the launch instructions. Branch `seq/b-language`; commit on it only.

**Touches**: `packages/seq/core/src/lang/**` (new), its `_test/`. Not in scope: queries and locations (lane C), the
protocol and store (lane D), resolving, rendering.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`,
`vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/`, `packages/seq/*/package.json` (ask in the report for
the entry exports and `_test` nulls you need), anything in `.seq/` and in `packages/core/`. Ask in the report for
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
