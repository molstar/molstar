# Lane S0b: `values.coarser`'s default to 1

**Goal**: the author's Gate 0 answer (item 10, 2026-11-18): `values.coarser` defaults to **1**, not 2. With the direct
path, level 1 at 100k positions costs about 2 ms, and 1 keeps the picture exact (a bin no wider than a device pixel).
- `packages/seq/core/src/render/config.ts`: the default, and the doc comment's wording of the default;
- the core `_spec` (`packages/seq/core/_spec/README.md`, Render, Configuration and the `values.coarser` bullet): the
  default and its reason (above), the factor still configurable;
- every test that names or depends on the default (search for `coarser`, and the committed draw-list records that a
  values layer at bars density writes: `render/_test/records/`); a record that changes is regenerated and the change
  explained in the report;
- the 3M-bars heavy bounds (`render/_test/frame-time.test.ts`) and the allocation tier still hold at the new default:
  rerun them and report the numbers.

**Model**: fast (sonnet).

**Base**: `packages-seq` at the commit named in the launch instructions. Branch `seq/s0b-coarser`; commit on it only.

**Touches**: `packages/seq/core/src/render/**` and its `_test/`; the core `_spec`. Nothing else.

**Hot files, never staged**: `pnpm-lock.yaml`, every `package.json` and `tsconfig*.json`, `biome.json`,
`vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/`, anything in `.scaffolds/seq/`.

**Headers**: a modified file keeps its header, year 2026, `@author David Sehnal <david.sehnal@gmail.com>`. No mention
of `packages/sequences` or `.scaffolds/seq/`.

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/seq/`
- `pnpm test:heavy packages/seq/`
- `pnpm -C packages/seq/core test:perf` (run in the background; about 10 minutes)
- `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`; `pnpm -C packages/seq/core typecheck:test`
- `pnpm lint`; `pnpm format:check`

**Stop when** green, or a check fails twice after a fix, or a heavy bound fails at the new default (report the numbers:
that is a decision for the author, not a bound to loosen).

**Report**: the commits; each check with its exit code and test count; the 3M-bars frame numbers (median, p95, max)
and the layout bytes at the new default; any record that changed, and why.
