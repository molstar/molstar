# Lane TA: the themes' layout allocation

**Goal**: the allocation tier's "p53 on its gene, whole protein" pan and zoom went from 42.6 KB a frame (after SL,
`2457f3ab9`) to 72.4 KB (after TH, `1772436b9`) against a 64 KB bound (`packages/seq/core/src/render/_test/frame.quick.ts`;
run with `pnpm -C packages/seq/core test:perf`). Find what TH's themes allocate per layout in that scenario (likely the
channel pass over the rows, the per-item field tables of `ItemRows`, theme or LUT objects rebuilt per layout, closures,
spreads) and remove it: compile once per theme and tokens' version, reuse buffers, read fields in place. Back under the
bound with room, near the SL number; every other scenario no worse; draw and hover still 0 bytes warm. The bound does not
move.

Not in scope: anything else. Touch only `packages/seq/core/src/{render,theme}/**` where the allocation is, and their
tests. Lanes RB and TV run beside you in `render/` (layouts, rings, text): keep your edits to the themes' channel and
item-reading code, small and local.

**Model**: fast (sonnet). **Base**: `packages-seq` at `63189b52e`: `git checkout -B seq/ta-theme-allocation 63189b52e`, check
`git merge-base --is-ancestor 63189b52e HEAD`.

**Test**: the tier itself is the test (show it failing on the base and passing after: name the command); if a quick-lane
case can hold the fix (a buffer reused, a compile cached), add it.

**Checks**: `pnpm test:unit packages/seq/`; `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`;
`pnpm -C packages/seq/core typecheck:test`; `pnpm lint`; `pnpm format:check`; `pnpm -C packages/seq/core test:perf` (the
whole tier, every scenario's number before and after). No heavy run needed; never one while `pgrep -fl "vitest.*heavy"`
shows one.

**Headers**: a modified file keeps its header, its year range 2026, `@author David Sehnal <david.sehnal@gmail.com>` added
if missing. No mention of `packages/sequences` or `.scaffolds/seq/` in code. **Hot files**: never stage `package.json`,
`tsconfig*`, the lockfile, `.v6/changelogs/`, or anything in `.scaffolds/seq/`.

**Report** (under 300 words): the cause, the commits, the tier's numbers before and after for every scenario, each check
with its exit code and test count.
