# Lane AT: the rows and pairs layers' accessible items, tested

**Goal** (AY's leftover, `followups.md`): tests for what the rows and pairs layers name to assistive technology and in
the exported SVG's accessible layer (core `_spec`, the accessibility section): a row by its label (a summary row by its
summary label), a pair by its arcs and by its dots (`render/pairs.ts`, `pairItems`), on a track, a ring and a grid; each
test fails if the layer's `items` returned nothing or the wrong label (show it by breaking the hook locally, never
committed). Fixtures: the existing `rows-fixtures.ts` and `pairs-fixtures.ts`; one new test file,
`packages/seq/core/src/render/_test/accessible-items.test.ts`, quick (under 0.3 s).

**Tests only**: no change outside the new file. A test that finds a real bug: stop and report it (file:line, the
scenario), do not fix. RK works in `core/src/render/rows.ts` beside you: if a rows test depends on code RK is changing,
note it.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/at-items <that
commit>`.

**Checks**: `pnpm test:unit packages/seq/core` (and its Duration), `pnpm -C packages/seq/core typecheck:test`, `pnpm
lint`, `pnpm format:check`.

**Headers**, **hot files** and **stop when**: as every lane's (see [qt-quick-lane.md](qt-quick-lane.md)).

**Report**: commits; checks; the file's time; how each test was shown to fail; any bug found.
