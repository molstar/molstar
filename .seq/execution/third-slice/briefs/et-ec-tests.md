# Lane ET: tests for the edge cases built but untested

**Goal**: EO found ECs whose behavior is built but has no test (`ec-outcomes-draft.md`, "not met"). Write one quick
test each, citing the EC in its name, in **new test files** under `packages/seq/core/src/**/_test/`:

- **EC-135**: the ESPript scheme's `gaps: count` against `gaps: ignore` on a column with a gap
  (`core/src/lang/theme-language.ts:138`, `core/src/theme/schemes.ts:275`).
- **EC-116**: a row of letters shows text for all its columns or none; then dots, then bars as density grows (the
  steps in `core/src/render/text.ts:62` and `properties.test.ts`'s table).
- **EC-115**: a stacked blocks section zoomed about a point keeps the block under the pointer (the row's case is in
  `blocks-chips-motion.test.ts`).
- **EC-99**: a text section's ligand chips flow in a line with the words and wrap with them (a small fixture of your
  own; do not import from `packages/seq/examples`).
- **EC-195**: a streamed sequence's default numbering (1-based labels) is known from the spec without its data.
- **EC-138**: the hover outline holds 3:1 against what it outlines, in both themes.
- **EC-63**: two pair kinds at one cell are sliced, each slice its own hit.

Read each EC's text in `.seq/design/edge-cases.md` first. A test that finds the behavior missing or different from the
EC: stop on that one and report it (file:line, the scenario), do not change library code; go on with the others.

**Tests only**: no change outside your new files. RK works in `core/src/render/{rows,tracks,grid,layer-type}.ts` and
the viewer; CS in examples and storybook; AT in `core/src/render/_test/accessible-items.test.ts`; avoid those files.
Each test file quick (under 0.3 s); the core folder's Duration reported before and after.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/et-ec-tests
<that commit>`.

**Checks**: `pnpm test:unit packages/seq/core`, `pnpm -C packages/seq/core typecheck:test`, `pnpm lint`, `pnpm
format:check`.

**Headers**, **hot files** and **stop when**: as every lane's (see [qt-quick-lane.md](qt-quick-lane.md)).

**Report**: commits; checks; each EC's test (file and name) or what it found; the folder's Duration before and after.
