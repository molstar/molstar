# Lane FX8: the navigator window's look and its readout, as decided

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx-8`. Model:
Sonnet (a decided design, specified). Follow `fx-fix-lane.md` (a commit per part, a test failing on its parent, a look
in light and dark).

**Build** navigators.md's new section "The window's look: one control row, one stroke, one readout place" (amended
2026-11-22, all its Decided items yes), and close these `fixes.md` items: "The navigator's window and zoom strip look
inconsistent and rough"; "The navigator's readout: one fixed place outside, and its text lagging a gesture" (both
parts: the place, and the text equal to the window drawn, in the same frame). Read the section's Consequences for the
files; the readout's old placement (`viewer/src/readout-place.ts`, FX5's) is removed as the section says.

**Beside you**: RG's branch (not merged) changes `core/src/data/*`, `render/{pairs*,grid*,matrix,layer-data,layer-type,
tracks*,renderer,built-ins,config,index,issues,region-hover}.ts`, `render/layout/{rings,warp}.ts`, `resolve/*`,
`view/view.ts`, `viewer/src/describe*.ts`, both `_spec/README.md`s (yours: edit only the navigator sections, small);
FX7 may touch `render/rail.ts` and the rail's mount (a follower's rail over the cyclic origin). Touch none of these
files; if the design needs one, stop that part and say so.

**Tests**: the section's Defaults as twin tests (two rows over a genome for chr6 and chr13 alike; the fit as one quiet
frame swept by a drag; one stroke weight in the draw list; the readout's box above the part, right-aligned, never
intersecting the window or the pointer, at rest and mid-drag; its text equal to the window's range after a drag and a
click-to-jump); the Configuration rows in the viewer `_spec` table (the config-keys check); both perf tiers (a warm
navigator frame allocates nothing).

**Where**: a worktree on `seq/fx-8`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files never
staged (`.seq/`, `.v6/changelogs/`, `packages/seq/TODO.md`, `package.json`s, `tsconfig*.json`): list edits. Stop any
Storybook server you start.

**Checks** (with numbers): `pnpm test:unit packages/seq/` (under 2 s a touched folder), `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`, each `typecheck:test`, the storybook typecheck, `pnpm lint`,
`pnpm format:check`, both perf tiers, `pnpm -s -C packages/seq/core size`, the `.seq/_checks` (features-planned with
`--plan .seq/execution/fifth-slice/plan.md`). No `pnpm test:heavy`.

**Report**: under 250 words: commits; checks; story ids for the look (BED cytobands, a ring that switches, the rows'
vertical overview, a genome row); what was not built and why; the `.v6/changelogs/` text.
