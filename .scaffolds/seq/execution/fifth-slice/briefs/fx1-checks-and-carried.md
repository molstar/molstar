# Lane FX1: the fold-back's checks and the carried items

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx-1`. Model:
fast. The fix lane's rules ([fx-fix-lane.md](fx-fix-lane.md)): one commit an item, a quick-lane test failing on the
parent for each fix, every modified file's header updated (`@author David Sehnal <david.sehnal@gmail.com>`, the year
range to 2026), hot files never staged (`package.json`s, the lockfile, `biome.json`, `.scaffolds/seq/`, `.v6/changelogs/`,
`storybook/packages.ts`, `packages/seq/TODO.md`): list what you need in the report. The checks under `.scaffolds/seq/_checks/`
are hot too: write them in your worktree, commit them on your branch, and the orchestrator reviews them like code (the
one exception to the hot-file rule for this lane, the author's fold-back B6, B7).

## Items, in order

1. **B6, `.scaffolds/seq/_checks/defaults-designed.ts`**: a design in `.scaffolds/seq/design/` whose `## Consequences` names a layer, a
   layout or a part has a `## Defaults` heading (§defaults-designed). Applies to the designs promoted from the fifth
   slice on (genome-coordinates, region-reads, genome-formats, across-sections, ring-extension) by a list in the check
   (older designs exempt until they are touched); `--worklist` prints what fails. Match the other checks' style
   (`.scaffolds/seq/_checks/ec-ids.ts`, `closed-work.ts`: read the git index, exit code, a short message).
2. **B7, `.scaffolds/seq/_checks/features-planned.ts`** (or in `ec-ids.ts`): a `features.md` row whose status says "planned" or
   "to come" for a lane `plan.md` marks merged (`- [x] **<LANE>**`) fails. Read the fifth slice's plan by default
   (`--plan <path>`).
3. **The services' structure section at the fit** (`packages/seq/TODO.md`): on `Examples/Services/Composing Services`
   the 2OCJ chain A section has its own ruler counting 1–219, unaligned with P04637's 94–312. At the least its ruler
   by author numbering (96–289 observed; the unobserved ends in the unobserved gray, ux-foundations section 1); better,
   the B-factor values drawn in P04637's section through the SIFTS map (as letters through a path). Pick the smaller
   change that reads under P04637's residues; say which in the report.
4. **The flake watch list** (`packages/seq/TODO.md`): Prism's `commands.test.ts` and `collapse.test.tsx` (failed once
   each in combined runs). Run the full `pnpm test:unit packages/seq/ packages/ui packages/prism` 6 times; a second
   failure is fixed at its cause (as FX8 did: no retries, no longer waits), else note "6 of 6" in the report.
5. **A story that reads a bigWig** (M4's gap): p53's conservation signal over the gene's region from the bigWig runtime
   (`runtimes/bigwig`, opt-in; its synthetic twin fixture), as `Examples/p53/Signal From A BigWig` or in the gene by
   region story — so a bigWig pan can be measured at M5. With its headless twin.
6. **`../fixes.md`'s open items**: the navigator's position readout (padding, a backdrop, one line with an ellipsis, a
   space after the colon; every navigator form, light and dark) and **letters capped at the label size**
   (`letters.maxPx`, ux-foundations section 2: cap heights equal, about 0.9× the label size for the mono face, regular
   weight by default; every letters path, canvas and SVG; core's Configuration table row; a test that a letter never
   exceeds the label size at the closest zoom).

## Checks before the report

`pnpm test:unit packages/seq/`, `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`,
`pnpm -s -C packages/seq/{core,viewer,examples} typecheck:test`, `pnpm -C packages/seq/storybook typecheck`, `pnpm lint`,
`pnpm format:check`, `node .scaffolds/seq/_checks/{ec-ids,closed-work,file-size}.ts` and your two new checks,
`pnpm -s -C packages/seq/core size` (bounds 230 / 830 / 1330 / 940 kB), and `pnpm -C packages/seq/{core,storybook}
test:perf` (item 6 is on a frame's path, B5). The look: a Storybook in your worktree on port 6016, light and dark, by
clicking and moving the pointer; describe what you saw — the orchestrator looks again.

**Report**: under 300 words: commits per item; checks; tests failing on the parent; public surface added; what is left
and why; the `.v6/changelogs/` entry's text.
