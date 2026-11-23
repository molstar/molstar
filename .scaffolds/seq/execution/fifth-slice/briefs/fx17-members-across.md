# Lane FX17: a member selected on a categorical grid lights its matrix row; the dot plot's cells

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx17`. Model:
Sonnet. TT is rebuilding the hover card as providers and owns `lang/`; LA4 is in maps, bridges and `lift.ts`: stay out
of `describe*.ts`, `hover-card.ts`, `lang/`, bridges and maps.

**Item** (`fixes.md`, Open, "A sample selected on the kinship grid does not light its row in the genotype matrix"; read
it whole): in Features › Across sections › Kinship beside the genotypes, a click on a member's label selects it (the
grid's row and column light), but the matrix on the same rows system shows nothing; across-sections.md says a member
selected on the grid lights its row in a matrix layer of the same system with no map. Also: the dot plot (Two genomes as
a dot plot) draws no lines between block-pair cells.

**Task**: (1) a selection (and a hover, through HC's contract) of a member of a rows system lights that member's row in
every layer of that system in view (a matrix's row band, an alignment's row), canvas and SVG; the twin asserts the
matrix row's emphasis primitive after the label click (fails on its parent). (2) the grid over two sections of blocks
draws the block-pair cells' boundaries as a hairline (the grid's own boundary token, as the single-section grid draws
its blocks' gaps), never over data; a test on the drawn lines' positions at the blocks' edges.

**Where**: a worktree on `seq/fx17`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); a commit per part;
hot files never staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.seq/`, `.v6/changelogs/`,
`storybook/packages.ts`, `packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`). No language change (the matrix's
`member_labels` is TT's design's separate change: not here). Anything in `@molstar/core` or `@molstar/io` stops the lane.
Sizes soft: report growth; never edit bounds.

**Checks**: `pnpm test:unit packages/seq/ cli/seqx-render/`, `npx tsc -b packages/seq/core packages/seq/viewer
packages/seq/examples`, each `typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf tiers
to completion, `pnpm -s -C packages/seq/core size` (read the rows), `pnpm test:medium packages/seq/` (HC's sweep), the
`.seq/_checks`. A fix's test fails on its parent.

**Report**: under 200 words: what was fixed, commits, checks with numbers, the `.v6/changelogs/` text.
